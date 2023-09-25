! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This tint_cross_section module

!> The temperature interpolation cross_section type and related functions
module micm_radXfer_tint_cross_section_type

  use micm_radXfer_abs_cross_section_type,     only : abs_cross_section_t
  use musica_constants,                        only : musica_dk, musica_ik

  implicit none

  private
  public :: tint_cross_section_t

  type cross_section_parms_t
    real(musica_dk), allocatable :: temperature(:)
    real(musica_dk), allocatable :: deltaT(:)
    real(musica_dk), allocatable :: array(:,:)
  end type cross_section_parms_t

  !> Calculator for tint_cross_section
  type, extends(abs_cross_section_t) :: tint_cross_section_t
    !> The cross section array
    type(cross_section_parms_t), allocatable :: cross_section_parms(:)
  contains
    !> Initialize the cross section
    procedure :: initialize
    !> Calculate the cross section
    procedure :: calculate => run
    !> clean up
    final     :: finalize
  end type tint_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize tint_cross_section_t object
  subroutine initialize( this, config, gridWareHouse, ProfileWareHouse )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use netcdf_util,                     only : netcdf_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg
    use micm_grid_warehouse,             only : grid_warehouse_t
    use micm_1d_grid,                    only : base_grid_t
    use micm_Profile_warehouse,     only : Profile_warehouse_t

    !> base cross section type
    class(tint_cross_section_t), intent(inout) :: this
    !> cross section configuration object
    type(config_t), intent(inout) :: config
    !> The warehouses
    type(grid_warehouse_t), intent(inout)         :: gridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse

!   local variables
    integer(musica_ik), parameter :: iONE = 1_musica_ik
    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'tint cross section initialize: '
    character(len=*), parameter :: Hdr = 'cross_section_'

    integer(musica_ik) :: retcode
    integer(musica_ik) :: parmNdx, fileNdx, Ndxl, Ndxu
    integer(musica_ik) :: nParms, nTemps
    real(musica_dk)              :: tmp
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_parameter(:)
    logical :: found, monopos
    character(len=:), allocatable :: msg
    type(netcdf_t), allocatable :: netcdf_obj
    type(string_t), allocatable :: netcdfFiles(:)
    class(base_grid_t), pointer :: lambdaGrid
    type(string_t)     :: Handle

    write(*,*) Iam,'entering'
    !> Get model wavelength grids
    Handle = 'Photolysis, wavelength'
    lambdaGrid => gridWareHouse%get_grid( Handle )

    !> Get cross section netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found=found )

has_netcdf_file: &
    if( found ) then
      allocate( this%cross_section_parms(size(netcdfFiles)) )
file_loop: &
      do fileNdx = 1,size(this%cross_section_parms)
        allocate( netcdf_obj )
    !> read netcdf cross section parameters
        call netcdf_obj%read_netcdf_file( filespec=netcdfFiles(fileNdx)%to_char(), Hdr=Hdr )
        nParms = size(netcdf_obj%parameters,dim=2)
    !> must have at least one parameter
        if( nParms < 2 ) then
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  array must have 2 or more parameters'
          call die_msg( 400000002, msg )
        endif
        associate( Xsection => this%cross_section_parms(fileNdx) )
    !> interpolation temperatures must be in netcdf file
        if( allocated(netcdf_obj%temperature) ) then
          Xsection%temperature = netcdf_obj%temperature
          nTemps = size( Xsection%temperature )
    !> must have two or more interpolation temperatures
          if( nTemps < 2 ) then
            write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  temperature array has < 2 entries'
            call die_msg( 400000003, msg )
          elseif( nTemps < nParms ) then
            write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  temperature array < number parameters'
            call die_msg( 400000004, msg )
          endif
          Xsection%deltaT = Xsection%temperature(2:nParms) - Xsection%temperature(1:nParms-1)
          monopos = all( Xsection%deltaT > rZERO )        
          if( .not. monopos ) then
            if( any( Xsection%deltaT > rZERO ) ) then
              write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  temperature array not monotonic'
              call die_msg( 400000005, msg )
            endif
            do Ndxl = 1,nParms/2
              Ndxu = nParms - Ndxl + 1
              tmp = Xsection%temperature(Ndxl)
              Xsection%temperature(Ndxl) = Xsection%temperature(Ndxu)
              Xsection%temperature(Ndxu) = tmp
              data_parameter = netcdf_obj%parameters(:,Ndxl)
              netcdf_obj%parameters(:,Ndxl) = netcdf_obj%parameters(:,Ndxu)
              netcdf_obj%parameters(:,Ndxu) = data_parameter
            enddo
            Xsection%deltaT = Xsection%temperature(2:nParms) - Xsection%temperature(1:nParms-1)
          endif
        else
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),' does not have interpolation temperatures'
          call die_msg( 400000006, msg )
        endif

    !> interpolate from data to model wavelength grid
        if( allocated(netcdf_obj%wavelength) ) then
          if( .not. allocated(Xsection%array) ) then
            allocate(Xsection%array(lambdaGrid%ncells_,nParms))
          endif
          do parmNdx = 1,nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,parmNdx)
            call this%addpnts( config, data_lambda, data_parameter )
            call inter2(xto=lambdaGrid%edge_, &
                        yto=Xsection%array(:,parmNdx), &
                        xfrom=data_lambda, &
                        yfrom=data_parameter,ierr=retcode)
          enddo
        else
          Xsection%array = netcdf_obj%parameters
        endif
        end associate
        deallocate( netcdf_obj )
      enddo file_loop
    else has_netcdf_file
      write(msg,*) Iam//'must have at least one netcdf input file'
      call die_msg( 400000008, msg )
    endif has_netcdf_file

    write(*,*) Iam,'exiting'

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, gridWareHouse, ProfileWareHouse ) result( cross_section )

    use micm_grid_warehouse,             only : grid_warehouse_t
    use micm_1d_grid,                    only : base_grid_t
    use micm_Profile_warehouse,     only : Profile_warehouse_t
    use micm_Profile,               only : base_profile_t
    use musica_string,                   only : string_t

    !> Arguments
    !> base cross section
    class(tint_cross_section_t), intent(in)       :: this
    !> The warehouses
    type(grid_warehouse_t), intent(inout)         :: gridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse
    !> Calculated cross section
    real(kind=musica_dk), allocatable             :: cross_section(:,:)

    !> Local variables
    character(len=*), parameter :: Iam = 'radXfer tint cross section calculate: '
    integer(musica_ik) :: nTemp
    integer(musica_ik) :: fileNdx, tNdx, layer
    real(musica_dk)    :: Tadj, Tstar
    real(musica_dk), allocatable  :: wrkCrossSection(:,:)
    class(base_grid_t), pointer :: zGrid
    class(base_grid_t), pointer :: lambdaGrid
    class(base_profile_t), pointer :: mdlTemperature
    type(string_t)     :: Handle

    write(*,*) Iam,'entering'

    Handle = 'Vertical Z'
    zGrid => gridWareHouse%get_grid( Handle )
    Handle = 'Photolysis, wavelength'
    lambdaGrid => gridWareHouse%get_grid( Handle )
    Handle = 'Temperature'
    mdlTemperature => ProfileWareHouse%get_Profile( Handle )

    allocate( wrkCrossSection(lambdaGrid%ncells_,zGrid%ncells_) )
    wrkCrossSection = 0.0_musica_dk

    do fileNdx = 1,size(this%cross_section_parms)
      associate( Temp => this%cross_section_parms(fileNdx)%temperature, wrkXsect => this%cross_section_parms(fileNdx) )
      nTemp = size( Temp )
      do layer = 1,zGrid%ncells_
        Tadj   = min( max( mdlTemperature%mid_val_(layer),Temp(1) ),Temp(nTemp) )
        do tNdx = 2,nTemp 
          if( Tadj <= Temp(tNdx) ) then
            exit
          endif
        enddo
        tNdx = tNdx - 1
        Tstar = (Tadj - Temp(tNdx))/wrkXsect%deltaT(tNdx)
        wrkCrossSection(:,layer) = wrkCrossSection(:,layer) + wrkXsect%array(:,tNdx) &
                                 + Tstar * (wrkXsect%array(:,tNdx+1) - wrkXsect%array(:,tNdx))
      enddo
      end associate
    enddo

    cross_section = transpose( wrkCrossSection )

    write(*,*) Iam,'exiting'

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the cross section type
   subroutine finalize( this )

   type(tint_cross_section_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'tint cross section finalize: '
   integer(musica_ik) :: ndx

   write(*,*) Iam,'entering'

   if( allocated(this%cross_section_parms) ) then
     do ndx = 1,size(this%cross_section_parms)
       associate( Xsection => this%cross_section_parms(ndx) )
       if( allocated(Xsection%array) ) then
         deallocate(Xsection%array)
       endif
       if( allocated(Xsection%temperature) ) then
         deallocate(Xsection%temperature)
       endif
       if( allocated(Xsection%deltaT) ) then
         deallocate(Xsection%deltaT)
       endif
       end associate
     enddo
     deallocate(this%cross_section_parms)
   endif

   write(*,*) Iam,'exiting'
   
   end subroutine finalize

end module micm_radXfer_tint_cross_section_type
