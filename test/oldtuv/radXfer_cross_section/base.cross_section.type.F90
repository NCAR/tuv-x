! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This base_cross_section module

!> The base_cross_section type and related functions
module micm_radXfer_base_cross_section_type

  use micm_radXfer_abs_cross_section_type, only : abs_cross_section_t
  use musica_constants,                    only : musica_dk, musica_ik

  implicit none

  private
  public :: base_cross_section_t

  type cross_section_parms_t
    real(musica_dk), allocatable :: temperature(:)
    real(musica_dk), allocatable :: array(:,:)
  end type cross_section_parms_t

  !> Calculator for base_cross_section
  type, extends(abs_cross_section_t) :: base_cross_section_t
    !> The cross section array
    type(cross_section_parms_t), allocatable :: cross_section_parms(:)
  contains
    !> Initialize the cross section
    procedure :: initialize
    !> Calculate the cross section
    procedure :: calculate => run
    !> clean up
    final     :: finalize
  end type base_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize base_cross_section_t object
  subroutine initialize( this, config, gridWareHouse, ProfileWareHouse )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use netcdf_util,                     only : netcdf_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg
    use micm_grid_warehouse,         only : grid_warehouse_t
    use micm_1d_grid,                only : base_grid_t
    use micm_Profile_warehouse, only : Profile_warehouse_t
    use micm_Profile,           only : base_profile_t

    !> base cross section type
    class(base_cross_section_t), intent(inout) :: this
    !> cross section configuration object
    type(config_t), intent(inout) :: config
    !> The warehouses
    type(grid_warehouse_t), intent(inout)         :: gridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse

!   local variables
    integer(musica_ik), parameter :: iONE = 1_musica_ik
    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'base cross section initialize: '
    character(len=*), parameter :: Hdr = 'cross_section_'

    integer(musica_ik) :: retcode
    integer(musica_ik) :: parmNdx, fileNdx
    integer(musica_ik) :: nParms
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_parameter(:)
    logical :: found
    character(len=:), allocatable :: msg
    type(netcdf_t), allocatable :: netcdf_obj
    type(string_t)              :: Handle
    type(string_t), allocatable :: netcdfFiles(:)
    class(base_grid_t), pointer :: lambdaGrid

    write(*,*) Iam,'entering'
    !> Get model wavelength grids
    Handle = 'Photolysis, wavelength'
    lambdaGrid => gridWareHouse%get_grid( Handle )

    !> get cross section netcdf filespec
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
        if( nParms < 1 ) then
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  parameters array has < 1 column'
          call die_msg( 400000002, msg )
        endif

    !> interpolate from data to model wavelength grid
        if( allocated(netcdf_obj%wavelength) ) then
          if( .not. allocated(this%cross_section_parms(fileNdx)%array) ) then
            allocate(this%cross_section_parms(fileNdx)%array(lambdaGrid%ncells_,nParms))
          endif
          do parmNdx = 1,nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,parmNdx)
            call this%addpnts( config, data_lambda, data_parameter )
            call inter2(xto=lambdaGrid%edge_, &
                        yto=this%cross_section_parms(fileNdx)%array(:,parmNdx), &
                        xfrom=data_lambda, &
                        yfrom=data_parameter,ierr=retcode)
          enddo
        else
          this%cross_section_parms(fileNdx)%array = netcdf_obj%parameters
        endif
        if( allocated(netcdf_obj%temperature) ) then
          this%cross_section_parms(fileNdx)%temperature = netcdf_obj%temperature
        endif
        deallocate( netcdf_obj )
      enddo file_loop
    endif has_netcdf_file

    write(*,*) Iam,'exiting'

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, gridWareHouse, ProfileWareHouse ) result( cross_section )

    use micm_grid_warehouse,         only : grid_warehouse_t
    use micm_1d_grid,                only : base_grid_t
    use micm_Profile_warehouse, only : Profile_warehouse_t
    use musica_string,               only : string_t

    !> base cross section
    class(base_cross_section_t), intent(in)    :: this
    !> The warehouses
    type(grid_warehouse_t), intent(inout)         :: gridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse
    !> Calculated cross section
    real(kind=musica_dk), allocatable          :: cross_section(:,:)

    !> Local variables
    integer :: colndx
    character(len=*), parameter :: Iam = 'radXfer base cross section calculate: '
    class(base_grid_t), pointer :: zGrid
    type(string_t)                :: Handle
    real(musica_dk), allocatable  :: wrkCrossSection(:,:)

    write(*,*) Iam,'entering'

    Handle = 'Vertical Z'
    zGrid => gridWareHouse%get_grid( Handle )

    allocate( wrkCrossSection(size(this%cross_section_parms(1)%array,dim=1),zGrid%ncells_) )

    !> Just copy the lambda interpolated array
    do colndx = 1,zGrid%ncells_
      wrkCrossSection(:,colndx) = this%cross_section_parms(1)%array(:,1)
    enddo
 
    cross_section = transpose( wrkCrossSection )

    write(*,*) Iam,'exiting'

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the cross section type
   subroutine finalize( this )

   type(base_cross_section_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'base cross section finalize: '
   integer(musica_ik) :: ndx

   write(*,*) Iam,'entering'

   if( allocated(this%cross_section_parms) ) then
     do ndx = 1,size(this%cross_section_parms)
       if( allocated(this%cross_section_parms(ndx)%array ) ) then
         deallocate(this%cross_section_parms(ndx)%array )
       endif
       if( allocated(this%cross_section_parms(ndx)%temperature ) ) then
         deallocate(this%cross_section_parms(ndx)%temperature )
       endif
     enddo
     deallocate(this%cross_section_parms)
   endif

   write(*,*) Iam,'exiting'
   
   end subroutine finalize

end module micm_radXfer_base_cross_section_type
