! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This tint quantum yield module

!> The tint quantum yield type and related functions
module micm_tint_quantum_yield_type

  use micm_abs_quantum_yield_type,     only : abs_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: tint_quantum_yield_t

  integer(musica_ik), parameter :: iONE = 1_musica_ik
  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  type quantum_yield_t
    real(musica_dk), allocatable :: temperature(:)
    real(musica_dk), allocatable :: deltaT(:)
    real(musica_dk), allocatable :: array(:,:)
  end type quantum_yield_t

  !> Calculator for tint quantum yield
  type, extends(abs_quantum_yield_t) :: tint_quantum_yield_t
    type(quantum_yield_t), allocatable :: quantum_yield(:)
  contains
    !> Initialize the quantum yield
    procedure :: initialize
    !> Calculate the quantum yield
    procedure :: calculate => run
    !> clean up
    final     :: finalize
  end type tint_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize tint quantum yield_t object
  subroutine initialize( this, config, mdlLambdaEdge )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use netcdf_util,                     only : netcdf_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg

    !> quantum yield configuration data
    type(config_t), intent(inout) :: config
    !> New tint quantum yield calculator
    class(tint_quantum_yield_t), intent(inout) :: this
    real(musica_dk), intent(in)                :: mdlLambdaEdge(:)

!   local variables
    character(len=*), parameter :: Iam = 'tint quantum yield constructor: '
    character(len=*), parameter :: Hdr = 'quantum_yield_'

    integer(musica_ik) :: retcode, nmdlLambda
    integer(musica_ik) :: nTemps, nParms
    integer(musica_ik) :: parmNdx, fileNdx, Ndxl, Ndxu
    real(musica_dk)    :: tmp
    real(musica_dk)    :: quantum_yield_constant
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_parameter(:)
    logical :: found, monopos
    character(len=:), allocatable :: msg
    type(netcdf_t), allocatable :: netcdf_obj
    type(string_t), allocatable :: netcdfFiles(:)

    write(*,*) Iam,'entering'

    !> set model wavelength array
    this%mdl_lambda_edge = mdlLambdaEdge
    nmdlLambda = size( this%mdl_lambda_edge )
    this%mdl_lambda_center = .5_musica_dk*(this%mdl_lambda_edge(2:) + this%mdl_lambda_edge(1:nmdlLambda-1))

    !> get quantum yield netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found=found )
has_netcdf_file: &
    if( found ) then
      allocate( this%quantum_yield(size(netcdfFiles)) )
file_loop: &
      do fileNdx = 1,size(netcdfFiles) 
        allocate( netcdf_obj )
    !> read netcdf file quantum yield data
        call netcdf_obj%read_netcdf_file( filespec=netcdfFiles(fileNdx)%to_char(), Hdr=Hdr )
        nParms = size(netcdf_obj%parameters,dim=2)
        if( nParms < 2 ) then
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  array must have 2 or more parameters'
          call die_msg( 400000002, msg )
        endif
        associate( Qyield => this%quantum_yield(fileNdx) )
    !> interpolation temperatures must be in netcdf file
        if( allocated(netcdf_obj%temperature) ) then
          Qyield%temperature = netcdf_obj%temperature
          nTemps = size( Qyield%temperature )
    !> must have two or more interpolation temperatures
          if( nTemps < 2 ) then
            write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  temperature array has < 2 entries'
            call die_msg( 400000003, msg )
          elseif( nTemps < nParms ) then
            write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  temperature array < number parameters'
            call die_msg( 400000004, msg )
          endif
          Qyield%deltaT = Qyield%temperature(2:nParms) - Qyield%temperature(1:nParms-1)
          monopos = all( Qyield%deltaT > rZERO )        
          if( .not. monopos ) then
            if( any( Qyield%deltaT > rZERO ) ) then
              write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  temperature array not monotonic'
              call die_msg( 400000005, msg )
            endif
            do Ndxl = 1,nParms/2
              Ndxu = nParms - Ndxl + 1
              tmp = Qyield%temperature(Ndxl)
              Qyield%temperature(Ndxl) = Qyield%temperature(Ndxu)
              Qyield%temperature(Ndxu) = tmp
              data_parameter = netcdf_obj%parameters(:,Ndxl)
              netcdf_obj%parameters(:,Ndxl) = netcdf_obj%parameters(:,Ndxu)
              netcdf_obj%parameters(:,Ndxu) = data_parameter
            enddo
            Qyield%deltaT = Qyield%temperature(2:nParms) - Qyield%temperature(1:nParms-1)
          endif
        else
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),' does not have interpolation temperatures'
          call die_msg( 400000006, msg )
        endif
    !> interpolate from data to model wavelength grid
        if( allocated(netcdf_obj%wavelength) ) then
          if( .not. allocated(this%quantum_yield(fileNdx)%array) ) then
            allocate(this%quantum_yield(fileNdx)%array(nmdlLambda-1,nParms))
          endif
          do parmNdx = 1,nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,parmNdx)
            call this%addpnts( config, data_lambda, data_parameter )
            call inter2(xto=this%mdl_lambda_edge, &
                        yto=this%quantum_yield(fileNdx)%array(:,parmNdx), &
                        xfrom=data_lambda, &
                        yfrom=data_parameter,ierr=retcode)
          enddo
        else
          this%quantum_yield(fileNdx)%array = netcdf_obj%parameters
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

  !> Calculate the photorate quantum yield for the environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> tint quantum yield
    class(tint_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                    :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'tint quantum yield calculate: '
    integer(musica_ik) :: nTemp
    integer(musica_ik) :: fileNdx, tNdx
    real(musica_dk)    :: Tadj, Tstar

    write(*,*) Iam,'entering'

    quantum_yield = 0.0_musica_dk
    do fileNdx = 1,size(this%quantum_yield)
      associate( Temp => this%quantum_yield(fileNdx)%temperature, wrkQyield => this%quantum_yield(fileNdx) )
      nTemp = size( Temp )
      Tadj   = min( max( environment%temperature,Temp(1) ),Temp(nTemp) )
      do tNdx = 2,nTemp 
        if( Tadj <= Temp(tNdx) ) then
          exit
        endif
      enddo
      tNdx = tNdx - 1
      Tstar = (Tadj - Temp(tNdx))/wrkQyield%deltaT(tNdx)
      quantum_yield = quantum_yield + wrkQyield%array(:,tNdx) &
                    + Tstar * (wrkQyield%array(:,tNdx+1) - wrkQyield%array(:,tNdx))
      end associate
    enddo

    write(*,*) Iam,'exiting'

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the quantum yield type

   subroutine finalize( this )

   type(tint_quantum_yield_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'tint quantum yield finalize: '
   integer(musica_ik) :: ndx

   write(*,*) Iam,'entering'
   if( allocated(this%quantum_yield) ) then
     do ndx = 1,size(this%quantum_yield)
       associate( Qyield => this%quantum_yield(ndx) )
       if( allocated(Qyield%array) ) then
         deallocate(Qyield%array)
       endif
       if( allocated(Qyield%temperature) ) then
         deallocate(Qyield%temperature)
       endif
       if( allocated(Qyield%deltaT) ) then
         deallocate(Qyield%deltaT)
       endif
       end associate
     enddo
     deallocate(this%quantum_yield)
   endif
   if( allocated(this%mdl_lambda_edge) ) then
     deallocate(this%mdl_lambda_edge)
   endif
   if( allocated(this%mdl_lambda_center) ) then
     deallocate(this%mdl_lambda_center)
   endif

   write(*,*) Iam,'exiting'

   end subroutine finalize

end module micm_tint_quantum_yield_type
