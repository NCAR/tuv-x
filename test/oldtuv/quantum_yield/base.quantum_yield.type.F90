! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This base quantum yield module

!> The base quantum yield type and related functions
module micm_base_quantum_yield_type

  use micm_abs_quantum_yield_type,     only : abs_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: base_quantum_yield_t

  integer(musica_ik), parameter :: iONE = 1_musica_ik
  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  type quantum_yield_t
    real(musica_dk), allocatable :: temperature(:)
    real(musica_dk), allocatable :: array(:,:)
  end type quantum_yield_t

  !> Calculator for base quantum yield
  type, extends(abs_quantum_yield_t) :: base_quantum_yield_t
    type(quantum_yield_t), allocatable :: quantum_yield(:)
  contains
    !> Initialize the quantum yield
    procedure :: initialize
    !> Calculate the quantum yield
    procedure :: calculate => run
    !> clean up
    final     :: finalize
  end type base_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize base quantum yield_t object
  subroutine initialize( this, config, mdlLambdaEdge )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use netcdf_util,                     only : netcdf_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg

    !> quantum yield configuration data
    type(config_t), intent(inout) :: config
    !> New base quantum yield calculator
    class(base_quantum_yield_t), intent(inout) :: this
    real(musica_dk), intent(in)                :: mdlLambdaEdge(:)

!   local variables
    character(len=*), parameter :: Iam = 'base quantum yield constructor: '
    character(len=*), parameter :: Hdr = 'quantum_yield_'

    integer(musica_ik) :: retcode, nmdlLambda
    integer(musica_ik) :: colNdx, fileNdx, nCols
    real(musica_dk)    :: quantum_yield_constant
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_parameter(:)
    logical :: found
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
        nCols = size(netcdf_obj%parameters,dim=2)
        if( nCols < 1 ) then
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  parameters array has < 1 column'
          call die_msg( 400000002, msg )
        endif
    !> interpolate from data to model wavelength grid
        if( allocated(netcdf_obj%wavelength) ) then
          if( .not. allocated(this%quantum_yield(fileNdx)%array) ) then
            allocate(this%quantum_yield(fileNdx)%array(nmdlLambda-1,nCols))
          endif
          do colNdx = 1,nCols
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,colNdx)
            call this%addpnts( config, data_lambda, data_parameter )
            call inter2(xto=this%mdl_lambda_edge, &
                        yto=this%quantum_yield(fileNdx)%array(:,colNdx), &
                        xfrom=data_lambda, &
                        yfrom=data_parameter,ierr=retcode)
          enddo
        else
          this%quantum_yield(fileNdx)%array = netcdf_obj%parameters
        endif
        if( allocated(netcdf_obj%temperature) ) then
          this%quantum_yield(fileNdx)%temperature = netcdf_obj%temperature
        endif
        deallocate( netcdf_obj )
      enddo file_loop
    else has_netcdf_file
    !> check for quantum yield constant
      call config%get( 'quantum yield constant', quantum_yield_constant, Iam, found=found )
      if( found ) then
        allocate( this%quantum_yield(1) )
        allocate(this%quantum_yield(1)%array(nmdlLambda-1,1))
        this%quantum_yield(1)%array(:,1) = quantum_yield_constant
      endif
    endif has_netcdf_file

    write(*,*) Iam,'exiting'
   
  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for the environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(base_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                    :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'base quantum yield calculate: '

    write(*,*) Iam,'entering'

    quantum_yield = this%quantum_yield(1)%array(:,1)

    write(*,*) Iam,'exiting'

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the quantum yield type
   subroutine finalize( this )

   type(base_quantum_yield_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'base quantum yield finalize: '
   integer(musica_dk) :: ndx

   write(*,*) Iam,'entering'

   if( allocated(this%quantum_yield) ) then
     do ndx = 1,size(this%quantum_yield)
       if( allocated(this%quantum_yield(ndx)%array ) ) then
         deallocate(this%quantum_yield(ndx)%array )
       endif
       if( allocated(this%quantum_yield(ndx)%temperature ) ) then
         deallocate(this%quantum_yield(ndx)%temperature )
       endif
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

end module micm_base_quantum_yield_type
