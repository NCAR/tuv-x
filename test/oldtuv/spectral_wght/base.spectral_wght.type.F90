! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This base_spectral_wght module

!> The spectral_wght type and related functions
module micm_base_spectral_wght_type

  use micm_abs_spectral_wght_type,     only : abs_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: base_spectral_wght_t

  type spectral_wght_t
    real(musica_dk), allocatable :: temperature(:)
    real(musica_dk), allocatable :: array(:,:)
  end type spectral_wght_t

  !> Calculator for base_spectral_wght
  type, extends(abs_spectral_wght_t) :: base_spectral_wght_t
    !> The spectral wght array
    type(spectral_wght_t), allocatable :: spectral_wght(:)
  contains
    !> Initialize the spectral wght
    procedure :: initialize
    !> Calculate the spectral wght
    procedure :: calculate => run
    !> clean up
    final     :: finalize
  end type base_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize base_spectral_wght_t object
  subroutine initialize( this, config, mdlLambdaEdge )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use netcdf_util,                     only : netcdf_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg

    !> base spectral wght type
    class(base_spectral_wght_t), intent(inout) :: this
    real(musica_dk), intent(in)   :: mdlLambdaEdge(:)
    !> spectral wght configuration object
    type(config_t), intent(inout) :: config

!   local variables
    integer(musica_ik), parameter :: iONE = 1_musica_ik
    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'base spectral wght section initialize: '
    character(len=*), parameter :: Hdr = 'spectral_weight_'

    integer(musica_ik) :: retcode, nmdlLambda
    integer(musica_ik) :: parmNdx, fileNdx
    integer(musica_ik) :: nParms
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
    !> get spectral wght netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found=found )

has_netcdf_file: &
    if( found ) then
      allocate( this%spectral_wght(size(netcdfFiles)) )
file_loop: &
      do fileNdx = 1,size(this%spectral_wght)
        allocate( netcdf_obj )
    !> read netcdf spectral wght parameters
        call netcdf_obj%read_netcdf_file( filespec=netcdfFiles(fileNdx)%to_char(), Hdr=Hdr )
        nParms = size(netcdf_obj%parameters,dim=2)
        if( nParms < 1 ) then
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  parameters array has < 1 column'
          call die_msg( 400000002, msg )
        endif

    !> interpolate from data to model wavelength grid
        if( allocated(netcdf_obj%wavelength) ) then
          if( .not. allocated(this%spectral_wght(fileNdx)%array) ) then
            allocate(this%spectral_wght(fileNdx)%array(nmdlLambda-1,nParms))
          endif
          do parmNdx = 1,nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,parmNdx)
            call this%addpnts( config, data_lambda, data_parameter )
            call inter2(xto=this%mdl_lambda_edge, &
                        yto=this%spectral_wght(fileNdx)%array(:,parmNdx), &
                        xfrom=data_lambda, &
                        yfrom=data_parameter,ierr=retcode)
          enddo
        else
          this%spectral_wght(fileNdx)%array = netcdf_obj%parameters
        endif
        if( allocated(netcdf_obj%temperature) ) then
          this%spectral_wght(fileNdx)%temperature = netcdf_obj%temperature
        endif
        deallocate( netcdf_obj )
      enddo file_loop
    endif has_netcdf_file

    write(*,*) Iam,'exiting'

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral wght for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> base spectral wght
    class(base_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral wght
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'base spectral wght calculate: '

    write(*,*) Iam,'entering'

    spectral_wght = this%spectral_wght(1)%array(:,1)

    write(*,*) Iam,'exiting'

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> finalize the spectral wght type
   subroutine finalize( this )

   type(base_spectral_wght_t), intent(inout) :: this

   character(len=*), parameter :: Iam = 'base spectral wght finalize: '
   integer(musica_ik) :: ndx

   write(*,*) Iam,'entering'
   if( allocated(this%spectral_wght) ) then
     do ndx = 1,size(this%spectral_wght)
       if( allocated(this%spectral_wght(ndx)%array ) ) then
         deallocate(this%spectral_wght(ndx)%array )
       endif
       if( allocated(this%spectral_wght(ndx)%temperature ) ) then
         deallocate(this%spectral_wght(ndx)%temperature )
       endif
     enddo
     deallocate(this%spectral_wght)
   endif
   if( allocated(this%mdl_lambda_edge) ) then
     deallocate(this%mdl_lambda_edge)
   endif
   if( allocated(this%mdl_lambda_center) ) then
     deallocate(this%mdl_lambda_center)
   endif
   write(*,*) Iam,'exiting'
   
   end subroutine finalize

end module micm_base_spectral_wght_type
