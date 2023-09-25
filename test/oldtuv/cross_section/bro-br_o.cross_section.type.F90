! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This bro+hv->br_o cross_section module

!> The bro+hv->br+o_cross_section type and related functions
module micm_bro_br_o_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: bro_br_o_cross_section_t

  !> Calculator for base_cross_section
  type, extends(base_cross_section_t) :: bro_br_o_cross_section_t
  contains
    !> Initialize the cross section
    procedure :: initialize
  end type bro_br_o_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize bro_br_o_cross_section_t object
  subroutine initialize( this, config, mdlLambdaEdge )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use netcdf_util,                     only : netcdf_t
    use photo_utils,                     only : inter4
    use musica_assert,                   only : die_msg

    !> base cross section type
    class(bro_br_o_cross_section_t), intent(inout) :: this
    !> cross section configuration object
    type(config_t), intent(inout) :: config
    real(musica_dk), intent(in)   :: mdlLambdaEdge(:)

!   local variables
    integer(musica_ik), parameter :: iONE = 1_musica_ik
    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'bro->br+o cross section initialize: '
    character(len=*), parameter :: Hdr = 'cross_section_'

    integer(musica_ik) :: retcode, nmdlLambda
    integer(musica_ik) :: colNdx, fileNdx, nCols
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
    !> get cross section netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found=found )

has_netcdf_file: &
    if( found ) then
      allocate( this%cross_section(size(netcdfFiles)) )
file_loop: &
      do fileNdx = 1,size(this%cross_section)
        allocate( netcdf_obj )
    !> read netcdf cross section parameters
        call netcdf_obj%read_netcdf_file( filespec=netcdfFiles(fileNdx)%to_char(), Hdr=Hdr )
        nCols = size(netcdf_obj%parameters,dim=2)
        if( nCols < 1 ) then
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  parameters array has < 1 column'
          call die_msg( 400000002, msg )
        endif

    !> interpolate from data to model wavelength grid
        if( allocated(netcdf_obj%wavelength) ) then
          if( .not. allocated(this%cross_section(fileNdx)%array) ) then
            allocate(this%cross_section(fileNdx)%array(nmdlLambda-1,nCols))
          endif
          do colNdx = 1,nCols
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,colNdx)
            call inter4(xto=this%mdl_lambda_edge, &
                        yto=this%cross_section(fileNdx)%array(:,colNdx), &
                        xfrom=data_lambda, &
                        yfrom=data_parameter,Foldin=1)
          enddo
        else
          this%cross_section(fileNdx)%array = netcdf_obj%parameters
        endif
        if( allocated(netcdf_obj%temperature) ) then
          this%cross_section(fileNdx)%temperature = netcdf_obj%temperature
        endif
        deallocate( netcdf_obj )
      enddo file_loop
    endif has_netcdf_file

    write(*,*) Iam,'exiting'

  end subroutine initialize

end module micm_bro_br_o_cross_section_type
