! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This rono2 cross_section module

!> The rono2_cross_section type and related functions
module micm_rono2_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: rono2_cross_section_t

  !> Calculator for rono2 cross section
  type, extends(base_cross_section_t) :: rono2_cross_section_t
  contains
    !> Initialize the cross section
    procedure :: initialize
    !> Calculate the cross section
    procedure :: calculate => run
  end type rono2_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize base_cross_section_t object
  subroutine initialize( this, config, mdlLambdaEdge )

    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use netcdf_util,                     only : netcdf_t
    use photo_utils,                     only : inter2
    use musica_assert,                   only : die_msg

    !> base cross section type
    class(rono2_cross_section_t), intent(inout) :: this
    real(musica_dk), intent(in)   :: mdlLambdaEdge(:)
    !> cross section configuration object
    type(config_t), intent(inout) :: config

!   local variables
    integer(musica_ik), parameter :: iONE = 1_musica_ik
    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'rono2 cross section initialize: '
    character(len=*), parameter :: Hdr = 'cross_section_'

    integer(musica_ik) :: retcode, nmdlLambda
    integer(musica_ik) :: parmNdx, fileNdx, nParms
    real(musica_dk), allocatable :: data_lambda(:)
    real(musica_dk), allocatable :: data_parameter(:)
    logical :: found
    character(len=:), allocatable :: msg
    character(len=:), allocatable :: addpntKey
    type(netcdf_t), allocatable :: netcdf_obj
    type(string_t), allocatable :: netcdfFiles(:)
    type(config_t)              :: tmp_config
    type(string_t)              :: addpntVal

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
        nParms = size(netcdf_obj%parameters,dim=2)
        if( nParms < 1 ) then
          write(msg,*) Iam//'File: ',trim(netcdfFiles(fileNdx)%to_char()),'  parameters array has < 1 column'
          call die_msg( 400000002, msg )
        endif

    !> interpolate from data to model wavelength grid
        if( allocated(netcdf_obj%wavelength) ) then
          if( .not. allocated(this%cross_section(fileNdx)%array) ) then
            allocate(this%cross_section(fileNdx)%array(nmdlLambda-1,nParms))
          endif
          do parmNdx = 1,nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters(:,parmNdx)
            if( parmNdx == 1 ) then
              call this%addpnts( config, data_lambda, data_parameter )
            elseif( parmNdx == 2 ) then
              tmp_config = config
              addpntKey = 'lower extrapolation'
              addpntVal = 'boundary'
              call tmp_config%add( addpntKey, addpntVal, Iam )
              addpntKey = 'upper extrapolation'
              addpntVal = 'boundary'
              call tmp_config%add( addpntKey, addpntVal, Iam )
              call this%addpnts( tmp_config, data_lambda, data_parameter )
            endif
            call inter2(xto=this%mdl_lambda_edge, &
                        yto=this%cross_section(fileNdx)%array(:,parmNdx), &
                        xfrom=data_lambda, &
                        yfrom=data_parameter,ierr=retcode)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use musica_constants,                only : musica_dk

    !> rono2 cross section
    class(rono2_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)               :: environment
    !> Calculated cross section
    real(kind=musica_dk)                           :: cross_section(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'rono2 cross section calculate: '
    real(musica_dk) :: Temp

    write(*,*) Iam,'entering'

    Temp = environment%temperature - 298._musica_dk 
    cross_section = this%cross_section(1)%array(:,1)*exp( this%cross_section(1)%array(:,2)*Temp )

    write(*,*) Iam,'exiting'

  end function run

end module micm_rono2_cross_section_type
