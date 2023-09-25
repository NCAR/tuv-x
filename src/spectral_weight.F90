! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_spectral_weight
  ! The spectral_weight type and related functions

  use musica_constants,       only : dk => musica_dk
  use musica_config,          only : config_t
  use tuvx_grid_warehouse,    only : grid_warehouse_ptr

  implicit none

  private
  public :: spectral_weight_t, spectral_weight_parms_t, spectral_weight_ptr,  &
            base_constructor

  type spectral_weight_parms_t
    real(dk), allocatable :: temperature(:) ! Temperature grid [K]
    real(dk), allocatable :: array(:,:)     ! Spectral weight parameters (wavelength, parameter type)
  end type spectral_weight_parms_t

  type :: spectral_weight_t
    ! Calculator for spectral weights
    type(spectral_weight_parms_t), allocatable :: spectral_weight_parms(:)
    ! Wavelength grid pointer
    type(grid_warehouse_ptr) :: wavelength_grid_
  contains
    procedure :: calculate => run
    procedure, private :: add_points
    final     :: finalize
    procedure :: pack_size
    procedure :: mpi_pack
    procedure :: mpi_unpack
  end type spectral_weight_t

  !> Pointer type for building sets of dose rate constants
  type :: spectral_weight_ptr
    class(spectral_weight_t), pointer :: val_ => null( )
  end type spectral_weight_ptr

  interface spectral_weight_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( new_spectral_weight )
    ! Instantiate the base spectral weight type

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_t),  pointer       :: new_spectral_weight  ! New :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    type(config_t),            intent(inout) :: config ! Spectral weight configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`


    ! Local variables
    type(string_t) :: required_keys(1), optional_keys(4)

    required_keys(1) = 'type'
    optional_keys(1) = 'netcdf files'
    optional_keys(2) = 'lower extrapolation'
    optional_keys(3) = 'upper extrapolation'
    optional_keys(4) = 'name'
    call assert_msg( 124969901,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "base spectral wght." )

    allocate( new_spectral_weight )
    call base_constructor( new_spectral_weight, config, grid_warehouse,       &
                           profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine base_constructor( this, config, grid_warehouse,                  &
      profile_warehouse )
      ! The base constructor, used by any class that doesn't need special
      ! initialization

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use tuvx_interpolate,              only : interpolator_conserving_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_t),  pointer       :: this   ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    type(config_t),            intent(inout) :: config ! Spectral weight configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! Local variables
    real(dk), parameter    :: rZERO = 0.0_dk
    real(dk), parameter    :: rONE  = 1.0_dk
    character(len=*), parameter    :: Iam = 'base spectral weight initialize: '
    character(len=*), parameter    :: Hdr = 'spectral_weight_'

    integer :: parmNdx, fileNdx
    integer :: nParms
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    logical :: found
    type(netcdf_t), allocatable :: netcdf_obj
    type(string_t), allocatable :: netcdfFiles(:)
    class(grid_t), pointer      :: lambdaGrid
    type(interpolator_conserving_t) :: interpolator

    ! Get model wavelength grid
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    ! Get spectral wght netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found = found )

    has_netcdf_file: if( found ) then
      allocate( this%spectral_weight_parms( size( netcdfFiles ) ) )
      file_loop: do fileNdx = 1, size( this%spectral_weight_parms )
        allocate( netcdf_obj )
        ! Read netcdf spectral wght parameters
        call netcdf_obj%read_netcdf_file(                                     &
                             file_path = netcdfFiles( fileNdx )%to_char( ),   &
                             variable_name = Hdr )
        nParms = size( netcdf_obj%parameters, dim = 2 )
        call assert_msg( 222975865, nParms > 0, Iam//'File: ' //              &
                         trim(netcdfFiles( fileNdx )%to_char( ) ) //          &
                         '  parameters array has < 1 column' )

        ! Interpolate from data to model wavelength grid
        if( allocated( netcdf_obj%wavelength ) ) then
          if( .not. allocated( this%spectral_weight_parms( fileNdx )%array ) )&
              then
            allocate( this%spectral_weight_parms( fileNdx                     &
                                       )%array( lambdaGrid%ncells_, nParms ) )
          endif
          do parmNdx = 1, nParms
            data_lambda    = netcdf_obj%wavelength
            data_parameter = netcdf_obj%parameters( :, parmNdx )
            call this%add_points( config, data_lambda, data_parameter )
            this%spectral_weight_parms( fileNdx )%array( :, parmNdx ) =       &
                interpolator%interpolate( x_target = lambdaGrid%edge_,        &
                                          x_source = data_lambda,             &
                                          y_source = data_parameter )
          enddo
        else
          this%spectral_weight_parms( fileNdx )%array = netcdf_obj%parameters
        endif
        if( allocated( netcdf_obj%temperature ) ) then
          this%spectral_weight_parms( fileNdx )%temperature =                 &
              netcdf_obj%temperature
        endif
        deallocate( netcdf_obj )
      enddo file_loop
    endif has_netcdf_file

    if( associated( lambdaGrid ) ) deallocate( lambdaGrid )

  end subroutine base_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function run( this, grid_warehouse, profile_warehouse )                     &
      result( spectral_weight )
    ! Calculate the spectral wght for a given set of environmental conditions

    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(spectral_weight_t),  intent(in)     :: this ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: spectral_weight(:) ! The calculated spectral weights (wavelength) [unitless]

    ! Local variables
    character(len=*), parameter :: Iam = 'spectral weight calculate: '

    spectral_weight = this%spectral_weight_parms(1)%array( :, 1 )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_points( this, config, data_lambda, data_parameter )
    ! Adds points to the spectral wght grid based on configuration data

    use musica_string,                 only : string_t
    use musica_assert,                 only : assert_msg, die_msg
    use tuvx_util,                     only : add_point

    class(spectral_weight_t), intent(in)    :: this ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    type(config_t),           intent(inout) :: config ! Spectral weight configuration data
    real(dk), allocatable,    intent(inout) :: data_lambda(:) ! Wavelength grid
    real(dk), allocatable,    intent(inout) :: data_parameter(:) ! Parameters (wavelength)

    ! Local variables
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk
    real(dk), parameter :: deltax = 1.e-5_dk
    character(len=*), parameter :: Iam = 'spectral_weight; add_points'

    integer  :: nRows
    real(dk) :: lowerLambda, upperLambda
    real(dk) :: addpnt_val_lower, addpnt_val_upper
    type(string_t)  :: addpnt_type
    type(config_t)  :: extrap_config
    logical         :: found
    type(string_t) :: required_keys(1), optional_keys(1)

    required_keys(1) = "type"
    optional_keys(1) = "value"

    ! add endpoints to data arrays; first the lower bound
    nRows = size( data_lambda )
    lowerLambda = data_lambda(1) ; upperLambda = data_lambda( nRows )

    addpnt_val_lower = rZERO
    call config%get( 'lower extrapolation', extrap_config, Iam, found = found )
    if( found ) then
      call assert_msg( 671608111,                                             &
                       extrap_config%validate( required_keys, optional_keys ),&
                       "Bad format for extrapolation in base spectral wght." )
      call extrap_config%get( "type", addpnt_type, Iam )
      if( addpnt_type == "boundary" ) then
        addpnt_val_lower = data_parameter(1)
      elseif( addpnt_type == "constant" ) then
        call extrap_config%get( "value", addpnt_val_lower, Iam )
      else
        call die_msg( 316405972,                                              &
                      "Bad extrapolation type: '"//addpnt_type//"'" )
      endif
    endif

    ! add endpoints to data arrays; now the upper bound
    addpnt_val_upper = rZERO
    call config%get( 'upper extrapolation', extrap_config, Iam, found = found )
    if( found ) then
      call assert_msg( 918590096,                                             &
                       extrap_config%validate( required_keys, optional_keys ),&
                       "Bad format for extrapolation in base spectral wght." )
      call extrap_config%get( "type", addpnt_type, Iam )
      if( addpnt_type == "boundary" ) then
        addpnt_val_upper = data_parameter( nRows )
      elseif( addpnt_type == "constant" ) then
        call extrap_config%get( "value", addpnt_val_upper, Iam )
      else
        call die_msg( 302970880,                                              &
                      "Bad extrapolation type: '"//addpnt_type//"'" )
      endif
    endif

    call add_point( x = data_lambda, y = data_parameter,                      &
                    xnew = ( rONE - deltax ) * lowerLambda,                   &
                    ynew = addpnt_val_lower )
    call add_point( x = data_lambda, y = data_parameter, xnew = rZERO,        &
                    ynew = addpnt_val_lower )
    call add_point( x = data_lambda, y = data_parameter,                      &
                    xnew = ( rONE + deltax ) * upperLambda,                   &
                    ynew = addpnt_val_upper )
    call add_point( x = data_lambda, y = data_parameter, xnew = 1.e38_dk,     &
                    ynew = addpnt_val_upper )

  end subroutine add_points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size(this, comm)
    ! Calculate the size, in bytes, needed to pack the spectral_weight_parms
    ! spectral_weight_parms is a variable length array of two more variable
    ! length arrays. The two sub-arrays may or may not exist.
    ! The pack size will represent the total size of the subarrays plus
    ! one byte per subarray that indicates if it exists or not

    use musica_mpi, only : musica_mpi_pack_size

    class(spectral_weight_t), intent(in) :: this ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    integer,           intent(in) :: comm        ! MPI communicator

    integer :: ndx
    logical :: is_allocated

#ifdef MUSICA_USE_MPI
    is_allocated = allocated( this%spectral_weight_parms )

    pack_size = musica_mpi_pack_size(is_allocated, comm)

    if (is_allocated) then
      pack_size = pack_size + musica_mpi_pack_size(                           &
        size( this%spectral_weight_parms ), comm )

      do ndx = 1,size( this%spectral_weight_parms )

        pack_size = pack_size + musica_mpi_pack_size(                         &
          this%spectral_weight_parms( ndx )%array, comm )

        pack_size = pack_size + musica_mpi_pack_size(                         &
          this%spectral_weight_parms( ndx )%temperature, comm )
      enddo
    endif
    pack_size = pack_size +                                                   &
                this%wavelength_grid_%pack_size( comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Pack the given value to the buffer, advancing position

    use musica_mpi, only : musica_mpi_pack
    use musica_assert,                 only : assert

    class(spectral_weight_t), intent(in) :: this ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    character, intent(inout) :: buffer(:) ! Memory buffer
    integer,   intent(inout) :: position  ! Current buffer position
    integer,   intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ndx
    logical :: is_allocated
    is_allocated = allocated( this%spectral_weight_parms )

    prev_position = position

    call musica_mpi_pack( buffer, position, is_allocated, comm )

    if (is_allocated) then
      ! first pack the number of arrays we have
      call musica_mpi_pack( buffer, position,                                 &
                            size( this%spectral_weight_parms ), comm )

      do ndx = 1,size( this%spectral_weight_parms )
        call musica_mpi_pack( buffer, position,                               &
          this%spectral_weight_parms( ndx )%array, comm )
        call musica_mpi_pack( buffer, position,                               &
          this%spectral_weight_parms( ndx )%temperature, comm )
      enddo
    endif
    call this%wavelength_grid_%mpi_pack( buffer, position, comm )
    call assert(243454577, &
         position - prev_position <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpack the given value from the buffer, advancing position

    use musica_mpi, only : musica_mpi_unpack
    use musica_assert,                 only : assert

    class(spectral_weight_t), intent(out) :: this ! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ndx, n_params
    logical :: is_allocated

    prev_position = position

    call musica_mpi_unpack( buffer, position, is_allocated, comm )

    if (is_allocated) then
      call musica_mpi_unpack( buffer, position, n_params, comm )

      allocate( this%spectral_weight_parms( n_params) )

      do ndx = 1, n_params
        call musica_mpi_unpack( buffer, position,                             &
          this%spectral_weight_parms( ndx )%array, comm )
        call musica_mpi_unpack( buffer, position,                             &
          this%spectral_weight_parms( ndx )%temperature, comm )
      enddo
    endif
    call this%wavelength_grid_%mpi_unpack( buffer, position, comm )
    call assert(567335412, &
         position - prev_position <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Finalizes the spectral wght type

    type(spectral_weight_t), intent(inout) :: this! This :f:type:`~tuvx_spectral_weight/spectral_weight_t`

    ! Local variables
    integer :: ndx

    if( allocated( this%spectral_weight_parms ) ) then
      do ndx = 1,size( this%spectral_weight_parms )
        if( allocated( this%spectral_weight_parms( ndx )%array ) ) then
          deallocate( this%spectral_weight_parms( ndx )%array )
        endif
        if( allocated( this%spectral_weight_parms( ndx )%temperature ) ) then
          deallocate( this%spectral_weight_parms( ndx )%temperature )
        endif
      enddo
      deallocate( this%spectral_weight_parms )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spectral_weight
