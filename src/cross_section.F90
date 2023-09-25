! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section
! The base cross section type and related functions

  use musica_constants,                only : dk => musica_dk
  use tuvx_grid_warehouse,             only : grid_warehouse_ptr
  use tuvx_profile_warehouse,          only : profile_warehouse_ptr

  implicit none

  private
  public :: cross_section_t, cross_section_parms_t, base_constructor,         &
            cross_section_ptr

  type cross_section_parms_t
    ! local working type for holding cross section parameters
    real(dk), allocatable :: temperature(:) ! Temperature grid [K]
    real(dk), allocatable :: deltaT(:)      ! Temperature difference between grid sections [K]
    real(dk), allocatable :: array(:,:)     ! Cross section parameters (wavelength, parameter type)
  contains
    ! Returns the number of bytes needed to pack the parameters onto a buffer
    procedure :: pack_size => parms_pack_size
    ! Packs the parameters onto a character buffer
    procedure :: mpi_pack => parms_mpi_pack
    ! Unpacks parameters from a character buffer
    procedure :: mpi_unpack => parms_mpi_unpack
  end type cross_section_parms_t

  type override_t
    ! local working type for specifying constant values for specific ranges
    integer :: min_wavelength_index_
    integer :: max_wavelength_index_
    real(kind=dk) :: value_
  contains
    procedure :: apply
    ! Returns the number of bytes needed to pack the override onto a buffer
    procedure :: pack_size => override_pack_size
    ! Packs the override onto a character buffer
    procedure :: mpi_pack => override_mpi_pack
    ! Unpacks override from a character buffer
    procedure :: mpi_unpack => override_mpi_unpack
  end type override_t

  interface override_t
    module procedure :: override_constructor
  end interface override_t

  type cross_section_t
    ! Calculator for cross_section
    ! The cross section array

    ! Cross section parameter sets
    type(cross_section_parms_t), allocatable :: cross_section_parms(:)
    ! Height grid pointer
    type(grid_warehouse_ptr) :: height_grid_
    ! Wavelength grid pointer
    type(grid_warehouse_ptr) :: wavelength_grid_
    ! Temperature profile pointer
    type(profile_warehouse_ptr) :: temperature_profile_
    ! Override values for specific bands
    type(override_t), allocatable :: overrides_(:)
  contains
    !> Calculate the cross section
    procedure :: calculate
    !> Add points to the cross section grid based on configuration data
    procedure :: add_points
    ! Returns the number of bytes needed to pack the cross section onto a
    ! buffer
    procedure :: pack_size
    ! Packs the cross section onto a character buffer
    procedure :: mpi_pack
    ! Unpacks the cross section from a character buffer into the object
    procedure :: mpi_unpack
    ! Processes a NetCDF input file
    procedure :: process_file
  end type cross_section_t

  interface cross_section_t
    module procedure :: constructor
  end interface

  type cross_section_ptr
    ! Pointer type for building sets of photo rate constants
    class(cross_section_t), pointer :: val_ => null( )
  end type cross_section_ptr

  real(dk), parameter    :: rZERO = 0.0_dk
  real(dk), parameter    :: rONE  = 1.0_dk

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Create an instance of the base cross section type

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(cross_section_t),    pointer       :: this ! Base :f:type:`~tuvx_cross_section/cross_section_t` type
    type(config_t),            intent(inout) :: config ! Cross section configuration object
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    type(string_t) :: required_keys(1), optional_keys(5)

    required_keys(1) = "type"
    optional_keys(1) = "netcdf files"
    optional_keys(2) = "name"
    optional_keys(3) = "merge data"
    optional_keys(4) = "override bands"
    optional_keys(5) = "apply O2 bands"
    call assert_msg( 124969900,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "base cross section." )
    allocate( this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine base_constructor( this, config, grid_warehouse,               &
      profile_warehouse )
    ! Initialize cross_section_t objects

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(cross_section_t),    pointer       :: this ! A :f:type:`~tuvx_cross_section/cross_section_t`
    type(config_t),            intent(inout) :: config ! Cross section configuration object
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    !   local variables
    character(len=*), parameter   :: Iam = 'base cross section initialize'
    integer :: i_param, i_file, i_override
    logical :: found
    type(config_t) :: netcdf_files, netcdf_file, overrides, override
    class(iterator_t), pointer :: iter
    logical :: merge_data
    class(grid_t), pointer :: wavelengths

    ! Get grid and profile pointers
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    this%temperature_profile_ = profile_warehouse%get_ptr( "temperature", "K" )

    ! get cross section netcdf filespec
    call config%get( 'netcdf files', netcdf_files, Iam, found = found )
    if( found ) then
      iter => netcdf_files%get_iterator( )
      call config%get( 'merge data', merge_data, Iam, default = .false. )
      allocate( this%cross_section_parms( netcdf_files%number_of_children( ) ) )
      i_file = 0
      do while( iter%next( ) )
        call netcdf_files%get( iter, netcdf_file, Iam )
        i_file = i_file + 1
        if( merge_data ) i_file = 1
        call this%process_file( netcdf_file, grid_warehouse,                  &
                                this%cross_section_parms( i_file ) )
      enddo
      deallocate( iter )
    end if

    ! get values to overlay for specific bands
    call config%get( "override bands", overrides, Iam, found = found )
    if( found ) then
      wavelengths => grid_warehouse%get_grid( this%wavelength_grid_ )
      iter => overrides%get_iterator( )
      allocate( this%overrides_( overrides%number_of_children( ) ) )
      i_override = 0
      do while( iter%next( ) )
        call overrides%get( iter, override, Iam )
        i_override = i_override + 1
        this%overrides_( i_override ) = override_t( override, wavelengths )
      end do
      deallocate( iter )
      deallocate( wavelengths )
    end if

  end subroutine base_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine process_file( this, config, grid_warehouse, parameters )
    ! Processes a NetCDF input file for the cross section

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_interpolate,              only : interpolator_t
    use tuvx_netcdf,                   only : netcdf_t

    class(cross_section_t),      intent(inout) :: this
    type(config_t),              intent(inout) :: config
    type(grid_warehouse_t),      intent(inout) :: grid_warehouse
    type(cross_section_parms_t), intent(inout) :: parameters

    character(len=*), parameter :: Iam = 'base cross section file reader'
    character(len=*), parameter :: hdr = 'cross_section_'
    class(grid_t), pointer :: wavelength_grid
    type(string_t) :: file_path
    type(netcdf_t), allocatable :: netcdf_obj
    type(config_t) :: interpolator_config
    class(interpolator_t), pointer :: interpolator
    integer :: n_params, i_param
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    real(dk), allocatable :: work_parameter(:)
    logical :: found
    type(string_t) :: required_keys(1), optional_keys(5)
    real(dk) :: bound
    integer :: i_low, i_high

    required_keys(1) = "file path"
    optional_keys(1) = "interpolator"
    optional_keys(2) = "lower extrapolation"
    optional_keys(3) = "upper extrapolation"
    optional_keys(4) = "zero below"
    optional_keys(5) = "zero above"
    call assert_msg( 538565562,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "base cross section NetCDF file." )

    wavelength_grid => grid_warehouse%get_grid( this%wavelength_grid_ )

    ! read netcdf cross section parameters
    allocate( netcdf_obj )
    call config%get( "file path", file_path, Iam )
    call netcdf_obj%read_netcdf_file( file_path = file_path%to_char( ),       &
                                      variable_name = hdr )
    n_params = size( netcdf_obj%parameters, dim = 2 )
    call assert_msg( 681779091, n_params >= 1,                                &
                     'File: '//file_path//' contains no parameters' )

    call config%get( "interpolator", interpolator_config, Iam, found = found )
    if( .not. found ) then
      call interpolator_config%empty( )
      call interpolator_config%add( "type", "conserving", Iam )
    end if
    interpolator => interpolator_t( interpolator_config )

    ! trim data if specified
    i_low = 1
    i_high = wavelength_grid%ncells_
    call config%get( "zero below", bound, Iam, found = found )
    if( found ) then
      do while( wavelength_grid%edge_( i_low ) < bound )
        i_low = i_low + 1
        call assert_msg( 395014916, i_low < i_high,                           &
                         "Invalid lower bound specified for file "//file_path )
      end do
    end if
    call config%get( "zero above", bound, Iam, found = found )
    if( found ) then
      do while( wavelength_grid%edge_( i_high ) > bound )
        i_high = i_high - 1
        call assert_msg( 434154076, i_high > i_low,                           &
                         "Invalid upper bound specified for file "//file_path )
      end do
    end if

    ! interpolate from data to model wavelength grid
    if( allocated( netcdf_obj%wavelength ) ) then
      if( .not. allocated( parameters%array ) )  then
        allocate( parameters%array( wavelength_grid%ncells_, n_params ) )
        parameters%array(:,:) = 0.0_dk
      endif
      call assert_msg( 622609278, size( parameters%array, 2 ) == n_params,    &
                       "Cannot merge data with different shapes" )
      allocate( work_parameter( wavelength_grid%ncells_ ) )
      do i_param = 1, n_params
        data_lambda    = netcdf_obj%wavelength
        data_parameter = netcdf_obj%parameters( :, i_param )
        call this%add_points( config, data_lambda, data_parameter )
        work_parameter(:) =                                                   &
            interpolator%interpolate( x_target = wavelength_grid%edge_,       &
                                      x_source = data_lambda,                 &
                                      y_source = data_parameter )
        parameters%array( i_low:i_high, i_param ) =                           &
                                  parameters%array( i_low:i_high, i_param ) + &
                                  work_parameter( i_low:i_high )
      enddo
    else
      parameters%array = netcdf_obj%parameters
    endif
    if( allocated( netcdf_obj%temperature ) ) then
      parameters%temperature = netcdf_obj%temperature
    endif

    deallocate( interpolator    )
    deallocate( netcdf_obj      )
    deallocate( wavelength_grid )

  end subroutine process_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculate the cross section for a given set of environmental conditions

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    real(kind=dk), allocatable               :: cross_section(:,:) ! Calculated cross section
    class(cross_section_t),    intent(in)    :: this ! A :f:type:`~tuvx_cross_section/cross_section_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,         intent(in)    :: at_mid_point ! Flag indicating that the cross section data is at grid mid-points. If omitted or false, data is assumed to be at interfaces

    !> Local variables
    integer :: colndx
    integer :: nzdim, i_override
    character(len=*), parameter :: Iam =                                      &
        'radXfer base cross section calculate: '
    class(grid_t), pointer     :: zGrid
    real(dk),      allocatable :: wrkCrossSection(:,:)

    zGrid => grid_warehouse%get_grid( this%height_grid_ )

    nzdim = zGrid%ncells_ + 1
    if( present( at_mid_point ) ) then
      if( at_mid_point ) then
        nzdim = nzdim - 1
      endif
    endif

    allocate( wrkCrossSection( size( this%cross_section_parms(1)%array,       &
                                                          dim = 1 ), nzdim ) )

    !> Just copy the lambda interpolated array
    do colndx = 1, nzdim
      wrkCrossSection( :, colndx ) = this%cross_section_parms(1)%array(:,1)
      if( allocated( this%overrides_ ) ) then
        do i_override = 1, size( this%overrides_ )
          call this%overrides_( i_override )%apply(                           &
                                                wrkCrossSection( :, colndx ) )
        end do
      end if
    enddo

    cross_section = transpose( wrkCrossSection )

    deallocate( zGrid )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_points( this, config, data_lambda, data_parameter )
    ! Adds points to the cross section grid based on configuration data

    use musica_assert,                 only : assert_msg, die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_util,                     only : add_point

    class(cross_section_t), intent(in)    :: this ! This :f:type:`~tuvx_cross_section/cross_section_t`
    type(config_t),         intent(inout) :: config ! The configuration used to build this object
    real(dk), allocatable,  intent(inout) :: data_lambda(:) ! Wavelength grid
    real(dk), allocatable,  intent(inout) :: data_parameter(:) ! Parameters (wavelength)

    real(dk), parameter :: deltax = 1.e-5_dk
    character(len=*), parameter :: Iam = 'cross_section; addpnts: '

    integer  :: nRows
    real(dk) :: lowerLambda, upperLambda
    real(dk) :: addpnt_val_lower, addpnt_val_upper
    type(string_t)  :: addpnt_type
    type(config_t)  :: extrap_config
    logical         :: found
    type(string_t) :: required_keys(1), optional_keys(1)

    required_keys(1) = "type"
    optional_keys(1) = "value"

    nRows = size( data_lambda )
    lowerLambda = data_lambda(1) ; upperLambda = data_lambda( nRows )

    ! add endpoints to data arrays; first the lower bound
    addpnt_val_lower = rZERO
    call config%get( 'lower extrapolation', extrap_config, Iam, found = found )
    if( found ) then
      call assert_msg( 671608110,                                             &
                       extrap_config%validate( required_keys, optional_keys ),&
                       "Bad format for extrapolation in base cross section." )
      call extrap_config%get( "type", addpnt_type, Iam )
      if( addpnt_type == "boundary" ) then
        addpnt_val_lower = data_parameter(1)
      elseif( addpnt_type == "constant" ) then
        call extrap_config%get( "value", addpnt_val_lower, Iam )
      else
        call die_msg( 316405971,                                              &
                      "Bad extrapolation type: '"//addpnt_type//"'" )
      endif
    endif

    !> add endpoints to data arrays; now the upper bound
    addpnt_val_upper = rZERO
    call config%get( 'upper extrapolation', extrap_config, Iam, found = found )
    if( found ) then
      call assert_msg( 918590095,                                             &
                       extrap_config%validate( required_keys, optional_keys ),&
                       "Bad format for extrapolation in base cross section." )
      call extrap_config%get( "type", addpnt_type, Iam )
      if( addpnt_type == "boundary" ) then
        addpnt_val_upper = data_parameter( nRows )
      elseif( addpnt_type == "constant" ) then
        call extrap_config%get( "value", addpnt_val_upper, Iam )
      else
        call die_msg( 302970879,                                              &
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

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the cross section onto a
    ! buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(cross_section_t), intent(in) :: this ! cross section to be packed
    integer,                intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_elem

    pack_size =                                                               &
        musica_mpi_pack_size( allocated( this%cross_section_parms ), comm )
    if( allocated( this%cross_section_parms ) ) then
      pack_size = pack_size +                                                 &
                musica_mpi_pack_size( size( this%cross_section_parms ), comm )
      do i_elem = 1, size( this%cross_section_parms )
        pack_size = pack_size +                                               &
                    this%cross_section_parms( i_elem )%pack_size( comm )
      end do
    end if
    pack_size = pack_size +                                                   &
                this%height_grid_%pack_size( comm ) +                         &
                this%wavelength_grid_%pack_size( comm ) +                     &
                this%temperature_profile_%pack_size( comm ) +                 &
           musica_mpi_pack_size( allocated( this%overrides_ ), comm )
    if( allocated( this%overrides_ ) ) then
      pack_size = pack_size +                                                 &
                  musica_mpi_pack_size( size( this%overrides_ ), comm )
      do i_elem = 1, size( this%overrides_ )
        pack_size = pack_size +                                               &
                    this%overrides_( i_elem )%pack_size( comm )
      end do
    end if
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the cross section onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(cross_section_t), intent(in)    :: this      ! cross section to be packed
    character,              intent(inout) :: buffer(:) ! memory buffer
    integer,                intent(inout) :: position  ! current buffer position
    integer,                intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem

    prev_pos = position
    call musica_mpi_pack( buffer, position,                                   &
                          allocated( this%cross_section_parms ), comm )
    if( allocated( this%cross_section_parms ) ) then
      call musica_mpi_pack( buffer, position,                                 &
                            size( this%cross_section_parms ), comm )
      do i_elem = 1, size( this%cross_section_parms )
        call this%cross_section_parms( i_elem )%mpi_pack( buffer, position,   &
                                                           comm )
      end do
    end if
    call this%height_grid_%mpi_pack(         buffer, position, comm )
    call this%wavelength_grid_%mpi_pack(     buffer, position, comm )
    call this%temperature_profile_%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, allocated( this%overrides_ ),     &
                          comm )
    if( allocated( this%overrides_ ) ) then
      call musica_mpi_pack( buffer, position, size( this%overrides_ ), comm )
      do i_elem = 1, size( this%overrides_ )
        call this%overrides_( i_elem )%mpi_pack( buffer, position, comm )
      end do
    end if
    call assert( 345613473, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks the cross section from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(cross_section_t), intent(out)   :: this      ! cross section to be unpacked
    character,              intent(inout) :: buffer(:) ! memory buffer
    integer,                intent(inout) :: position  ! current buffer position
    integer,                intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_elem, n_elem
    logical :: alloced

    prev_pos = position
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( allocated( this%cross_section_parms ) )                               &
        deallocate( this%cross_section_parms )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elem, comm )
      allocate( this%cross_section_parms( n_elem ) )
      do i_elem = 1, n_elem
        call this%cross_section_parms( i_elem )%mpi_unpack( buffer, position,&
                                                            comm )
      end do
    end if
    call this%height_grid_%mpi_unpack(         buffer, position, comm )
    call this%wavelength_grid_%mpi_unpack(     buffer, position, comm )
    call this%temperature_profile_%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( allocated( this%overrides_ ) ) deallocate( this%overrides_ )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_elem, comm )
      allocate( this%overrides_( n_elem ) )
      do i_elem = 1, n_elem
        call this%overrides_( i_elem )%mpi_unpack( buffer, position, comm )
      end do
    end if
    call assert( 764657896, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function parms_pack_size( this, comm ) result( pack_size )
    ! Returns the size of a character buffer required to pack the parameters

    use musica_mpi,                    only : musica_mpi_pack_size

    class(cross_section_parms_t), intent(in) :: this ! parameters to be packed
    integer,                      intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%temperature, comm ) +              &
                musica_mpi_pack_size( this%deltaT,      comm ) +              &
                musica_mpi_pack_size( this%array,       comm )
#else
    pack_size = 0
#endif

  end function parms_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parms_mpi_pack( this, buffer, position, comm )
    ! Packs the parameters onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(cross_section_parms_t), intent(in)    :: this      ! parameters to be packed
    character,                    intent(inout) :: buffer(:) ! memory buffer
    integer,                      intent(inout) :: position  ! current buffer position
    integer,                      intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%temperature, comm )
    call musica_mpi_pack( buffer, position, this%deltaT,      comm )
    call musica_mpi_pack( buffer, position, this%array,       comm )
    call assert( 841935272, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine parms_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine parms_mpi_unpack( this, buffer, position, comm )
    ! Unpacks parameters from a character buffer into the object

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(cross_section_parms_t), intent(out)   :: this      ! parameters to be unpacked
    character,                    intent(inout) :: buffer(:) ! memory buffer
    integer,                      intent(inout) :: position  ! current buffer position
    integer,                      intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%temperature, comm )
    call musica_mpi_unpack( buffer, position, this%deltaT,      comm )
    call musica_mpi_unpack( buffer, position, this%array,       comm )
    call assert( 363886174, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine parms_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function override_constructor( config, wavelengths ) result( this )
    ! Constructor for override_t objects

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_la_sr_bands,              only : get_band_min_index,             &
                                              get_band_max_index

    type(override_t)               :: this
    type(config_t),  intent(inout) :: config
    class(grid_t),   intent(in)    :: wavelengths

    character(len=*), parameter :: my_name =                                  &
        "cross section band override constructor"
    type(string_t) :: type_name
    type(string_t) :: required_keys(2), optional_keys(0)

    required_keys(1) = "band"
    required_keys(2) = "value"
    call assert_msg( 145671194,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration for cross section band averride" )
    call config%get( "band", type_name, my_name )
    this%min_wavelength_index_ = get_band_min_index( type_name, wavelengths )
    this%max_wavelength_index_ = get_band_max_index( type_name, wavelengths )
    call config%get( "value", this%value_, my_name )

  end function override_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine apply( this, cross_section )
    ! Overwrites cross section values for the specified band

    class(override_t), intent(in)    :: this
    real(kind=dk),           intent(inout) :: cross_section(:)

    cross_section( this%min_wavelength_index_ : this%max_wavelength_index_ )  &
        = this%value_

  end subroutine apply

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function override_pack_size( this, comm ) result( pack_size )
    ! Returns the size of a character buffer required to pack the override

    use musica_mpi,                    only : musica_mpi_pack_size

    class(override_t), intent(in) :: this ! override to be packed
    integer,           intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%min_wavelength_index_, comm ) +    &
                musica_mpi_pack_size( this%max_wavelength_index_, comm ) +    &
                musica_mpi_pack_size( this%value_,                comm )
#else
    pack_size = 0
#endif

  end function override_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine override_mpi_pack( this, buffer, position, comm )
    ! Packs the override onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(override_t), intent(in)    :: this      ! override to be packed
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%min_wavelength_index_, comm )
    call musica_mpi_pack( buffer, position, this%max_wavelength_index_, comm )
    call musica_mpi_pack( buffer, position, this%value_,                comm )
    call assert( 408243741, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine override_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine override_mpi_unpack( this, buffer, position, comm )
    ! Unpacks override from a character buffer into the object

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(override_t), intent(out)   :: this      ! override to be unpacked
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%min_wavelength_index_,comm )
    call musica_mpi_unpack( buffer, position, this%max_wavelength_index_,comm )
    call musica_mpi_unpack( buffer, position, this%value_,               comm )
    call assert( 687360217, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine override_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section
