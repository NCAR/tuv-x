! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_from_csv_file
  ! Profile from csv file type. See
  ! :ref:`configuration-profiles-from-csv` for more information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t

  implicit none

  public :: profile_from_csv_file_t

  type, extends(profile_t) :: profile_from_csv_file_t
  contains
    final     :: finalize
  end type profile_from_csv_file_t

  interface profile_from_csv_file_t
    module procedure constructor
  end interface profile_from_csv_file_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result( this )
    ! Initialize profile

    use musica_assert,                 only : assert_msg, die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_interpolate

    type(profile_from_csv_file_t), pointer       :: this ! This f:type:`~tuvx_profile_from_csv_file/profile_from_csv_file_t`
    type(config_t),                intent(inout) :: config ! A profile config
    type(grid_warehouse_t),        intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    ! local variables
    character(len=*), parameter :: Iam = 'From_csv_file profile initialize: '

    integer, parameter :: Ok = 0
    integer, parameter :: inUnit = 20
    real(dk), parameter    :: km2cm = 1.e5_dk
    class(grid_t), pointer :: grid
    type(config_t) :: grid_config
    type(string_t) :: grid_name, grid_units

    integer  :: istat, k
    real(dk) :: zd, Value, unit_conv
    real(dk) :: exo_layer_dens, accum
    logical  :: found
    character(len=132) :: InputLine
    type(string_t)     :: Filespec, Interpolator
    real(dk), allocatable :: zdata(:)
    real(dk), allocatable :: profile(:)
    class(interpolator_t), pointer :: theInterpolator
    type(string_t) :: required_keys(5), optional_keys(2)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "file path"
    required_keys(4) = "grid"
    required_keys(5) = "name"
    optional_keys(1) = "interpolator"
    optional_keys(2) = "scale height"

    call assert_msg( 789059716,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "profile from csv file." )

    allocate( this )

    ! Get the configuration settings
    call config%get( 'file path', Filespec, Iam )
    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )
    call config%get( 'interpolator', Interpolator, Iam, default = 'linear' )
    call config%get( 'scale height', this%hscale_, Iam, default = 0._dk )
    call config%get( 'grid', grid_config, Iam )
    call grid_config%get( 'name', grid_name, Iam )
    call grid_config%get( 'units', grid_units, Iam )

    ! locate the grid
    grid => grid_warehouse%get_grid( grid_name, grid_units )
    this%ncells_ = grid%ncells_

    ! Does input grid file exist?
    inquire( file = Filespec%to_char( ), exist = found )
    if( .not. found ) then
      call die_msg( 560768215, "File " // Filespec%to_char( ) // " not found" )
    endif

    open( unit = inUnit, file = Filespec%to_char( ), iostat = istat )
    if( istat /= Ok ) then
      call die_msg( 560768231, "Error opening " // Filespec%to_char( ) )
    endif

    ! Skip the header
    do
      read( inUnit, '(a)', iostat = istat ) InputLine
      if( istat /= Ok ) then
        call die_msg( 560768227, "Error reading " // Filespec%to_char( ) )
      elseif( verify( InputLine(1:1),'#!$%*' ) /= 0 ) then
        exit
      endif
    enddo

    allocate( profile(0) )
    allocate( zdata(0) )

    ! Read the data
    do
      read( InputLine, *, iostat = istat ) zd, Value
      if( istat /= Ok ) then
        call die_msg( 560768229, "Invalid data format in " // &
          Filespec%to_char( ) )
      endif
      profile = [ profile, Value ]
      zdata = [ zdata, zd ]
      read( inUnit, '(a)', iostat = istat) InputLine
      if( istat /= Ok ) then
        exit
      endif
    enddo

    close( unit = inUnit )

    ! assign actual interpolator for this profile
    select case( Interpolator%to_char( ) )
      case( 'linear' )
        allocate( interpolator_linear_t :: theInterpolator )
      case( 'conserving' )
        allocate( interpolator_conserving_t :: theInterpolator )
      case( 'fractional source' )
        allocate( interpolator_fractional_source_t :: theInterpolator )
      case( 'fractional target' )
        allocate( interpolator_fractional_target_t :: theInterpolator )
      case default
        call die_msg( 560768275, "interpolator " // Interpolator%to_char( )   &
          // " not a valid selection" )
    end select

    this%edge_val_ = theInterpolator%interpolate( grid%edge_, zdata, profile )

    this%mid_val_ = .5_dk * ( this%edge_val_( 1 : this%ncells_ ) +            &
                              this%edge_val_( 2 : this%ncells_ + 1 ) )

    this%delta_val_  = ( this%edge_val_( 2 : this%ncells_ + 1 ) -             &
                         this%edge_val_( 1 : this%ncells_) )

    ! This can be removed as part of issue #139 for adopting SI units
    if( grid_name == "height" .and. grid_units == "km" ) then
      unit_conv = km2cm
    else
      unit_conv = 1.0_dk
    end if

    this%layer_dens_ = this%mid_val_ * grid%delta_ * unit_conv

    this%layer_dens_( this%ncells_ ) = this%layer_dens_( this%ncells_ )

    exo_layer_dens = this%edge_val_( this%ncells_ + 1 ) * this%hscale_ *      &
                     unit_conv

    this%exo_layer_dens_ = [ this%layer_dens_, exo_layer_dens ]

    this%layer_dens_( this%ncells_ ) = this%layer_dens_( this%ncells_ ) +     &
      exo_layer_dens

    allocate( this%burden_dens_( grid%ncells_ ) )
    accum = this%layer_dens_( grid%ncells_ )
    this%burden_dens_( grid%ncells_ ) = this%layer_dens_( this%ncells_ )
    do k = grid%ncells_ - 1, 1, -1
      accum = accum + this%layer_dens_( k )
      this%burden_dens_( k ) = accum
    enddo

    deallocate( grid )
    deallocate( theInterpolator )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleanup the memory used by this object

    type(profile_from_csv_file_t), intent(inout) :: this ! This f:type:`~tuvx_profile_from_csv_file/profile_from_csv_file_t`


    if( allocated( this%edge_val_ ) ) then
      deallocate( this%edge_val_ )
    endif
    if( allocated( this%mid_val_ ) ) then
      deallocate( this%mid_val_ )
    endif
    if( allocated( this%delta_val_ ) ) then
      deallocate( this%delta_val_ )
    endif
    if( allocated( this%layer_dens_ ) ) then
      deallocate( this%layer_dens_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_from_csv_file
