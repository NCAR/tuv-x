! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_o2
  ! O2 profile from csv file type

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t

  implicit none

  public :: profile_o2_t

  type, extends(profile_t) :: profile_o2_t
  contains
    final     :: finalize
  end type profile_o2_t

  interface profile_o2_t
    module procedure constructor
  end interface profile_o2_t

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

    type(profile_o2_t),     pointer       :: this ! This f:type:`~tuvx_profile_o2/profile_o2_t`
    type(config_t),         intent(inout) :: config ! A profile config
    type(grid_warehouse_t), intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    ! local variables
    integer, parameter :: Ok = 0
    integer, parameter :: inUnit = 20
    real(dk), parameter    :: o2Vmr = .2095_dk
    real(dk), parameter    :: km2cm = 1.e5_dk

    character(len=*), parameter :: Iam = 'From_csv_file profile initialize: '

    integer  :: istat, nData
    real(dk) :: zd, Value
    real(dk) :: exo_layer_dens
    logical  :: found
    character(len=132) :: InputLine
    type(string_t)     :: Filespec, Interpolator
    real(dk), allocatable :: zdata(:)
    real(dk), allocatable :: profile(:)
    real(dk), allocatable :: airlog(:)
    class(grid_t), pointer :: zGrid
    class(interpolator_t), pointer :: theInterpolator
    type(string_t) :: required_keys(4), optional_keys(2)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "file path"
    required_keys(4) = "name"
    optional_keys(1) = "interpolator"
    optional_keys(2) = "scale height"

    call assert_msg( 355866380,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "O2 profile." )

    allocate( this )

    ! Get the configuration settings
    call config%get( 'file path', Filespec, Iam )
    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )
    call config%get( 'interpolator', Interpolator, Iam, default = 'linear' )
    call config%get( 'scale height', this%hscale_, Iam, default = 8.01_dk )

    ! Does input grid file exist?
    inquire( file=Filespec%to_char( ), exist=found )
    if( .not. found ) then
      call die_msg( 560768215, "File " // Filespec%to_char( ) // " not found" )
    endif

    open(unit=inUnit,file=Filespec%to_char( ),iostat=istat)
    if( istat /= Ok ) then
      call die_msg( 560768231, "Error opening " // Filespec%to_char( ) )
    endif

    ! Skip the header
    do
      read(inUnit,'(a)',iostat=istat) InputLine
      if( istat /= Ok ) then
        call die_msg( 560768227, "Error reading " // Filespec%to_char( ) )
        exit
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
        call die_msg( 560768229, "Invalid data format in " &
          // Filespec%to_char( ) )
      endif
      profile = [ profile, Value ]
      zdata = [ zdata, zd ]
      read( inUnit, '(a)', iostat = istat ) InputLine
      if( istat /= Ok ) then
        exit
      endif
    enddo

    close( unit = inUnit )

    zGrid => grid_warehouse%get_grid( "height", "km" )
    this%ncells_ = zGrid%ncells_

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

    nData = size( zdata )
    zdata( nData ) = zdata( nData ) + .001_dk
    airlog = log( profile )
    this%edge_val_ = theInterpolator%interpolate( zGrid%edge_, zdata, airlog )
    this%edge_val_ = exp( this%edge_val_ )
    this%edge_val_ = o2Vmr * this%edge_val_

    this%mid_val_ = .5_dk * ( this%edge_val_( 1 : this%ncells_ ) +            &
                              this%edge_val_( 2 : this%ncells_ + 1 ) )

    this%delta_val_ = ( this%edge_val_( 2 : this%ncells_ + 1 ) -              &
                        this%edge_val_( 1 : this%ncells_ ) )

    this%layer_dens_ = zGrid%delta_ *                                         &
      sqrt( this%edge_val_( 1 : this%ncells_ ) ) *                            &
      sqrt( this%edge_val_( 2 : this%ncells_ + 1 ) ) * km2cm

    exo_layer_dens = this%edge_val_( this%ncells_ + 1 ) * this%hscale_ * km2cm

    this%exo_layer_dens_ = [ this%layer_dens_, exo_layer_dens ]

    this%layer_dens_( this%ncells_ ) = this%layer_dens_( this%ncells_ ) +     &
                                       exo_layer_dens

    deallocate( zGrid )
    deallocate( theInterpolator )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleanup the memory used by this object

    type(profile_o2_t), intent(inout) :: this ! This f:type:`~tuvx_profile_o2/profile_o2_t`

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
    if( allocated( this%exo_layer_dens_ ) ) then
      deallocate( this%exo_layer_dens_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_o2
