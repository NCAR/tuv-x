! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_o3
  ! O3 profile type

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t

  implicit none

  public :: profile_o3_t

  type, extends(profile_t) :: profile_o3_t
  contains
    final     :: finalize
  end type profile_o3_t

  !> Constructor
  interface profile_o3_t
    module procedure constructor
  end interface profile_o3_t

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

    type(profile_o3_t),     pointer       :: this ! This f:type:`~tuvx_profile_o3/profile_o3_t`
    type(config_t),         intent(inout) :: config ! A profile config
    type(grid_warehouse_t), intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    ! local variables
    integer, parameter :: Ok = 0
    integer, parameter :: inUnit = 20
    real(dk), parameter    :: km2cm = 1.e5_dk
    character(len=*), parameter :: Iam = &
      'O3_from_csv_file profile initialize: '

    integer :: istat, m
    real(dk)    :: zd, Value, zTop
    real(dk)    :: rfact
    real(dk)    :: Scale2DU, InputDU, O3ScaleFactor
    logical :: found
    character(len=132) :: InputLine
    type(string_t)     :: Filespec, Interpolator
    real(dk), allocatable :: zdata(:)
    real(dk), allocatable :: profile(:)
    class(grid_t), pointer :: zGrid
    class(interpolator_t), pointer :: theInterpolator
    type(string_t) :: required_keys(4), optional_keys(3)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "file path"
    required_keys(4) = "name"
    optional_keys(1) = "interpolator"
    optional_keys(2) = "scale height"
    optional_keys(3) = "reference column"

    call assert_msg( 495794327,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "O3 profile." )

    allocate( this )

    ! Get the configuration settings
    call config%get( 'file path', Filespec, Iam )
    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )
    call config%get( 'interpolator', Interpolator, Iam, default = 'linear' )
    call config%get( 'scale height', this%hscale_, Iam, default = 4.5_dk )
    call config%get( 'reference column', Scale2DU, Iam, default = 300._dk )

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
      elseif( verify( InputLine(1:1), '#!$%*' ) /= 0 ) then
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

    ! Set o3 concentration if data ztop < mdl top
    ztop = zGrid%edge_( this%ncells_ + 1 )
    if( this%hscale_ /= 0.0_dk ) then
      rfact = exp( -1.0_dk/this%hscale_ )
    else
      rfact = 0.0_dk
    endif
    m = size( profile )
    do
      zdata   = [ zdata, zdata( m ) + 1.0_dk ]
      profile = [ profile, profile( m ) * rfact ]
      m = m + 1
      if( zdata( m ) > zTop ) then
        exit
      endif
    enddo

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

    this%edge_val_ = theInterpolator%interpolate( zGrid%edge_, zdata, profile )

    this%mid_val_ = .5_dk * ( this%edge_val_( 1 : this%ncells_ ) +            &
                              this%edge_val_( 2 : this%ncells_ + 1 ) )

    this%delta_val_ = this%edge_val_( 2 : this%ncells_ + 1 ) -                &
                      this%edge_val_( 1 : this%ncells_ )

    this%layer_dens_ = zGrid%delta_ * this%mid_val_ * km2cm

    this%layer_dens_( this%ncells_ ) = this%layer_dens_( this%ncells_ ) +     &
        this%edge_val_( this%ncells_ + 1 ) * this%hscale_ * km2cm

    if( Scale2DU /= 1.0_dk ) then
      InputDU = sum( this%layer_dens_ ) / 2.687e16_dk
      O3ScaleFactor = Scale2DU / InputDU

      if( O3ScaleFactor /= 1.0_dk ) then
        this%edge_val_ = O3ScaleFactor * this%edge_val_

        this%mid_val_ = .5_dk * ( this%edge_val_( 1 : this%ncells_ ) +        &
                                  this%edge_val_( 2 : this%ncells_ + 1 ) )

        this%delta_val_ = this%edge_val_( 2 : this%ncells_ + 1 ) -            &
                          this%edge_val_( 1 : this%ncells_ )

        this%layer_dens_ = zGrid%delta_ * this%mid_val_ * km2cm

        this%layer_dens_( this%ncells_ ) = this%layer_dens_( this%ncells_ ) + &
            this%edge_val_( this%ncells_ + 1 ) * this%hscale_ * km2cm

      endif
    endif

    deallocate( zGrid )
    deallocate( theInterpolator )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleanup the memory used by this object

    type(profile_o3_t), intent(inout) :: this ! This f:type:`~tuvx_profile_o3/profile_o3_t`

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

end module tuvx_profile_o3
