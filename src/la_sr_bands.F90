! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_la_sr_bands
  ! Calculator of properties of the Lyman-Alpha and Shuman-Runge bands

  use musica_constants,                only : dk => musica_dk
  use tuvx_grid_warehouse,             only : grid_warehouse_ptr
  use tuvx_profile_warehouse,          only : profile_warehouse_ptr

  implicit none

  private
  public :: la_sr_bands_t, get_band_min_index, get_band_max_index

  integer,  parameter :: nPoly = 20      ! order of the Chebyshev polynomials
  real(dk), parameter :: kLowerLimit = 38.0_dk ! Lower bound of Chebyshev polynomial
  real(dk), parameter :: kUpperLimit = 56.0_dk ! Upper bound of Chebyshev polynomial
  integer,  parameter :: kla = 2         ! dimension of the Lymann-Alpha wavelength grid
  integer,  parameter :: nla = kla - 1   ! number of Lymann-Alpha wavelength bins
  real(dk), parameter :: wlla(kla) = (/ 121.4_dk, 121.9_dk /) ! Lymann-Alpha wavelength grid [nm]
  integer,  parameter :: ksrb = 18       ! dimension of the Schumann-Runge wavelength grid
  integer,  parameter :: nsrb = ksrb - 1 ! number of Schumann-Runge wavelength bins
  real(dk), parameter :: wlsrb(ksrb) =                                        &
  (/ 175.4_dk, 177.0_dk, 178.6_dk, 180.2_dk, 181.8_dk, 183.5_dk, 185.2_dk,    &
     186.9_dk, 188.7_dk, 190.5_dk, 192.3_dk, 194.2_dk, 196.1_dk, 198.0_dk,    &
     200.0_dk, 202.0_dk, 204.1_dk, 206.2_dk/) ! Schumann-Runge wavelength grid [nm]

  integer,  parameter :: iONE = 1
  real(dk), parameter :: rZERO = 0.0_dk
  real(dk), parameter :: rONE  = 1.0_dk
  real(dk), parameter :: rTWO  = 2.0_dk
  real(dk), parameter :: rTEN  = 10.0_dk
  real(dk), parameter :: precis  = 1.e-7_dk
  real(dk), parameter :: largest = 1.e36_dk

  type :: la_sr_bands_t
    ! Calculator of properties of the Lyman-Alpha and Shuman-Runge bands
    integer  :: ila           ! TUV-x photolysis wavelength index where the Lymann-Alpha band starts
    integer  :: isrb          ! TUV-x photolysis wavelength index where the Schumann-Runge band starts
    logical  :: has_la        ! .true. if TUV-x photolysis spectrum includes the Lymann-Alpha band
    logical  :: has_srb       ! .true. if TUV-x photolysis spectrum includes the Schumann-Runge band
    logical  :: has_la_srb    ! .true. if has_la OR has_srb are .true.
    logical  :: do_scaled_O2_ ! .true. if the O2 profile should be scaled from total air
    real(dk) :: O2_scale_factor_ = 0.0_dk ! fraction of air that is O2 by volume
    real(dk) :: AC( nPoly, nsrb ) ! Chebyshev polynomial coefficients
    real(dk) :: BC( nPoly, nsrb ) ! Chebyshev polynomial coefficients
    type(grid_warehouse_ptr) :: height_grid_
    type(grid_warehouse_ptr) :: wavelength_grid_
    type(profile_warehouse_ptr) :: temperature_profile_
    type(profile_warehouse_ptr) :: O2_profile_
  contains
    procedure :: optical_depth => la_srb_OD
    procedure :: cross_section => la_srb_xs
    ! Returns the number of bytes required to pack the calculator onto a buffer
    procedure :: pack_size
    ! Pack the calculator onto a character buffer
    procedure :: mpi_pack
    ! Unpack a calculator from a character buffer
    procedure :: mpi_unpack
    procedure, private :: lymana_OD
    procedure, private :: lymana_xs
    procedure, private :: schum_OD
    procedure, private :: schum_xs
    procedure, private :: init_srb_xs
    procedure, private :: effxs
    procedure, private :: calc_params
    procedure, private :: chebyshev_evaluation
    final :: finalize
  end type la_sr_bands_t

  interface la_sr_bands_t
    ! la_sr_bands_t object constructor
    module procedure constructor
  end interface la_sr_bands_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Checks that the user wavelength grid, WL(IW), is compatible
    ! with the wavelengths for the parameterizations of the Lyman-alpha and
    ! SRB.
    ! Also computes and saves corresponding grid indices (ILA, ISRB)

    use musica_assert,                 only : assert_msg, die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t, to_char
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(la_sr_bands_t),       pointer       :: this
    type(config_t),            intent(inout) :: config            ! Lyman-Alpha Shumann Runge configuration
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse    ! Grid warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! profile warehouse

    character(len=*), parameter :: Iam = 'la_srb initialize: '

    integer :: iw, nw
    type(string_t)         :: file_path
    class(grid_t), pointer :: lambdaGrid
    type(config_t)         :: o2_config

    allocate( this )

    this%height_grid_         = grid_warehouse%get_ptr( "height", "km" )
    this%wavelength_grid_     = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%temperature_profile_ = profile_warehouse%get_ptr( "temperature", "K" )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    ! Are la and srb grids fully "inside" the model grid?
    this%has_la  = lambdaGrid%edge_( iONE ) <= wlla( iONE ) .and.             &
        lambdaGrid%edge_( lambdaGrid%ncells_ + iONE ) >= wlla( kla )
    this%has_srb = lambdaGrid%edge_( iONE ) <= wlsrb( iONE ) .and.            &
        lambdaGrid%edge_( lambdaGrid%ncells_ + iONE ) >= wlsrb( ksrb )
    this%has_la_srb = this%has_la .or. this%has_srb

    has_la_srb: if( this%has_la_srb ) then
      nw = lambdaGrid%ncells_ + iONE
      if( this%has_la ) then
        ! locate Lyman-alpha wavelengths on grid
        this%ila = 0
        do iw = iONE, nw
          if( abs( lambdaGrid%edge_( iw ) - wlla( iONE ) ) < rTEN * precis )  &
              then
            this%ila = iw
            exit
          endif
        enddo
        ! check Lyman-alpha wavelength grid
        call assert_msg( 592167903, this%ila .ne. 0,                          &
                         'Lyman alpha grid mis-match. '//                     &
                         'For wavelengths below 205.8 nm, only the '//        &
                         'pre-specified wavelength grid is permitted.' )
        do iw = 2, nla + iONE
          call assert_msg( 236653044,                                         &
                           abs( lambdaGrid%edge_( this%ila + iw - iONE )      &
                                - wlla( iw ) ) <= rTEN * precis,              &
                           "Lymann-Alpha grid mismatch. Expected "//          &
                           trim( to_char( wlla( iw ) ) )//" but got "//       &
                           trim( to_char( lambdaGrid%edge_( this%ila + iw     &
                                                            - iONE ) ) ) )
        enddo
      endif
      if( this%has_srb ) then
        ! locate Schumann-Runge wavelengths on grid
        this%isrb = 0
        do iw = iONE, nw
          if( abs( lambdaGrid%edge_( iw ) - wlsrb( iONE ) ) < rTEN * precis ) &
              then
            this%isrb = iw
            exit
          endif
        enddo
        ! check Schumann-Runge wavelength grid
        call assert_msg( 469773239,  this%isrb .ne. 0,                        &
                         'Schumann-Runge grid mis-match. '//                  &
                         'For wavelengths below 205.8 nm, only the '//        &
                         'pre-specified wavelength grid is permitted.' )
        do iw = 2, nsrb + iONE
          call assert_msg( 947796740,                                         &
                           abs( lambdaGrid%edge_( this%isrb + iw - iONE )     &
                                - wlsrb( iw ) ) <= rTEN * precis,             &
                           "Shumann-Runge grid mismatch. Expected "//         &
                           trim( to_char( wlsrb( iw ) ) )//" but got "//      &
                           trim( to_char( lambdaGrid%edge_( this%isrb + iw    &
                                                            - iONE ) ) ) )
        enddo
        ! Loads Chebyshev polynomial Coeff.
        call config%get( 'cross section parameters file', file_path, Iam )
        call this%init_srb_xs( file_path )
      endif

      ! Determine how to handle O2 profile
      call config%get( 'O2 estimate', o2_config, Iam,                         &
                       found = this%do_scaled_O2_ )
      if( this%do_scaled_O2_ ) then
        call o2_config%get( 'scale factor', this%O2_scale_factor_, Iam )
      else
        this%O2_profile_ = profile_warehouse%get_ptr( "O2", "molecule cm-3" )
      end if

    endif has_la_srb

    deallocate( lambdaGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine la_srb_OD( this, grid_warehouse, profile_warehouse,              &
      air_vertical_column, air_slant_column, o2_optical_depth,                &
      spherical_geometry )
    ! Computes equivalent optical depths for O2 absorption, and O2 effective
    ! absorption cross sections, parameterized in the Lyman-alpha and SR bands

    use musica_assert,                 only : assert_msg
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t

    class(la_sr_bands_t),      intent(inout) :: this
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(spherical_geometry_t),intent(inout) :: spherical_geometry

    real(dk), intent(in)    :: air_vertical_column(:), air_slant_column(:)
    real(dk), intent(inout) :: o2_optical_depth(:,:)

    ! Local variables
    character(len=*), parameter :: Iam = 'la_srb OD: '

    integer :: nz ! number of specified altitude levels in the working grid
    integer :: nzm1
    integer :: iw
    real(dk)    :: secchi(size(air_slant_column))
    real(dk)    :: o2vcol(size(air_vertical_column))
    real(dk)    :: o2scol(size(air_slant_column))
    class(grid_t),    pointer :: zGrid ! specified altitude working grid [km]
    class(grid_t),    pointer :: lambdaGrid
    class(profile_t), pointer :: temperature, O2_profile

    ! Lyman-alpha variables
    ! O2 optical depth and equivalent cross section in the Lyman-alpha region
    ! Wavelengths for Lyman alpha and SRB parameterizations:
    real(dk) :: o2_optical_depth_la( size( air_vertical_column ), nla )

    ! grid on which Koppers' parameterization is defined
    ! O2 optical depth and equivalent cross section on Koppers' grid
     real(dk) :: o2_optical_depth_k( size( air_vertical_column ), nsrb )

    has_la_srb: if( this%has_la_srb ) then

      ! get specific grids and vertical profiles
      zGrid => grid_warehouse%get_grid( this%height_grid_ )
      lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
      temperature => profile_warehouse%get_profile( this%temperature_profile_ )

      nzm1 = zGrid%ncells_
      nz   = nzm1 +  iONE

      call assert_msg( 543202773, size( air_slant_column, dim=1 ) == nz,      &
                     'invalid dimension size of the air slant column.' //     &
                     'The slant air column must be one larger than the ' //   &
                     'number of cells in the vertical grid')

      ! O2 slant column
      if( this%do_scaled_O2_ ) then
        o2scol(:) = this%O2_scale_factor_ * air_slant_column(:)
      else
        O2_profile => profile_warehouse%get_profile( this%O2_profile_ )
        call spherical_geometry%air_mass( O2_profile%exo_layer_dens_, o2vcol, &
                                          o2scol )
        deallocate( O2_profile )
      end if

      ! Effective secant of solar zenith angle.
      ! Use 2.0 if no direct sun (value for isotropic radiation)
      ! For nz, use value at nz-1
      where( air_slant_column( 1 : nzm1 ) <= .1_dk * largest )
        secchi( 1 : nzm1 ) = air_slant_column( 1 : nzm1 )                     &
                             / air_vertical_column( 1 : nzm1 )
      elsewhere
        secchi( 1 : nzm1 ) = rTWO
      endwhere
      secchi( nz ) = secchi( nzm1 )

      ! Lyman-Alpha parameterization, output values of O2 optical depth
      ! and O2 effective (equivalent) cross section
      if( this%has_la ) then
        call this%lymana_OD( o2scol, secchi, o2_optical_depth_la )
        do iw = this%ila, this%ila + nla - iONE
          o2_optical_depth( :, iw ) =                                         &
              o2_optical_depth_la( :, iw - this%ila + iONE )
        enddo
      endif

      ! Koppers' parameterization of the SR bands, output values of O2
      ! optical depth and O2 equivalent cross section
      if( this%has_srb ) then
        call this%schum_OD( o2scol, temperature%edge_val_, secchi,            &
                            o2_optical_depth_k )
        do iw = this%isrb, this%isrb + nsrb - iONE
          o2_optical_depth( :, iw ) =                                         &
              o2_optical_depth_k( :, iw - this%isrb + iONE )
        enddo
      endif

      deallocate( zGrid )
      deallocate( lambdaGrid )
      deallocate( temperature )

    endif has_la_srb

  end subroutine la_srb_OD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine la_srb_xs( this, grid_warehouse, profile_warehouse,              &
      air_vertical_column, air_slant_column, o2_cross_section,                &
      spherical_geometry )
    ! Computes equivalent optical depths for O2 absorption, and O2 effective
    ! absorption cross sections, parameterized in the Lyman-alpha and SR bands

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t

    class(la_sr_bands_t),      intent(inout) :: this
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse
    type(spherical_geometry_t),intent(inout) :: spherical_geometry

    real(dk), intent(in)    :: air_vertical_column(:), air_slant_column(:)
    real(dk), intent(inout) :: o2_cross_section(:,:)

    ! Local variables
    character(len=*), parameter :: Iam = 'la_srb xs: '

    integer :: nz, nzm1, iw
    real(dk)    :: secchi(size(air_slant_column))
    real(dk)    :: o2vcol(size(air_vertical_column))
    real(dk)    :: o2scol(size(air_slant_column))
    class(grid_t), pointer :: zGrid
    class(grid_t), pointer :: lambdaGrid
    class(profile_t), pointer :: temperature, O2_profile

    ! Lyman-alpha variables
    ! O2 optical depth and equivalent cross section in the Lyman-alpha
    ! region.
    ! Wavelengths for Lyman alpha and SRB parameterizations:
    real(dk) :: o2_cross_section_la( size( air_slant_column ), nla )

    ! grid on which Koppers' parameterization is defined
    ! O2 optical depth and equivalent cross section on Koppers' grid
    real(dk) :: o2_cross_section_k( size( air_slant_column ), nsrb )

    has_la_srb: if( this%has_la_srb ) then

      ! get specific grids and vertical profiles
      zGrid => grid_warehouse%get_grid( this%height_grid_ )
      lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
      temperature => profile_warehouse%get_profile( this%temperature_profile_ )

      nzm1 = zGrid%ncells_
      nz   = nzm1 +  iONE

      ! O2 slant column
      if( this%do_scaled_O2_ ) then
        o2scol(:) = this%O2_scale_factor_ * air_slant_column(:)
      else
        O2_profile => profile_warehouse%get_profile( this%O2_profile_ )
        call spherical_geometry%air_mass( O2_profile%exo_layer_dens_, o2vcol, &
                                          o2scol )
        deallocate( O2_profile )
      end if

      ! Effective secant of solar zenith angle.
      ! Use 2.0 if no direct sun (value for isotropic radiation)
      ! For nz, use value at nz-1
      where( air_slant_column( 1 : nzm1 ) <= .1_dk * largest )
        secchi( 1 : nzm1 ) = air_slant_column( 1 : nzm1 )                     &
                             / air_vertical_column( 1 : nzm1 )
      elsewhere
        secchi( 1 : nzm1 ) = rTWO
      endwhere
      secchi( nz ) = secchi( nzm1 )

      ! Lyman-Alpha parameterization, output values of O2 optical depth
      ! and O2 effective (equivalent) cross section
      if( this%has_la ) then
        call this%lymana_xs( o2scol,  secchi,  o2_cross_section_la )
        do iw = this%ila, this%ila + nla - iONE
          o2_cross_section( :, iw ) =                                         &
              o2_cross_section_la( :, iw - this%ila + iONE )
        enddo
      endif

      ! Koppers' parameterization of the SR bands, output values of O2
      ! optical depth and O2 equivalent cross section
      if( this%has_srb ) then
        call this%schum_xs( o2scol, temperature%edge_val_, secchi,            &
                            o2_cross_section_k )
        do iw = this%isrb, this%isrb + nsrb - iONE
          o2_cross_section( :, iw ) =                                         &
              o2_cross_section_k( :, iw - this%isrb + iONE )
        enddo
      endif

      deallocate( zGrid )
      deallocate( lambdaGrid )
      deallocate( temperature )

    endif has_la_srb

  end subroutine la_srb_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the calculator onto a
    ! buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(la_sr_bands_t), intent(in) :: this ! calculator to be packed
    integer,              intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    real(dk), allocatable :: ac(:,:), bc(:,:)

    ac = this%AC
    bc = this%BC
    pack_size = musica_mpi_pack_size( this%ila,              comm ) +         &
                musica_mpi_pack_size( this%isrb,             comm ) +         &
                musica_mpi_pack_size( this%has_la,           comm ) +         &
                musica_mpi_pack_size( this%has_srb,          comm ) +         &
                musica_mpi_pack_size( this%has_la_srb,       comm ) +         &
                musica_mpi_pack_size( this%do_scaled_O2_,    comm ) +         &
                musica_mpi_pack_size( this%O2_scale_factor_, comm ) +         &
                musica_mpi_pack_size( ac,                    comm ) +         &
                musica_mpi_pack_size( bc,                    comm ) +         &
                this%height_grid_%pack_size(                 comm ) +         &
                this%wavelength_grid_%pack_size(             comm ) +         &
                this%temperature_profile_%pack_size(         comm ) +         &
                this%O2_profile_%pack_size(                  comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the calculator onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(la_sr_bands_t), intent(in)    :: this      ! calculator to be packed
    character,            intent(inout) :: buffer(:) ! memory buffer
    integer,              intent(inout) :: position  ! current buffer position
    integer,              intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    real(dk), allocatable :: ac(:,:), bc(:,:)

    prev_pos = position
    ac = this%AC
    bc = this%BC
    call musica_mpi_pack( buffer, position, this%ila,              comm )
    call musica_mpi_pack( buffer, position, this%isrb,             comm )
    call musica_mpi_pack( buffer, position, this%has_la,           comm )
    call musica_mpi_pack( buffer, position, this%has_srb,          comm )
    call musica_mpi_pack( buffer, position, this%has_la_srb,       comm )
    call musica_mpi_pack( buffer, position, this%do_scaled_O2_,    comm )
    call musica_mpi_pack( buffer, position, this%O2_scale_factor_, comm )
    call musica_mpi_pack( buffer, position, ac,                    comm )
    call musica_mpi_pack( buffer, position, bc,                    comm )
    call this%height_grid_%mpi_pack(         buffer, position, comm )
    call this%wavelength_grid_%mpi_pack(     buffer, position, comm )
    call this%temperature_profile_%mpi_pack( buffer, position, comm )
    call this%O2_profile_%mpi_pack(          buffer, position, comm )
    call assert( 708555791, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks the calculator onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(la_sr_bands_t), intent(out)   :: this      ! calculator to be unpacked
    character,            intent(inout) :: buffer(:) ! memory buffer
    integer,              intent(inout) :: position  ! current buffer position
    integer,              intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    real(dk), allocatable :: ac(:,:), bc(:,:)

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%ila,              comm )
    call musica_mpi_unpack( buffer, position, this%isrb,             comm )
    call musica_mpi_unpack( buffer, position, this%has_la,           comm )
    call musica_mpi_unpack( buffer, position, this%has_srb,          comm )
    call musica_mpi_unpack( buffer, position, this%has_la_srb,       comm )
    call musica_mpi_unpack( buffer, position, this%do_scaled_O2_,    comm )
    call musica_mpi_unpack( buffer, position, this%O2_scale_factor_, comm )
    call musica_mpi_unpack( buffer, position, ac,                    comm )
    call musica_mpi_unpack( buffer, position, bc,                    comm )
    call this%height_grid_%mpi_unpack(         buffer, position, comm )
    call this%wavelength_grid_%mpi_unpack(     buffer, position, comm )
    call this%temperature_profile_%mpi_unpack( buffer, position, comm )
    call this%O2_profile_%mpi_unpack(          buffer, position, comm )
    this%ac = ac
    this%bc = bc
    call assert( 683591141, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lymana_OD( this,o2col,secchi,o2_optical_depth_la )
    ! Calculate the effective absorption cross section of O2 in the
    ! Lyman-Alpha bands and an effective O2 optical depth at all altitudes.
    ! Parameterized after:  Chabrillat, S., and G. Kockarts, Simple
    ! parameterization of the absorption of the solar Lyman-Alpha line,
    ! Geophysical Research Letters, Vol.24, No.21, pp 2659-2662, 1997.

    class(la_sr_bands_t), intent(inout)  :: this
    real(dk), intent(in) :: o2col(:) ! Slant overhead O2 column (molec/cc) at each specified altitude
    real(dk), intent(in) :: secchi(:)

    real(dk), intent(out) :: o2_optical_depth_la(:,:)

    ! Local variables
    real(dk), parameter :: xsmin     = 1.e-20_dk
    real(dk), parameter :: tiny_val  = 1.e-100_dk
    real(dk), parameter :: exp_lim   = 100.e8_dk
    real(dk), parameter :: large_od  = 1000._dk
    real(dk), parameter :: b(3) =                                             &
        (/  6.8431e-01_dk, 2.29841e-01_dk,  8.65412e-02_dk /)
    real(dk), parameter :: c(3) =                                             &
        (/ 8.22114e-21_dk, 1.77556e-20_dk,  8.22112e-21_dk /)
    real(dk), parameter :: d(3) =                                             &
        (/  6.0073e-21_dk, 4.28569e-21_dk,  1.28059e-20_dk /)
    real(dk), parameter :: e(3) =                                             &
        (/ 8.21666e-21_dk, 1.63296e-20_dk,  4.85121e-17_dk /)

    integer :: nz
    integer :: iz
    real(dk) :: coldens
    real(dk) :: sigma(3), tau(3)
    real(dk) :: rm( size( o2col ) )

    ! calculate reduction factors at every altitude
    rm(:)  = rZERO
    nz = size( o2col )
    do iz = iONE, nz
      coldens = o2col( iz )
      sigma   = c * coldens
      where( sigma < exp_lim )
        tau = exp( -sigma )
      elsewhere
        tau = rZERO
      endwhere
      rm( iz ) = dot_product( tau, b )
    enddo

    ! calculate effective O2 optical depth
    do iz = iONE, nz - iONE
      if( rm( iz ) > tiny_val ) then
        if( rm( iz + iONE ) > rZERO ) then
          o2_optical_depth_la( iz, iONE ) =                                   &
              log( rm( iz + 1 ) ) / secchi( iz + 1 )                          &
              - log( rm( iz ) ) / secchi( iz )
        else
          o2_optical_depth_la( iz, iONE ) = large_od
        endif
      else
        o2_optical_depth_la( iz, iONE ) = large_od
      endif
    enddo

  end subroutine lymana_OD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lymana_xs( this,o2col,secchi,o2_cross_section_la )
    ! Calculates the effective absorption cross section of O2 in the
    ! Lyman-Alpha bands and an effective O2 optical depth at all altitudes.
    ! Parameterized after:  Chabrillat, S., and G. Kockarts, Simple
    ! parameterization of the absorption of the solar Lyman-Alpha line,
    ! Geophysical Research Letters, Vol.24, No.21, pp 2659-2662, 1997.

    class(la_sr_bands_t), intent(inout)  :: this
    real(dk), intent(in) :: o2col(:) ! Slant overhead O2 column (molec/cc) at each specified altitude
    real(dk), intent(in) :: secchi(:)

    real(dk), intent(out) :: o2_cross_section_la(:,:)

    ! Local variables
    real(dk), parameter :: xsmin     = 1.e-20_dk
    real(dk), parameter :: tiny_val  = 1.e-100_dk
    real(dk), parameter :: exp_lim   = 100.e8_dk
    real(dk), parameter :: large_od  = 1000._dk
    real(dk), parameter :: b(3) =                                             &
        (/  6.8431e-01_dk, 2.29841e-01_dk,  8.65412e-02_dk /)
    real(dk), parameter :: c(3) =                                             &
        (/ 8.22114e-21_dk, 1.77556e-20_dk,  8.22112e-21_dk /)
    real(dk), parameter :: d(3) =                                             &
        (/  6.0073e-21_dk, 4.28569e-21_dk,  1.28059e-20_dk /)
    real(dk), parameter :: e(3) =                                             &
        (/ 8.21666e-21_dk, 1.63296e-20_dk,  4.85121e-17_dk /)

    integer :: nz
    integer :: iz
    real(dk) :: coldens
    real(dk) :: sigma(3), tau(3)
    real(dk) :: rm( size( o2col ) ), ro2( size( o2col ) )

    ! calculate reduction factors at every altitude
    rm(:)  = rZERO
    ro2(:) = rZERO
    nz = size( o2col )
    do iz = iONE, nz
      coldens = o2col( iz )
      sigma   = c * coldens
      where( sigma < exp_lim )
        tau = exp( -sigma )
      elsewhere
        tau = rZERO
      endwhere
      rm( iz ) = dot_product( tau, b )

      sigma   = e * coldens
      where( sigma < exp_lim )
        tau = exp( -sigma )
      elsewhere
        tau = rZERO
      endwhere
      ro2( iz ) = dot_product( tau, d )
    enddo

    ! calculate effective O2 cross sections
    do iz = iONE, nz - iONE
      if( rm( iz ) > tiny_val ) then
        if( ro2( iz ) > tiny_val ) then
          o2_cross_section_la( iz, iONE ) = ro2( iz ) / rm( iz )
        else
          o2_cross_section_la( iz, iONE ) = xsmin
        endif
      else
        o2_cross_section_la( iz, iONE ) = xsmin
      endif
    enddo

    ! do top layer separately
    if( rm( nz ) > tiny_val ) then
      o2_cross_section_la( nz, 1 ) = ro2( nz ) / rm( nz )
    else
      o2_cross_section_la( nz, 1 ) = xsmin
    endif

  end subroutine lymana_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine schum_OD( this, o2col, tlev, secchi, o2_optical_depth_k )
    ! Calculate the equivalent absorption cross section of O2 in the SR bands.
    ! The algorithm is based on parameterization of G.A. Koppers, and
    ! D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]
    ! Final values do include effects from the Herzberg continuum.

    class(la_sr_bands_t), intent(inout)  :: this
    real(dk), intent(in) :: o2col(:)  ! Slant overhead O2 column (molec/cc) at each specified altitude
    real(dk), intent(in) :: tlev(:)   ! Temperature at each level
    real(dk), intent(in) :: secchi(:) ! Ratio of slant to vertical O2 columns
    real(dk), intent(out) :: o2_optical_depth_k(:,:)

    ! Local variables
    real(dk), parameter :: colmin = exp( 38._dk )
    real(dk), parameter :: xslod(nsrb) =                                      &
    (/ 6.2180730E-21_dk, 5.8473627E-22_dk, 5.6996334E-22_dk,                  &
       4.5627094E-22_dk, 1.7668250E-22_dk, 1.1178808E-22_dk,                  &
       1.2040544E-22_dk, 4.0994668E-23_dk, 1.8450616E-23_dk,                  &
       1.5639540E-23_dk, 8.7961075E-24_dk, 7.6475608E-24_dk,                  &
       7.6260556E-24_dk, 7.5565696E-24_dk, 7.6334338E-24_dk,                  &
       7.4371992E-24_dk, 7.3642966E-24_dk /)

    integer :: nz
    integer :: k, kp1, ktop, ktop1, kbot, lambdaNdx

    real(dk) :: X, NORM
    real(dk) :: o2col1( size( o2col ) )
    real(dk), allocatable :: o2_cross_section_k( :, : )

    allocate( o2_cross_section_k( size( o2col ), nsrb ) )

    ! Initialize cross sections to values at large optical depth
    do lambdaNdx = iONE, nsrb
      o2_cross_section_k( :, lambdaNdx ) = xslod( lambdaNdx )
    enddo

    ! Calculate cross sections
    ! Set smallest O2col = exp(38.) molec cm-2
    ! to stay in range of parameterization
    ! given by Koppers et al. at top of atm.
    nz = size( o2col )
    ktop = nz
    kbot = 0
    NORM = rONE/ real( nz - iONE, dk )
    do k = iONE, nz
      o2col1( k ) = max( o2col( k ), colmin )
      x  = log( o2col1( k ) )
      if( x < 38.0_dk ) then
        ktop1 = k - iONE
        ktop  = min( ktop1, ktop )
      else if( x > 56.0_dk ) then
        kbot = k
      else
        o2_cross_section_k( k, : nsrb ) = this%effxs( x, tlev( k ) )
      endif
    enddo

    ! Fill in cross section where X is out of range
    ! by repeating edge table values.
    ! Do not allow kbot = nz to avoid division by zero in no light case.
    if( kbot == nz) then
      kbot = nz - iONE
    endif

    do lambdaNdx = iONE, nsrb
      o2_cross_section_k( iONE : kbot, lambdaNdx ) =                          &
          o2_cross_section_k( kbot + iONE, lambdaNdx )
      o2_cross_section_k( ktop + iONE : nz, lambdaNdx ) =                     &
          o2_cross_section_k( ktop, lambdaNdx )
    enddo

    !  Calculate incremental optical depth
    do lambdaNdx = iONE, nsrb
      do k = iONE, nz - iONE
        kp1 = k + iONE

        !... calculate an optical depth weighted by density
        ! put in mean value estimate, if in shade
        if( abs( rONE - o2col1( kp1 ) / o2col1( k ) ) <= rTWO * precis ) then
          o2_optical_depth_k( k, lambdaNdx ) =                                &
              o2_cross_section_k( kp1, lambdaNdx ) * o2col1( kp1 ) * NORM
        else
          o2_optical_depth_k( k, lambdaNdx ) =                                &
              abs( ( o2_cross_section_k( kp1, lambdaNdx ) * o2col1( kp1 )     &
                     - o2_cross_section_k( k, lambdaNdx ) * o2col1( k ) )     &
                   / ( rONE + log( o2_cross_section_k( kp1, lambdaNdx )       &
                                   / o2_cross_section_k( k, lambdaNdx ) )     &
                              / log( o2col1( kp1 ) / o2col1( k ) ) ) )

          !... change to vertical optical depth
          o2_optical_depth_k( k, lambdaNdx ) =                                &
              rTWO * o2_optical_depth_k( k, lambdaNdx )                       &
              / ( secchi( k ) + secchi( kp1 ) )
        endif
      enddo
    enddo

  end subroutine schum_OD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine schum_xs( this, o2col, tlev, secchi, o2_cross_section_k )
    ! Calculate the equivalent absorption cross section of O2 in the SR bands.
    ! The algorithm is based on parameterization of G.A. Koppers, and
    ! D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]
    ! Final values do include effects from the Herzberg continuum.

    class(la_sr_bands_t), intent(inout)  :: this
    real(dk), intent(in) :: o2col(:)  ! Slant overhead O2 column (molec/cc) at each specified altitude
    real(dk), intent(in) :: tlev(:)   ! Temperature at each level
    real(dk), intent(in) :: secchi(:) ! Ratio of slant to vertical O2 columns

    real(dk), intent(out) :: o2_cross_section_k(:,:)

    ! Local variables
    real(dk), parameter :: colmin = exp( 38._dk )
    real(dk), parameter :: xslod( nsrb ) =                                    &
    (/ 6.2180730E-21_dk, 5.8473627E-22_dk, 5.6996334E-22_dk,                  &
       4.5627094E-22_dk, 1.7668250E-22_dk, 1.1178808E-22_dk,                  &
       1.2040544E-22_dk, 4.0994668E-23_dk, 1.8450616E-23_dk,                  &
       1.5639540E-23_dk, 8.7961075E-24_dk, 7.6475608E-24_dk,                  &
       7.6260556E-24_dk, 7.5565696E-24_dk, 7.6334338E-24_dk,                  &
       7.4371992E-24_dk, 7.3642966E-24_dk /)

    integer :: nz
    integer :: k, ktop, ktop1, kbot, lambdaNdx

    real(dk) :: X
    real(dk) :: o2col1( size( o2col ) )

    ! Initialize cross sections to values at large optical depth
    do lambdaNdx = iONE, nsrb
      o2_cross_section_k( :, lambdaNdx ) = xslod( lambdaNdx )
    enddo

    ! Calculate cross sections
    ! Set smallest O2col = exp(38.) molec cm-2
    ! to stay in range of parameterization
    ! given by Koppers et al. at top of atm.
    nz = size( o2col )
    ktop = nz
    kbot = 0
    do k = iONE, nz
      o2col1( k ) = max( o2col( k ), colmin )
      x = log( o2col1( k ) )
      if( x < 38.0_dk ) then
        ktop1 = k - 1
        ktop  = min( ktop1, ktop )
      else if ( x > 56.0_dk ) then
        kbot = k
      else
        o2_cross_section_k( k, : nsrb ) = this%effxs( x, tlev( k ) )
      endif
    enddo

    ! Fill in cross section where X is out of range
    ! by repeating edge table values
    ! Do not allow kbot = nz to avoid division by zero in
    ! no light case.
    if( kbot == nz ) then
    kbot = nz - iONE
    endif

    do lambdaNdx = iONE, nsrb
      o2_cross_section_k( iONE : kbot, lambdaNdx ) =                          &
          o2_cross_section_k( kbot + iONE, lambdaNdx )
      o2_cross_section_k( ktop + iONE : nz, lambdaNdx ) =                     &
          o2_cross_section_k( ktop, lambdaNdx )
    enddo

  end subroutine schum_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_srb_xs( this, file_path )
    !	polynomial coeffs necessary to calculate O2 effective cross-sections

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t

    class(la_sr_bands_t), intent(inout)  :: this
    type(string_t),       intent(in)     :: file_path ! Path to cross section parameters file

    integer :: in_lun ! file unit number
    integer :: i, j, istat

    in_lun = 11

    open( unit = in_lun, file = file_path%to_char( ), form = 'formatted',     &
          iostat = istat )
    call assert_msg( 245202775, istat == 0,                                   &
                     'failed to open ' // file_path )

    read( in_lun, 901 )
    do I = iONE, nPoly
      read( in_lun, 903, iostat = istat ) ( this%AC( I, J ), J = 1, nsrb )
      call assert_msg( 176835602, istat == 0,                                 &
                       'failed to read ' // file_path )
    enddo
    read( in_lun, 901 )
    do I = iONE, nPoly
      read( in_lun, 903, iostat = istat ) ( this%BC( I, J ), J = 1, nsrb )
      call assert_msg( 685853075, istat == 0,                                 &
                       'failed to read ' // file_path )
    enddo

    901  format( / )
    903  format( 17( E23.14, 1x ) )

    CLOSE( in_lun )

  end subroutine init_srb_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function effxs( this, X, T ) result( xs )
    ! Evaluates the effective cross section
    ! of O2 in the Schumann-Runge bands using parameterization
    ! of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
    ! 68-79, 1996]
    !
    ! method:
    ! ln(xs) = A(X)[T-220]+B(X)
    ! X = log of slant column of O2
    ! A,B calculated from Chebyshev polynomial coeffs
    ! AC and BC using NR routine chebev.  Assume interval
    ! is 38<ln(NO2)<56.
    !
    ! Revision History:
    !
    ! drm 2/97  initial coding

    class(la_sr_bands_t), intent(inout)  :: this
    real(dk), intent(in)  :: T
    real(dk), intent(in)  :: X
    real(dk), allocatable :: XS(:)

    ! Local variables
    real(dk), parameter :: T0 = 220._dk
    real(dk)            :: A( nsrb ), B( nsrb )

    call this%calc_params( X, A, B )

    XS = exp( A * ( T - T0 ) + B )

  end function effxs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine calc_params( this, X, A, B )
    ! Calculates coefficients (A,B), used in calculating the effective 
    ! cross section, for 17 wavelength intervals as a function of log O2 
    ! column density (X), Wavelength intervals are defined in WMO1985

    class(la_sr_bands_t), intent(inout)  :: this
    real(dk), intent(in)  :: X
    real(dk), intent(out) :: A(:)
    real(dk), intent(out) :: B(:)

    integer :: I

    ! call Chebyshev Evaluation routine to calc A and B from
    !	set of 20 coeficients for each wavelength
    do I = 1,size( A )
      A(I) = this%chebyshev_evaluation( this%AC( 1, I ), X )
      B(I) = this%chebyshev_evaluation( this%BC( 1, I ), X )
    enddo

  end subroutine calc_params

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function chebyshev_evaluation( this, coefficients, x ) result( ret_val )
    ! Calculates the value of the Chebyshev polynomial at
    ! y = (x - lower_limit/2 - upper_limit/2)/((upper_limit-lower_limit)/2)
    ! to approximate f(x)
    use musica_assert,                 only : assert_msg

    class(la_sr_bands_t), intent(in) :: this
    real(dk),             intent(in) :: coefficients( nPoly )
    real(dk),             intent(in) :: x

    real(dk) :: ret_val

    integer  :: i_coeff
    real(dk) :: di, di1, dtemp, y, y2

    call assert_msg( 560585657, x >= kLowerLimit .and. x <= kUpperLimit,     &
                     "x out-of-bounds for Chebyshev polynomial" )

    di  = rZERO
    di1 = rZERO
    y  = ( rTWO * x - ( kLowerLimit + kUpperLimit ) )                         &
         / ( kUpperLimit - kLowerLimit )
    y2 = rTWO * y
    do i_coeff = nPoly, 2, -1
      dtemp = di
      di  = y2 * di - di1 + coefficients( i_coeff )
      di1 = dtemp
    enddo
    ret_val = y * di - di1 + 0.5_dk * coefficients( 1 )

  end function chebyshev_evaluation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )

    type(la_sr_bands_t), intent(inout) :: this

    ! nothing to do for now

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function get_band_min_index( band_name, wavelengths ) result( idx )
    ! Returns the minimum wavelength bin index that includes the specified
    ! band

    use musica_assert,                 only : die_msg, almost_equal
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    type(string_t), intent(in) :: band_name
    class(grid_t),  intent(in) :: wavelengths

    integer :: i_wl
    real(kind=dk) :: search_val

    select case( band_name%to_char( ) )
    case( "lyman-alpha" )
      search_val = wlla(1)
    case( "schumann-runge" )
      search_val = wlsrb(1)
    case( "schumann-runge continuum" )
      search_val = wlla(kla)
    case default
      call die_msg( 943741887, "Unknown wavelength band '"//band_name//"'" )
    end select
    idx = 0
    do i_wl = 1, wavelengths%size( ) + 1
      if( almost_equal( wavelengths%edge_( i_wl ), search_val ) ) then
        idx = i_wl
        return
      end if
    end do
    call die_msg( 651190319, "Wavelength grid does not map to band '"//       &
                             band_name//"'" )

  end function get_band_min_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function get_band_max_index( band_name, wavelengths ) result( idx )
    ! Returns the maximum wavelength bin index that includes the specified
    ! band

    use musica_assert,                 only : die_msg, almost_equal
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t

    type(string_t), intent(in) :: band_name
    class(grid_t),  intent(in) :: wavelengths

    integer :: i_wl
    real(kind=dk) :: search_val

    select case( band_name%to_char( ) )
    case( "lyman-alpha" )
      search_val = wlla(kla)
    case( "schumann-runge" )
      search_val = wlsrb(ksrb)
    case( "schumann-runge continuum" )
      search_val = wlsrb(1)
    case default
      call die_msg( 943741887, "Unknown wavelength band '"//band_name//"'" )
    end select
    idx = 0
    do i_wl = 1, wavelengths%size( ) + 1
      if( almost_equal( wavelengths%edge_( i_wl ), search_val ) ) then
        idx = i_wl - 1
        return
      end if
    end do
    call die_msg( 316052437, "Wavelength grid does not map to band '"//       &
                             band_name//"'" )

  end function get_band_max_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_la_sr_bands
