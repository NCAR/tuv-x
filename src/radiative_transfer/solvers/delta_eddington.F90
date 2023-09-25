! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_solver_delta_eddington

  use tuvx_solver,                    only : solver_t, radiation_field_t
  use musica_constants,               only : dk => musica_dk
  use tuvx_constants,                 only : pi
  use tuvx_grid_warehouse,            only : grid_warehouse_ptr
  use tuvx_profile_warehouse,         only : profile_warehouse_ptr

  implicit none

  private
  public :: solver_delta_eddington_t

  type, extends(solver_t) :: solver_delta_eddington_t
     ! Radiative flux calculator that applies the delta-Eddington Approximation.
     !
     ! Solves two-stream equations for multiple layers. These routines are based
     ! on equations from: Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.
     ! It contains 9 two-stream methods to choose from. A pseudo-spherical
     ! correction has also been added.
     !
     ! The original delta-Eddington paper is:
     ! Joseph and Wiscombe, J. Atmos. Sci., 33, 2453-2459, 1976
     type(grid_warehouse_ptr) :: height_grid_
     type(grid_warehouse_ptr) :: wavelength_grid_
     type(profile_warehouse_ptr) :: surface_albedo_profile_
  contains
    procedure :: update_radiation_field
    procedure :: pack_size
    procedure :: mpi_pack
    procedure :: mpi_unpack
  end type solver_delta_eddington_t

  interface solver_delta_eddington_t
    procedure :: constructor
  end interface solver_delta_eddington_t

  real(dk), parameter :: rZERO = 0.0_dk
  real(dk), parameter :: rONE  = 1.0_dk
  real(dk), parameter :: rTWO  = 2.0_dk
  real(dk), parameter :: d2r   = pi/180._dk

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructs a delta Eddington solver
  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( solver )

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(config_t),                  intent(inout) :: config
    class(solver_delta_eddington_t), pointer       :: solver
    type(grid_warehouse_t),          intent(in)    :: grid_warehouse
    type(profile_warehouse_t),       intent(in)    :: profile_warehouse

    type(string_t) :: required_keys(1), optional_keys(0)

    required_keys(1) = "type"

    call assert_msg( 657111982,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration format for delta Eddington solver" )

    allocate( solver )
    solver%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    solver%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    solver%surface_albedo_profile_ =                                          &
        profile_warehouse%get_ptr( "surface albedo", "none" )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function update_radiation_field( this, solar_zenith_angle, n_layers,        &
      spherical_geometry, grid_warehouse, profile_warehouse,                  &
      radiator_warehouse ) result( radiation_field )

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_linear_algebra_linpack,   only : linear_algebra_linpack_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_radiator,                 only : radiator_state_t
    use tuvx_radiator_warehouse,       only : radiator_warehouse_t
    use tuvx_solver,                   only : slant_optical_depth
    use tuvx_spherical_geometry,       only : spherical_geometry_t

    class(solver_delta_eddington_t), intent(inout) :: this ! Delta-Eddington solver

    integer,                    intent(in)    :: n_layers  ! number of vertical layers
    real(dk),                   intent(in)    :: solar_zenith_angle ! solar zenith angle [degrees]
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse
    type(profile_warehouse_t),  intent(inout) :: profile_warehouse
    type(radiator_warehouse_t), intent(inout) :: radiator_warehouse
    type(spherical_geometry_t), intent(inout) :: spherical_geometry

    type(radiation_field_t),   pointer       :: radiation_field

    ! Local variables
    character(len=*), parameter :: Iam = 'Update radiation field: '
    real(dk) :: mu
    real(dk) :: tausla( 0 : n_layers ), tauc( 0 : n_layers )
    real(dk) :: mu2( 0 : n_layers )

    ! internal coefficients and matrix
    integer     :: row
    real(dk)    :: lam( n_layers ), taun( n_layers ), bgam( n_layers )
    real(dk)    :: e1( n_layers ), e2( n_layers )
    real(dk)    :: e3( n_layers ), e4( n_layers )
    real(dk)    :: cup( n_layers ), cdn( n_layers )
    real(dk)    :: cuptn( n_layers ), cdntn( n_layers )
    real(dk)    :: mu1( n_layers )
    real(dk)    :: a( 2 * n_layers ), b( 2 * n_layers ), d( 2 * n_layers )
    real(dk)    :: e( 2 * n_layers ), y( 2 * n_layers )

    real(dk) :: pifs, fdn0, surfem, tempg
    real(dk) :: f, g, om
    real(dk) :: gam1, gam2, gam3, gam4
    real(dk) :: gi(n_layers), omi(n_layers)

    integer     :: mrows, lev
    integer     :: i, j
    real(dk) :: expon, expon0, expon1, divisr, temp, up, dn
    real(dk) :: ssfc

    ! Linear algebra package, radiation field type
    type(linear_algebra_linpack_t) :: linpack

    ! Local variables
    real(dk), parameter                  :: largest = 1.e36_dk
    real(dk), parameter                  :: kfloor = rONE/largest
    real(dk), parameter                  :: precis = 1.e-7_dk
    real(dk), parameter                  :: eps    = 1.e-3_dk

    integer                              :: nlambda, lambdaNdx
    type(radiator_state_t)               :: atmRadiatorState
    class(grid_t),    pointer            :: zGrid
    class(grid_t),    pointer            :: lambdaGrid
    class(profile_t), pointer            :: surfaceAlbedo

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    surfaceAlbedo =>                                                          &
        profile_warehouse%get_profile( this%surface_albedo_profile_ )

    nlambda = lambdaGrid%ncells_
    radiation_field => radiation_field_t( n_layers + 1, nlambda )

    allocate( atmRadiatorState%layer_G_( n_layers, nlambda, 1 ) )
    ! Create cumulative state from all radiators
    call radiator_warehouse%accumulate_states( atmRadiatorState )

    ! MU = cosine of solar zenith angle
    ! RSFC = surface albedo
    ! TAUU =  unscaled optical depth of each layer
    ! OMU  =  unscaled single scattering albedo
    ! GU   =  unscaled asymmetry factor
    ! N_LAYERS = number of layers in the atmosphere
    ! N_LEVELS = nlayer + 1 = number of levels

    mu = cos( solar_zenith_angle*d2r )
    associate( nid  => spherical_geometry%nid_,                               &
               dsdh => spherical_geometry%dsdh_ )

    wavelength_loop: do lambdaNdx = 1, nlambda
      associate( rsfc => surfaceAlbedo%mid_val_( lambdaNdx ),                 &
             tauu => atmRadiatorState%layer_OD_( n_layers:1:-1, lambdaNdx ),  &
             omu  => atmRadiatorState%layer_SSA_( n_layers:1:-1, lambdaNdx ), &
             gu   => atmRadiatorState%layer_G_( n_layers:1:-1, lambdaNdx, 1 ) )

      ! initial conditions:  pi*solar flux = 1;  diffuse incidence = 0
      pifs = rONE
      fdn0 = rZERO
      ! emission at surface (for night light pollution, set pifs = 0, surfem = 1.)
      surfem = rZERO
      !************* compute coefficients for each layer:
      ! GAM1 - GAM4 = 2-stream coefficients, different for different approximations
      ! expON0 = calculation of e when TAU is zero
      ! expON1 = calculation of e when TAU is TAUN
      ! CUP and CDN = calculation when TAU is zero
      ! CUPTN and CDNTN = calc. when TAU is TAUN
      ! DIVISR = prevents division by zero
      tauc   = rZERO
      tausla = rZERO
      mu2    = rONE / sqrt( largest )
      ! delta-scaling. Has to be done for delta-Eddington approximation,
      ! delta discrete ordinate, Practical Improved Flux Method, delta function,
      ! and Hybrid modified Eddington-delta function methods approximations

      do i = 1, n_layers
        f         = gu( i ) * gu( i )
        gi( i )   = ( gu( i ) - f ) / ( rONE - f )
        omi( i )  = ( rONE - f ) * omu( i ) / ( rONE - omu( i ) * f )
        taun( i ) = ( rONE - omu( i ) * f ) * tauu( i )
      end do

      ! calculate slant optical depth at the top of the atmosphere when
      ! the solar zenith angle is > 90 degrees.
      if( mu < rZERO ) then
        tausla(0) = slant_optical_depth( 0, nid(0), dsdh(0,:), taun )
      end if

      layer_loop: do i = 1, n_layers

        g  = gi( i )
        om = omi( i )
        tauc( i ) = tauc( i - 1 ) + taun( i )

        ! stay away from 1 by precision.  For g, also stay away from -1

        tempg = min( abs( g ), rONE - precis )
        g = sign( tempg, g )
        om = min( om, rONE - precis )

        ! calculate slant optical depth
        tausla( i ) = slant_optical_depth( i, nid( i ), dsdh( i, : ), taun )

        if( nid( i ) >= 0 ) then
          if( tausla( i ) == tausla( i - 1 ) ) then
            mu2( i ) = sqrt( largest )
          else
            mu2( i ) = ( tauc( i ) - tauc( i - 1 ) )                          &
                       / ( tausla( i ) - tausla( i - 1 ) )
            mu2( i ) = sign( max( abs( mu2( i ) ), rONE / sqrt( largest ) ),  &
                             mu2( i ) )
          end if
        end if

        !** the following gamma equations are from pg 16,289, Table 1
        !** save mu1 for each approx. for use in converting irradiance to actinic flux
        ! Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):

        gam1 =   ( 7._dk - om * ( 4._dk + 3._dk * g ) ) / 4._dk
        gam2 = - ( rONE - om * ( 4._dk - 3._dk * g ) ) / 4._dk
        gam3 =   ( rTWO - 3._dk * g * mu ) / 4._dk
        gam4 =   rONE - gam3
        mu1( i ) = 0.5_dk

        lam( i ) = sqrt( gam1 * gam1 - gam2 * gam2 )

        if( gam2 /= rZERO) then
          bgam( i ) = ( gam1 - lam( i ) ) / gam2
        else
          bgam( i ) = rZERO
        endif

        expon = exp( - lam( i ) * taun( i ) )

        ! e1 - e4 = pg 16,292 equation 44

        e1( i ) = rONE + bgam( i ) * expon
        e2( i ) = rONE - bgam( i ) * expon
        e3( i ) = bgam( i ) + expon
        e4( i ) = bgam( i ) - expon

        ! the following sets up for the C equations 23, and 24
        ! found on page 16,290
        ! prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
        ! which is approx equiv to shifting MU by 0.5*EPS* (MU)**3

        expon0 = exp( -tausla( i - 1 ) )
        expon1 = exp( -tausla( i ) )

        divisr = lam( i ) * lam( i ) - rONE / ( mu2( i ) * mu2( i ) )
        temp   = max( eps, abs( divisr ) )
        divisr = sign( temp, divisr )

        up = om * pifs * ( ( gam1 - rONE / mu2( i ) ) * gam3 + gam4 * gam2 ) &
             / divisr
        dn = om * pifs * ( ( gam1 + rONE / mu2( i ) ) * gam4 + gam2 * gam3 ) &
             / divisr

        ! cup and cdn are when tau is equal to zero
        ! cuptn and cdntn are when tau is equal to taun

        cup(i) = up*expon0
        cdn(i) = dn*expon0
        cuptn(i) = up*expon1
        cdntn(i) = dn*expon1

      enddo layer_loop

      !**************** set up matrix ******
      ! ssfc = pg 16,292 equation 37  where pi Fs is one (unity).

      ssfc = rsfc * mu * exp( -tausla( n_layers ) ) * pifs + surfem

      ! MROWS = the number of rows in the matrix

      mrows = 2 * n_layers

      ! the following are from pg 16,292  equations 39 - 43.
      ! set up first row of matrix:

      a(1) = rZERO
      b(1) = e1(1)
      d(1) = -e2(1)
      e(1) = fdn0 - cdn(1)

      ! set up odd rows 3 thru (MROWS - 1):

      i = 0
      do row = 3, mrows - 1, 2
         i = i + 1
         a( row ) = e2( i ) * e3( i ) - e4( i ) * e1( i )
         b( row ) = e1( i ) * e1( i + 1 ) - e3( i ) * e3( i + 1 )
         d( row ) = e3( i ) * e4( i + 1 ) - e1( i ) * e2( i + 1 )
         e( row ) = e3( i ) * ( cup( i + 1 ) - cuptn( i ) )                   &
                    + e1( i ) * ( cdntn( i ) - cdn( i + 1 ) )
      enddo

      ! set up even rows 2 thru (MROWS - 2):

      i = 0
      do row = 2, mrows - 2, 2
         i = i + 1
         a( row ) = e2( i + 1 ) * e1( i ) - e3( i ) * e4( i + 1 )
         b( row ) = e2( i ) * e2( i + 1 ) - e4( i ) * e4( i + 1 )
         d( row ) = e1( i + 1 ) * e4( i + 1 ) - e2( i + 1) * e3( i + 1 )
         e( row ) = ( cup( i + 1 ) - cuptn( i ) ) * e2( i + 1 )               &
                    - ( cdn( i + 1 ) - cdntn( i ) ) * e4( i + 1 )
      enddo

      ! set up last row of matrix at MROWS:

      a( mrows ) = e1( n_layers ) - rsfc * e3( n_layers )
      b( mrows ) = e2( n_layers ) - rsfc * e4( n_layers )
      d( mrows ) = rZERO
      e( mrows ) = ssfc - cuptn( n_layers ) + rsfc * cdntn( n_layers )

      ! solve tri-diagonal system:

      y = linpack%tridiag( a, b, d, e )

      !*** unfold solution of matrix, compute output fluxes:
      ! the following equations are from pg 16,291  equations 31 & 32

      associate( edr => radiation_field%edr_( :, lambdaNdx ),                 &
                 eup => radiation_field%eup_( :, lambdaNdx ),                 &
                 edn => radiation_field%edn_( :, lambdaNdx ),                 &
                 fdr => radiation_field%fdr_( :, lambdaNdx ),                 &
                 fup => radiation_field%fup_( :, lambdaNdx ),                 &
                 fdn => radiation_field%fdn_( :, lambdaNdx ) )
      fdr(1) = pifs * exp( -tausla(0) )
      edr(1) = mu * fdr(1)
      edn(1) = fdn0
      eup(1) =  y(1) * e3(1) - y(2) * e4(1) + cup(1)
      fdn(1) = edn(1) / mu1(1)
      fup(1) = eup(1) / mu1(1)

      j   = 1
      row = 1
      do lev = 2, n_layers + 1
         fdr( lev ) = pifs * exp( -tausla( lev - 1 ) )
         edr( lev ) = mu * fdr( lev )
         edn( lev ) = y( row ) * e3( j ) + y( row + 1 ) * e4( j ) + cdntn( j )
         eup( lev ) = y( row ) * e1( j ) + y( row + 1 ) * e2( j ) + cuptn( j )
         fdn( lev ) = edn( lev ) / mu1( j )
         fup( lev ) = eup( lev ) / mu1( j )

         row = row + 2
         j = j + 1
      enddo
      ! transform from top-down to buttom-up
      fdr = fdr( n_layers + 1 : 1 : -1 )
      fup = fup( n_layers + 1 : 1 : -1 )
      fdn = fdn( n_layers + 1 : 1 : -1 )
      edr = edr( n_layers + 1 : 1 : -1 )
      eup = eup( n_layers + 1 : 1 : -1 )
      edn = edn( n_layers + 1 : 1 : -1 )

      end associate

    end associate

    enddo wavelength_loop

    end associate

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( surfaceAlbedo )

  end function update_radiation_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the solver

    use musica_mpi,                    only : musica_mpi_pack_size

    class(solver_delta_eddington_t), intent(in) :: this ! solver to be packed
    integer,                         intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = this%height_grid_%pack_size(            comm ) +              &
                this%wavelength_grid_%pack_size(        comm ) +              &
                this%surface_albedo_profile_%pack_size( comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the solver onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(solver_delta_eddington_t), intent(in)    :: this      ! solver to be packed
    character,                       intent(inout) :: buffer(:) ! memory buffer
    integer,                         intent(inout) :: position  ! current buffer position
    integer,                         intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    prev_pos = position

    call this%height_grid_%mpi_pack(            buffer, position, comm )
    call this%wavelength_grid_%mpi_pack(        buffer, position, comm )
    call this%surface_albedo_profile_%mpi_pack( buffer, position, comm )
    call assert( 485414316, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a solver from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(solver_delta_eddington_t), intent(out)   :: this      ! solver to be packed
    character,                       intent(inout) :: buffer(:) ! memory buffer
    integer,                         intent(inout) :: position  ! current buffer position
    integer,                         intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    prev_pos = position

    call this%height_grid_%mpi_unpack(            buffer, position, comm )
    call this%wavelength_grid_%mpi_unpack(        buffer, position, comm )
    call this%surface_albedo_profile_%mpi_unpack( buffer, position, comm )
    call assert( 764530792, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_solver_delta_eddington
