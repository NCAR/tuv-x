!> Driver for cross section / quantum yield tests based on Doug's version of TUV
module tuv_doug


  ! Parameters from params file

  ! output file unit number
  integer, parameter :: kout = 53
  ! input file unit number
  integer, parameter :: kin = 12

  ! altitude grid size
  integer, parameter :: kz = 151
  ! wavelength grid size
  integer, parameter :: kw = 700

  ! nubmer of wavelength-dependent weighting functions
  integer, parameter :: ks = 60
  ! Number of wavelength and altitude-dependent weighting functions
  integer, parameter :: kj = 150

  ! delta for adding points at beginning or end of data grids
  real, parameter :: deltax = 1.E-4

  ! Pi
  real, parameter ::  pi = 3.1415926535898
  ! radius of the earth:
  real, parameter :: radius = 6.371E+3
  ! largest number of the machine:
  real, parameter :: largest = 1.E+36
  ! small positive number
  real, parameter :: pzero = +10./largest
  ! small negative number
  real, parameter :: nzero = -10./largest
  !  machine precision
  real, parameter :: precis = 1.e-7

  ! TUV grids

  ! altitude grid
  integer :: nz
  real :: z(kz)

  ! wavelength grid
  integer :: nw
  real :: wl(kw), wc(kw), wu(kw)

contains

!==============================================================================

  ! Initialize and return grids
  function get_grids( ) result( grids )

    use musica_assert,                 only : die
    use musica_constants,              only : dk => musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_from_host,           only : grid_from_host_t, grid_updater_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t

    class(grid_warehouse_t), pointer :: grids

    character(len=80) :: pathname
    class(grid_t), pointer :: grid
    type(grid_updater_t) :: updater
    real, allocatable :: wl_edges(:), wl_mids(:)

    ! used in set gridz
    integer :: izout
    real :: z1, zout, zaird, ztemp

    grids => grid_warehouse_t( )

    ! initialize the altitude grid
    z1 = 0.0 ! surface elevation (km)
    call gridz(z1,zout,izout,nz,z)
    grid => grid_from_host_t( "height", "km", nz-1 )
    call grids%add( grid )
    updater = grids%get_updater( grid )
    call updater%update( edges = real( z(:nz), kind=dk ) )
    deallocate( grid )

    ! initialize wavelength grid
    pathname = 'tuv_doug/INPUT/GRIDS/waccm2_ref_101.grid'
    call gridw(nw,wl,wc,wu,pathname)
    grid => grid_from_host_t( "wavelength", "nm", nw )
    call grids%add( grid )
    allocate( wl_edges(nw+1) )
    wl_edges(:nw) = wl(:nw)
    wl_edges(nw+1) = wl(nw) + ( wl(nw) - wl(nw-1) )
    allocate( wl_mids(nw) )
    wl_mids(:nw-1) = wc(:nw-1)
    wl_mids(nw) = wl_edges(nw) + ( wl_edges(nw+1) - wl_edges(nw) ) * 0.5_dk
    updater = grids%get_updater( grid )
    call updater%update( edges = real( wl_edges(:nw+1), kind=dk ),            &
                         mid_points = real( wl_mids(:nw), kind=dk ) )
    deallocate( grid )

  end function get_grids

!==============================================================================

  ! Calculate the product of the cross section and quantum yield
  subroutine calculate( label, temperature, air_density, xsqy )

    use musica_assert,                 only : die

    ! Label for the photolysis reaction
    character(len=*),  intent(in) :: label
    ! Product of the cross section and quantum yield (height, wavelength) [cm2]
    real, allocatable, intent(out) :: xsqy(:,:)
    ! Temperature at each vertical level [K]
    real,              intent(in)  :: temperature(:)
    ! Air density at each vertical level [molecule cm-3]
    real,              intent(in)  :: air_density(:)

    real :: l_xsqy(kj,kz,kw)
    character(len=60) :: all_labels(kj)
    character(len=80) :: pn = "tuv_doug/INPUT/XSQY/"
    integer :: j

    j = 0

    allocate( xsqy(nz,nw) )
    select case( label )
      case( "H2O + hv -> H + OH" )
        call XSQY_H2O(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "H2O + hv -> H2 + O(1D)" )
        call XSQY_H2O(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(2,:nz,:nw)
      case( "H2O + hv -> 2H + O(3P)" )
        call XSQY_H2O(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(3,:nz,:nw)
      case( "N2O5 + hv -> NO3 + NO2" )
        call XSQY_N2O5(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "N2O5 + hv -> NO3 + NO + O" )
        call XSQY_N2O5(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(2,:nz,:nw)
      case( "CH2Br2 + hv -> 2Br" )
        call XSQY_CH2BR2(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "BRO + hv -> Br + O" )
        call XSQY_BRO(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "Cl2O2 + hv -> Cl + ClOO" )
        call XSQY_CL2O2(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "ClO + hv -> Cl + O" )
        call XSQY_CLO(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case ( "HNO3 + hv -> OH + NO2" )
        call XSQY_HNO3(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case ( "CF2Cl2 + hv -> 2Cl" )
        call XSQY_CF2CL2(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case ( "CFC113 + hv -> 3Cl" )
        call XSQY_CFC113(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case ( "CFC114 + hv -> 2Cl" )
        call XSQY_CFC114(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "CFC115 + hv -> Cl" )
        call XSQY_CFC115(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case ( "CFCl3 + hv -> 3Cl" )
        call XSQY_CFCL3(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "CH3Br + hv -> Br" )
        call XSQY_CH3BR(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case ( "CHBr3 + hv -> 3Br" )
        call XSQY_CHBR3(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "H1301 + hv -> Br" )
        call XSQY_H1301(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "H2402 + hv -> 2Br" )
        call XSQY_H2402(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "HCFC22 + hv -> Cl" )
        call XSQY_HCFC22(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "HCFC141b + hv -> 2Cl" )
        call XSQY_HCFC141b(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "HCFC142b + hv -> Cl" )
        call XSQY_HCFC142b(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "BrONO2 + hv -> Br + NO3" )
        call XSQY_BRONO2(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "BrONO2 + hv -> BrO + NO2" )
        call XSQY_BRONO2(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(2,:nz,:nw)
      case( "HO2NO2 + hv -> OH + NO3" )
        call XSQY_HO2NO2(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "HO2NO2 + hv -> HO2 + NO2" )
        call XSQY_HO2NO2(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(2,:nz,:nw)
      case( "CH3Cl + hv -> Cl" )
        call XSQY_CH3CL(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "SO2 + hv -> SO + O" )
        call XSQY_SO2(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case( "CH3COCH3 + hv -> CH3CO3 + CH3O2" )
        call xsqy_acetone(nw,wl,wc,nz,temperature,air_density,j,l_xsqy,all_labels,pn)
        xsqy(:,:) = l_xsqy(1,:nz,:nw)
      case default
        call die( 946669022 )
    end select

  end subroutine calculate

!==============================================================================

end module tuv_doug
