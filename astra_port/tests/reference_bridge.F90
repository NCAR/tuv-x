! Test-only ISO C bindings to the original TUV-x implementation.
! No reference source is linked into the production C library.
module reference_bridge
  use iso_c_binding
  implicit none
contains
  subroutine ref_interpolate(method, fold, nx, x, y, nt, target, out) bind(c)
    use tuvx_interpolate
    integer(c_int), value :: method, fold, nx, nt
    real(c_double), intent(in) :: x(nx), y(nx), target(nt)
    real(c_double), intent(out) :: out(*)
    class(interpolator_t), allocatable :: interp
    real(c_double), allocatable :: result(:)
    select case(method)
    case(0)
      allocate(interpolator_linear_t :: interp)
    case(1)
      allocate(interpolator_conserving_t :: interp)
    case(2)
      allocate(interpolator_fractional_source_t :: interp)
    case(3)
      allocate(interpolator_fractional_target_t :: interp)
    case default
      error stop 'unknown interpolation method'
    end select
    interp%fold_in_ = fold /= 0
    result = interp%interpolate(target, x, y, 'C roundoff comparison')
    out(1:size(result)) = result
  end subroutine

  subroutine ref_tridiagonal(n, a, b, c, rhs, out) bind(c)
    use tuvx_linear_algebra_linpack
    integer(c_int), value :: n
    real(c_double), intent(in) :: a(n), b(n), c(n), rhs(n)
    real(c_double), intent(out) :: out(n)
    type(linear_algebra_linpack_t) :: algebra
    out = algebra%tridiag(a, b, c, rhs)
  end subroutine

  subroutine ref_solar(year, day, hour, lat, lon, distance, elevation) bind(c)
    use tuvx_profile_utils
    integer(c_int), value :: year, day
    real(c_double), value :: hour, lat, lon
    real(c_double), intent(out) :: distance, elevation
    distance = earth_sun_distance(year, day, hour)
    elevation = solar_zenith_angle(year, day, hour, lat, lon)
  end subroutine

  subroutine ref_geometry(n, z, zen, air, nid, paths, vertical, slant) bind(c)
    use tuvx_grid
    use tuvx_grid_from_host
    use tuvx_grid_warehouse
    use tuvx_spherical_geometry
    integer(c_int), value :: n
    real(c_double), value :: zen
    real(c_double), intent(in) :: z(n), air(n)
    integer(c_int), intent(out) :: nid(n)
    real(c_double), intent(out) :: paths(n-1,n), vertical(n), slant(n)
    class(grid_warehouse_t), pointer :: grids
    class(grid_t), pointer :: grid
    type(spherical_geometry_t), pointer :: geometry
    integer :: i
    grids => grid_warehouse_t()
    grid => grid_from_host_t('height', 'km', n-1)
    grid%edge_ = z
    grid%mid_ = (z(:n-1) + z(2:)) * 0.5_c_double
    grid%delta_ = z(2:) - z(:n-1)
    call grids%add(grid)
    deallocate(grid)
    geometry => spherical_geometry_t(grids)
    call geometry%set_parameters(zen, grids)
    nid = geometry%nid_
    paths = transpose(geometry%dsdh_)
    ! The original air_mass leaves vertical(n) unspecified.
    vertical = 0.0_c_double
    call geometry%air_mass(air, vertical, slant)
    deallocate(geometry, grids)
  end subroutine

  subroutine ref_delta(n, nw, nr, z, zen, surface, tau, ssa, asym, is_air, mixed, radiation) bind(c)
    use tuvx_grid
    use tuvx_grid_from_host
    use tuvx_grid_warehouse
    use tuvx_profile
    use tuvx_profile_from_host
    use tuvx_profile_warehouse
    use tuvx_radiator
    use tuvx_radiator_from_host
    use tuvx_radiator_warehouse
    use tuvx_solver
    use tuvx_solver_delta_eddington
    use tuvx_spherical_geometry
    integer(c_int), value :: n, nw, nr
    real(c_double), value :: zen
    real(c_double), intent(in) :: z(n+1), surface(nw), tau(n,nw,nr), ssa(n,nw,nr), asym(n,nw,nr)
    integer(c_int), intent(in) :: is_air(nr)
    real(c_double), intent(out) :: mixed(n,nw,3), radiation(n+1,nw,6)
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(radiator_warehouse_t), pointer :: radiators
    class(grid_t), pointer :: heights, wavelengths
    class(profile_t), pointer :: albedo
    class(radiator_t), pointer :: radiator
    type(radiator_state_t) :: state
    type(solver_delta_eddington_t) :: solver
    type(spherical_geometry_t), pointer :: geometry
    type(radiation_field_t), pointer :: field
    integer :: j
    character(len=32) :: name
    grids => grid_warehouse_t()
    profiles => profile_warehouse_t()
    radiators => radiator_warehouse_t()
    heights => grid_from_host_t('height', 'km', n)
    heights%edge_ = z
    heights%mid_ = (z(:n)+z(2:))*0.5_c_double
    heights%delta_ = z(2:)-z(:n)
    wavelengths => grid_from_host_t('wavelength', 'nm', nw)
    wavelengths%edge_ = [(real(j,c_double), j=1,nw+1)]
    wavelengths%mid_ = [(real(j,c_double)+0.5_c_double, j=1,nw)]
    wavelengths%delta_ = 1.0_c_double
    call grids%add(heights)
    call grids%add(wavelengths)
    albedo => profile_from_host_t('surface albedo', 'none', nw)
    albedo%mid_val_ = surface
    call profiles%add(albedo)
    do j=1,nr
      write(name,'(I0)') j
      radiator => radiator_from_host_t(trim(name), heights, wavelengths)
      radiator%state_%layer_OD_ = tau(:,:,j)
      radiator%state_%layer_SSA_ = ssa(:,:,j)
      radiator%state_%layer_G_(:,:,1) = asym(:,:,j)
      radiator%is_air_ = is_air(j) /= 0
      call radiators%add(radiator)
      deallocate(radiator)
    enddo
    allocate(state%layer_G_(n,nw,1))
    call radiators%accumulate_states(state)
    mixed(:,:,1) = state%layer_OD_
    mixed(:,:,2) = state%layer_SSA_
    mixed(:,:,3) = state%layer_G_(:,:,1)
    solver%height_grid_ = grids%get_ptr('height','km')
    solver%wavelength_grid_ = grids%get_ptr('wavelength','nm')
    solver%surface_albedo_profile_ = profiles%get_ptr('surface albedo','none')
    geometry => spherical_geometry_t(grids)
    call geometry%set_parameters(zen, grids)
    field => solver%update_radiation_field(zen,n,geometry,grids,profiles,radiators)
    radiation(:,:,1) = field%fdr_
    radiation(:,:,2) = field%fdn_
    radiation(:,:,3) = field%fup_
    radiation(:,:,4) = field%edr_
    radiation(:,:,5) = field%edn_
    radiation(:,:,6) = field%eup_
    deallocate(field,geometry,heights,wavelengths,albedo,grids,profiles,radiators)
  end subroutine
end module
