! Test-only binding to the original complete discrete-ordinate solver.
module reference_discrete
  use iso_c_binding
  implicit none
contains
  subroutine ref_discrete(n, nw, nr, ns, z, zen, surface, tau, ssa, asym, is_air, mixed, radiation) bind(c)
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
    use tuvx_solver_discrete_ordinate
    use tuvx_discrete_ordinate_util, only : solver_constants_t
    use tuvx_spherical_geometry
    integer(c_int), value :: n, nw, nr, ns
    real(c_double), value :: zen
    real(c_double), intent(in) :: z(n+1), surface(nw), tau(n,nw,nr), ssa(n,nw,nr), asym(n,nw,nr)
    integer(c_int), intent(in) :: is_air(nr)
    real(c_double), intent(out) :: mixed(n,nw,ns+2), radiation(n+1,nw,6)
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(radiator_warehouse_t), pointer :: radiators
    class(grid_t), pointer :: heights, wavelengths
    class(profile_t), pointer :: albedo
    class(radiator_t), pointer :: radiator
    type(radiator_state_t) :: state
    type(solver_discrete_ordinate_t) :: solver
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
    allocate(state%layer_G_(n,nw,ns))
    call radiators%accumulate_states(state)
    mixed(:,:,1) = state%layer_OD_
    mixed(:,:,2) = state%layer_SSA_
    mixed(:,:,3:ns+2) = state%layer_G_
    solver%n_streams_ = ns
    solver%solver_constants_ = solver_constants_t()
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
