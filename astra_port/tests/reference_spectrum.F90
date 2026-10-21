! Full upstream constructors and evaluators; this is test-only code.
module reference_spectrum
  use iso_c_binding
  use musica_config
  use tuvx_grid
  use tuvx_grid_from_host
  use tuvx_grid_warehouse
  use tuvx_profile
  use tuvx_profile_from_host
  use tuvx_profile_warehouse
  implicit none
contains
  subroutine ref_spectrum(role, length, text, nw, nz, midpoint, edges, mids, temp, air, out) bind(c)
    use tuvx_cross_section
    use tuvx_cross_section_factory
    use tuvx_quantum_yield
    use tuvx_quantum_yield_factory
    use tuvx_spectral_weight
    use tuvx_spectral_weight_factory
    integer(c_int), value :: role, length, nw, nz, midpoint
    character(c_char), intent(in) :: text(length)
    real(c_double), intent(in) :: edges(nw+1), mids(nw), temp(nz), air(nz)
    real(c_double), intent(out) :: out(*)
    type(config_t) :: config
    character(:), allocatable :: input
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(grid_t), pointer :: grid
    class(profile_t), pointer :: profile
    class(cross_section_t), pointer :: cs
    class(quantum_yield_t), pointer :: qy
    class(spectral_weight_t), pointer :: weight
    real(c_double), allocatable :: result(:,:), weights(:)
    integer :: i, n
    allocate(character(length) :: input)
    do i=1,length
      input(i:i)=text(i)
    end do
    config = input
    n=nz-1+midpoint
    grids => grid_warehouse_t()
    profiles => profile_warehouse_t()
    grid => grid_from_host_t('height','km',n)
    call grids%add(grid)
    deallocate(grid)
    grid => grid_from_host_t('wavelength','nm',nw)
    grid%edge_=edges
    grid%mid_=mids
    grid%delta_=edges(2:)-edges(:nw)
    call grids%add(grid)
    deallocate(grid)
    profile => profile_from_host_t('temperature','K',n)
    if (midpoint == 0) then
      profile%edge_val_=temp
    else
      profile%mid_val_=temp
      profile%edge_val_=temp(1)
    end if
    call profiles%add(profile)
    deallocate(profile)
    profile => profile_from_host_t('air','molecule cm-3',n)
    if (midpoint == 0) then
      profile%edge_val_=air
    else
      profile%mid_val_=air
      profile%edge_val_=air(1)
    end if
    call profiles%add(profile)
    deallocate(profile)
    select case(role)
    case(0)
      cs => cross_section_builder(config,grids,profiles)
      result=cs%calculate(grids,profiles,midpoint/=0)
      out(1:nw*nz)=reshape(result(:nz,:),[nw*nz])
      deallocate(cs)
    case(1)
      qy => quantum_yield_builder(config,grids,profiles)
      result=qy%calculate(grids,profiles)
      out(1:nw*nz)=reshape(result(:nz,:),[nw*nz])
      deallocate(qy)
    case(2)
      weight => spectral_weight_builder(config,grids,profiles)
      weights=weight%calculate(grids,profiles)
      out(1:nw)=weights
      deallocate(weight)
    end select
    deallocate(grids,profiles)
  end subroutine
end module
