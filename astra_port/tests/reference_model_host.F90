module reference_model_host
  use iso_c_binding
  use musica_config
  use tuvx_core
  use tuvx_grid_from_host, only: grid_from_host_t, grid_updater_t
  use tuvx_profile_from_host, only: profile_from_host_t, profile_updater_t
  use tuvx_radiator_from_host, only: radiator_from_host_t, radiator_updater_t
  use tuvx_grid_warehouse
  use tuvx_profile_warehouse
  use tuvx_radiator_warehouse
  use reference_model, only: model
  implicit none
  type(grid_updater_t) :: height_update
  type(profile_updater_t) :: air_update,temperature_update
  type(radiator_updater_t) :: cloud_update
contains
  subroutine ref_model_host_create(length,text,nz,nw,z,w,air,temperature,columns,hscale) bind(c)
    integer(c_int), value :: length,nz,nw
    character(c_char), intent(in) :: text(length)
    real(c_double), intent(in) :: z(nz),w(nw+1),air(nz),temperature(nz),columns(nz-1)
    real(c_double), value :: hscale
    character(:), allocatable :: input
    type(config_t) :: config
    type(grid_from_host_t), pointer :: height,wavelength
    type(profile_from_host_t), pointer :: air_profile,temperature_profile
    type(radiator_from_host_t), pointer :: cloud
    type(grid_updater_t) :: wave_update
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(radiator_warehouse_t), pointer :: radiators
    integer :: i
    if(associated(model)) deallocate(model)
    allocate(character(length) :: input)
    do i=1,length
      input(i:i)=text(i)
    end do
    config=input
    height=>grid_from_host_t('height','km',nz-1)
    height_update=grid_updater_t(height)
    call height_update%update(edges=z)
    wavelength=>grid_from_host_t('wavelength','nm',nw)
    wave_update=grid_updater_t(wavelength)
    call wave_update%update(edges=w)
    grids=>grid_warehouse_t()
    call grids%add(height)
    call grids%add(wavelength)
    air_profile=>profile_from_host_t('air','molecule cm-3',nz-1)
    air_update=profile_updater_t(air_profile)
    call air_update%update(edge_values=air,layer_densities=columns,scale_height=hscale)
    temperature_profile=>profile_from_host_t('temperature','K',nz-1)
    temperature_update=profile_updater_t(temperature_profile)
    call temperature_update%update(edge_values=temperature)
    profiles=>profile_warehouse_t()
    call profiles%add(air_profile)
    call profiles%add(temperature_profile)
    cloud=>radiator_from_host_t('host cloud',height,wavelength)
    radiators=>radiator_warehouse_t()
    call radiators%add(cloud)
    model=>core_t(config,grids,profiles,radiators)
    height_update=model%get_updater(height)
    air_update=model%get_updater(air_profile)
    temperature_update=model%get_updater(temperature_profile)
    cloud_update=model%get_updater(cloud)
    deallocate(height,wavelength,air_profile,temperature_profile,cloud,grids,profiles,radiators)
  end subroutine
  subroutine ref_model_host_update(nz,nw,mode,z,mid,air,am,temperature,tm,columns,top,tau,ssa,g) bind(c)
    integer(c_int), value :: nz,nw,mode
    real(c_double), value :: top
    real(c_double), intent(in) :: z(nz),mid(nz-1),air(nz),am(nz-1),temperature(nz),tm(nz-1),columns(nz-1)
    real(c_double), intent(in) :: tau(nz-1,nw),ssa(nz-1,nw),g(nz-1,nw)
    if(mode==0) then
      call height_update%update(edges=z)
      call air_update%update(edge_values=air,layer_densities=columns,scale_height=top)
      call temperature_update%update(edge_values=temperature)
      call cloud_update%update(tau)
    else if(mode==1) then
      call height_update%update(edges=z,mid_points=mid)
      call air_update%update(edge_values=air,mid_point_values=am,layer_densities=columns,exo_density=top)
      call temperature_update%update(edge_values=temperature,mid_point_values=tm)
      call cloud_update%update(tau,ssa,g)
    else
      call height_update%update(edges=z)
      call air_update%update(edge_values=air)
      call temperature_update%update(edge_values=temperature)
      call cloud_update%update(tau,asymmetry_factors=g)
    end if
  end subroutine
end module
