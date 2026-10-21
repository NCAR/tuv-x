module reference_optics
  use iso_c_binding
  use musica_config
  use tuvx_grid
  use tuvx_grid_from_host
  use tuvx_grid_warehouse
  use tuvx_profile
  use tuvx_profile_from_host
  use tuvx_profile_warehouse
  use tuvx_radiator
  implicit none
contains
  subroutine setup(n,nw,z,w,grids,profiles)
    integer, intent(in) :: n,nw
    real(c_double), intent(in) :: z(n+1),w(nw+1)
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(grid_t), pointer :: grid
    grids => grid_warehouse_t()
    profiles => profile_warehouse_t()
    grid => grid_from_host_t('height','km',n)
    grid%edge_=z
    grid%mid_=.5_c_double*(z(:n)+z(2:))
    grid%delta_=z(2:)-z(:n)
    call grids%add(grid)
    deallocate(grid)
    grid => grid_from_host_t('wavelength','nm',nw)
    grid%edge_=w
    grid%mid_=.5_c_double*(w(:nw)+w(2:))
    grid%delta_=w(2:)-w(:nw)
    call grids%add(grid)
    deallocate(grid)
  end subroutine

  subroutine ref_aerosol(n,nw,z,w,nd,depth,ssa,g,tau550,output) bind(c)
    use tuvx_radiator_aerosol
    integer(c_int), value :: n,nw,nd
    real(c_double), intent(in) :: z(n+1),w(nw+1),depth(nd)
    real(c_double), value :: ssa,g,tau550
    real(c_double), intent(out) :: output(n,nw,3)
    type(config_t) :: config
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(radiator_t), pointer :: radiator
    call setup(n,nw,z,w,grids,profiles)
    call config%add('name','aerosol','reference')
    call config%add('type','aerosol','reference')
    call config%add('optical depths',depth,'reference')
    call config%add('single scattering albedo',ssa,'reference')
    call config%add('asymmetry factor',g,'reference')
    call config%add('550 nm optical depth',tau550,'reference')
    radiator => radiator_aerosol_t(config,grids,profiles)
    output(:,:,1)=radiator%state_%layer_OD_
    output(:,:,2)=radiator%state_%layer_SSA_
    output(:,:,3)=radiator%state_%layer_G_(:,:,1)
    deallocate(radiator,grids,profiles)
  end subroutine

  subroutine ref_gas(n,nw,z,w,columns,cross,is_air,output) bind(c)
    use tuvx_cross_section
    use tuvx_cross_section_warehouse
    integer(c_int), value :: n,nw,is_air
    real(c_double), intent(in) :: z(n+1),w(nw+1),columns(n),cross(nw)
    real(c_double), intent(out) :: output(n,nw,3)
    type(config_t) :: config
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(profile_t), pointer :: p
    type(cross_section_warehouse_t) :: warehouse
    class(cross_section_t), pointer :: spectrum
    class(radiator_t), pointer :: radiator
    call setup(n,nw,z,w,grids,profiles)
    p => profile_from_host_t('profile','molecule cm-3',n)
    p%layer_dens_=columns
    call profiles%add(p)
    deallocate(p)
    allocate(spectrum)
    spectrum%height_grid_=grids%get_ptr('height','km')
    spectrum%wavelength_grid_=grids%get_ptr('wavelength','nm')
    allocate(spectrum%cross_section_parms(1))
    allocate(spectrum%cross_section_parms(1)%array(nw,1))
    spectrum%cross_section_parms(1)%array(:,1)=cross
    allocate(warehouse%handles_(1),warehouse%cross_sections_(1))
    warehouse%handles_(1)='cross'
    warehouse%cross_sections_(1)%val_=>spectrum
    call config%add('name','gas','reference')
    call config%add('type','base','reference')
    call config%add('vertical profile','profile','reference')
    call config%add('vertical profile units','molecule cm-3','reference')
    call config%add('cross section','cross','reference')
    call config%add('treat as air',is_air/=0,'reference')
    radiator => radiator_t(config,grids,profiles,warehouse)
    call radiator%update_state(grids,profiles,warehouse)
    output(:,:,1)=radiator%state_%layer_OD_
    output(:,:,2)=radiator%state_%layer_SSA_
    output(:,:,3)=radiator%state_%layer_G_(:,:,1)
    deallocate(radiator,grids,profiles)
  end subroutine

  subroutine ref_radiator_dataset(output) bind(c)
    use tuvx_radiator_from_netcdf_file
    real(c_double), intent(out) :: output(3,4,3)
    type(config_t) :: config
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(radiator_t), pointer :: radiator
    call setup(3,4,[0._c_double,10._c_double,20._c_double,30._c_double], &
        [550._c_double,600._c_double,650._c_double,700._c_double,720._c_double],grids,profiles)
    call config%add('name','foo','reference')
    call config%add('type','from netcdf file','reference')
    call config%add('file path',TT_REFERENCE_DATA//'/../test/data/radiator.nc','reference')
    radiator => radiator_from_netcdf_file_t(config,grids,profiles)
    output(:,:,1)=radiator%state_%layer_OD_
    output(:,:,2)=radiator%state_%layer_SSA_
    output(:,:,3)=radiator%state_%layer_G_(:,:,1)
    deallocate(radiator,grids,profiles)
  end subroutine
end module
