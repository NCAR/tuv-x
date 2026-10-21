! Exercise original rate constructors and integration routines.
module reference_rates
  use iso_c_binding
  implicit none
contains
  subroutine ref_rates(kind,nw,nz,edges,mids,temp,air,solar,field,weights,threshold,scale,out) bind(c)
    use reference_yields, only: setup
    use musica_config
    use tuvx_grid_warehouse
    use tuvx_profile_warehouse
    use tuvx_profile
    use tuvx_profile_from_host
    use tuvx_photolysis_rates
    use tuvx_heating_rates
    use tuvx_dose_rates
    use tuvx_spectral_weight
    use tuvx_solver
    use tuvx_spherical_geometry
    use tuvx_la_sr_bands
    integer(c_int), value :: kind,nw,nz
    real(c_double), intent(in) :: edges(nw+1),mids(nw),temp(nz),air(nz),solar(nw)
    real(c_double), intent(in) :: field(nz,nw,6),weights(nw)
    real(c_double), value :: threshold,scale
    real(c_double), intent(out) :: out(nz,3)
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    class(profile_t), pointer :: profile
    class(photolysis_rates_t), pointer :: photo
    type(heating_rates_t), pointer :: heat
    type(dose_rates_t) :: dose
    type(radiation_field_t) :: radiation
    type(spherical_geometry_t) :: geometry
    type(la_sr_bands_t) :: bands
    type(config_t) :: config,reactions,reaction,cross,yield,heating
    call setup(nw,nz,edges,mids,temp,air,grids,profiles)
    profile => profile_from_host_t('extraterrestrial flux','photon cm-2 s-1',nw)
    profile%mid_val_=solar
    call profiles%add(profile)
    deallocate(profile)
    call cross%add('type','air','reference')
    select case(kind)
    case(0)
      call yield%add('type','NO3-_(aq)+hv->NO2(aq)+O-','reference')
    case(1)
      call yield%add('type','HO2','reference')
    case(2)
      call yield%add('type','O3+hv->O2+O(1D)','reference')
    end select
    call heating%add('energy term',threshold,'reference')
    call reaction%add('name','test reaction','reference')
    call reaction%add('cross section',cross,'reference')
    call reaction%add('quantum yield',yield,'reference')
    call reaction%add('heating',heating,'reference')
    call reaction%add('scaling factor',scale,'reference')
    call reactions%empty()
    call config%empty()
    call reactions%add('test reaction',reaction,'reference')
    call config%add('reactions',reactions,'reference')
    photo => photolysis_rates_t(config,grids,profiles)
    heat => heating_rates_t(config,grids,profiles)
    dose%height_grid_=grids%get_ptr('height','km')
    dose%wavelength_grid_=grids%get_ptr('wavelength','nm')
    dose%etfl_profile_=profiles%get_ptr('extraterrestrial flux','photon cm-2 s-1')
    dose%enable_diagnostics_=.false.
    allocate(dose%spectral_weights_(1))
    allocate(spectral_weight_t :: dose%spectral_weights_(1)%val_)
    associate(weight => dose%spectral_weights_(1)%val_)
      allocate(weight%spectral_weight_parms(1))
      allocate(weight%spectral_weight_parms(1)%array(nw,1))
      weight%spectral_weight_parms(1)%array(:,1)=weights
    end associate
    radiation%fdr_=field(:,:,1)
    radiation%fdn_=field(:,:,2)
    radiation%fup_=field(:,:,3)
    radiation%edr_=field(:,:,4)
    radiation%edn_=field(:,:,5)
    radiation%eup_=field(:,:,6)
    call photo%get(bands,geometry,grids,profiles,radiation,out(:,1:1),'')
    call heat%get(bands,geometry,grids,profiles,radiation,out(:,2:2))
    call dose%get(grids,profiles,radiation,out(:,3:3),'')
    deallocate(photo,heat,grids,profiles)
  end subroutine
end module
