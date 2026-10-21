! Original constructors, including their data parsing and integration paths.
module reference_profiles
  use iso_c_binding
  use musica_config
  use musica_string
  use tuvx_grid
  use tuvx_grid_from_host
  use tuvx_grid_warehouse
  use tuvx_profile
  implicit none
contains
  subroutine ref_grid(kind, n, edges, first, last, delta, path, length, output) bind(c)
    use tuvx_grid_equal_delta
    use tuvx_grid_from_config
    use tuvx_grid_from_csv_file
    integer(c_int), value :: kind,n,length
    real(c_double), intent(in) :: edges(n+1)
    real(c_double), value :: first,last,delta
    character(c_char), intent(in) :: path(length)
    real(c_double), intent(out) :: output(n+1,3)
    character(len=length) :: filename
    type(config_t) :: config
    class(grid_t), pointer :: grid
    integer :: i
    call config%add('name','height','reference')
    call config%add('units','km','reference')
    select case(kind)
    case(0)
      call config%add('type','from config file','reference')
      call config%add('values',edges,'reference')
      grid => grid_from_config_t(config)
    case(1)
      call config%add('type','equal interval','reference')
      call config%add('begins at',first,'reference')
      call config%add('ends at',last,'reference')
      call config%add('cell delta',delta,'reference')
      grid => grid_equal_delta_t(config)
    case(2)
      do i=1,length
        filename(i:i)=path(i)
      enddo
      call config%add('type','from csv file','reference')
      call config%add('file path',TT_REFERENCE_DATA//'/'//filename,'reference')
      grid => grid_from_csv_file_t(config)
    end select
    if (grid%ncells_ /= n) error stop 'grid size differs'
    output=0
    output(:,1)=grid%edge_
    output(:n,2)=grid%mid_
    output(:n,3)=grid%delta_
    deallocate(grid)
  end subroutine

  subroutine setup_grid(n, edges, name, units, grids)
    integer, intent(in) :: n
    real(c_double), intent(in) :: edges(n+1)
    character(*), intent(in) :: name,units
    class(grid_warehouse_t), pointer :: grids
    class(grid_t), pointer :: grid
    grids => grid_warehouse_t()
    grid => grid_from_host_t(name, units, n)
    grid%edge_=edges
    grid%mid_=.5_c_double*(edges(:n)+edges(2:))
    grid%delta_=edges(2:)-edges(:n)
    call grids%add(grid)
    deallocate(grid)
  end subroutine

  subroutine copy_profile(p,n,out)
    class(profile_t), intent(in) :: p
    integer, intent(in) :: n
    real(c_double), intent(out) :: out(n+1,6)
    out=0
    if (allocated(p%edge_val_)) out(:,1)=p%edge_val_
    if (allocated(p%mid_val_)) out(:n,2)=p%mid_val_
    if (allocated(p%delta_val_)) out(:n,3)=p%delta_val_
    if (allocated(p%layer_dens_)) out(:n,4)=p%layer_dens_
    if (allocated(p%exo_layer_dens_)) out(:,5)=p%exo_layer_dens_
    if (allocated(p%burden_dens_)) out(:n,6)=p%burden_dens_
  end subroutine

  subroutine ref_profile(kind,n,edges,values,hscale,column,convert,output) bind(c)
    use tuvx_profile_from_config
    use tuvx_profile_from_csv_file
    use tuvx_profile_air
    use tuvx_profile_o2
    use tuvx_profile_o3
    integer(c_int), value :: kind,n
    real(c_double), intent(in) :: edges(n+1),values(n+1)
    real(c_double), value :: hscale,column,convert
    real(c_double), intent(out) :: output(n+1,6)
    class(grid_warehouse_t), pointer :: grids
    class(profile_t), pointer :: profile
    type(config_t) :: config,grid
    type(string_t) :: grid_name,grid_units
    grid_name='height'
    grid_units='km'
    if (convert == 1) grid_units='none'
    call setup_grid(n,edges,grid_name%to_char(),grid_units%to_char(),grids)
    call config%add('name','profile','reference')
    call config%add('units','molecule cm-3','reference')
    if (kind /= 1) call config%add('scale height',hscale,'reference')
    select case(kind)
    case(0,4)
      call grid%add('name',grid_name,'reference')
      call grid%add('units',grid_units,'reference')
      call config%add('grid',grid,'reference')
      if (kind == 0) then
        call config%add('type','from config file','reference')
        call config%add('values',values,'reference')
        profile => profile_from_config_t(config,grids)
      else
        call config%add('type','from csv file','reference')
        call config%add('file path',TT_REFERENCE_DATA//'/profiles/atmosphere/ussa.temp','reference')
        profile => profile_from_csv_file_t(config,grids)
      endif
    case(1,2)
      call config%add('file path',TT_REFERENCE_DATA//'/profiles/atmosphere/ussa.dens','reference')
      if(kind==1) then
        call config%add('type','air','reference')
        profile => profile_air_t(config,grids)
      else
        call config%add('type','O2','reference')
        profile => profile_o2_t(config,grids)
      endif
    case(3)
      call config%add('type','O3','reference')
      call config%add('file path',TT_REFERENCE_DATA//'/profiles/atmosphere/ussa.ozone','reference')
      call config%add('reference column',column,'reference')
      profile => profile_o3_t(config,grids)
    end select
    call copy_profile(profile,n,output)
    deallocate(profile,grids)
  end subroutine

  subroutine ref_solar_profiles(n,edges,year,month,day,zone,lat,lon,zenith,distance) bind(c)
    use tuvx_profile_solar_zenith_angle
    use tuvx_profile_earth_sun_distance
    integer(c_int), value :: n,year,month,day
    real(c_double), intent(in) :: edges(n+1)
    real(c_double), value :: zone,lat,lon
    real(c_double), intent(out) :: zenith(n+1,6),distance(n+1,6)
    type(config_t) :: c
    class(grid_warehouse_t), pointer :: grids
    class(profile_t), pointer :: p
    call setup_grid(n,edges,'time','hours',grids)
    call c%add('name','Earth-Sun distance','reference')
    call c%add('type','Earth-Sun distance','reference')
    call c%add('units','AU','reference')
    call c%add('year',year,'reference')
    call c%add('month',month,'reference')
    call c%add('day',day,'reference')
    call c%add('time zone',zone,'reference')
    p => profile_earth_sun_distance_t(c,grids)
    call copy_profile(p,n,distance)
    deallocate(p)
    call c%add('name','solar zenith angle','reference')
    call c%add('type','solar zenith angle','reference')
    call c%add('units','degrees','reference')
    call c%add('latitude',lat,'reference')
    call c%add('longitude',lon,'reference')
    p => profile_solar_zenith_angle_t(c,grids)
    call copy_profile(p,n,zenith)
    deallocate(p,grids)
  end subroutine

  subroutine ref_solar_flux(n,edges,nfiles,methods,output) bind(c)
    use tuvx_profile_extraterrestrial_flux
    integer(c_int), value :: n,nfiles
    integer(c_int), intent(in) :: methods(nfiles)
    real(c_double), intent(in) :: edges(n+1)
    real(c_double), intent(out) :: output(n)
    type(string_t) :: files(4),interpolators(nfiles)
    type(config_t) :: config
    class(grid_warehouse_t), pointer :: grids
    class(profile_t), pointer :: profile
    integer :: i
    files(1)=TT_REFERENCE_DATA//'/profiles/solar/susim_hi.flx'
    files(2)=TT_REFERENCE_DATA//'/profiles/solar/atlas3_1994_317_a.dat'
    files(3)=TT_REFERENCE_DATA//'/profiles/solar/sao2010.solref.converted'
    files(4)=TT_REFERENCE_DATA//'/profiles/solar/neckel.flx'
    do i=1,nfiles
      select case(methods(i))
      case(0)
        interpolators(i)='linear'
      case(1)
        interpolators(i)='conserving'
      case(2)
        interpolators(i)='fractional source'
      case(3)
        interpolators(i)='fractional target'
      end select
    enddo
    call setup_grid(n,edges,'wavelength','nm',grids)
    call config%add('name','extraterrestrial flux','reference')
    call config%add('type','extraterrestrial flux','reference')
    call config%add('units','photon cm-2 s-1','reference')
    call config%add('file path',files(:nfiles),'reference')
    call config%add('interpolator',interpolators,'reference')
    profile => profile_extraterrestrial_flux_t(config,grids)
    output=profile%mid_val_
    deallocate(profile,grids)
  end subroutine
  subroutine ref_host_profile(n,mode,edges,mids,columns,hscale,exo,output) bind(c)
    use tuvx_profile_from_host
    integer(c_int), value :: n,mode
    real(c_double), intent(in) :: edges(n+1),mids(n),columns(n)
    real(c_double), value :: hscale,exo
    real(c_double), intent(out) :: output(n+1,6)
    type(profile_from_host_t), pointer :: p
    type(profile_updater_t) :: updater
    p => profile_from_host_t('host','units',n)
    updater=profile_updater_t(p)
    select case(mode)
    case(0)
      call updater%update(edge_values=edges)
    case(1)
      call updater%update(edge_values=edges,mid_point_values=mids,layer_densities=columns)
    case(2)
      call updater%update(edge_values=edges,layer_densities=columns,scale_height=hscale)
    case(3)
      call updater%update(edge_values=edges,layer_densities=columns,exo_density=exo)
    case(4)
      call updater%update(edge_values=edges,layer_densities=columns,scale_height=hscale)
      call updater%update(edge_values=edges,mid_point_values=mids)
    end select
    call copy_profile(p,n,output)
    deallocate(p)
  end subroutine

  subroutine ref_host_grid(n,has_mids,edges,mids,output) bind(c)
    integer(c_int), value :: n,has_mids
    real(c_double), intent(in) :: edges(n+1),mids(n)
    real(c_double), intent(out) :: output(n+1,3)
    type(grid_from_host_t), pointer :: grid
    type(grid_updater_t) :: updater
    grid => grid_from_host_t('host','units',n)
    updater=grid_updater_t(grid)
    if(has_mids==0) then
      call updater%update(edges=edges)
    else
      call updater%update(edges=edges,mid_points=mids)
    endif
    output=0
    output(:,1)=grid%edge_
    output(:n,2)=grid%mid_
    output(:n,3)=grid%delta_
    deallocate(grid)
  end subroutine
end module
