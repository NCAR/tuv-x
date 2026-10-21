! Original O2 constructors and calculations, including original coefficient parsing.
module reference_oxygen
  use iso_c_binding
  implicit none
contains
  subroutine ref_oxygen(nw,nz,edges,temp,vertical,slant,scale,depth,cross,coefficients) bind(c)
    use reference_yields, only: setup
    use musica_config
    use tuvx_la_sr_bands
    use tuvx_grid_warehouse
    use tuvx_profile_warehouse
    use tuvx_spherical_geometry
    integer(c_int), value :: nw,nz
    real(c_double), intent(in) :: edges(nw+1),temp(nz),vertical(nz-1),slant(nz)
    real(c_double), value :: scale
    real(c_double), intent(inout) :: depth(nz-1,nw),cross(nz,nw)
    real(c_double), intent(out) :: coefficients(20,17,2)
    real(c_double) :: mids(nw)
    type(config_t) :: config,estimate
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    type(spherical_geometry_t) :: geometry
    type(la_sr_bands_t), pointer :: bands
    mids=.5_c_double*(edges(:nw)+edges(2:))
    call setup(nw,nz,edges,mids,temp,temp,grids,profiles)
    call config%add('cross section parameters file', &
                   TT_REFERENCE_DATA//'/cross_sections/O2_parameters.txt','reference')
    call estimate%add('scale factor',scale,'reference')
    call config%add('O2 estimate',estimate,'reference')
    bands => la_sr_bands_t(config,grids,profiles)
    call bands%optical_depth(grids,profiles,vertical,slant,depth,geometry)
    call bands%cross_section(grids,profiles,vertical,slant,cross,geometry)
    coefficients(:,:,1)=bands%AC
    coefficients(:,:,2)=bands%BC
    deallocate(bands,grids,profiles)
  end subroutine
end module
