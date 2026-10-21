module reference_profile_config
  use iso_c_binding
  use musica_config
  use musica_string
  use tuvx_grid_warehouse
  use tuvx_profile
  use tuvx_profile_factory
  use reference_profiles, only: copy_profile
  implicit none
contains
  subroutine ref_profile_config(ng,grid_text,np,profile_text,n,fields,out) bind(c)
    integer(c_int), value :: ng,np,n
    character(c_char), intent(in) :: grid_text(ng),profile_text(np)
    integer(c_int), intent(out) :: fields(6)
    real(c_double), intent(out) :: out(n+1,6)
    character(:), allocatable :: input
    type(config_t) :: config
    type(string_t) :: kind
    class(grid_warehouse_t), pointer :: grids
    class(profile_t), pointer :: profile
    integer :: i
    allocate(character(ng) :: input)
    do i=1,ng
      input(i:i)=grid_text(i)
    end do
    config=input
    grids=>grid_warehouse_t(config)
    deallocate(input)
    allocate(character(np) :: input)
    do i=1,np
      input(i:i)=profile_text(i)
    end do
    config=input
    call config%get('type',kind,'reference')
    if(kind=='extraterrestrial flux') call config%add('enable diagnostics',.false.,'reference')
    profile=>profile_builder(config,grids)
    call copy_profile(profile,n,out)
    fields=0
    if(allocated(profile%edge_val_)) fields(1)=1
    if(allocated(profile%mid_val_)) fields(2)=1
    if(allocated(profile%delta_val_)) fields(3)=1
    if(allocated(profile%layer_dens_)) fields(4)=1
    if(allocated(profile%exo_layer_dens_)) fields(5)=1
    if(allocated(profile%burden_dens_)) fields(6)=1
    deallocate(profile,grids)
  end subroutine
end module
