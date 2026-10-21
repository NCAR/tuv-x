module reference_model
  use iso_c_binding
  use musica_config
  use musica_string
  use tuvx_core
  use tuvx_grid
  use tuvx_solver
  use tuvx_solver_from_host
  implicit none
  class(core_t), pointer :: model => null()
contains
  subroutine ref_model_create(length,text,shape) bind(c)
    integer(c_int), value :: length
    character(c_char), intent(in) :: text(length)
    integer(c_int), intent(out) :: shape(5)
    character(:), allocatable :: input
    type(config_t) :: config
    class(grid_t), pointer :: grid
    integer :: i
    if(associated(model)) deallocate(model)
    allocate(character(length) :: input)
    do i=1,length
      input(i:i)=text(i)
    end do
    config=input
    model=>core_t(config)
    grid=>model%get_grid('height','km')
    shape(1)=grid%size()+1
    deallocate(grid)
    grid=>model%get_grid('wavelength','nm')
    shape(2)=grid%size()
    deallocate(grid)
    shape(3)=model%number_of_photolysis_reactions()
    shape(4)=model%number_of_heating_rates()
    shape(5)=model%number_of_dose_rates()
  end subroutine
  subroutine ref_model_run(nz,nw,np,nh,nd,zenith,scale,mode,host,radiation,photo,heat,dose) bind(c)
    integer(c_int), value :: nz,nw,np,nh,nd,mode
    real(c_double), value :: zenith,scale
    real(c_double), intent(in) :: host(nz,nw,6)
    real(c_double), intent(out) :: radiation(nz,nw,6),photo(nz,np),heat(nz,nh),dose(nz,nd)
    type(radiation_field_updater_t) :: updater
    type(radiation_field_t) :: field
    logical :: found
    if(mode>0) then
      updater=model%get_radiation_field_updater(found)
      if(.not.found) error stop 'missing host updater'
      if(mode==1) then
        call updater%update(host(:,:,1),host(:,:,3),host(:,:,2))
      else
        call updater%update(host(:,:,1),host(:,:,3),host(:,:,2),host(:,:,4),host(:,:,6),host(:,:,5))
      end if
    end if
    call model%run(zenith,scale,photo,dose,heat)
    field=model%get_radiation_field()
    radiation(:,:,1)=field%fdr_
    radiation(:,:,2)=field%fdn_
    radiation(:,:,3)=field%fup_
    radiation(:,:,4)=field%edr_
    radiation(:,:,5)=field%edn_
    radiation(:,:,6)=field%eup_
  end subroutine
  subroutine ref_model_label(kind,index,length,text) bind(c)
    integer(c_int), value :: kind,index,length
    character(c_char), intent(out) :: text(length)
    type(string_t), allocatable :: labels(:)
    character(256) :: label
    integer :: i
    select case(kind)
    case(0)
      labels=model%photolysis_reaction_labels()
    case(1)
      labels=model%heating_rate_labels()
    case(2)
      labels=model%dose_rate_labels()
    end select
    label=labels(index+1)
    text=c_null_char
    do i=1,min(len_trim(label),length-1)
      text(i)=label(i:i)
    end do
  end subroutine
  subroutine ref_model_destroy() bind(c)
    if(associated(model)) deallocate(model)
    nullify(model)
  end subroutine
end module
