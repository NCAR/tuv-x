! Test bridge to the original, unmodified temperature parameterizations.
module reference_temperature
  use iso_c_binding
  use tuvx_temperature_parameterization
  use tuvx_temperature_parameterization_taylor_series
  use tuvx_temperature_parameterization_burkholder
  use tuvx_temperature_parameterization_harwood
  implicit none
  type, bind(c) :: temperature_range
    real(c_double) :: minimum, maximum, fixed_value
    integer(c_int) :: is_fixed
  end type
  type, bind(c) :: temperature_model
    integer(c_int) :: kind
    integer(c_size_t) :: first,count,terms
    type(c_ptr) :: a,b,powers,sigma
    real(c_double) :: base_temperature,base_wavelength
    integer(c_int) :: base_10,invert_temperature
    real(c_double) :: burkholder_a,burkholder_b
    integer(c_size_t) :: range_count
    type(c_ptr) :: ranges
  end type
contains
  subroutine make_parameter(p,n,wavelengths,param)
    type(temperature_model), intent(in) :: p
    integer, intent(in) :: n
    real(c_double), intent(in) :: wavelengths(n)
    class(temperature_parameterization_t), pointer :: param
    type(temperature_range), pointer :: ranges(:)
    real(c_double), pointer :: a(:),b(:),powers(:),sigma(:),coeffs(:,:)
    integer :: i,nc,nt,nr
    nc=int(p%count)
    nt=int(p%terms)
    select case(p%kind)
    case(0)
      allocate(temperature_parameterization_t :: param)
    case(1)
      allocate(temperature_parameterization_taylor_series_t :: param)
    case(2)
      allocate(temperature_parameterization_burkholder_t :: param)
    case(3)
      allocate(temperature_parameterization_harwood_t :: param)
    end select
    param%min_wavelength_index_=int(p%first)+1
    param%max_wavelength_index_=int(p%first+p%count)
    param%wavelengths_=wavelengths(param%min_wavelength_index_:param%max_wavelength_index_)
    param%base_temperature_=p%base_temperature
    param%base_wavelength_=p%base_wavelength
    param%is_base_10_=p%base_10/=0
    param%is_temperature_inverted_=p%invert_temperature/=0
    select type(param)
    type is(temperature_parameterization_t)
      call c_f_pointer(p%a,a,[nt])
      call c_f_pointer(p%b,b,[nt])
      call c_f_pointer(p%powers,powers,[nt])
      param%AA_=a
      param%BB_=b
      param%lp_=powers
    type is(temperature_parameterization_taylor_series_t)
      call c_f_pointer(p%a,coeffs,[nc,nt])
      call c_f_pointer(p%sigma,sigma,[nc])
      param%A_=transpose(coeffs)
      param%sigma_=sigma
    type is(temperature_parameterization_burkholder_t)
      call c_f_pointer(p%a,a,[nc])
      call c_f_pointer(p%b,b,[nc])
      param%AA_=a
      param%BB_=b
      param%A_=p%burkholder_a
      param%B_=p%burkholder_b
    type is(temperature_parameterization_harwood_t)
      call c_f_pointer(p%a,a,[nc])
      call c_f_pointer(p%b,b,[nc])
      param%AA_=a
      param%BB_=b
    end select
    nr=int(p%range_count)
    allocate(param%ranges_(max(1,nr)))
    if(nr>0)then
      call c_f_pointer(p%ranges,ranges,[nr])
      do i=1,nr
        param%ranges_(i)%min_temperature_=ranges(i)%minimum
        param%ranges_(i)%max_temperature_=ranges(i)%maximum
        param%ranges_(i)%fixed_temperature_=ranges(i)%fixed_value
        param%ranges_(i)%is_fixed_=ranges(i)%is_fixed/=0
      enddo
    endif
  end subroutine

  subroutine ref_temperature(p,t,n,wl,out) bind(c)
    type(temperature_model), intent(in) :: p
    real(c_double), value :: t
    integer(c_int), value :: n
    real(c_double), intent(in) :: wl(n)
    real(c_double), intent(inout) :: out(n)
    class(temperature_parameterization_t), pointer :: param
    call make_parameter(p,n,wl,param)
    call param%calculate(t,wl,out)
    deallocate(param)
  end subroutine

  subroutine ref_temperature_cross(p,n,wl,base,nw,nz,edges,mids,temp,out) bind(c)
    use tuvx_cross_section_temperature_based
    use reference_yields, only: setup
    use tuvx_grid_warehouse
    use tuvx_profile_warehouse
    type(temperature_model), intent(in) :: p
    integer(c_int), value :: n,nw,nz
    real(c_double), intent(in) :: wl(n),base(n),edges(nw+1),mids(nw),temp(nz)
    real(c_double), intent(out) :: out(nz,nw)
    type(cross_section_temperature_based_t) :: cross
    class(grid_warehouse_t), pointer :: grids
    class(profile_warehouse_t), pointer :: profiles
    call setup(nw,nz,edges,mids,temp,temp,grids,profiles)
    cross%height_grid_=grids%get_ptr('height','km')
    cross%wavelength_grid_=grids%get_ptr('wavelength','nm')
    cross%temperature_profile_=profiles%get_ptr('temperature','K')
    cross%raw_wavelengths_=wl
    cross%raw_data_=base
    call make_parameter(p,n,wl,cross%parameterization_)
    out=cross%calculate(grids,profiles)
    deallocate(grids,profiles)
  end subroutine
end module
