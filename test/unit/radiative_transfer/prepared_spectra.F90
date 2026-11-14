! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module custom_prepared_fixture
  use musica_constants, only : dk => musica_dk
  use tuvx_cross_section, only : cross_section_t
  use tuvx_grid, only : grid_t
  use tuvx_grid_warehouse, only : grid_warehouse_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  implicit none
  integer :: custom_calls = 0
  type, extends(cross_section_t) :: custom_cross_section_t
  contains
    procedure :: calculate => custom_calculate
  end type
contains
  function custom_calculate(this,grid_warehouse,profile_warehouse,at_mid_point) result(values)
    class(custom_cross_section_t), intent(in) :: this
    type(grid_warehouse_t), intent(inout) :: grid_warehouse
    type(profile_warehouse_t), intent(inout) :: profile_warehouse
    logical, optional, intent(in) :: at_mid_point
    real(dk), allocatable :: values(:,:)
    class(grid_t), pointer :: height,wl
    height => grid_warehouse%get_grid('height','km')
    wl => grid_warehouse%get_grid('wavelength','nm')
    allocate(values(height%ncells_+1,wl%ncells_))
    custom_calls = custom_calls+1
    values = real(custom_calls,dk)
    deallocate(height,wl)
  end function
end module

program test_prepared_spectra
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
  use, intrinsic :: ieee_exceptions, only : ieee_set_halting_mode, ieee_invalid, ieee_overflow, ieee_divide_by_zero, &
      ieee_set_flag, ieee_all
  use custom_prepared_fixture, only : custom_cross_section_t, custom_calls
  use musica_assert, only : assert
  use musica_config, only : config_t
  use musica_constants, only : dk => musica_dk
  use musica_mpi, only : musica_mpi_init, musica_mpi_finalize
  use musica_string, only : string_t
  use tuvx_grid_warehouse, only : grid_warehouse_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  use tuvx_profile_from_host, only : profile_from_host_t
  use tuvx_cross_section, only : cross_section_t
  use tuvx_quantum_yield, only : quantum_yield_t
  use tuvx_cross_section_factory, only : cross_section_builder
  use tuvx_quantum_yield_factory, only : quantum_yield_builder
  use tuvx_prepared_spectra, only : spectrum_cache_t, spectrum_plan_t
  implicit none
  integer, parameter :: nl = 40
  type(config_t) :: catalogue, grid_set, height_config, wl_config, spec_config
  type(config_t), allocatable :: cases(:)
  type(spectrum_cache_t) :: cache
  type(custom_cross_section_t), target :: custom
  type(spectrum_plan_t), pointer :: plan, duplicate, midpoint_plan
  class(grid_warehouse_t), pointer :: grids
  type(profile_warehouse_t) :: profiles
  type(profile_from_host_t), pointer :: temperature, air, unused_temperature, unused_air
  class(cross_section_t), pointer :: cs, cs_duplicate
  class(quantum_yield_t), pointer :: qy, qy_duplicate
  real(dk), allocatable :: expected(:,:)
  type(string_t) :: origin
  integer :: i,j,role,column,comparisons,variable,total
  real(dk), parameter :: edges(*) = [180.0_dk,200.0_dk,210.0_dk,218.0_dk,233.0_dk, &
      250.0_dk,273.0_dk,295.0_dk,298.0_dk,300.0_dk,320.0_dk,350.0_dk,400.0_dk]
  ! Temperatures beyond every clamp, down to and up to the prepared limits
  real(dk), parameter :: extremes(*) = [1.0_dk,60.0_dk,120.0_dk,150.0_dk,179.0_dk,263.0_dk, &
      420.0_dk,600.0_dk,1000.0_dk,3000.0_dk,10000.0_dk]
  integer, parameter :: columns = 8

  call musica_mpi_init()
  ! The constructor catalogue includes deliberately large unit-test table
  ! coefficients. Compare their IEEE results too, with bounds checks enabled.
  call ieee_set_halting_mode(ieee_invalid,.false.)
  call ieee_set_halting_mode(ieee_overflow,.false.)
  call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
  call catalogue%from_file('test/data/prepared_spectra.json')
  call catalogue%get('cases',cases,'prepared spectra test')
  height_config = '{"name":"height","type":"equal interval","units":"km",'// &
                   '"begins at":0,"ends at":40,"cell delta":1}'
  allocate(profiles%profiles_(4))
  temperature => profile_from_host_t('temperature','K',nl)
  air => profile_from_host_t('air','molecule cm-3',nl)
  profiles%profiles_(1)%val_ => temperature
  profiles%profiles_(2)%val_ => air
  unused_temperature => profile_from_host_t('temperature','degC',nl)
  unused_air => profile_from_host_t('air','kg m-3',nl)
  unused_temperature%edge_val_ = -100
  unused_temperature%mid_val_ = -100
  unused_air%edge_val_ = -100
  unused_air%mid_val_ = -100
  profiles%profiles_(3)%val_ => unused_temperature
  profiles%profiles_(4)%val_ => unused_air
  comparisons = 0
  variable = 0
  total = 0
  do i = 1, size(cases)
    call cases(i)%get('origin',origin,'prepared spectra test')
    call cases(i)%get('role',role,'prepared spectra test')
    call cases(i)%get('grid',wl_config,'prepared spectra test')
    call cases(i)%get('spectrum',spec_config,'prepared spectra test')
    call grid_set%empty()
    call grid_set%add('grids',[height_config,wl_config],'prepared spectra test')
    call grid_set%get('grids',wl_config,'prepared spectra test')
    grids => grid_warehouse_t(wl_config)
    temperature%edge_val_ = 260
    temperature%mid_val_ = 260
    air%edge_val_ = 2.45e19_dk
    air%mid_val_ = 2.45e19_dk
    if( role == 0 ) then
      cs => cross_section_builder(spec_config,grids,profiles)
      cs_duplicate => cross_section_builder(spec_config,grids,profiles)
    else
      qy => quantum_yield_builder(spec_config,grids,profiles)
      qy_duplicate => quantum_yield_builder(spec_config,grids,profiles)
    end if
    do column = 1, columns
      ! Interfaces and mid-points have different values, so a plan that read
      ! the wrong set would be detected. Column 7 uses the branch
      ! temperatures exactly; column 8 goes beyond every clamp.
      select case(column)
      case(7)
        temperature%edge_val_ = [(edges(mod(j*7+column,size(edges))+1),j=1,nl+1)]
        temperature%mid_val_ = [(edges(mod(j*5+column,size(edges))+1),j=1,nl)]
      case(8)
        temperature%edge_val_ = [(extremes(mod(j*7,size(extremes))+1),j=1,nl+1)]
        temperature%mid_val_ = [(extremes(mod(j*5+3,size(extremes))+1),j=1,nl)]
      case default
        temperature%edge_val_ = [(edges(mod(j*7+column,size(edges))+1)+0.123_dk,j=1,nl+1)]
        temperature%mid_val_ = [(edges(mod(j*5+column+2,size(edges))+1)+0.321_dk,j=1,nl)]
      end select
      air%edge_val_ = [(real(mod(j+column,7),dk)*.7e19_dk,j=1,nl+1)]
      air%mid_val_ = [(real(mod(j+column+3,7),dk)*.7e19_dk,j=1,nl)]
      call cache%begin_column(grids,profiles)
      call assert(826130016,cache%usable)
      if( role == 0 ) then
        plan => cache%cross_section(cs,grids,profiles)
        duplicate => cache%cross_section(cs_duplicate,grids,profiles)
        expected = cs%calculate(grids,profiles)
      else
        plan => cache%quantum_yield(qy,grids,profiles)
        duplicate => cache%quantum_yield(qy_duplicate,grids,profiles)
        expected = qy%calculate(grids,profiles)
      end if
      call assert(826130001,associated(plan,duplicate))
      if( .not. plan%guarded .and. all(ieee_is_finite(expected)) ) then
        print *, 'Constant classification guard failed:',i,origin%to_char()
        stop 1
      end if
      call compare(plan%materialize(),expected,i,column,origin%to_char())
      comparisons = comparisons+size(expected)
      if( role == 0 ) then
        ! Radiators use plans evaluated at layer mid-points.
        midpoint_plan => cache%cross_section(cs,grids,profiles,.true.)
        duplicate => cache%cross_section(cs_duplicate,grids,profiles,.true.)
        expected = cs%calculate(grids,profiles,.true.)
        call assert(826130017,associated(midpoint_plan,duplicate) .and. .not. associated(midpoint_plan,plan))
        if( .not. midpoint_plan%guarded .and. all(ieee_is_finite(expected)) ) then
          print *, 'Constant classification guard failed at mid-points:',i,origin%to_char()
          stop 1
        end if
        call compare(midpoint_plan%materialize(),expected,i,column,origin%to_char())
        comparisons = comparisons+size(expected)
      end if
    end do
    variable = variable+count(plan%varying)
    total = total+size(plan%varying)
    call cache%clear()
    if( role == 0 ) then
      deallocate(cs,cs_duplicate)
    else
      deallocate(qy,qy_duplicate)
    end if
    deallocate(grids)
  end do
  ! Unsupported dynamic types use the actual object, exactly once per request,
  ! including repeated requests within a column. No factory or clone is needed.
  wl_config = '{"name":"wavelength","type":"equal interval","units":"nm",'// &
               '"begins at":200,"ends at":400,"cell delta":10}'
  call grid_set%empty()
  call grid_set%add('grids',[height_config,wl_config],'custom prepared fixture')
  call grid_set%get('grids',wl_config,'custom prepared fixture')
  grids => grid_warehouse_t(wl_config)
  temperature%edge_val_ = 260
  temperature%mid_val_ = 260
  call cache%begin_column(grids,profiles)
  plan => cache%cross_section(custom,grids,profiles)
  call assert(826130006,custom_calls == 1 .and. all(plan%values == 1.0_dk))
  cache%epoch = huge(cache%epoch)
  call cache%begin_column(grids,profiles)
  call assert(826130015,cache%epoch == 1 .and. cache%usable)
  plan => cache%cross_section(custom,grids,profiles)
  call assert(826130007,custom_calls == 2 .and. all(plan%values == 2.0_dk))
  temperature%edge_val_(1) = ieee_value(0.0_dk,ieee_quiet_nan)
  call cache%begin_column(grids,profiles)
  call assert(826130008,.not. cache%usable)
  temperature%edge_val_ = 11000
  call cache%begin_column(grids,profiles)
  call assert(826130009,.not. cache%usable)
  call cache%clear()
  deallocate(grids)
  print *, 'Compared spectrum values:',comparisons,'cases:',size(cases)
  print *, 'Varying wavelengths:',variable,'of',total
  call ieee_set_flag(ieee_all,.false.)
  call musica_mpi_finalize()
contains
  subroutine compare(actual,reference,case_index,column,origin)
    real(dk), intent(in) :: actual(:,:),reference(:,:)
    integer, intent(in) :: case_index,column
    character(*), intent(in) :: origin
    real(dk), allocatable :: bound(:,:)
    integer :: worst(2)
    call assert(826130002,all(shape(actual) == shape(reference)))
    call assert(826130003,all(ieee_is_finite(actual) .eqv. ieee_is_finite(reference)))
    call assert(826130004,all(ieee_is_nan(actual) .eqv. ieee_is_nan(reference)))
    call assert(826130005,all(ieee_is_finite(reference) .or. ieee_is_nan(reference) .or. actual == reference))
    bound = spread(max(maxval(merge(abs(actual),0.0_dk,ieee_is_finite(actual)),dim=1), &
                       maxval(merge(abs(reference),0.0_dk,ieee_is_finite(reference)),dim=1), &
                       tiny(1.0_dk)),1,size(actual,1))
    if( any(ieee_is_finite(reference) .and. abs(actual-reference) > 64*epsilon(1.0_dk)*bound) ) then
      print *, 'Spectrum mismatch:',case_index,column,origin
      print *, 'Max normalized eps:',maxval(abs(actual-reference)/bound/epsilon(1.0_dk))
      worst = maxloc(merge(abs(actual-reference)/bound,0.0_dk,ieee_is_finite(reference)))
      print *, 'Height, wavelength:',worst,'actual, reference:',actual(worst(1),worst(2)),reference(worst(1),worst(2))
      stop 3
    end if
  end subroutine
end program
