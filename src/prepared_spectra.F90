! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Spectrum/rate plans adapt TurboTUV's spectrum_plan.c and model_plan.c.
! Prepared spectra evaluate the existing formulas only at varying wavelengths.
! Constant bins are proved from the kernel's branches/data, then guarded with
! generic evaluations. No science formula is replaced by a fitted surrogate.
module tuvx_prepared_spectra
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: ieee_exceptions, only : ieee_status_type, ieee_get_status, ieee_set_status, &
      ieee_set_halting_mode, ieee_invalid, ieee_overflow, ieee_divide_by_zero
  use musica_constants, only : dk => musica_dk
  use tuvx_constants, only : pi, gas_constant
  use musica_string, only : string_t
  use tuvx_util, only : add_point
  use tuvx_grid, only : grid_t
  use tuvx_grid_warehouse, only : grid_warehouse_t
  use tuvx_profile_warehouse, only : profile_warehouse_t
  use tuvx_cross_section, only : cross_section_t
  use tuvx_quantum_yield, only : quantum_yield_t
  use tuvx_cross_section_factory, only : cross_section_type_name
  use tuvx_quantum_yield_factory, only : quantum_yield_type_name
  use tuvx_temperature_parameterization, only : temperature_parameterization_t
  use tuvx_temperature_parameterization_taylor_series, only : temperature_parameterization_taylor_series_t
  use tuvx_temperature_parameterization_burkholder, only : temperature_parameterization_burkholder_t
  use tuvx_temperature_parameterization_harwood, only : temperature_parameterization_harwood_t
  use tuvx_cross_section_hcfc, only : cross_section_hcfc_t
  use tuvx_cross_section_o3_jpl06, only : cross_section_o3_jpl06_t
  use tuvx_cross_section_n2o5_no2_no3, only : cross_section_n2o5_no2_no3_t
  use tuvx_cross_section_ch2o, only : cross_section_ch2o_t
  use tuvx_cross_section_cfc11, only : cross_section_cfc11_t
  use tuvx_cross_section_nitroxy_acetone, only : cross_section_nitroxy_acetone_t
  use tuvx_cross_section_hobr_oh_br, only : cross_section_hobr_oh_br_t
  use tuvx_cross_section_hno3_oh_no2, only : cross_section_hno3_oh_no2_t
  use tuvx_cross_section_temperature_based, only : cross_section_temperature_based_t
  use tuvx_cross_section_chcl3, only : cross_section_chcl3_t
  use tuvx_cross_section_nitroxy_ethanol, only : cross_section_nitroxy_ethanol_t
  use tuvx_cross_section_oclo, only : cross_section_oclo_t
  use tuvx_cross_section_clono2, only : cross_section_clono2_t
  use tuvx_cross_section_rono2, only : cross_section_rono2_t
  use tuvx_cross_section_o3_tint, only : cross_section_o3_tint_t
  use tuvx_cross_section_rayliegh, only : cross_section_rayliegh_t
  use tuvx_cross_section_n2o_n2_o1d, only : cross_section_n2o_n2_o1d_t
  use tuvx_cross_section_no2_tint, only : cross_section_no2_tint_t
  use tuvx_cross_section_ccl4, only : cross_section_ccl4_t
  use tuvx_cross_section_cl2_cl_cl, only : cross_section_cl2_cl_cl_t
  use tuvx_cross_section_ch3ono2_ch3o_no2, only : cross_section_ch3ono2_ch3o_no2_t
  use tuvx_cross_section_tint, only : cross_section_tint_t
  use tuvx_cross_section_ch3coch3_ch3co_ch3, only : cross_section_ch3coch3_ch3co_ch3_t
  use tuvx_cross_section_t_butyl_nitrate, only : cross_section_t_butyl_nitrate_t
  use tuvx_cross_section_chbr3, only : cross_section_chbr3_t
  use tuvx_cross_section_h2o2_oh_oh, only : cross_section_h2o2_oh_oh_t
  use tuvx_quantum_yield_ch2chcho, only : quantum_yield_ch2chcho_t
  use tuvx_quantum_yield_mvk, only : quantum_yield_mvk_t
  use tuvx_quantum_yield_ch2o_h2_co, only : quantum_yield_ch2o_h2_co_t
  use tuvx_quantum_yield_o3_o2_o1d, only : quantum_yield_o3_o2_o1d_t
  use tuvx_quantum_yield_clo_cl_o3p, only : quantum_yield_clo_cl_o3p_t
  use tuvx_quantum_yield_o3_o2_o3p, only : quantum_yield_o3_o2_o3p_t
  use tuvx_quantum_yield_clo_cl_o1d, only : quantum_yield_clo_cl_o1d_t
  use tuvx_quantum_yield_c2h5cho_c2h5_hco, only : quantum_yield_c2h5cho_c2h5_hco_t
  use tuvx_quantum_yield_ch3cocho_ch3co_hco, only : quantum_yield_ch3cocho_ch3co_hco_t
  use tuvx_quantum_yield_taylor_series, only : quantum_yield_taylor_series_t
  use tuvx_quantum_yield_ho2_oh_o, only : quantum_yield_ho2_oh_o_t
  use tuvx_quantum_yield_ch3cho_ch3_hco, only : quantum_yield_ch3cho_ch3_hco_t
  use tuvx_quantum_yield_no2_tint, only : quantum_yield_no2_tint_t
  use tuvx_quantum_yield_no3m_aq, only : quantum_yield_no3m_aq_t
  use tuvx_quantum_yield_clono2_cl_no3, only : quantum_yield_clono2_cl_no3_t
  use tuvx_quantum_yield_h2so4_mills, only : quantum_yield_h2so4_mills_t
  use tuvx_quantum_yield_tint, only : quantum_yield_tint_t
  use tuvx_quantum_yield_ch3coch2ch3, only : quantum_yield_ch3coch2ch3_t
  use tuvx_quantum_yield_clono2_clo_no2, only : quantum_yield_clono2_clo_no2_t
  use tuvx_quantum_yield_ch3coch3_ch3co_ch3, only : quantum_yield_ch3coch3_ch3co_ch3_t
  implicit none
  private
  public :: spectrum_cache_t, spectrum_plan_t, rate_plan_t, prepare_rate, integrate_rate

  type interval_t
    integer :: first, last
    type(grid_warehouse_t) :: grids
    class(cross_section_t), allocatable :: cs
    class(quantum_yield_t), allocatable :: qy
    real(dk), allocatable :: raw_grid(:), raw_values(:), target(:), param_powers(:,:)
    integer, allocatable :: raw_positions(:), segment_start(:), segment_index(:)
    real(dk), allocatable :: segment_delta(:), segment_left(:), segment_right(:), segment_area(:), target_width(:)
  end type

  type spectrum_plan_t
    class(cross_section_t), pointer :: source_cs => null()
    class(quantum_yield_t), pointer :: source_qy => null()
    character(:), allocatable :: key
    logical :: midpoint = .false., supported = .true., guarded = .true.
    integer, allocatable :: table_index(:)
    integer :: epoch = -1, kernel = 0
    real(dk), allocatable :: coefficient(:,:), level(:,:), wavelength(:)
    real(dk), allocatable :: constant(:), values(:,:)
    logical, allocatable :: varying(:)
    integer, allocatable :: position(:)
    type(interval_t), allocatable :: intervals(:)
    real(dk), allocatable :: oxygen(:,:)
    integer :: oxygen_epoch = -1
  contains
    procedure :: materialize
  end type
  type plan_ptr_t
    type(spectrum_plan_t), pointer :: p => null()
  end type
  type alias_t
    class(cross_section_t), pointer :: cs => null()
    class(quantum_yield_t), pointer :: qy => null()
    type(spectrum_plan_t), pointer :: p => null()
    logical :: midpoint = .false.
  end type
  type spectrum_cache_t
    type(plan_ptr_t), allocatable :: plans(:)
    type(alias_t), allocatable :: aliases(:)
    integer :: epoch = 0, generation = 0
    logical :: usable = .false., enabled = .true.
    integer :: levels = 0
    real(dk), allocatable :: wavelengths(:), edges(:)
  contains
    procedure :: begin_column
    procedure :: cross_section => get_cross
    procedure :: quantum_yield => get_yield
    procedure :: clear
    final :: finalize_cache
  end type
  type rate_plan_t
    type(spectrum_plan_t), pointer :: cs => null(), qy => null()
    integer, allocatable :: indices(:)
    logical, allocatable :: constant_term(:)
    real(dk), allocatable :: product(:), energy(:)
    logical :: heating = .false., oxygen = .false.
    integer :: generation = -1
  end type
contains
  function materialize(this) result(values)
    class(spectrum_plan_t), intent(in) :: this
    real(dk), allocatable :: values(:,:)
    integer :: w
    values = spread(this%constant,1,size(this%values,1))
    do w = 1, size(this%varying)
      if( this%varying(w) ) values(:,w) = this%values(:,this%position(w))
    end do
  end function

  subroutine clear(this)
    class(spectrum_cache_t), intent(inout) :: this
    integer :: i
    if( allocated(this%aliases) ) deallocate(this%aliases)
    if( allocated(this%plans) ) then
      do i = 1, size(this%plans)
        deallocate(this%plans(i)%p)
      end do
      deallocate(this%plans)
    end if
    allocate(this%plans(0),this%aliases(0))
    this%generation = this%generation+1
  end subroutine
  subroutine finalize_cache(this)
    type(spectrum_cache_t), intent(inout) :: this
    call this%clear()
  end subroutine

  subroutine begin_column(this,grids,profiles)
    class(spectrum_cache_t), intent(inout) :: this
    type(grid_warehouse_t), intent(inout) :: grids
    type(profile_warehouse_t), intent(inout) :: profiles
    class(grid_t), pointer :: wl, height
    integer :: i
    logical :: changed, finite_grid
    wl => grids%get_grid('wavelength','nm')
    height => grids%get_grid('height','km')
    finite_grid = all(ieee_is_finite(wl%mid_)) .and. all(ieee_is_finite(wl%edge_))
    changed = .not. allocated(this%wavelengths)
    if( .not. changed ) changed = size(this%wavelengths) /= size(wl%mid_) .or. this%levels /= height%ncells_+1
    if( .not. changed ) then
      if( finite_grid .and. all(ieee_is_finite(this%wavelengths)) .and. all(ieee_is_finite(this%edges)) ) then
        changed = any(this%wavelengths /= wl%mid_) .or. any(this%edges /= wl%edge_)
      else
        changed = .true.
      end if
    end if
    if( changed ) then
      call this%clear()
      this%wavelengths = wl%mid_
      this%edges = wl%edge_
      this%levels = height%ncells_+1
    end if
    if( this%epoch == huge(this%epoch) ) then
      this%epoch = 0
      do i = 1, size(this%plans)
        this%plans(i)%p%epoch = -1
        this%plans(i)%p%oxygen_epoch = -1
      end do
    end if
    this%epoch = this%epoch+1
    this%usable = this%enabled .and. finite_grid
    do i = 1, size(profiles%profiles_)
      associate(p => profiles%profiles_(i)%val_)
        if( p%handle_ == 'temperature' .and. p%units_ == 'K' ) then
          if( allocated(p%edge_val_) ) this%usable = this%usable .and. valid_temperature(p%edge_val_)
          if( allocated(p%mid_val_) ) this%usable = this%usable .and. valid_temperature(p%mid_val_)
        else if( p%handle_ == 'air' .and. p%units_ == 'molecule cm-3' ) then
          if( allocated(p%edge_val_) ) this%usable = this%usable .and. valid_density(p%edge_val_)
          if( allocated(p%mid_val_) ) this%usable = this%usable .and. valid_density(p%mid_val_)
        end if
      end associate
    end do
    deallocate(wl,height)
  end subroutine

  logical function valid_temperature(v) result(valid)
    real(dk), intent(in) :: v(:)
    valid = .false.
    if( .not. all(ieee_is_finite(v)) ) return
    valid = all(v > 0) .and. all(v <= 1.e4_dk)
  end function
  logical function valid_density(v) result(valid)
    real(dk), intent(in) :: v(:)
    valid = .false.
    if( .not. all(ieee_is_finite(v)) ) return
    valid = all(v >= 0)
  end function
  logical function known_cross(s) result(known)
    class(cross_section_t), intent(in) :: s
    known = .true.
    select type(s)
    type is(cross_section_t)
    type is(cross_section_hcfc_t)
    type is(cross_section_o3_jpl06_t)
    type is(cross_section_n2o5_no2_no3_t)
    type is(cross_section_ch2o_t)
    type is(cross_section_cfc11_t)
    type is(cross_section_nitroxy_acetone_t)
    type is(cross_section_hobr_oh_br_t)
    type is(cross_section_hno3_oh_no2_t)
    type is(cross_section_temperature_based_t)
    type is(cross_section_chcl3_t)
    type is(cross_section_nitroxy_ethanol_t)
    type is(cross_section_oclo_t)
    type is(cross_section_clono2_t)
    type is(cross_section_rono2_t)
    type is(cross_section_o3_tint_t)
    type is(cross_section_rayliegh_t)
    type is(cross_section_n2o_n2_o1d_t)
    type is(cross_section_no2_tint_t)
    type is(cross_section_ccl4_t)
    type is(cross_section_cl2_cl_cl_t)
    type is(cross_section_ch3ono2_ch3o_no2_t)
    type is(cross_section_tint_t)
    type is(cross_section_ch3coch3_ch3co_ch3_t)
    type is(cross_section_t_butyl_nitrate_t)
    type is(cross_section_chbr3_t)
    type is(cross_section_h2o2_oh_oh_t)
    class default
      known = .false.
    end select
  end function

  logical function known_yield(s) result(known)
    class(quantum_yield_t), intent(in) :: s
    known = .true.
    select type(s)
    type is(quantum_yield_t)
    type is(quantum_yield_ch2chcho_t)
    type is(quantum_yield_mvk_t)
    type is(quantum_yield_ch2o_h2_co_t)
    type is(quantum_yield_o3_o2_o1d_t)
    type is(quantum_yield_clo_cl_o3p_t)
    type is(quantum_yield_o3_o2_o3p_t)
    type is(quantum_yield_clo_cl_o1d_t)
    type is(quantum_yield_c2h5cho_c2h5_hco_t)
    type is(quantum_yield_ch3cocho_ch3co_hco_t)
    type is(quantum_yield_taylor_series_t)
    type is(quantum_yield_ho2_oh_o_t)
    type is(quantum_yield_ch3cho_ch3_hco_t)
    type is(quantum_yield_no2_tint_t)
    type is(quantum_yield_no3m_aq_t)
    type is(quantum_yield_clono2_cl_no3_t)
    type is(quantum_yield_h2so4_mills_t)
    type is(quantum_yield_tint_t)
    type is(quantum_yield_ch3coch2ch3_t)
    type is(quantum_yield_clono2_clo_no2_t)
    type is(quantum_yield_ch3coch3_ch3co_ch3_t)
    class default
      known = .false.
    end select
  end function

  function key_cross(cs) result(key)
    class(cross_section_t), intent(in) :: cs
    character(:), allocatable :: key
    type(string_t) :: name
    real(dk), allocatable :: v(:)
    integer :: i
    name = cross_section_type_name(cs)
    v = [0.0_dk]
    if( allocated(cs%cross_section_parms) ) then
    v = [real(size(cs%cross_section_parms),dk)]
    do i = 1, size(cs%cross_section_parms)
      associate(p => cs%cross_section_parms(i))
        call append(v,p%array)
        call append1(v,p%temperature)
        call append1(v,p%deltaT)
      end associate
    end do
    end if
    v = [v,0.0_dk]
    if( allocated(cs%overrides_) ) then
      v(size(v)) = real(size(cs%overrides_),dk)
      do i = 1, size(cs%overrides_)
        v = [v,real(cs%overrides_(i)%min_wavelength_index_,dk), &
             real(cs%overrides_(i)%max_wavelength_index_,dk),cs%overrides_(i)%value_]
      end do
    end if
    select type(cs)
    type is(cross_section_o3_tint_t)
      v = [v,cs%v185,cs%v195,cs%v345]
    type is(cross_section_temperature_based_t)
      call append1(v,cs%raw_wavelengths_)
      call append1(v,cs%raw_data_)
      associate(p => cs%parameterization_)
        call append1(v,p%AA_)
        call append1(v,p%BB_)
        call append1(v,p%lp_)
        call append1(v,p%wavelengths_)
        v = [v,p%base_temperature_,p%base_wavelength_,real(p%min_wavelength_index_,dk), &
             real(p%max_wavelength_index_,dk),real(merge(1,0,p%is_base_10_),dk), &
             real(merge(1,0,p%is_temperature_inverted_),dk)]
        v = [v,real(size(p%ranges_),dk)]
        do i = 1, size(p%ranges_)
          v = [v,p%ranges_(i)%min_temperature_,p%ranges_(i)%max_temperature_, &
               p%ranges_(i)%fixed_temperature_,real(merge(1,0,p%ranges_(i)%is_fixed_),dk)]
        end do
        select type(p)
        type is(temperature_parameterization_taylor_series_t)
          call append1(v,p%sigma_)
          call append(v,p%A_)
          name = name//' Taylor'
        type is(temperature_parameterization_burkholder_t)
          v = [v,p%A_,p%B_]
          name = name//' Burkholder'
        type is(temperature_parameterization_harwood_t)
          name = name//' Harwood'
        end select
      end associate
    end select
    key = "cs:"//name%to_char()//achar(0)//bytes(v)
  end function
  function key_yield(qy) result(key)
    class(quantum_yield_t), intent(in) :: qy
    character(:), allocatable :: key
    type(string_t) :: name
    real(dk), allocatable :: v(:)
    integer :: i
    name = quantum_yield_type_name(qy)
    v = [0.0_dk]
    if( allocated(qy%quantum_yield_parms) ) then
      v = [real(size(qy%quantum_yield_parms),dk)]
      do i = 1, size(qy%quantum_yield_parms)
        call append(v,qy%quantum_yield_parms(i)%array)
        call append1(v,qy%quantum_yield_parms(i)%temperature)
      end do
    end if
    v = [v,0.0_dk]
    if( allocated(qy%overrides_) ) then
      v(size(v)) = real(size(qy%overrides_),dk)
      do i = 1, size(qy%overrides_)
        v = [v,real(qy%overrides_(i)%min_wavelength_index_,dk), &
             real(qy%overrides_(i)%max_wavelength_index_,dk),qy%overrides_(i)%value_]
      end do
    end if
    select type(qy)
    type is(quantum_yield_tint_t)
      v = [v,real(size(qy%parameters),dk)]
      do i = 1, size(qy%parameters)
        call append(v,qy%parameters(i)%array)
        call append1(v,qy%parameters(i)%temperature)
        call append1(v,qy%parameters(i)%deltaT)
      end do
    type is(quantum_yield_no2_tint_t)
      v = [v,real(size(qy%parameters),dk)]
      do i = 1, size(qy%parameters)
        call append(v,qy%parameters(i)%array)
        call append1(v,qy%parameters(i)%temperature)
        call append1(v,qy%parameters(i)%deltaT)
      end do
    type is(quantum_yield_ch3coch3_ch3co_ch3_t)
      v = [v,real(merge(1,0,qy%do_CO_),dk),real(merge(1,0,qy%do_CH3CO_),dk), &
           qy%minimum_temperature_,qy%maximum_temperature_,qy%low_wavelength_value_,qy%high_wavelength_value_]
    type is(quantum_yield_h2so4_mills_t)
      v = [v,real(size(qy%wavelength_indices_),dk),real(qy%wavelength_indices_,dk),qy%molecular_diameter_,qy%molecular_weight_]
      call append1(v,qy%collision_rate_)
    end select
    key = "qy:"//name%to_char()//achar(0)//bytes(v)
  end function
  function bytes(v) result(s)
    real(dk), intent(in) :: v(:)
    character(:), allocatable :: s
    allocate(character(len=size(v)*storage_size(v)/8)::s)
    s = transfer(v,s)
  end function
  subroutine append(v,a)
    real(dk), allocatable, intent(inout) :: v(:)
    real(dk), allocatable, intent(in) :: a(:,:)
    if( allocated(a) ) then
      v = [v,real(shape(a),dk),reshape(a,[size(a)])]
    else
      v = [v,-1.0_dk,-1.0_dk]
    end if
  end subroutine
  subroutine append1(v,a)
    real(dk), allocatable, intent(inout) :: v(:)
    real(dk), allocatable, intent(in) :: a(:)
    if( allocated(a) ) then
      v = [v,real(size(a),dk),a]
    else
      v = [v,-1.0_dk]
    end if
  end subroutine

  function get_cross(this,cs,grids,profiles,midpoint) result(p)
    class(spectrum_cache_t), intent(inout) :: this
    class(cross_section_t), target, intent(in) :: cs
    type(grid_warehouse_t), intent(inout) :: grids
    type(profile_warehouse_t), intent(inout) :: profiles
    logical, optional, intent(in) :: midpoint
    type(spectrum_plan_t), pointer :: p
    type(alias_t) :: alias
    type(plan_ptr_t) :: entry
    character(:), allocatable :: key
    logical :: mid, created
    integer :: i
    mid = .false.
    if( present(midpoint) ) mid = midpoint
    if( .not. allocated(this%plans) ) call this%begin_column(grids,profiles)
    do i = 1, size(this%aliases)
      if( associated(this%aliases(i)%cs,cs) .and. (this%aliases(i)%midpoint .eqv. mid) ) then
        p => this%aliases(i)%p
        call evaluate(p,this%epoch,grids,profiles)
        return
      end if
    end do
    if( known_cross(cs) ) then
      key = key_cross(cs)
    else
      key = 'custom cs:'//bytes([real(size(this%plans),dk)])
    end if
    created = .false.
    nullify(p)
    do i = 1, size(this%plans)
      if( len(this%plans(i)%p%key) /= len(key) ) cycle
      if( this%plans(i)%p%key /= key ) cycle
      if( this%plans(i)%p%midpoint .neqv. mid ) cycle
      p => this%plans(i)%p
      exit
    end do
    if( .not. associated(p) ) then
      created = .true.
      allocate(p)
      p%source_cs => cs
      p%key = key
      p%midpoint = mid
      call build(p,grids,profiles)
      entry%p => p
      this%plans = [this%plans,entry]
    end if
    alias%cs => cs
    alias%p => p
    alias%midpoint = mid
    this%aliases = [this%aliases,alias]
    if( .not. created ) call evaluate(p,this%epoch,grids,profiles)
  end function
  function get_yield(this,qy,grids,profiles) result(p)
    class(spectrum_cache_t), intent(inout) :: this
    class(quantum_yield_t), target, intent(in) :: qy
    type(grid_warehouse_t), intent(inout) :: grids
    type(profile_warehouse_t), intent(inout) :: profiles
    type(spectrum_plan_t), pointer :: p
    type(alias_t) :: alias
    type(plan_ptr_t) :: entry
    character(:), allocatable :: key
    logical :: created
    integer :: i
    if( .not. allocated(this%plans) ) call this%begin_column(grids,profiles)
    do i = 1, size(this%aliases)
      if( associated(this%aliases(i)%qy,qy) ) then
        p => this%aliases(i)%p
        call evaluate(p,this%epoch,grids,profiles)
        return
      end if
    end do
    if( known_yield(qy) ) then
      key = key_yield(qy)
    else
      key = 'custom qy:'//bytes([real(size(this%plans),dk)])
    end if
    created = .false.
    nullify(p)
    do i = 1, size(this%plans)
      if( len(this%plans(i)%p%key) /= len(key) ) cycle
      if( this%plans(i)%p%key /= key ) cycle
      p => this%plans(i)%p
      exit
    end do
    if( .not. associated(p) ) then
      created = .true.
      allocate(p)
      p%source_qy => qy
      p%key = key
      call build(p,grids,profiles)
      entry%p => p
      this%plans = [this%plans,entry]
    end if
    alias%qy => qy
    alias%p => p
    this%aliases = [this%aliases,alias]
    if( .not. created ) call evaluate(p,this%epoch,grids,profiles)
  end function

  subroutine build(p,grids,profiles)
    type(spectrum_plan_t), intent(inout) :: p
    type(grid_warehouse_t), intent(inout) :: grids
    type(profile_warehouse_t), intent(inout) :: profiles
    class(grid_t), pointer :: wl
    real(dk), allocatable :: reference(:,:)
    integer :: i, j, runs, first, last
    type(ieee_status_type) :: saved_status
    wl => grids%get_grid('wavelength','nm')
    allocate(p%varying(wl%ncells_))
    if( associated(p%source_cs) ) then
      reference = p%source_cs%calculate(grids,profiles,p%midpoint)
    else
      reference = p%source_qy%calculate(grids,profiles)
    end if
    ! Only the generic evaluation above is an evaluation of the caller's
    ! current column. Configuration-only preparation/probes preserve its FP
    ! environment, even for tables that overflow in other temperature ranges.
    call ieee_get_status(saved_status)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    if( associated(p%source_cs) ) then
      call mask_cross(p%source_cs,wl,p%varying,p%supported)
    else
      call mask_yield(p%source_qy,wl,p%varying,p%supported)
    end if
    p%constant = reference(1,:)
    p%values = reference
    call prepare_kernel(p,wl)
    call guard(p,grids,profiles)
    allocate(p%position(size(p%varying)))
    p%position = 0
    j = 0
    do i = 1, size(p%varying)
      if( .not. p%varying(i) ) cycle
      j = j+1
      p%position(i) = j
    end do
    p%values = reference(:,pack([(i,i=1,size(p%varying))],p%varying))
    runs = 0
    i = 1
    do while(i <= size(p%varying))
      if( .not. p%varying(i) ) then
        i = i+1
        cycle
      end if
      runs = runs+1
      do while(i < size(p%varying))
        if( .not. p%varying(i+1) ) exit
        i = i+1
      end do
      i = i+1
    end do
    allocate(p%intervals(runs))
    i = 1
    do j = 1, runs
      do while(.not. p%varying(i))
        i = i+1
      end do
      first = i
      do while(i < size(p%varying))
        if( .not. p%varying(i+1) ) exit
        i = i+1
      end do
      last = i
      p%intervals(j)%first = first
      p%intervals(j)%last = last
      if( p%supported .and. p%kernel == 0 ) call make_interval(p,p%intervals(j),grids,first,last)
      i = i+1
    end do
    call ieee_set_status(saved_status)
    deallocate(wl)
  end subroutine

  subroutine guard(p,grids,profiles)
    use, intrinsic :: ieee_exceptions, only : ieee_status_type, ieee_get_status, ieee_set_status, &
        ieee_set_halting_mode, ieee_invalid, ieee_overflow, ieee_divide_by_zero
    type(spectrum_plan_t), intent(inout) :: p
    type(grid_warehouse_t), intent(inout) :: grids
    type(profile_warehouse_t), intent(in) :: profiles
    type(profile_warehouse_t) :: checks
    real(dk), allocatable :: reference(:,:)
    integer :: i, k, r, w
    type(ieee_status_type) :: saved_status
    if( all(p%varying) ) return
    ! Probe columns must not alter the caller's flags or halting modes. Some
    ! synthetic fixture coefficients overflow outside their intended range.
    call ieee_get_status(saved_status)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    allocate(checks%profiles_(size(profiles%profiles_)))
    do i = 1, size(checks%profiles_)
      allocate(checks%profiles_(i)%val_,source=profiles%profiles_(i)%val_)
    end do
reference_loop: do r = 1, 3
      do i = 1, size(checks%profiles_)
        associate(q => checks%profiles_(i)%val_)
          if( q%handle_ == 'temperature' .and. q%units_ == 'K' ) then
            q%edge_val_ = [(180.0_dk+real(mod(k*31+r*71,221),dk),k=1,size(q%edge_val_))]
            q%mid_val_ = [(180.0_dk+real(mod(k*31+r*71,221),dk),k=1,size(q%mid_val_))]
          else if( q%handle_ == 'air' .and. q%units_ == 'molecule cm-3' ) then
            q%edge_val_ = [(real(mod(k+r,5),dk)*1.e20_dk,k=1,size(q%edge_val_))]
            q%mid_val_ = [(real(mod(k+r,5),dk)*1.e20_dk,k=1,size(q%mid_val_))]
          end if
        end associate
      end do
      if( associated(p%source_cs) ) then
        reference = p%source_cs%calculate(grids,checks,p%midpoint)
      else
        reference = p%source_qy%calculate(grids,checks)
      end if
      do w = 1, size(p%varying)
        if( p%varying(w) ) cycle
        if( any(reference(:,w) /= p%constant(w)) ) then
          ! A failed guard never permits stale or incorrect constants to run.
          p%guarded = .false.
          p%supported = .false.
          p%varying = .true.
          exit reference_loop
        end if
      end do
    end do reference_loop
    call ieee_set_status(saved_status)
  end subroutine

  subroutine make_interval(p,part,grids,first,last)
    type(spectrum_plan_t), intent(in) :: p
    type(interval_t), intent(inout) :: part
    type(grid_warehouse_t), intent(in) :: grids
    integer, intent(in) :: first,last
    integer :: i,n,k
    integer, allocatable :: selected(:)
    allocate(part%grids%grids_(size(grids%grids_)))
    do i = 1, size(grids%grids_)
      allocate(part%grids%grids_(i)%val_,source=grids%grids_(i)%val_)
      associate(w => part%grids%grids_(i)%val_)
        if( w%handle_ /= 'wavelength' .or. w%units_ /= 'nm' ) cycle
        w%ncells_ = last-first+1
        w%mid_ = w%mid_(first:last)
        w%edge_ = w%edge_(first:last+1)
        w%delta_ = w%delta_(first:last)
      end associate
    end do
    if( associated(p%source_cs) ) then
      allocate(part%cs,source=p%source_cs)
      select type(c => part%cs)
      type is(cross_section_temperature_based_t)
        ! Intrinsic assignment copies pointers; give this clone its own parameterization.
        nullify(c%parameterization_)
        select type(s => p%source_cs)
        type is(cross_section_temperature_based_t)
          allocate(c%parameterization_,source=s%parameterization_)
        end select
        call prepare_temperature(part,c,grids,first,last)
        return
      type is(cross_section_hcfc_t)
        return ! These rows are polynomial coefficients, not wavelength data.
      end select
      if( .not. allocated(part%cs%cross_section_parms) ) return
      do i = 1, size(part%cs%cross_section_parms)
        if( .not. allocated(part%cs%cross_section_parms(i)%array) ) cycle
        part%cs%cross_section_parms(i)%array = part%cs%cross_section_parms(i)%array(first:last,:)
      end do
    else
      allocate(part%qy,source=p%source_qy)
      if( allocated(part%qy%quantum_yield_parms) ) then
        do i = 1, size(part%qy%quantum_yield_parms)
          if( .not. allocated(part%qy%quantum_yield_parms(i)%array) ) cycle
          part%qy%quantum_yield_parms(i)%array = part%qy%quantum_yield_parms(i)%array(first:last,:)
        end do
      end if
      if( allocated(part%qy%overrides_) ) then
        part%qy%overrides_ = pack(part%qy%overrides_, &
            part%qy%overrides_%min_wavelength_index_ <= last .and. &
            part%qy%overrides_%max_wavelength_index_ >= first)
        part%qy%overrides_%min_wavelength_index_ = max(part%qy%overrides_%min_wavelength_index_,first)-first+1
        part%qy%overrides_%max_wavelength_index_ = min(part%qy%overrides_%max_wavelength_index_,last)-first+1
      end if
      select type(q => part%qy)
      type is(quantum_yield_tint_t)
        do i = 1, size(q%parameters)
          q%parameters(i)%array = q%parameters(i)%array(first:last,:)
        end do
      type is(quantum_yield_no2_tint_t)
        do i = 1, size(q%parameters)
          q%parameters(i)%array = q%parameters(i)%array(first:last,:)
        end do
      type is(quantum_yield_h2so4_mills_t)
        selected = pack([(i,i=1,size(q%wavelength_indices_))], &
                         q%wavelength_indices_ >= first .and. q%wavelength_indices_ <= last)
        q%collision_rate_ = q%collision_rate_(selected)
        q%wavelength_indices_ = q%wavelength_indices_(selected)-first+1
      end select
    end if
  end subroutine

  subroutine evaluate(p,epoch,grids,profiles)
    type(spectrum_plan_t), intent(inout) :: p
    integer, intent(in) :: epoch
    type(grid_warehouse_t), intent(inout) :: grids
    type(profile_warehouse_t), intent(inout) :: profiles
    integer :: i
    if( p%epoch == epoch .and. p%supported ) return
    if( .not. p%supported ) then
      if( associated(p%source_cs) ) then
        p%values = p%source_cs%calculate(grids,profiles,p%midpoint)
      else
        p%values = p%source_qy%calculate(grids,profiles)
      end if
    else if( evaluate_kernel(p,profiles) ) then
      ! Wavelength constants and level-only expressions are already prepared.
    else
      do i = 1, size(p%intervals)
        associate(g => p%intervals(i))
          if( allocated(g%raw_grid) ) then
            call evaluate_temperature(g,profiles,p%midpoint,p%values(:,p%position(g%first):p%position(g%last)))
          else if( associated(p%source_cs) ) then
            p%values(:,p%position(g%first):p%position(g%last)) = g%cs%calculate(g%grids,profiles,p%midpoint)
          else
            p%values(:,p%position(g%first):p%position(g%last)) = g%qy%calculate(g%grids,profiles)
          end if
        end associate
      end do
    end if
    p%epoch = epoch
  end subroutine

  subroutine prepare_temperature(g,cs,grids,first,last)
    type(interval_t), intent(inout) :: g
    type(cross_section_temperature_based_t), intent(in) :: cs
    type(grid_warehouse_t), intent(in) :: grids
    integer, intent(in) :: first,last
    class(grid_t), pointer :: wl
    integer :: i,k,n,nb
    real(dk) :: a1,a2
    real(dk), parameter :: delta = real(1.e-5,dk)
    g%raw_grid = cs%raw_wavelengths_
    g%raw_values = cs%raw_data_
    call add_point(g%raw_grid,g%raw_values,(1-delta)*g%raw_grid(1),0.0_dk)
    call add_point(g%raw_grid,g%raw_values,0.0_dk,0.0_dk)
    call add_point(g%raw_grid,g%raw_values,(1+delta)*g%raw_grid(size(g%raw_grid)),0.0_dk)
    call add_point(g%raw_grid,g%raw_values,1.e38_dk,0.0_dk)
    allocate(g%raw_positions(size(cs%raw_data_)))
    do i = 1, size(g%raw_positions)
      g%raw_positions(i) = findloc(g%raw_grid,cs%raw_wavelengths_(i),dim=1)
    end do
    select type(p => cs%parameterization_)
    type is(temperature_parameterization_t)
      allocate(g%param_powers(p%max_wavelength_index_-p%min_wavelength_index_+1,size(p%lp_)))
      do i = 1, size(p%lp_)
        g%param_powers(:,i) = (cs%raw_wavelengths_(p%min_wavelength_index_:p%max_wavelength_index_) &
                              -p%base_wavelength_)**p%lp_(i)
      end do
    end select
    wl => grids%get_grid('wavelength','nm')
    g%target = wl%edge_(first:last+1)
    deallocate(wl)
    nb = size(g%target)-1
    allocate(g%segment_start(nb+1))
    allocate(g%segment_index(0),g%segment_delta(0),g%segment_left(0),g%segment_right(0),g%segment_area(0))
    g%target_width = g%target(2:)-g%target(:nb)
    n = 1
    do i = 1, nb
      g%segment_start(i) = n
      do k = 1, size(g%raw_grid)-1
        if( g%raw_grid(k+1) <= g%target(i) .or. g%raw_grid(k) >= g%target(i+1) ) cycle
        a1 = max(g%raw_grid(k),g%target(i))
        a2 = min(g%raw_grid(k+1),g%target(i+1))
        g%segment_index = [g%segment_index,k]
        g%segment_delta = [g%segment_delta,g%raw_grid(k+1)-g%raw_grid(k)]
        g%segment_left = [g%segment_left,a1-g%raw_grid(k)]
        g%segment_right = [g%segment_right,a2-g%raw_grid(k)]
        g%segment_area = [g%segment_area,.5_dk*(a2-a1)]
        n = n+1
      end do
    end do
    g%segment_start(nb+1) = n
  end subroutine

  subroutine evaluate_temperature(g,profiles,midpoint,values)
    type(interval_t), intent(inout) :: g
    type(profile_warehouse_t), intent(in) :: profiles
    logical, intent(in) :: midpoint
    real(dk), intent(out) :: values(:,:)
    real(dk), pointer :: temperatures(:)
    real(dk), allocatable :: unique_temperatures(:), raw(:,:), keys(:,:), unique_keys(:,:), interpolated(:,:)
    integer, allocatable :: level_index(:)
    real(dk) :: slope,b1,b2,area
    integer :: i,r,j,k,n,nranges,u,w,segment
    nullify(temperatures)
    do i = 1, size(profiles%profiles_)
      associate(p => profiles%profiles_(i)%val_)
        if( p%handle_ /= 'temperature' .or. p%units_ /= 'K' ) cycle
        if( midpoint ) then
          temperatures => p%mid_val_
        else
          temperatures => p%edge_val_
        end if
      end associate
    end do
    select type(cs => g%cs)
    type is(cross_section_temperature_based_t)
      associate(p => cs%parameterization_)
        nranges = size(p%ranges_)
        n = size(temperatures)
        allocate(keys(2*nranges,n),unique_keys(2*nranges,n),level_index(n),unique_temperatures(n))
        keys = 0
        do i = 1, n
          do r = 1, nranges
            if( temperatures(i) < p%ranges_(r)%min_temperature_ .or. &
                temperatures(i) > p%ranges_(r)%max_temperature_ ) cycle
            keys(2*r-1,i) = 1
            keys(2*r,i) = temperatures(i)
            if( p%ranges_(r)%is_fixed_ ) keys(2*r,i) = p%ranges_(r)%fixed_temperature_
          end do
        end do
        ! Equal effective temperatures in every active range have identical
        ! raw spectra, including overlapping/fixed ranges and out-of-range rows.
        u = 0
        do i = 1, n
          do j = 1, u
            if( all(keys(:,i) == unique_keys(:,j)) ) exit
          end do
          if( j > u ) then
            u = j
            unique_keys(:,u) = keys(:,i)
            unique_temperatures(u) = temperatures(i)
          end if
          level_index(i) = j
        end do
        raw = spread(cs%raw_data_,2,u)
        if( allocated(g%param_powers) ) then
          call p%calculate_column(unique_temperatures(:u),cs%raw_wavelengths_,raw,g%param_powers)
        else
          call p%calculate_column(unique_temperatures(:u),cs%raw_wavelengths_,raw)
        end if
      end associate
      allocate(interpolated(u,size(values,2)))
      do j = 1, u
        g%raw_values(g%raw_positions) = raw(:,j)
        do w = 1, size(values,2)
          area = 0
          do segment = g%segment_start(w),g%segment_start(w+1)-1
            k = g%segment_index(segment)
            slope = (g%raw_values(k+1)-g%raw_values(k))/g%segment_delta(segment)
            b1 = g%raw_values(k)+slope*g%segment_left(segment)
            b2 = g%raw_values(k)+slope*g%segment_right(segment)
            area = area+g%segment_area(segment)*(b2+b1)
          end do
          interpolated(j,w) = area/g%target_width(w)
        end do
      end do
      do i = 1, size(values,1)
        values(i,:) = interpolated(level_index(i),:)
      end do
    end select
  end subroutine

  subroutine prepare_kernel(p,w)
    type(spectrum_plan_t), intent(inout) :: p
    class(grid_t), intent(in) :: w
    integer :: nw,nc,k,j
    real(dk) :: x
    nw = w%ncells_
    nc = 4
    p%wavelength = w%mid_
    if( associated(p%source_cs) ) then
      select type(c => p%source_cs)
      type is(cross_section_ccl4_t)
        p%kernel = 1
      type is(cross_section_chcl3_t)
        p%kernel = 2
      type is(cross_section_cfc11_t)
        p%kernel = 3
      type is(cross_section_ch3ono2_ch3o_no2_t)
        p%kernel = 4
      type is(cross_section_hno3_oh_no2_t)
        if( .not. p%midpoint ) p%kernel = 4
      type is(cross_section_rono2_t)
        p%kernel = 4
      type is(cross_section_chbr3_t)
        p%kernel = 5
      type is(cross_section_cl2_cl_cl_t)
        if( .not. p%midpoint ) p%kernel = 6
      type is(cross_section_h2o2_oh_oh_t)
        if( .not. p%midpoint ) p%kernel = 7
      type is(cross_section_n2o_n2_o1d_t)
        p%kernel = 8
      type is(cross_section_ch3coch3_ch3co_ch3_t)
        p%kernel = 14
      type is(cross_section_ch2o_t)
        p%kernel = 15
      type is(cross_section_clono2_t)
        p%kernel = 16
      type is(cross_section_oclo_t)
        p%kernel = 17
      type is(cross_section_n2o5_no2_no3_t)
        p%kernel = 10
      type is(cross_section_o3_jpl06_t)
        p%kernel = 11
      type is(cross_section_tint_t)
        p%kernel = 12
        nc = max(nc,2*size(c%cross_section_parms))
      type is(cross_section_no2_tint_t)
        p%kernel = 12
        nc = max(nc,2*size(c%cross_section_parms))
      type is(cross_section_o3_tint_t)
        p%kernel = 13
        nc = max(nc,2*size(c%cross_section_parms))
        allocate(p%table_index(nw))
        p%table_index = 1
        where(w%edge_(:nw) < c%v345(1)) p%table_index = 2
        where(w%edge_(:nw) < c%v195(1)) p%table_index = 4
        where(w%edge_(:nw) < c%v185(1)) p%table_index = 3
      type is(cross_section_hcfc_t)
        p%kernel = 9
        nc = max(nc,size(c%cross_section_parms(1)%array,1))
      end select
    else
      select type(q => p%source_qy)
      type is(quantum_yield_o3_o2_o1d_t)
        p%kernel = 101
      type is(quantum_yield_o3_o2_o3p_t)
        p%kernel = 102
      type is(quantum_yield_mvk_t)
        p%kernel = 103
      type is(quantum_yield_ch2o_h2_co_t)
        p%kernel = 105
      type is(quantum_yield_ch3cho_ch3_hco_t)
        p%kernel = 106
      type is(quantum_yield_c2h5cho_c2h5_hco_t)
        p%kernel = 107
      type is(quantum_yield_ch3cocho_ch3co_hco_t)
        p%kernel = 108
      type is(quantum_yield_no3m_aq_t)
        p%kernel = 109
      type is(quantum_yield_ch2chcho_t)
        p%kernel = 110
      type is(quantum_yield_ch3coch2ch3_t)
        p%kernel = 111
      type is(quantum_yield_h2so4_mills_t)
        p%kernel = 112
      type is(quantum_yield_tint_t)
        p%kernel = 113
        nc = max(nc,2*size(q%parameters))
      type is(quantum_yield_no2_tint_t)
        p%kernel = 114
        nc = max(nc,2*size(q%parameters))
      type is(quantum_yield_ch3coch3_ch3co_ch3_t)
        p%kernel = 104
      end select
    end if
    if( p%kernel == 0 .or. .not. any(p%varying) ) return
    allocate(p%coefficient(nw,nc),p%level(size(p%values,1),max(12,nc)))
    p%coefficient = 0
    if( p%kernel == 6 ) then
      ! Retain the reference's contiguous wavelength evaluation of log. Using
      ! scalar log here can amplify a vector/scalar libm rounding difference
      ! in the very small tails of the chlorine spectrum.
      p%coefficient(:,1) = log(329.5_dk/w%mid_)**2
      p%coefficient(:,2) = log(406.5_dk/w%mid_)**2
    end if
    do k = 1, nw
      if( .not. p%varying(k) ) cycle
      x = w%mid_(k)
      associate(a => p%coefficient(k,:))
        select case(p%kernel)
        case(1)
          a(1) = 1.0739_dk+x*(-1.6275e-2_dk+x*(8.8141e-5_dk+x*(-1.9811e-7_dk+1.5022e-10_dk*x)))
        case(2)
          a(1) = 3.7973_dk+x*(-7.0913e-2_dk+x*(4.9397e-4_dk+x*(-1.5226e-6_dk+x*1.7555e-9_dk)))
        case(3)
          a(1) = x-184.9_dk
        case(5)
          a(1) = .06183_dk-.000241_dk*x
          a(2) = 2.376_dk+.14757_dk*x
        case(7)
          a(1) = ((((((1.5534675e-13_dk*x-2.652014e-10_dk)*x+1.6878206e-7_dk)*x-4.035101e-5_dk)*x &
                    -4.4589016e-3_dk)*x+4.535649_dk)*x-9.2170972e2_dk)*x+6.4761e4_dk
          a(2) = (((-1.0924e-7_dk*x-3.0493e-5_dk)*x+1.1522e-1_dk)*x-5.1351e1_dk)*x+6.8123e3_dk
        case(8)
          a(1) = (((2.520672e-7_dk*x-1.777846e-4_dk)*x+4.301146e-2_dk)*x-4.071805_dk)*x+68.21023_dk
          a(2) = exp(((-1.881058e-5_dk*x+1.111572e-2_dk)*x-2.116255_dk)*x+123.4014_dk)
        case(9)
          block
            integer :: j
            do j = 1, size(p%source_cs%cross_section_parms(1)%array,1)
              a(j) = (x-206.214_dk)**(j-1)
            end do
          end block
        case(10)
          a(1) = 1000.0_dk*p%source_cs%cross_section_parms(2)%array(k,1)
        case(11)
          a(1) = p%source_cs%cross_section_parms(2)%array(k,1)
          a(2) = (p%source_cs%cross_section_parms(3)%array(k,1)-a(1))/(298.0_dk-218.0_dk)
          a(3) = p%source_cs%cross_section_parms(3)%array(k,1)
        case(101,102)
          a(1) = exp(-((304.225_dk-x)/5.576_dk)**4)
          a(2) = exp(-((314.957_dk-x)/6.601_dk)**2)
          a(3) = exp(-((310.737_dk-x)/2.187_dk)**2)
        case(103)
          a(1) = exp(-.055_dk*(x-308.0_dk))
        case(105)
          block
            real(dk) :: q1,q2,qt
            q1 = p%source_qy%quantum_yield_parms(1)%array(k,1)
            q2 = p%source_qy%quantum_yield_parms(1)%array(k,2)
            qt = 1-q1
            a(1) = (1-(q1+q2))/(2.45e19_dk*q2*qt)
            a(2) = .05_dk*(x-329.0_dk)
            a(3) = 1/qt
          end block
        case(106)
          block
            real(dk) :: q1,q2,q
            q1 = p%source_qy%quantum_yield_parms(1)%array(k,2)
            q2 = 1-p%source_qy%quantum_yield_parms(1)%array(k,1)
            q = 0
            if( q1 > 0 ) q = q2/q1-1
            a(1) = q1*(1+q)
            a(2) = q
          end block
        case(107)
          a(1) = 1/p%source_qy%quantum_yield_parms(1)%array(k,1)-1
        case(108)
          a(1) = max(min(1-(x-380.0_dk)/60.0_dk,1.0_dk),0.0_dk)
          a(2) = 1.36e8_dk*exp(-8793.0_dk/x)
        case(112)
          select type(q => p%source_qy)
          type is(quantum_yield_h2so4_mills_t)
            do j = 1, size(q%wavelength_indices_)
              if( q%wavelength_indices_(j) == k ) a(1) = q%collision_rate_(j)
            end do
          end select
        case(104)
          a(1) = x-248.0_dk
          a(2) = 1.e7_dk/x
          a(3) = a(2)-33113.0_dk
          a(4) = a(2)-30488.0_dk
        end select
      end associate
    end do
  end subroutine

  logical function evaluate_kernel(p,profiles) result(handled)
    type(spectrum_plan_t), intent(inout) :: p
    type(profile_warehouse_t), intent(in) :: profiles
    real(dk), pointer :: temperature(:),air(:)
    real(dk) :: t,bb,alpha,chi,q2,qfac1,qfac2,t300,arg,fco,fac,c1,c2,c3,c4,ozone,value
    integer :: i,j,k,nz,table,first_table
    real(dk) :: fraction
    handled = p%kernel /= 0
    if( .not. handled ) return
    if( .not. any(p%varying) ) return
    nullify(temperature,air)
    do i = 1, size(profiles%profiles_)
      associate(q => profiles%profiles_(i)%val_)
        if( q%handle_ == 'temperature' .and. q%units_ == 'K' ) then
          if( p%midpoint ) then
            temperature => q%mid_val_
          else
            temperature => q%edge_val_
          end if
        else if( q%handle_ == 'air' .and. q%units_ == 'molecule cm-3' ) then
          air => q%edge_val_
        end if
      end associate
    end do
    nz = size(p%values,1)
    do i = 1, nz
      select case(p%kernel)
      case(1,2)
        p%level(i,1) = min(max(temperature(i),210.0_dk),300.0_dk)-295.0_dk
      case(3)
        p%level(i,1) = 1.e-4_dk*(temperature(i)-298.0_dk)
      case(4)
        p%level(i,1) = temperature(i)-298.0_dk
      case(5)
        p%level(i,1) = 273.0_dk-temperature(i)
      case(6)
        bb = exp(402.7_dk/temperature(i))
        alpha = (bb*bb-1)/(bb*bb+1)
        p%level(i,1) = -99.0_dk*alpha
        p%level(i,2) = -91.5_dk*alpha
        p%level(i,3) = 1.e-20_dk*sqrt(alpha)
      case(7)
        t = min(max(temperature(i),200.0_dk),400.0_dk)
        p%level(i,1) = 1/(1+exp(-1265.0_dk/t))
        p%level(i,2) = 1-p%level(i,1)
      case(8)
        p%level(i,1) = max(194.0_dk,min(temperature(i),320.0_dk))-300.0_dk
      case(9)
        t = min(295.0_dk,max(203.0_dk,temperature(i)))-p%source_cs%cross_section_parms(1)%temperature(1)
        associate(a => p%source_cs%cross_section_parms(1)%array)
          do j = 1, size(a,1)
            p%level(i,j) = a(j,1)+t*(a(j,2)+t*a(j,3))
          end do
        end associate
      case(14)
        p%level(i,1) = max(235.0_dk,min(298.0_dk,temperature(i)))
      case(15)
        p%level(i,1) = temperature(i)-298.0_dk
      case(16)
        p%level(i,1) = temperature(i)-296.0_dk
      case(17)
        associate(q => p%source_cs%cross_section_parms)
          if( temperature(i) <= q(1)%temperature(1) ) then
            p%level(i,1) = -1
          else if( temperature(i) >= q(size(q))%temperature(1) ) then
            p%level(i,1) = -2
          else
            do j = 2, size(q)
              if( q(j)%temperature(1) > temperature(i) ) exit
            end do
            j = j-1
            p%level(i,1) = real(j,dk)
            p%level(i,2) = (temperature(i)-q(j)%temperature(1))/(q(j+1)%temperature(1)-q(j)%temperature(1))
          end if
        end associate
      case(10)
        t = max(233.0_dk,min(temperature(i),300.0_dk))
        p%level(i,1) = 300.0_dk-t
        p%level(i,2) = 300.0_dk*t
      case(105)
        p%level(i,1) = (300.0_dk-temperature(i))/80.0_dk
      case(106)
        p%level(i,1) = air(i)/2.465e19_dk
      case(107)
        p%level(i,1) = air(i)/2.45e19_dk
      case(108)
        p%level(i,1) = air(i)*760.0_dk/2.456e19_dk
      case(109)
        p%level(i,1) = exp(-2400.0_dk/temperature(i)+3.6_dk)
      case(110)
        p%level(i,1) = .004_dk
        if( air(i) <= 2.6e19_dk ) then
          if( air(i) <= 8.e17_dk ) then
            t = .086_dk+1.613e-17_dk*8.e17_dk
          else
            t = .086_dk+1.613e-17_dk*air(i)
          end if
          p%level(i,1) = .004_dk+1/t
        end if
      case(111)
        t = 760.0_dk*air(i)/2.69e19_dk
        p%level(i,1) = min(1/(.96_dk+2.22e-3_dk*t),1.0_dk)
      case(112)
        select type(q => p%source_qy)
        type is(quantum_yield_h2so4_mills_t)
          p%level(i,1) = sqrt(16.0_dk*gas_constant*temperature(i)/(pi*q%molecular_weight_)) &
                         *(pi*q%molecular_diameter_**2*air(i)*1.e6_dk)
        end select
      case(101,102)
        q2 = exp(-825.518_dk/(.695_dk*temperature(i)))
        qfac1 = 1/(1+q2)
        qfac2 = q2/(1+q2)
        t300 = temperature(i)/300.0_dk
        p%level(i,1) = .8036_dk*qfac1
        p%level(i,2) = 8.9061_dk*t300*t300*qfac2
        p%level(i,3) = .1192_dk*t300**1.5_dk
      case(103)
        p%level(i,1) = 5.5_dk+9.2e-19_dk*air(i)
      case(104)
        select type(q => p%source_qy)
        type is(quantum_yield_ch3coch3_ch3co_ch3_t)
          t = max(q%minimum_temperature_,min(q%maximum_temperature_,temperature(i)))/295.0_dk
        end select
        p%level(i,1) = .350_dk*t**(-1.28_dk)
        p%level(i,2) = .068_dk*t**(-2.65_dk)
        p%level(i,3) = 1.600e-19_dk*t**(-2.38_dk)
        p%level(i,4) = .55e-3_dk*t**(-3.19_dk)
        p%level(i,5) = 1.62e-17_dk*t**(-10.03_dk)
        p%level(i,6) = 1.79e-3_dk*t**(-1.364_dk)
        p%level(i,7) = 26.29_dk*t**(-6.59_dk)
        p%level(i,8) = 5.72e-7_dk*t**(-2.93_dk)
        p%level(i,9) = 30006.0_dk*t**(-.064_dk)
        p%level(i,10) = 1.67e-15_dk*t**(-7.25_dk)
        p%level(i,11) = 2.08e-3_dk*t**(-1.16_dk)
      end select
    end do
    ! Table brackets and fractions depend only on height, not wavelength.
    if( p%kernel == 12 .or. p%kernel == 13 ) then
      first_table = 1
      if( p%kernel == 13 ) first_table = 2
      do table = first_table, size(p%source_cs%cross_section_parms)
        associate(q => p%source_cs%cross_section_parms(table))
          do i = 1, nz
            call bracket_value(temperature(i),q%temperature,q%deltaT,.true.,j,fraction)
            p%level(i,2*table-1) = real(j,dk)
            p%level(i,2*table) = fraction
          end do
        end associate
      end do
    else if( p%kernel == 113 .or. p%kernel == 114 ) then
      select type(q => p%source_qy)
      type is(quantum_yield_tint_t)
        do table = 1, size(q%parameters)
          do i = 1, nz
            call bracket_value(temperature(i),q%parameters(table)%temperature,q%parameters(table)%deltaT, &
                               .true.,j,fraction)
            p%level(i,2*table-1) = real(j,dk)
            p%level(i,2*table) = fraction
          end do
        end do
      type is(quantum_yield_no2_tint_t)
        do table = 1, size(q%parameters)
          do i = 1, nz
            call bracket_value(temperature(i),q%parameters(table)%temperature,q%parameters(table)%deltaT, &
                               .false.,j,fraction)
            p%level(i,2*table-1) = real(j,dk)
            p%level(i,2*table) = fraction
          end do
        end do
      end select
    end if
    do k = 1, size(p%varying)
      if( .not. p%varying(k) ) cycle
      do i = 1, nz
        associate(a => p%coefficient(k,:),z => p%level(i,:))
          select case(p%kernel)
          case(1,2)
            value = p%source_cs%cross_section_parms(1)%array(k,1)*10.0_dk**(a(1)*z(1))
          case(3)
            value = p%source_cs%cross_section_parms(1)%array(k,1)*exp(z(1)*a(1))
          case(4)
            value = p%source_cs%cross_section_parms(1)%array(k,1)* &
                      exp(p%source_cs%cross_section_parms(1)%array(k,2)*z(1))
          case(5)
            value = p%source_cs%cross_section_parms(1)%array(k,1)
            if( temperature(i) > 210 .and. temperature(i) < 300 ) value = exp(a(1)*z(1)-a(2))
          case(6)
            value = z(3)*(27.3_dk*exp(z(1)*a(1))+.932_dk*exp(z(2)*a(2)))
          case(7)
            value = (z(1)*a(1)+z(2)*a(2))*1.e-21_dk
          case(8)
            value = exp(a(1)+z(1)*a(2))
          case(9)
            value = 0
            do j = 1, size(p%source_cs%cross_section_parms(1)%array,1)
              value = value+z(j)*a(j)
            end do
            value = exp(value)
          case(14)
            associate(q => p%source_cs%cross_section_parms(1)%array)
              value = q(k,1)*(1.0_dk+z(1)*(q(k,2)+z(1)*(q(k,3)+z(1)*q(k,4))))
            end associate
          case(15)
            associate(q => p%source_cs%cross_section_parms(1)%array)
              value = q(k,1)+q(k,2)*z(1)
            end associate
          case(16)
            associate(q => p%source_cs%cross_section_parms(1)%array)
              value = q(k,1)*(1.0_dk+z(1)*(q(k,2)+z(1)*q(k,3)))
            end associate
          case(17)
            j = int(z(1))
            associate(q => p%source_cs%cross_section_parms)
              if( j == -1 ) then
                value = q(1)%array(k,1)
              else if( j == -2 ) then
                value = q(size(q))%array(k,1)
              else
                value = q(j)%array(k,1)+z(2)*(q(j+1)%array(k,1)-q(j)%array(k,1))
              end if
            end associate
          case(10)
            value = p%source_cs%cross_section_parms(1)%array(k,1)*10.0_dk**(a(1)*z(1)/z(2))
          case(11)
            if( temperature(i) < 218 ) then
              value = a(1)
            else if( temperature(i) <= 298 ) then
              value = a(1)+a(2)*(temperature(i)-218.0_dk)
            else
              value = a(3)
            end if
          case(12,13)
            value = 0
            first_table = 1
            table = size(p%source_cs%cross_section_parms)
            if( p%kernel == 13 ) then
              first_table = p%table_index(k)
              table = first_table
            end if
            do j = first_table, table
              block
                integer :: index
                index = int(z(2*j-1))
                associate(q => p%source_cs%cross_section_parms(j)%array)
                  value = value+q(k,index)+z(2*j)*(q(k,index+1)-q(k,index))
                end associate
              end block
            end do
          case(101,102)
            value = .0765_dk+z(1)*a(1)+z(2)*a(2)+z(3)*a(3)
            if( p%kernel == 102 ) value = 1-value
          case(103)
            value = min(a(1)/z(1),1.0_dk)
          case(105)
            t = a(1)*(1+a(2)*z(1))
            value = 1/(a(3)+t*air(i))
          case(106)
            value = min(1.0_dk,max(0.0_dk,a(1)/(1+a(2)*z(1))))
          case(107)
            value = min(1.0_dk,1/(1+a(1)*z(1)))
          case(108)
            value = a(1)/(a(1)+a(2)*z(1))
          case(109,110,111)
            value = z(1)
          case(112)
            value = a(1)/(a(1)+z(1))
            if( i == nz ) value = 1
          case(113,114)
            value = 0
            select type(q => p%source_qy)
            type is(quantum_yield_tint_t)
              do table = 1, size(q%parameters)
                j = int(z(2*table-1))
                value = value+q%parameters(table)%array(k,j)+z(2*table) &
                        *(q%parameters(table)%array(k,j+1)-q%parameters(table)%array(k,j))
              end do
            type is(quantum_yield_no2_tint_t)
              do table = 1, size(q%parameters)
                j = int(z(2*table-1))
                value = value+q%parameters(table)%array(k,j)+z(2*table) &
                        *(q%parameters(table)%array(k,j+1)-q%parameters(table)%array(k,j))
              end do
              value = max(value,0.0_dk)
            end select
          case(104)
            arg = z(2)*a(1)
            if( arg > 80 ) then
              c1 = 5.e34_dk
            else
              c1 = exp(arg)*z(1)/(1-z(1))
            end if
            fco = 1/(1+c1)
            if( p%wavelength(k) < 302 ) then
              c1 = z(3)*exp(-z(4)*a(3))
              fac = (1-fco)/(1+c1*air(i))
            else
              c2 = z(5)*exp(-z(6)*a(4))
              c3 = z(7)*exp(-z(8)*(a(2)-z(9))**2)
              c4 = z(10)*exp(-z(11)*a(4))
              fac = (1-fco)*(1+c3+c4*air(i))/((1+c3+c2*air(i))*(1+c4*air(i)))
            end if
            value = 0
            select type(q => p%source_qy)
            type is(quantum_yield_ch3coch3_ch3co_ch3_t)
              if( q%do_CO_ ) value = value+fco
              if( q%do_CH3CO_ ) value = value+fac
            end select
          end select
          p%values(i,p%position(k)) = value
        end associate
      end do
    end do
  end function

  subroutine bracket_value(temperature,axis,delta,clamp,index,fraction)
    real(dk), intent(in) :: temperature,axis(:),delta(:)
    logical, intent(in) :: clamp
    integer, intent(out) :: index
    real(dk), intent(out) :: fraction
    real(dk) :: t
    t = temperature
    if( clamp ) t = min(max(t,axis(1)),axis(size(axis)))
    index = 1
    do while(index < size(axis)-1)
      if( t <= axis(index+1) ) exit
      index = index+1
    end do
    fraction = (t-axis(index))/delta(index)
  end subroutine

  subroutine mask_cross(cs,w,v,supported)
    class(cross_section_t), intent(in) :: cs
    class(grid_t), intent(in) :: w
    logical, intent(out) :: v(:), supported
    real(dk), allocatable :: coeff(:)
    integer :: i,k,j,lo,hi,table
    v = .true.
    supported = .true.
    select type(cs)
    type is(cross_section_t)
      v = .false.
    type is(cross_section_rayliegh_t)
      v = .false.
    type is(cross_section_nitroxy_acetone_t)
      v = .false.
    type is(cross_section_nitroxy_ethanol_t)
      v = .false.
    type is(cross_section_hobr_oh_br_t)
      v = .false.
    type is(cross_section_t_butyl_nitrate_t)
      v = .false.
    type is(cross_section_ccl4_t)
      v = w%mid_ > 194 .and. w%mid_ < 250
      coeff = 1.0739_dk+w%mid_*(-1.6275e-2_dk+w%mid_*(8.8141e-5_dk+ &
                   w%mid_*(-1.9811e-7_dk+1.5022e-10_dk*w%mid_)))
      v = v .and. .not.(cs%cross_section_parms(1)%array(:,1) == 0 .and. abs(coeff)*85 < 300)
    type is(cross_section_chcl3_t)
      v = w%mid_ > 190 .and. w%mid_ < 240
      coeff = 3.7973_dk+w%mid_*(-7.0913e-2_dk+w%mid_*(4.9397e-4_dk+ &
                   w%mid_*(-1.5226e-6_dk+w%mid_*1.7555e-9_dk)))
      v = v .and. .not.(cs%cross_section_parms(1)%array(:,1) == 0 .and. abs(coeff)*85 < 300)
    type is(cross_section_cfc11_t)
      v = .not.(cs%cross_section_parms(1)%array(:,1) == 0 .and. &
                abs(w%mid_-184.9_dk)*1.e-4_dk*9702 < 700)
    type is(cross_section_ch3ono2_ch3o_no2_t)
      call exponential_mask(cs,v)
    type is(cross_section_hno3_oh_no2_t)
      call exponential_mask(cs,v)
    type is(cross_section_rono2_t)
      call exponential_mask(cs,v)
    type is(cross_section_chbr3_t)
      v = w%mid_ > 290 .and. w%mid_ < 340
    type is(cross_section_h2o2_oh_oh_t)
      v = w%edge_(:w%ncells_) >= 260 .and. w%edge_(:w%ncells_) < 350
    type is(cross_section_hcfc_t)
      v = w%mid_ >= 190 .and. w%mid_ <= cs%cross_section_parms(1)%temperature(2)
    type is(cross_section_n2o_n2_o1d_t)
      v = w%mid_ >= 173 .and. w%mid_ <= 240
    type is(cross_section_n2o5_no2_no3_t)
      coeff = 1000*cs%cross_section_parms(2)%array(:,1)
      v = .not.(coeff == 0 .or. (cs%cross_section_parms(1)%array(:,1) == 0 .and. abs(coeff)/1043 < 300))
    type is(cross_section_no2_tint_t)
      call table_mask(cs,v)
    type is(cross_section_tint_t)
      call table_mask(cs,v)
    type is(cross_section_o3_jpl06_t)
      v = w%edge_(:w%ncells_) >= 196.078_dk .and. w%edge_(2:) <= 342.5_dk
    type is(cross_section_o3_tint_t)
      v = .false.
      do i = 1, w%ncells_
        if( w%edge_(i) < cs%v185(1) ) then
          table = 3
        else if( w%edge_(i) < cs%v195(1) ) then
          table = 4
        else if( w%edge_(i) < cs%v345(1) ) then
          table = 2
        else
          cycle
        end if
        v(i) = any(cs%cross_section_parms(table)%array(i,:) /= cs%cross_section_parms(table)%array(i,1))
      end do
    type is(cross_section_temperature_based_t)
      v = .false.
      associate(p => cs%parameterization_)
        lo = p%min_wavelength_index_
        hi = p%max_wavelength_index_
        do i = 1, w%ncells_
          do k = max(1,lo-1), min(size(cs%raw_wavelengths_)-1,hi)
            if( cs%raw_wavelengths_(k+1) > w%edge_(i) .and. cs%raw_wavelengths_(k) < w%edge_(i+1) ) v(i)=.true.
          end do
          ! The extrapolation segments also touch parameterized end points.
          if( lo == 1 .and. w%edge_(i) < cs%raw_wavelengths_(1) .and. &
              w%edge_(i+1) > (1.0_dk-real(1.e-5,dk))*cs%raw_wavelengths_(1) ) v(i)=.true.
          if( hi == size(cs%raw_wavelengths_) .and. w%edge_(i+1) > cs%raw_wavelengths_(hi) .and. &
              w%edge_(i) < (1.0_dk+real(1.e-5,dk))*cs%raw_wavelengths_(hi) ) v(i)=.true.
        end do
      end associate
    type is(cross_section_ch3coch3_ch3co_ch3_t)
    type is(cross_section_ch2o_t)
    type is(cross_section_clono2_t)
    type is(cross_section_cl2_cl_cl_t)
    type is(cross_section_oclo_t)
    class default
      supported = .false.
    end select
  end subroutine
  subroutine exponential_mask(cs,v)
    class(cross_section_t), intent(in) :: cs
    logical, intent(out) :: v(:)
    v = .not.(cs%cross_section_parms(1)%array(:,2) == 0 .or. &
        (cs%cross_section_parms(1)%array(:,1) == 0 .and. abs(cs%cross_section_parms(1)%array(:,2))*9702 < 700))
  end subroutine
  subroutine table_mask(cs,v)
    class(cross_section_t), intent(in) :: cs
    logical, intent(out) :: v(:)
    integer :: i,k
    v = .false.
    do k = 1, size(cs%cross_section_parms)
      do i = 1, size(v)
        v(i) = v(i) .or. any(cs%cross_section_parms(k)%array(i,:) /= cs%cross_section_parms(k)%array(i,1))
      end do
    end do
  end subroutine
  subroutine mask_yield(q,w,v,supported)
    class(quantum_yield_t), intent(in) :: q
    class(grid_t), intent(in) :: w
    logical, intent(out) :: v(:), supported
    integer :: i,k
    real(dk), allocatable :: a(:),b(:)
    v = .true.
    supported = .true.
    select type(q)
    type is(quantum_yield_t)
      v = .false.
    type is(quantum_yield_taylor_series_t)
      v = .false.
    type is(quantum_yield_ho2_oh_o_t)
      v = .false.
    type is(quantum_yield_clo_cl_o1d_t)
      v = .false.
    type is(quantum_yield_clo_cl_o3p_t)
      v = .false.
    type is(quantum_yield_clono2_cl_no3_t)
      v = .false.
    type is(quantum_yield_clono2_clo_no2_t)
      v = .false.
    type is(quantum_yield_tint_t)
      v = .false.
      do k = 1, size(q%parameters)
        do i = 1, size(v)
          v(i) = v(i) .or. any(q%parameters(k)%array(i,:) /= q%parameters(k)%array(i,1))
        end do
      end do
    type is(quantum_yield_no2_tint_t)
      v = .false.
      do k = 1, size(q%parameters)
        do i = 1, size(v)
          v(i) = v(i) .or. any(q%parameters(k)%array(i,:) /= q%parameters(k)%array(i,1))
        end do
      end do
    type is(quantum_yield_o3_o2_o1d_t)
      v = w%mid_ > 305 .and. w%mid_ <= 328
    type is(quantum_yield_o3_o2_o3p_t)
      v = w%mid_ > 305 .and. w%mid_ <= 328
    type is(quantum_yield_ch2o_h2_co_t)
      v = w%mid_ >= 330 .and. w%mid_ < 360 .and. q%quantum_yield_parms(1)%array(:,2) > 0
    type is(quantum_yield_ch3cho_ch3_hco_t)
      a = q%quantum_yield_parms(1)%array(:,2)
      b = 0*a
      where(a > 0) b = (1-q%quantum_yield_parms(1)%array(:,1))/a-1
      v = b /= 0
    type is(quantum_yield_c2h5cho_c2h5_hco_t)
      v = q%quantum_yield_parms(1)%array(:,1) >= 10.0_dk/1.e36_dk
      where(q%quantum_yield_parms(1)%array(:,1) == 1) v = .false.
    type is(quantum_yield_ch3coch3_ch3co_ch3_t)
      v = w%mid_ >= 279 .and. w%mid_ <= 327
    type is(quantum_yield_ch3cocho_ch3co_hco_t)
      v = w%mid_ >= 380 .and. w%mid_ < 440
    type is(quantum_yield_h2so4_mills_t)
      v = .false.
      v(q%wavelength_indices_) = .true.
    type is(quantum_yield_no3m_aq_t)
    type is(quantum_yield_ch2chcho_t)
    type is(quantum_yield_ch3coch2ch3_t)
    type is(quantum_yield_mvk_t)
    class default
      supported = .false.
    end select
  end subroutine

  subroutine prepare_rate(r,cs,qy,generation,oxygen,energy,oxygen_bands)
    type(rate_plan_t), intent(inout) :: r
    type(spectrum_plan_t), pointer, intent(in) :: cs,qy
    integer, intent(in) :: generation
    logical, intent(in) :: oxygen
    real(dk), optional, intent(in) :: energy(:)
    logical, optional, intent(in) :: oxygen_bands(:)
    integer :: i,n
    logical, allocatable :: keep(:)
    r%cs => cs
    r%qy => qy
    r%generation = generation
    r%oxygen = oxygen
    r%heating = present(energy)
    n = size(cs%constant)
    r%constant_term = .not.(cs%varying .or. qy%varying .or. oxygen)
    if( present(oxygen_bands) .and. oxygen ) &
      r%constant_term = .not.(cs%varying .or. qy%varying .or. oxygen_bands)
    r%product = cs%constant*qy%constant
    if( present(energy) ) then
      r%energy = energy
      r%product = energy*r%product
    end if
    keep = .not.(r%constant_term .and. r%product == 0)
    r%indices = pack([(i,i=1,n)],keep)
  end subroutine
  subroutine integrate_rate(r,flux,scale,result,corrected)
    type(rate_plan_t), intent(in) :: r
    real(dk), intent(in) :: flux(:,:),scale
    real(dk), intent(out) :: result(:)
    real(dk), optional, intent(in) :: corrected(:,:)
    real(dk) :: term,constant
    integer :: t,w,z
    result = 0
    do t = 1, size(r%indices)
      w = r%indices(t)
      if( r%constant_term(w) ) then
        constant = r%product(w)
        do z = 1, size(result)
          result(z) = result(z)+flux(z,w)*constant
        end do
      else if( present(corrected) ) then
        if( .not. r%qy%varying(w) ) then
          constant = r%qy%constant(w)
          do z = 1, size(result)
            term = corrected(z,w)*constant
            if( r%heating ) term = r%energy(w)*term
            result(z) = result(z)+flux(z,w)*term
          end do
        else
          do z = 1, size(result)
            term = corrected(z,w)*r%qy%values(z,r%qy%position(w))
            if( r%heating ) term = r%energy(w)*term
            result(z) = result(z)+flux(z,w)*term
          end do
        end if
      else if( .not. r%cs%varying(w) ) then
        constant = r%cs%constant(w)
        do z = 1, size(result)
          term = constant*r%qy%values(z,r%qy%position(w))
          if( r%heating ) term = r%energy(w)*term
          result(z) = result(z)+flux(z,w)*term
        end do
      else if( .not. r%qy%varying(w) ) then
        constant = r%qy%constant(w)
        do z = 1, size(result)
          term = r%cs%values(z,r%cs%position(w))*constant
          if( r%heating ) term = r%energy(w)*term
          result(z) = result(z)+flux(z,w)*term
        end do
      else
        do z = 1, size(result)
          term = r%cs%values(z,r%cs%position(w))*r%qy%values(z,r%qy%position(w))
          if( r%heating ) term = r%energy(w)*term
          result(z) = result(z)+flux(z,w)*term
        end do
      end if
    end do
    result = result*scale
  end subroutine
end module tuvx_prepared_spectra
