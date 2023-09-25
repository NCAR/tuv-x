! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_factory
  ! The quantum_yield_factory module
  ! Builder of quantum yield calculators

  use tuvx_quantum_yield,              only : quantum_yield_t
  use tuvx_quantum_yield_tint,         only : quantum_yield_tint_t
  use tuvx_quantum_yield_no2_tint,     only : quantum_yield_no2_tint_t
  use tuvx_quantum_yield_o3_o2_o1d,    only : quantum_yield_o3_o2_o1d_t
  use tuvx_quantum_yield_o3_o2_o3p,    only : quantum_yield_o3_o2_o3p_t
  use tuvx_quantum_yield_ho2_oh_o,     only : quantum_yield_ho2_oh_o_t
  use tuvx_quantum_yield_no3m_aq,      only : quantum_yield_no3m_aq_t
  use tuvx_quantum_yield_ch2o_h2_co,   only : quantum_yield_ch2o_h2_co_t
  use tuvx_quantum_yield_ch3cho_ch3_hco,                                      &
    only : quantum_yield_ch3cho_ch3_hco_t
  use tuvx_quantum_yield_c2h5cho_c2h5_hco,                                    &
    only : quantum_yield_c2h5cho_c2h5_hco_t
  use tuvx_quantum_yield_ch2chcho,     only : quantum_yield_ch2chcho_t
  use tuvx_quantum_yield_mvk,          only : quantum_yield_mvk_t
  use tuvx_quantum_yield_ch3coch3_ch3co_ch3,                                  &
    only : quantum_yield_ch3coch3_ch3co_ch3_t
  use tuvx_quantum_yield_ch3coch2ch3,  only : quantum_yield_ch3coch2ch3_t
  use tuvx_quantum_yield_ch3cocho_ch3co_hco,                                  &
    only : quantum_yield_ch3cocho_ch3co_hco_t
  use tuvx_quantum_yield_clo_cl_o1d,   only : quantum_yield_clo_cl_o1d_t
  use tuvx_quantum_yield_clo_cl_o3p,   only : quantum_yield_clo_cl_o3p_t
  use tuvx_quantum_yield_clono2_cl_no3,only : quantum_yield_clono2_cl_no3_t
  use tuvx_quantum_yield_clono2_clo_no2,                                      &
    only : quantum_yield_clono2_clo_no2_t

  implicit none

  private
  public :: quantum_yield_builder, quantum_yield_type_name,                   &
            quantum_yield_allocate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function quantum_yield_builder( config, grid_warehouse, profile_warehouse ) &
      result( quantum_yield )
    ! Builds quantum yield calculators based off of the configuration

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_t), pointer :: quantum_yield ! New :f:type:`~tuvx_quantum_yield/quantum_yield_t` calculator
    type(config_t),            intent(inout) :: config ! Quantum yield configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    type(string_t) :: quantum_yield_type
    character(len=*), parameter :: Iam = 'quantum yield builder'

    quantum_yield => null()
    call config%get( 'type', quantum_yield_type, Iam )

    select case( quantum_yield_type%to_char() )
      case( 'base' )
        quantum_yield => quantum_yield_t( config, grid_warehouse,             &
                                          profile_warehouse )
      case( 'tint' )
        quantum_yield => quantum_yield_tint_t( config, grid_warehouse,        &
                                               profile_warehouse )
      case( 'NO2 tint' )
        quantum_yield => quantum_yield_no2_tint_t( config, grid_warehouse,    &
                                                   profile_warehouse )
      case( 'O3+hv->O2+O(1D)' )
        quantum_yield => quantum_yield_o3_o2_o1d_t( config, grid_warehouse,   &
                                                    profile_warehouse )
      case( 'O3+hv->O2+O(3P)' )
        quantum_yield => quantum_yield_o3_o2_o3p_t( config, grid_warehouse,   &
                                                    profile_warehouse )
      case( 'HO2' )
        quantum_yield => quantum_yield_ho2_oh_o_t( config, grid_warehouse,    &
                                                   profile_warehouse )
      case( 'NO3-_(aq)+hv->NO2(aq)+O-' )
        quantum_yield => quantum_yield_no3m_aq_t( config, grid_warehouse,     &
                                                  profile_warehouse )
      case( 'CH2O' )
        quantum_yield => quantum_yield_ch2o_h2_co_t( config, grid_warehouse,  &
                                                     profile_warehouse )
      case( 'CH3CHO+hv->CH3+HCO' )
        quantum_yield => quantum_yield_ch3cho_ch3_hco_t( config,              &
                                                         grid_warehouse,      &
                                                         profile_warehouse )
      case( 'C2H5CHO' )
        quantum_yield => quantum_yield_c2h5cho_c2h5_hco_t( config,            &
                                                           grid_warehouse,    &
                                                           profile_warehouse )
      case( 'CH2CHCHO+hv->Products' )
        quantum_yield => quantum_yield_ch2chcho_t( config, grid_warehouse,    &
                                                   profile_warehouse )
      case( 'MVK+hv->Products' )
        quantum_yield => quantum_yield_mvk_t( config, grid_warehouse,         &
                                              profile_warehouse )
      case( 'CH3COCH3+hv->CH3CO+CH3' )
        quantum_yield => quantum_yield_ch3coch3_ch3co_ch3_t( config,          &
                                                            grid_warehouse,   &
                                                            profile_warehouse )
      case( 'CH3COCH2CH3+hv->CH3CO+CH2CH3' )
        quantum_yield => quantum_yield_ch3coch2ch3_t( config, grid_warehouse, &
                                                      profile_warehouse )
      case( 'CH3COCHO+hv->CH3CO+HCO' )
        quantum_yield => quantum_yield_ch3cocho_ch3co_hco_t( config,          &
                                                            grid_warehouse,   &
                                                            profile_warehouse )
      case( 'ClO+hv->Cl+O(1D)' )
        quantum_yield => quantum_yield_clo_cl_o1d_t( config, grid_warehouse,  &
                                                     profile_warehouse )
      case( 'ClO+hv->Cl+O(3P)' )
        quantum_yield => quantum_yield_clo_cl_o3p_t( config, grid_warehouse,  &
                                                     profile_warehouse )
      case( 'ClONO2+hv->Cl+NO3' )
        quantum_yield => quantum_yield_clono2_cl_no3_t( config,               &
                                                        grid_warehouse,       &
                                                        profile_warehouse )
      case( 'ClONO2+hv->ClO+NO2' )
        quantum_yield => quantum_yield_clono2_clo_no2_t( config,              &
                                                         grid_warehouse,      &
                                                         profile_warehouse )
      case default
        call die_msg( 450768214, "Invalid quantum yield type: '"//            &
                                 quantum_yield_type%to_char( )//"'" )
    end select

  end function quantum_yield_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function quantum_yield_type_name( quantum_yield )            &
      result( name )
    ! Returns the type of a quantum yield as a string

    use musica_assert,                 only : die
    use musica_string,                 only : string_t

    class(quantum_yield_t), intent(in) :: quantum_yield ! quantum yield to return type for

    select type( quantum_yield )
      type is( quantum_yield_t )
        name = "quantum_yield_t"
      type is( quantum_yield_tint_t )
        name = "quantum_yield_tint_t"
      type is( quantum_yield_no2_tint_t )
        name = "quantum_yield_no2_tint_t"
      type is( quantum_yield_o3_o2_o1d_t )
        name = "quantum_yield_o3_o2_o1d_t"
      type is( quantum_yield_o3_o2_o3p_t )
        name = "quantum_yield_o3_o2_o3p_t"
      type is( quantum_yield_ho2_oh_o_t )
        name = "quantum_yield_ho2_oh_o_t"
      type is( quantum_yield_no3m_aq_t )
        name = "quantum_yield_no3m_aq_t"
      type is( quantum_yield_ch2o_h2_co_t )
        name = "quantum_yield_ch2o_h2_co_t"
      type is( quantum_yield_ch3cho_ch3_hco_t )
        name = "quantum_yield_ch3cho_ch3_hco_t"
      type is( quantum_yield_c2h5cho_c2h5_hco_t )
        name = "quantum_yield_c2h5cho_c2h5_hco_t"
      type is( quantum_yield_ch2chcho_t )
        name = "quantum_yield_ch2chcho_t"
      type is( quantum_yield_mvk_t )
        name = "quantum_yield_mvk_t"
      type is( quantum_yield_ch3coch3_ch3co_ch3_t )
        name = "quantum_yield_ch3coch3_ch3co_ch3_t"
      type is( quantum_yield_ch3coch2ch3_t )
        name = "quantum_yield_ch3coch2ch3_t"
      type is( quantum_yield_ch3cocho_ch3co_hco_t )
        name = "quantum_yield_ch3cocho_ch3co_hco_t"
      type is( quantum_yield_clo_cl_o1d_t )
        name = "quantum_yield_clo_cl_o1d_t"
      type is( quantum_yield_clo_cl_o3p_t )
        name = "quantum_yield_clo_cl_o3p_t"
      type is( quantum_yield_clono2_cl_no3_t )
        name = "quantum_yield_clono2_cl_no3_t"
      type is( quantum_yield_clono2_clo_no2_t )
        name = "quantum_yield_clono2_clo_no2_t"
    class default
      call die( 853572483 )
    end select

  end function quantum_yield_type_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function quantum_yield_allocate( type_name ) result( quantum_yield )
    ! Allocates a quantum yield pointer as the base or a subclass by type

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t

    type(string_t),         intent(in) :: type_name     ! Name of the type to allocate
    class(quantum_yield_t), pointer    :: quantum_yield ! Allocated quantum yield

    select case( type_name%to_char( ) )
      case( 'quantum_yield_t' )
        allocate( quantum_yield_t :: quantum_yield )
      case( 'quantum_yield_tint_t' )
        allocate( quantum_yield_tint_t :: quantum_yield )
      case( 'quantum_yield_no2_tint_t' )
        allocate( quantum_yield_no2_tint_t :: quantum_yield )
      case( 'quantum_yield_o3_o2_o1d_t' )
        allocate( quantum_yield_o3_o2_o1d_t :: quantum_yield )
      case( 'quantum_yield_o3_o2_o3p_t' )
        allocate( quantum_yield_o3_o2_o3p_t :: quantum_yield )
      case( 'quantum_yield_ho2_oh_o_t' )
        allocate( quantum_yield_ho2_oh_o_t :: quantum_yield )
      case( 'quantum_yield_no3m_aq_t' )
        allocate( quantum_yield_no3m_aq_t :: quantum_yield )
      case( 'quantum_yield_ch2o_h2_co_t' )
        allocate( quantum_yield_ch2o_h2_co_t :: quantum_yield )
      case( 'quantum_yield_ch3cho_ch3_hco_t' )
        allocate( quantum_yield_ch3cho_ch3_hco_t :: quantum_yield )
      case( 'quantum_yield_c2h5cho_c2h5_hco_t' )
        allocate( quantum_yield_c2h5cho_c2h5_hco_t :: quantum_yield )
      case( 'quantum_yield_ch2chcho_t' )
        allocate( quantum_yield_ch2chcho_t :: quantum_yield )
      case( 'quantum_yield_mvk_t' )
        allocate( quantum_yield_mvk_t :: quantum_yield )
      case( 'quantum_yield_ch3coch3_ch3co_ch3_t' )
        allocate( quantum_yield_ch3coch3_ch3co_ch3_t :: quantum_yield )
      case( 'quantum_yield_ch3coch2ch3_t' )
        allocate( quantum_yield_ch3coch2ch3_t :: quantum_yield )
      case( 'quantum_yield_ch3cocho_ch3co_hco_t' )
        allocate( quantum_yield_ch3cocho_ch3co_hco_t :: quantum_yield )
      case( 'quantum_yield_clo_cl_o1d_t' )
        allocate( quantum_yield_clo_cl_o1d_t :: quantum_yield )
      case( 'quantum_yield_clo_cl_o3p_t' )
        allocate( quantum_yield_clo_cl_o3p_t :: quantum_yield )
      case( 'quantum_yield_clono2_cl_no3_t' )
        allocate( quantum_yield_clono2_cl_no3_t :: quantum_yield )
      case( 'quantum_yield_clono2_clo_no2_t' )
        allocate( quantum_yield_clono2_clo_no2_t :: quantum_yield )
    case default
      call die_msg( 894617177, "Invalid quantum yield type: '"//type_name//"'" )
    end select

  end function quantum_yield_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_factory
