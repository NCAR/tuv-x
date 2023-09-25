! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_factory
! Builder of cross section calculators for use by 
! :f:type:`~tuvx_cross_section_warehouse/cross_section_warehouse_t`.

  use tuvx_cross_section,              only : cross_section_t
  use tuvx_cross_section_ccl4,         only : cross_section_ccl4_t
  use tuvx_cross_section_cfc11,        only : cross_section_cfc11_t
  use tuvx_cross_section_ch3coch3_ch3co_ch3,                                  &
    only : cross_section_ch3coch3_ch3co_ch3_t
  use tuvx_cross_section_chbr3,        only : cross_section_chbr3_t
  use tuvx_cross_section_chcl3,        only : cross_section_chcl3_t
  use tuvx_cross_section_ch2o,         only : cross_section_ch2o_t
  use tuvx_cross_section_ch3ono2_ch3o_no2,                                    &
    only : cross_section_ch3ono2_ch3o_no2_t
  use tuvx_cross_section_cl2_cl_cl,    only : cross_section_cl2_cl_cl_t
  use tuvx_cross_section_clono2,       only : cross_section_clono2_t
  use tuvx_cross_section_h2o2_oh_oh,   only : cross_section_h2o2_oh_oh_t
  use tuvx_cross_section_hno3_oh_no2,  only : cross_section_hno3_oh_no2_t
  use tuvx_cross_section_hcfc,         only : cross_section_hcfc_t
  use tuvx_cross_section_hobr_oh_br,   only : cross_section_hobr_oh_br_t
  use tuvx_cross_section_n2o_n2_o1d,   only : cross_section_n2o_n2_o1d_t
  use tuvx_cross_section_n2o5_no2_no3, only : cross_section_n2o5_no2_no3_t
  use tuvx_cross_section_nitroxy_acetone,                                    &
    only : cross_section_nitroxy_acetone_t
  use tuvx_cross_section_nitroxy_ethanol,                                    &
    only : cross_section_nitroxy_ethanol_t
  use tuvx_cross_section_no2_tint,     only : cross_section_no2_tint_t
  use tuvx_cross_section_o3_jpl06,     only : cross_section_o3_jpl06_t
  use tuvx_cross_section_o3_tint,      only : cross_section_o3_tint_t
  use tuvx_cross_section_oclo,         only : cross_section_oclo_t
  use tuvx_cross_section_rayliegh,     only : cross_section_rayliegh_t
  use tuvx_cross_section_rono2,        only : cross_section_rono2_t
  use tuvx_cross_section_temperature_based,                                  &
    only : cross_section_temperature_based_t
  use tuvx_cross_section_t_butyl_nitrate,                                    &
    only : cross_section_t_butyl_nitrate_t
  use tuvx_cross_section_tint,         only : cross_section_tint_t

  implicit none

  private
  public :: cross_section_builder, cross_section_type_name,                   &
            cross_section_allocate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function cross_section_builder( config, grid_warehouse, profile_warehouse ) &
      result( new_cross_section )
      ! Build cross sections from a configuration

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(cross_section_t),    pointer       :: new_cross_section ! New :f:type:`~tuvx_cross_section/cross_section_t`
    type(config_t),            intent(inout) :: config ! Cross section configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    type(string_t) :: cross_section_type
    character(len=*), parameter :: Iam = 'cross section builder'

    new_cross_section => null( )
    call config%get( 'type', cross_section_type, Iam )

    select case( cross_section_type%to_char() )
      case( 'air' )
        new_cross_section => cross_section_rayliegh_t( config, grid_warehouse,&
                                                           profile_warehouse )
      case( 'base' )
        new_cross_section => cross_section_t( config, grid_warehouse,         &
                                                           profile_warehouse )
      case( 'CCl4+hv->Products' )
        new_cross_section => cross_section_ccl4_t( config, grid_warehouse,    &
                                                           profile_warehouse )
      case( 'CCl3F+hv->Products' )
        new_cross_section => cross_section_cfc11_t( config, grid_warehouse,   &
                                                           profile_warehouse )
      case( 'CHCl3+hv->Products' )
        new_cross_section => cross_section_chcl3_t( config, grid_warehouse,   &
                                                           profile_warehouse )
      case( 'CH2(OH)CH2(ONO2)+hv->CH2(OH)CH2(O.)+NO2' )
        new_cross_section => cross_section_nitroxy_ethanol_t( config,         &
                                           grid_warehouse, profile_warehouse )
      case( 'CH2O' )
        new_cross_section => cross_section_ch2o_t( config, grid_warehouse,    &
                                                           profile_warehouse )
      case( 'CH3COCH2(ONO2)+hv->CH3COCH2(O.)+NO2' )
        new_cross_section => cross_section_nitroxy_acetone_t( config,         &
                                           grid_warehouse, profile_warehouse )
      case( 'CH3COCH3+hv->CH3CO+CH3' )
        new_cross_section => cross_section_ch3coch3_ch3co_ch3_t( config,      &
                                           grid_warehouse, profile_warehouse )
      case( 'CH3ONO2+hv->CH3O+NO2' )
        new_cross_section => cross_section_ch3ono2_ch3o_no2_t( config,        &
                                           grid_warehouse, profile_warehouse )
      case( 'CHBr3+hv->Products' )
        new_cross_section => cross_section_chbr3_t( config, grid_warehouse,   &
                                                           profile_warehouse )
      case( 'Cl2+hv->Cl+Cl' )
        new_cross_section => cross_section_cl2_cl_cl_t( config,               &
                                           grid_warehouse, profile_warehouse )
      case( 'ClONO2' )
        new_cross_section => cross_section_clono2_t( config, grid_warehouse,  &
                                                           profile_warehouse )
      case( 'H2O2+hv->OH+OH' )
        new_cross_section => cross_section_h2o2_oh_oh_t( config,              &
                                           grid_warehouse, profile_warehouse )
      case( 'HCFC+hv->Products' )
        new_cross_section => cross_section_hcfc_t( config, grid_warehouse,    &
                                                           profile_warehouse )
      case( 'HNO3+hv->OH+NO2' )
        new_cross_section => cross_section_hno3_oh_no2_t( config,             &
                                           grid_warehouse, profile_warehouse )
      case( 'HOBr+hv->OH+Br' )
        new_cross_section => cross_section_hobr_oh_br_t( config,              &
                                           grid_warehouse, profile_warehouse )
      case( 'N2O+hv->N2+O(1D)' )
        new_cross_section => cross_section_n2o_n2_o1d_t( config,              &
                                           grid_warehouse, profile_warehouse )
      case( 'N2O5+hv->NO2+NO3' )
        new_cross_section => cross_section_n2o5_no2_no3_t( config,            &
                                           grid_warehouse, profile_warehouse )
      case( 'NO2 tint' )
        new_cross_section => cross_section_no2_tint_t( config, grid_warehouse,&
                                                           profile_warehouse )
      case( 'O3 JPL06' )
        new_cross_section => cross_section_o3_jpl06_t( config, grid_warehouse,&
                                                           profile_warehouse )
      case( 'O3' )
        new_cross_section => cross_section_o3_tint_t( config, grid_warehouse, &
                                                           profile_warehouse )
      case( 'OClO+hv->Products' )
        new_cross_section => cross_section_oclo_t( config, grid_warehouse,    &
                                                           profile_warehouse )
      case( 'RONO2' )
        new_cross_section => cross_section_rono2_t( config, grid_warehouse,   &
                                                           profile_warehouse )
      case( 'SO2' )
        new_cross_section => cross_section_t( config, grid_warehouse,         &
                                                           profile_warehouse )
      case( 'temperature based' )
        new_cross_section => cross_section_temperature_based_t(               &
                                   config, grid_warehouse, profile_warehouse )
      case( 't_butyl_nitrate+hv->Products' )
        new_cross_section => cross_section_t_butyl_nitrate_t ( config,        &
                                           grid_warehouse, profile_warehouse )
      case( 'tint' )
        new_cross_section => cross_section_tint_t( config, grid_warehouse,    &
                                                           profile_warehouse )
      case default
        call die_msg( 450768214, "Invalid cross section type: '"//            &
                                 cross_section_type%to_char( )//"'" )
    end select

  end function cross_section_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function cross_section_type_name( cross_section )            &
      result( name )
    ! Returns the type of a cross section as a string

    use musica_assert,                 only : die
    use musica_string,                 only : string_t

    class(cross_section_t), intent(in) :: cross_section ! Cross section to return type for

    select type( cross_section )
      type is( cross_section_t )
        name = "cross_section_t"
      type is( cross_section_rayliegh_t )
        name = "cross_section_rayliegh_t"
      type is( cross_section_ccl4_t )
        name = "cross_section_ccl4_t"
      type is( cross_section_cfc11_t )
        name = "cross_section_cfc11_t"
      type is( cross_section_chcl3_t )
        name = "cross_section_chcl3_t"
      type is( cross_section_nitroxy_ethanol_t )
        name = "cross_section_nitroxy_ethanol_t"
      type is( cross_section_ch2o_t )
        name = "cross_section_ch2o_t"
      type is( cross_section_nitroxy_acetone_t )
        name = "cross_section_nitroxy_acetone_t"
      type is( cross_section_ch3coch3_ch3co_ch3_t )
        name = "cross_section_ch3coch3_ch3co_ch3_t"
      type is( cross_section_ch3ono2_ch3o_no2_t )
        name = "cross_section_ch3ono2_ch3o_no2_t"
      type is( cross_section_chbr3_t )
        name = "cross_section_chbr3_t"
      type is( cross_section_cl2_cl_cl_t )
        name = "cross_section_cl2_cl_cl_t"
      type is( cross_section_clono2_t )
        name = "cross_section_clono2_t"
      type is( cross_section_h2o2_oh_oh_t )
        name = "cross_section_h2o2_oh_oh_t"
      type is( cross_section_hcfc_t )
        name = "cross_section_hcfc_t"
      type is( cross_section_hno3_oh_no2_t )
        name = "cross_section_hno3_oh_no2_t"
      type is( cross_section_hobr_oh_br_t )
        name = "cross_section_hobr_oh_br_t"
      type is( cross_section_n2o_n2_o1d_t )
        name = "cross_section_n2o_n2_o1d_t"
      type is( cross_section_n2o5_no2_no3_t )
        name = "cross_section_n2o5_no2_no3_t"
      type is( cross_section_no2_tint_t )
        name = "cross_section_no2_tint_t"
      type is( cross_section_o3_jpl06_t )
        name = "cross_section_o3_jpl06_t"
      type is( cross_section_o3_tint_t )
        name = "cross_section_o3_tint_t"
      type is( cross_section_oclo_t )
        name = "cross_section_oclo_t"
      type is( cross_section_rono2_t )
        name = "cross_section_rono2_t"
      type is( cross_section_temperature_based_t )
        name = "cross_section_temperature_based_t"
      type is( cross_section_t_butyl_nitrate_t )
        name = "cross_section_t_butyl_nitrate_t"
      type is( cross_section_tint_t )
        name = "cross_section_tint_t"
      class default
        call die( 692607640 )
    end select

  end function cross_section_type_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function cross_section_allocate( type_name ) result( cross_section )
    ! Allocates a cross section pointer as the base or a subclass by type

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t

    type(string_t),         intent(in) :: type_name     ! Name of the type to allocate
    class(cross_section_t), pointer    :: cross_section ! Allocated cross section

    select case( type_name%to_char( ) )
      case( 'cross_section_t' )
        allocate( cross_section_t :: cross_section )
      case( 'cross_section_rayliegh_t' )
        allocate( cross_section_rayliegh_t :: cross_section )
      case( 'cross_section_ccl4_t' )
        allocate( cross_section_ccl4_t :: cross_section )
      case( 'cross_section_cfc11_t' )
        allocate( cross_section_cfc11_t :: cross_section )
      case( 'cross_section_chcl3_t' )
        allocate( cross_section_chcl3_t :: cross_section )
      case( 'cross_section_nitroxy_ethanol_t' )
        allocate( cross_section_nitroxy_ethanol_t :: cross_section )
      case( 'cross_section_ch2o_t' )
        allocate( cross_section_ch2o_t :: cross_section )
      case( 'cross_section_nitroxy_acetone_t' )
        allocate( cross_section_nitroxy_acetone_t :: cross_section )
      case( 'cross_section_ch3coch3_ch3co_ch3_t' )
        allocate( cross_section_ch3coch3_ch3co_ch3_t :: cross_section )
      case( 'cross_section_ch3ono2_ch3o_no2_t' )
        allocate( cross_section_ch3ono2_ch3o_no2_t :: cross_section )
      case( 'cross_section_chbr3_t' )
        allocate( cross_section_chbr3_t :: cross_section )
      case( 'cross_section_cl2_cl_cl_t' )
        allocate( cross_section_cl2_cl_cl_t :: cross_section )
      case( 'cross_section_clono2_t' )
        allocate( cross_section_clono2_t :: cross_section )
      case( 'cross_section_h2o2_oh_oh_t' )
        allocate( cross_section_h2o2_oh_oh_t :: cross_section )
      case( 'cross_section_hcfc_t' )
        allocate( cross_section_hcfc_t :: cross_section )
      case( 'cross_section_hno3_oh_no2_t' )
        allocate( cross_section_hno3_oh_no2_t :: cross_section )
      case( 'cross_section_hobr_oh_br_t' )
        allocate( cross_section_hobr_oh_br_t :: cross_section )
      case( 'cross_section_n2o_n2_o1d_t' )
        allocate( cross_section_n2o_n2_o1d_t :: cross_section )
      case( 'cross_section_n2o5_no2_no3_t' )
        allocate( cross_section_n2o5_no2_no3_t :: cross_section )
      case( 'cross_section_no2_tint_t' )
        allocate( cross_section_no2_tint_t :: cross_section )
      case( 'cross_section_o3_jpl06_t' )
        allocate( cross_section_o3_jpl06_t :: cross_section )
      case( 'cross_section_o3_tint_t' )
        allocate( cross_section_o3_tint_t :: cross_section )
      case( 'cross_section_oclo_t' )
        allocate( cross_section_oclo_t :: cross_section )
      case( 'cross_section_rono2_t' )
        allocate( cross_section_rono2_t :: cross_section )
      case( 'cross_section_temperature_based_t' )
        allocate( cross_section_temperature_based_t :: cross_section )
      case( 'cross_section_t_butyl_nitrate_t' )
        allocate( cross_section_t_butyl_nitrate_t :: cross_section )
      case( 'cross_section_tint_t' )
        allocate( cross_section_tint_t :: cross_section )
      case default
        call die_msg( 828010698, "Invalid cross section type: '"//            &
                                 type_name//"'" )
    end select

  end function cross_section_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_factory
