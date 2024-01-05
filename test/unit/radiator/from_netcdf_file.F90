! Copyright (C) 2023 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests the radiator_from_netcdf_file_t type
program test_radiator_from_netcdf_file

  use musica_constants,                only : dk => musica_dk

  implicit none

  call test_radiator_from_netcdf_file_t()

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine test_radiator_from_netcdf_file_t( )

    use musica_assert,                 only : assert
    use musica_config,                 only : config_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t
    use tuvx_radiator,                 only : radiator_t
    use tuvx_radiator_warehouse,       only : radiator_warehouse_t
    use tuvx_cross_section_warehouse,  only : cross_section_warehouse_t

    character(len=*), parameter :: Iam = "radiator_from_netcdf_file_t test"
    type(config_t) :: config, grid_config, radiator_config
    type(grid_warehouse_t), pointer :: grids
    type(profile_warehouse_t), pointer :: profiles
    type(cross_section_warehouse_t), pointer :: cross_sections
    type(radiator_warehouse_t), pointer :: radiators
    class(radiator_t), pointer :: radiator

    call config%from_file( "test/data/radiators.from_netcdf_file.config.json" )
    call config%get( "grids", grid_config, Iam )
    call config%get( "radiators", radiator_config, Iam )

    grids => grid_warehouse_t( grid_config )
    allocate( profiles )
    allocate( cross_sections )
    radiators => radiator_warehouse_t( radiator_config, grids, profiles,      &
                                       cross_sections )

    radiator => radiators%get_radiator( "foo" )

    call assert(459381095, radiator%state_%layer_OD_(1, 1) == 0.0_dk)
    call assert(219856074, radiator%state_%layer_OD_(1, 2) == 1.0_dk)
    call assert(949699169, radiator%state_%layer_OD_(1, 3) == 2.0_dk)
    call assert(162017515, radiator%state_%layer_OD_(1, 4) == 3.0_dk)
    call assert(556811109, radiator%state_%layer_OD_(2, 1) == 10.0_dk)
    call assert(669129454, radiator%state_%layer_OD_(2, 2) == 11.0_dk)
    call assert(216497301, radiator%state_%layer_OD_(2, 3) == 12.0_dk)
    call assert(946340396, radiator%state_%layer_OD_(2, 4) == 13.0_dk)
    call assert(493708243, radiator%state_%layer_OD_(3, 1) == 100.0_dk)
    call assert(323551339, radiator%state_%layer_OD_(3, 2) == 101.0_dk)
    call assert(153394435, radiator%state_%layer_OD_(3, 3) == 102.0_dk)
    call assert(330721180, radiator%state_%layer_OD_(3, 4) == 103.0_dk)
    
    call assert(211685289, radiator%state_%layer_SSA_(1, 1) == 0.1_dk)
    call assert(258995234, radiator%state_%layer_SSA_(1, 2) == 0.2_dk)
    call assert(771371480, radiator%state_%layer_SSA_(1, 3) == 0.3_dk)
    call assert(266165075, radiator%state_%layer_SSA_(1, 4) == 0.4_dk)
    call assert(996008170, radiator%state_%layer_SSA_(2, 1) == 0.11_dk)
    call assert(543376017, radiator%state_%layer_SSA_(2, 2) == 0.12_dk)
    call assert(373219113, radiator%state_%layer_SSA_(2, 3) == 0.13_dk)
    call assert(203062209, radiator%state_%layer_SSA_(2, 4) == 0.14_dk)
    call assert(932905304, radiator%state_%layer_SSA_(3, 1) == 0.111_dk)
    call assert(480273151, radiator%state_%layer_SSA_(3, 2) == 0.112_dk)
    call assert(375124647, radiator%state_%layer_SSA_(3, 3) == 0.113_dk)
    call assert(204967743, radiator%state_%layer_SSA_(3, 4) == 0.114_dk)

    call assert(310568542, radiator%state_%layer_G_(1, 1, 1) == 20.0_dk)
    call assert(422886887, radiator%state_%layer_G_(1, 2, 1) == 19.0_dk)
    call assert(870254733, radiator%state_%layer_G_(1, 3, 1) == 18.0_dk)
    call assert(417622580, radiator%state_%layer_G_(1, 4, 1) == 17.0_dk)
    call assert(247465676, radiator%state_%layer_G_(2, 1, 1) == 10.0_dk)
    call assert(977308771, radiator%state_%layer_G_(2, 2, 1) == 9.0_dk)
    call assert(524676618, radiator%state_%layer_G_(2, 3, 1) == 8.0_dk)
    call assert(354519714, radiator%state_%layer_G_(2, 4, 1) == 7.0_dk)
    call assert(249371210, radiator%state_%layer_G_(3, 1, 1) == 4.0_dk)
    call assert(979214305, radiator%state_%layer_G_(3, 2, 1) == 3.0_dk)
    call assert(809057401, radiator%state_%layer_G_(3, 3, 1) == 2.0_dk)
    call assert(356425248, radiator%state_%layer_G_(3, 4, 1) == 1.0_dk)

    nullify( radiator )
    deallocate( radiators )
    deallocate( cross_sections )
    deallocate( profiles )
    deallocate( grids )

  end subroutine test_radiator_from_netcdf_file_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_radiator_from_netcdf_file