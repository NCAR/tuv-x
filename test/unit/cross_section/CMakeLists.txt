################################################################################
# Test utilities

include(test_util)

################################################################################
# Photo-decomp tests

create_standard_test(NAME cross_section_warehouse SOURCES cross_section_warehouse.F90)
create_standard_test(NAME cross_section SOURCES base.F90)
create_standard_test(NAME cross_section_acetone_ch3co_ch3_t SOURCES acetone-ch3co_ch3_test.F90 )
create_standard_test(NAME cross_section_ccl4 SOURCES ccl4_test.F90 )
create_standard_test(NAME cross_section_cfc-11 SOURCES cfc-11_test.F90 )
create_standard_test(NAME cross_section_ch2o SOURCES ch2o_test.F90 )
create_standard_test(NAME cross_section_ch3ono2_ch3o_no2 SOURCES ch3ono2_ch3o_no2_test.F90 )
create_standard_test(NAME cross_section_chbr3 SOURCES chbr3_test.F90 )
create_standard_test(NAME cross_section_chcl3 SOURCES chcl3_test.F90 )
create_standard_test(NAME cross_section_cl2-cl_cl SOURCES cl2-cl_cl_test.F90 )
create_standard_test(NAME cross_section_clono2 SOURCES clono2_test.F90 )
create_standard_test(NAME cross_section_h2o2-oh_oh SOURCES h2o2-oh_oh_test.F90 )
create_standard_test(NAME cross_section_hcfc SOURCES hcfc_test.F90 )
create_standard_test(NAME cross_section_hno3-oh_no2 SOURCES hno3-oh_no2_test.F90 )
create_standard_test(NAME cross_section_hobr-oh_br SOURCES hobr-oh_br_test.F90 )
create_standard_test(NAME cross_section_n2o-n2_o1d SOURCES n2o-n2_o1d_test.F90 )
create_standard_test(NAME cross_section_n2o5-no2_no3 SOURCES n2o5-no2_no3_test.F90 )
create_standard_test(NAME cross_section_nitroxy_acetone SOURCES nitroxy_acetone_test.F90 )
create_standard_test(NAME cross_section_nitroxy_ethanol SOURCES nitroxy_ethanol_test.F90 )
create_standard_test(NAME cross_section_no2_tint SOURCES no2_tint_test.F90 )
create_standard_test(NAME cross_section_o3_tint SOURCES o3_tint_test.F90 )
create_standard_test(NAME cross_section_oclo SOURCES oclo_test.F90 )
create_standard_test(NAME cross_section_rayliegh SOURCES rayliegh_test.F90 )
create_standard_test(NAME cross_section_rono2 SOURCES rono2_test.F90 )
create_standard_test(NAME cross_section_t_butyl_nitrate SOURCES t_butyl_nitrate_test.F90 )
create_standard_test(NAME cross_section_tint SOURCES tint_test.F90 )

################################################################################
