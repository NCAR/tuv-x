################################################################################
# cross section and quantum yield calculators from Doug's TUV code

target_sources(tuv_doug
  PRIVATE
    XSQY_ACETONE.f
    XSQY_BRO.f
    XSQY_BRONO2.f
    XSQY_CF2CL2.f
    XSQY_CFC113.f
    XSQY_CFC114.f
    XSQY_CFC115.f
    XSQY_CFCL3.f
    XSQY_CH2BR2.f
    XSQY_CH3BR.f
    XSQY_CH3CL.f
    XSQY_CHBR3.f
    XSQY_CL2O2.f
    XSQY_CLO.f
    XSQY_H2O.f
    XSQY_H1301.f
    XSQY_H2402.f
    XSQY_HCFC22.f
    XSQY_HCFC141b.f
    XSQY_HCFC142b.f
    XSQY_HNO3.f
    XSQY_HO2NO2.f
    XSQY_N2O5.f
    XSQY_SO2.f
)

################################################################################
