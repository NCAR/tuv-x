################################################################################
# cross section and quantum yield calculators from Doug's TUV code

target_sources(tuv_doug
  PRIVATE
    XSQY_BRO.f
    XSQY_CH2BR2.f
    XSQY_CL2O2.f
    XSQY_CLO.f
    XSQY_H2O.f
    XSQY_HNO3.f
    XSQY_N2O5.f
)

################################################################################
