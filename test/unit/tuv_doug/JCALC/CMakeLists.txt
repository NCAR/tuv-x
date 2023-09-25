################################################################################
# cross section and quantum yield calculators from Doug's TUV code

target_sources(tuv_doug
  PRIVATE
    XSQY_CH2BR2.f
    XSQY_H2O.f
    XSQY_N2O5.f
)

################################################################################
