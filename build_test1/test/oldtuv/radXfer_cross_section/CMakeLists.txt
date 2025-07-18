################################################################################
# Radiation transfer source

set(RAD_SRC radXfer_cross_section/abstract.radXfer.cross_section.type.F90
            radXfer_cross_section/base.cross_section.type.F90
            radXfer_cross_section/o3.tint.cross_section.type.F90
            radXfer_cross_section/tint.cross_section.type.F90
            radXfer_cross_section/rayliegh.cross_section.type.F90
            radXfer_cross_section/radXfer_xsect_warehouse.F90
            radXfer_cross_section/radXfer_xsect_factory.F90
            PARENT_SCOPE)

################################################################################
