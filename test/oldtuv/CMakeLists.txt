################################################################################
# Copy test data

add_custom_target(copy_oldtuv_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/../data/old_tuv ${CMAKE_BINARY_DIR}/odat)

################################################################################
# Photo-decomp tool source

set(PHOTO_SRC abs_linalgebra.f90
              # abstract.cross_section.type.F90
              abstract_radXfer.f90
              delta_eddington.f90
              disord.f90
              disord_subs.f
              environment.F90
              functs.f
              grids.f
              la_srb.f
              lepoly.f
              linalgebra.f90
              netcdf_util.F90
              newlst.f
              numer.f
              odo3.f
              odrl.f
              orbit.f
              photo_kinetics.F90
              photo_utils.F90
              psndo.f
              quantum_yield_factory.F90
              qys.f
              # radXfer_xsect_factory.F90
              # radXfer_xsect_warehouse.F90
              rdetfl.f
              rdinp.f
              rdxs.f
              rtlink.f
              rtrans.f
              rxn.f
              savout.f
              setaer.f
              setalb.f
              setcld.f
              #              seth2o.f - Flagged as inactive in source comments
              setmtx.f
              setno2.f
              seto2.f
              setsnw.f
              setso2.f
              sphers.f
              solvec.f
              # spectral_wght_factory.F90
              # spectral_wght_warehouse.F90
              swbiol.f
              swchem.f
              swdom.f
              swphys.f
              terint.f
              tuv_params.f
              vpair.f
              vpo3.f
              vptmp.f
              waters.f
              #   wavelength.F90 - not used
              wrflut.f
              wshift.f
              zeroal.f
              diagout.f90
              )

file(GLOB LITTLE_F_SRC *.f)
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set_property(SOURCE ${LITTLE_F_SRC} APPEND PROPERTY COMPILE_OPTIONS "-fdefault-real-8")
  set_property(SOURCE ${LITTLE_F_SRC} APPEND PROPERTY COMPILE_OPTIONS "-fdefault-double-8")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set_property(SOURCE ${LITTLE_F_SRC} APPEND PROPERTY COMPILE_OPTIONS "-r8")
  set_property(SOURCE ${LITTLE_F_SRC} APPEND PROPERTY COMPILE_OPTIONS "-autodouble")
endif()

file(GLOB LITTLE_F90_SRC *.f90)
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set_property(SOURCE ${LITTLE_F90_SRC} APPEND PROPERTY COMPILE_OPTIONS "-fdefault-real-8")
  set_property(SOURCE ${LITTLE_F90_SRC} APPEND PROPERTY COMPILE_OPTIONS "-fdefault-double-8")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set_property(SOURCE ${LITTLE_F90_SRC} APPEND PROPERTY COMPILE_OPTIONS "-r8")
  set_property(SOURCE ${LITTLE_F90_SRC} APPEND PROPERTY COMPILE_OPTIONS "-autodouble")
endif()

add_subdirectory(cross_section)
add_subdirectory(grid)
add_subdirectory(quantum_yield)
add_subdirectory(radiator)
add_subdirectory(radXfer_cross_section)
add_subdirectory(spectral_wght)
add_subdirectory(util)
add_subdirectory(Profile)

add_library(oldphotolib ${PHOTO_SRC} ${GRID_SRC} ${PROFILE_SRC} ${CROSS_SRC} ${QY_SRC} ${RADIATOR_SRC} ${RAD_SRC} ${SW_SRC} ${UTIL_SRC})
target_link_libraries(oldphotolib musica::musicacore ${NETCDF_LIBS} ${JSON_LIB})
set_target_properties(oldphotolib PROPERTIES OUTPUT_NAME oldphoto)

add_executable(oldtuv tuv.f)
target_link_libraries(oldtuv oldphotolib)

################################################################################
