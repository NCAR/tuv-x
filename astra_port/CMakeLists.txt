cmake_minimum_required(VERSION 3.20)
project(TurboTUV VERSION 0.1.0 LANGUAGES C)
include(CTest)
include(GNUInstallDirs)

option(TURBOTUV_SANITIZE "Enable address and undefined-behavior sanitizers" OFF)
option(TURBOTUV_BUILD_CLI "Build the standalone POSIX runner" ON)
option(TURBOTUV_NETCDF_OUTPUT "Enable NetCDF output in the standalone runner" OFF)
option(TURBOTUV_BUILD_BENCHMARK "Build the column benchmark (not installed)" ON)
set(TURBOTUV_REFERENCE_BUILD "" CACHE PATH "Configured and built TUV-x reference directory")
set(TURBOTUV_REFERENCE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../TUV-x" CACHE PATH "TUV-x reference checkout")
set(TURBOTUV_REFERENCE_PREFIX "" CACHE PATH "Compiler-matched NetCDF Fortran installation")
set(TURBOTUV_NETCDF_PREFIX "" CACHE PATH "NetCDF C installation (reference tests only)")

add_library(turbotuv src/numerics.c src/geometry.c src/data.c src/radiation.c src/spectra.c
    src/yields.c src/weights.c src/temperature.c src/rates.c src/oxygen.c src/linalg.c
    src/discrete_eigen.c src/discrete_quadrature.c src/discrete.c
    src/grid.c src/profiles.c src/tables.c src/solar_flux.c src/optics.c
    src/config.c src/profile_config.c src/build.c src/spectrum.c src/spectrum_temperature.c
    src/model_build.c src/model_run.c src/model_host.c src/decimal.c
    src/spectrum_plan.c src/model_plan.c src/radiation_block.c)
add_library(TurboTUV::turbotuv ALIAS turbotuv)
target_compile_features(turbotuv PUBLIC c_std_11)
set_target_properties(turbotuv PROPERTIES C_EXTENSIONS OFF POSITION_INDEPENDENT_CODE ON)
target_include_directories(turbotuv PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_include_directories(turbotuv PRIVATE third_party/libyaml third_party/libyaml/include)
add_library(turbotuv_yaml OBJECT
    third_party/libyaml/src/api.c third_party/libyaml/src/loader.c
    third_party/libyaml/src/parser.c third_party/libyaml/src/reader.c
    third_party/libyaml/src/scanner.c)
set_target_properties(turbotuv_yaml PROPERTIES POSITION_INDEPENDENT_CODE ON C_VISIBILITY_PRESET hidden)
target_include_directories(turbotuv_yaml PRIVATE third_party/libyaml/include)
target_compile_definitions(turbotuv_yaml PRIVATE YAML_DECLARE_STATIC
    YAML_VERSION_MAJOR=0 YAML_VERSION_MINOR=2 YAML_VERSION_PATCH=5 YAML_VERSION_STRING="0.2.5")
if(MSVC)
    target_compile_options(turbotuv_yaml PRIVATE /FI${CMAKE_CURRENT_SOURCE_DIR}/third_party/libyaml/namespace.h)
else()
    target_compile_options(turbotuv_yaml PRIVATE -include ${CMAKE_CURRENT_SOURCE_DIR}/third_party/libyaml/namespace.h)
endif()
target_sources(turbotuv PRIVATE $<TARGET_OBJECTS:turbotuv_yaml>)
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(turbotuv PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off)
endif()
if(UNIX)
    target_link_libraries(turbotuv PUBLIC m)
endif()
if(TURBOTUV_SANITIZE)
    target_compile_options(turbotuv_yaml PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_compile_options(turbotuv PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(turbotuv PUBLIC -fsanitize=address,undefined)
endif()
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
if(TURBOTUV_BUILD_CLI)
    add_executable(turbotuv_cli app/main.c app/output.c)
    set_target_properties(turbotuv_cli PROPERTIES OUTPUT_NAME turbotuv C_EXTENSIONS OFF)
    target_link_libraries(turbotuv_cli PRIVATE turbotuv)
    target_compile_definitions(turbotuv_cli PRIVATE TT_VERSION="${PROJECT_VERSION}")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(turbotuv_cli PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off)
    endif()
    if(TURBOTUV_NETCDF_OUTPUT)
        find_package(PkgConfig REQUIRED)
        pkg_check_modules(NETCDF REQUIRED IMPORTED_TARGET netcdf)
        target_link_libraries(turbotuv_cli PRIVATE PkgConfig::NETCDF)
        target_compile_definitions(turbotuv_cli PRIVATE TT_NETCDF_OUTPUT)
    endif()
    install(TARGETS turbotuv_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
if(TURBOTUV_BUILD_BENCHMARK)
    add_executable(turbotuv_bench bench/turbotuv_bench.c)
    set_target_properties(turbotuv_bench PROPERTIES C_EXTENSIONS OFF)
    target_link_libraries(turbotuv_bench PRIVATE turbotuv)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(turbotuv_bench PRIVATE -Wall -Wextra -Wpedantic -ffp-contract=off)
    endif()
endif()
install(TARGETS turbotuv EXPORT TurboTUVTargets)
install(FILES include/turbotuv.h include/turbotuv_config.h include/turbotuv_model.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES third_party/libyaml/License DESTINATION ${CMAKE_INSTALL_DATADIR}/TurboTUV/licenses RENAME LibYAML.txt)
install(FILES LICENSE DESTINATION ${CMAKE_INSTALL_DATADIR}/TurboTUV/licenses RENAME TurboTUV.txt)
install(EXPORT TurboTUVTargets NAMESPACE TurboTUV:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/TurboTUV)
include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/TurboTUVConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/TurboTUVConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/TurboTUV)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/TurboTUVConfigVersion.cmake"
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/TurboTUVConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/TurboTUVConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/TurboTUV)
install(DIRECTORY examples/ DESTINATION ${CMAKE_INSTALL_DATADIR}/TurboTUV/examples
    FILES_MATCHING PATTERN "*.json" PATTERN "*.yml" PATTERN "*.c")
