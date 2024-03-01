######################################################################
# Utility source

target_sources(tuvx_object
  PRIVATE
    array.F90
    assert.F90
    constants.F90
    convert.F90
    data_type.F90
    datetime.F90
    iterator.F90
    io.F90
    map.F90
    mpi.F90
    string.F90
    target.F90
)

target_sources(tuvx_object
  PRIVATE
    config.F90
    config.cpp
    yaml_util.F90
)

add_subdirectory(io)

######################################################################
