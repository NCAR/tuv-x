######################################################################
# Utility source

target_sources(tuvx_object
  PRIVATE
    array.F90
    assert.F90
    config.F90
    config.cpp
    constants.F90
    iterator.F90
    io.F90
    map.F90
    mpi.F90
    string.F90
    yaml_util.F90
)

add_subdirectory(io)

######################################################################
