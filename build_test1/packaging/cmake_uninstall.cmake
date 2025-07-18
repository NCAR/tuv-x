# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake

if(NOT EXISTS "/home/runner/work/tuv-x/tuv-x/build_test1/install_manifest.txt")
  message(FATAL_ERROR "Cannot find install manifest: /home/runner/work/tuv-x/tuv-x/build_test1/install_manifest.txt")
endif()

file(READ "/home/runner/work/tuv-x/tuv-x/build_test1/install_manifest.txt" files)
string(REGEX REPLACE "\n" ";" files "${files}")
foreach(file ${files})
  message(STATUS "Uninstalling $ENV{DESTDIR}${file}")
  if(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    exec_program(
      "/usr/local/bin/cmake" ARGS "-E remove \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    if(NOT "${rm_retval}" STREQUAL 0)
      message(FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}")
    endif()
  else(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    message(STATUS "File $ENV{DESTDIR}${file} does not exist.")
  endif()
endforeach()

set(uninstall_dirs "tuvx-0.12.0;tuvx-0.12.0/cmake/tuvx-0.12.0")

foreach(dir IN LISTS uninstall_dirs)
  message(STATUS "Uninstalling ${dir}")
  exec_program(
    "/usr/local/bin/cmake" ARGS "-E remove_directory ${dir}"
    OUTPUT_VARIABLE rm_out
    RETURN_VALUE rm_retval
    )
  if(NOT "${rm_retval}" STREQUAL 0)
    message(FATAL_ERROR "Problem when removing ${dir}")
  endif()
endforeach()
