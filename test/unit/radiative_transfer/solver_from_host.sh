#!/bin/bash

# turn on command echoing
set -v
# move to the build directory, which holds the test data
cd ${0%/*}/../../..
# define a function for failure tests
failure_test () {
  local expected_failure=$1
  local output=$(./solver_from_host_failure $1 2>&1)
  local failure_code=$(echo $output | sed -n 's/.*ERROR (Musica-\([0-9][0-9]*\).*/\1/p')
  if ! [ "$failure_code" = "$expected_failure" ]; then
    echo "Expected failure $expected_failure"
    echo "Got output: $output"
    exit 1
  else
    local failure_code=$(cat error.json | sed -n 's/[[:space:]]*\"code\" : \"\([0-9][0-9]*\).*/\1/p')
    if ! [ "$failure_code" = "$expected_failure" ]; then
      echo "Expected failure $expected_failure in file 'error.json'"
      echo "Got: $(cat error.json)"
      rm -f error.json
      exit 1
    else
      rm -f error.json
      echo $output
    fi
  fi
}

# bad shape for a required actinic flux component
failure_test 254866372
# bad shape for an optional irradiance component
failure_test 785646610
# no from host solver, and the caller omitted the found flag
failure_test 921177509
# an update through an updater that has no solver behind it
failure_test 419530284

exit 0
