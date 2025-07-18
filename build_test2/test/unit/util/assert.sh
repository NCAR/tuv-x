#!/bin/bash

# turn on command echoing
set -v
# move to the directory this script is in
cd ${0%/*}
# define a function for failure tests
failure_test () {
  local expected_failure=$(echo $1 | sed -n 's/\([[:digit:]]\+\).*/\1/p')
  local output=$(../../../util_assert_failure $1 2>&1)
  local failure_code=$(echo $output | sed -n 's/[[:space:]]*ERROR (Musica-\([[:digit:]]\+\).*/\1/p')
  if ! [ "$failure_code" = "$expected_failure" ]; then
    echo "Expected failure $expected_failure"
    echo "Got output: $output"
    exit 1
  else
    local expected_failure=$(echo $1 | sed -n 's/\([[:digit:]]\+\)/\1/p')
    local failure_code=$(cat error.json | sed -n 's/[[:space:]]*\"code\" : \"\([[:digit:]]\+\).*/\1/p')
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

failure_test 903602145
failure_test 151700878

exit 0
