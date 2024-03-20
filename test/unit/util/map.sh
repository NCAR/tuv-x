#!/bin/bash

# turn on command echoing
set -v
# move to the directory this script is in
cd ${0%/*}
# define a function for failure tests
failure_test () {
  local expected_failure=$(echo $1 | sed -n 's/\([[:digit:]]\+\).*/\1/p')
  local output=$(../../../util_map_failure $1 2>&1)
  local failure_code=$(echo $output | sed -n 's/[[:space:]]*ERROR (Musica-\([[:digit:]]\+\).*/\1/p')
  if ! [ "$failure_code" = "$expected_failure" ]; then
    echo "Expected failure $expected_failure"
    echo "Got output: $output"
    exit 1
  else
    echo $output
  fi
}

failure_test 170733942
failure_test 764798475
failure_test 133386338
failure_test 956987954
failure_test 200274675
failure_test 240867074
failure_test 309595761
failure_test 122570601
failure_test 740547646
failure_test 548594113

exit 0
