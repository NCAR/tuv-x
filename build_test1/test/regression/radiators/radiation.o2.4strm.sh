#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v

tmpdir=$(mktemp -d)
basedir=$(pwd)
ln -s $tmpdir `basename $tmpdir`

cp -r odat $tmpdir/odat
cp oldtuv $tmpdir
ln -s $basedir/data $tmpdir/data
ln -s $basedir/test $tmpdir/test

cd $tmpdir

exec_oldtuv() {
  ./oldtuv DO_O2 < $basedir/test/regression/tuv_scenario_6.in
}
exec_newtuv() {
  $basedir/tuv-x $basedir/test/data/radiators.o2.4strm.config.json
}
exec_analysis() {
  python3 $basedir/tool/diagnostics/var.compare.py $basedir/test/regression/radiators/radiation.o2.compare.json
}

if ! exec_oldtuv; then
  echo FAIL - old TUV
  exit 1
fi

if ! exec_newtuv; then
  echo FAIL - new TUV
  exit 1
fi

if ! exec_analysis; then
  echo FAIL - analysis
  exit 1
fi

echo PASS
exit 0
