#!/bin/bash

# exit on error
set -e
# turn on command echoing
set -v
# copy photolysis rate configuration into template
sed -e '/#PHOTO_RATES/ {' -e 'r data/photolysis_rate_constants.json' -e 'd' -e '}' -i test/data/photorates.test.config.json

tmpdir=$(mktemp -d)
basedir=$(pwd)
ln -s $tmpdir `basename $tmpdir`

cp -r odat $tmpdir/odat
cp oldtuv $tmpdir
ln -s $basedir/data $tmpdir/data
ln -s $basedir/test $tmpdir/test

cd $tmpdir

exec_oldtuv() {
  ./oldtuv DO_RAYLEIGH DO_O2 DO_O3 DO_AEROSOLS DO_CLOUDS < $basedir/test/regression/tuv_scenario_2.in
}
exec_newtuv() {
  valgrind --error-exitcode=1 --trace-children=yes --leak-check=full --gen-suppressions=all --suppressions=test/valgrind.supp $basedir/tuv-x $basedir/test/data/photorates.test.config.json
}
exec_analysis() {
  python3 $basedir/test/regression/photolysis_rates/xsqy.compare.py $basedir/test/regression/photolysis_rates odat/OUTPUTS output
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
