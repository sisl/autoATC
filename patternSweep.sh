#!/bin/sh

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 3 ] || die "3 arguments required, $# provided"

script=$1
n0=$2
nend=$3

echo $n0   | grep -E -q '^[0-9]+$' || die "Numeric argument required, $n0 provided"
echo $nend | grep -E -q '^[0-9]+$' || die "Numeric argument required, $nend provided"

[ -f "$script" ] || die "$script does not exist"

echo Sweeping $script with nPhases in [$n0 - $nend] 

for n in `seq $n0 $nend`
do
  sed -i "s/const nPhases =/const nPhases = $n #/g" pattern.jl
  echo Running $script with nPhases = $n
  sleep 1 #precaution to let NFS sync...
  julia $script
done