#!/bin/sh

for n in `seq 1 5`
do
  sed -i "s/const nPhases =/const nPhases = $n #/g" pattern.jl
  echo Running with nPhases = $n
  julia runSims.jl
done