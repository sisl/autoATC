#!/bin/sh


for n in `seq 1 $2`
do
  sed -i "s/const nPhases =/const nPhases = $n #/g" pattern.jl
  echo Running $1 with nPhases = $n
  julia $1
done