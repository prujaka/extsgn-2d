#!/bin/bash
cd src

file1=parameters.f90
file2=methods.f90
main=main.f90
exefile=exe

gfortran -O3 -o $exefile $file1 $file2 $main
rm *mod*
mv $exefile ../$exefile
cd ..
./exe

gnuplot plot2d.gnu
