#!/bin/bash
cd src

file1=parameters.f90
file2=aux.f90
file3=model.f90
file4=methods.f90
main=main.f90
exefile=exe

gfortran -funderscoring -C -g -Wall -o $exefile $file1 $file2 $file3 $file4 $main
rm *.mod*
mv $exefile ../$exefile
cd ..
./exe

gnuplot plot2d.gnu
gnuplot plot1d.gnu
