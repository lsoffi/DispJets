#! /bin/sh

path="~/www/Plots/DispJets/GenStudies/"

root -l -b -q "Multiplot.C("\"${path}\"")"

rm *~
rm *_cpp*
rm *_C.d
