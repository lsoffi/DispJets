#! /bin/sh

path="../samples/"
outpath="~/www/Plots/DispJets/GenStudies/"

echo "Running samples"
root -l -b -q "Plotter.C("\"${path}\"","\"${outpath}\"")"
echo "Done plotter"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
