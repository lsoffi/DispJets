#! /bin/sh

path="~/www/DispJets/GenStudies/"

echo "Running samples"
root -l -b -q "Overlay.C("\"${path}\"")"
echo "Done overlay"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
