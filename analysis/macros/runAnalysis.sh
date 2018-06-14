#! /bin/sh

echo "Running analysis"
root -l -b -q "Analysis.C()"
echo "Done"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
