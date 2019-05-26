#! /bin/sh

#path="/eos/cms/store/group/phys_egamma/soffi/XXQQQQ-GEN-SIM/DisplJetAnalyzer/XXQQQQ_m50/"
path="../python/"
outpath="~/www/DispJets/GenStudies/"

echo "Running samples"
root -l -b -q "Plotter.C("\"${path}\"","\"${outpath}\"")"
echo "Done plotter"

echo "Cleaning"
rm *~
rm *_cpp*
rm *_C.d
