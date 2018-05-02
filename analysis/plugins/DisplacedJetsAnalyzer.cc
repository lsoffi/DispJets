// -*- C++ -*-
//
// Package:    DispJets/DisplacedJetsAnalyzer
// Class:      DisplacedJetsAnalyzer
// 
/**\class DisplacedJetsAnalyzer DisplacedJetsAnalyzer.cc DispJets/plugins/DisplacedJetsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Margaret Zientek
//         Created:  Mon, 30 Apr 2018 17:04:03 GMT
//
//

//---------------------------------------------------------------------------------------------------
// system include files
#include <memory>

// user include files
// generic
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
// specific
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
// root include files
#include "TTree.h"

//---------------------------------------------------------------------------------------------------
// ttree structure

struct tree_struc_{
  long unsigned int		event;
  int				run;
  int				lumi;
  int				ngenjets;
  std::vector<float>            genjet_pt;
  std::vector<float>            genjet_e;
  std::vector<float>            genjet_eta;
  std::vector<float>            genjet_phi;
};

//---------------------------------------------------------------------------------------------------
// class declaration

class DisplacedJetsAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DisplacedJetsAnalyzer(const edm::ParameterSet&);
      ~DisplacedJetsAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void initTreeStructure();

      // ----------member data ---------------------------
      edm::Service<TFileService> fs;
     
      bool verbose_;

      // tokens
      edm::EDGetTokenT<std::vector<reco::GenJet> > genjetToken_;

      // setup tree;
      TTree* tree;
      tree_struc_ tree_;
     
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//---------------------------------------------------------------------------------------------------
// ------------ constructor  ------------
DisplacedJetsAnalyzer::DisplacedJetsAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   verbose_		= iConfig.getUntrackedParameter<bool>("verbose");
   genjetToken_		= consumes<std::vector<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genjets"));


}


//---------------------------------------------------------------------------------------------------
// ------------ destructor  ------------
DisplacedJetsAnalyzer::~DisplacedJetsAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if (verbose_) std::cout << "Done" << std::endl;

}


//---------------------------------------------------------------------------------------------------
// ------------ method called for each event  ------------
void
DisplacedJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;

   edm::Handle<std::vector<reco::GenJet> > genjets;
   iEvent.getByToken(genjetToken_, genjets);

   // general event info 
   unsigned long int event = iEvent.id().event();   
   int run   = iEvent.id().run();
   int lumi  = iEvent.luminosityBlock();

   // genjet info
   int ngenjets = 0;
   std::vector<float> genjet_pt;
   std::vector<float> genjet_e;
   std::vector<float> genjet_eta;
   std::vector<float> genjet_phi;
   if (genjets.isValid()){ // make sure have genjet collection
     ngenjets = genjets->size();

     for (const auto & genjet_iter : *genjets){ // loop over genjets
 
       genjet_pt.push_back(genjet_iter.pt()); 
       genjet_e.push_back(genjet_iter.energy()); 
       genjet_eta.push_back(genjet_iter.eta()); 
       genjet_phi.push_back(genjet_iter.phi()); 
 
     }// end loop over genjets 
   }// end if genjets.isValid
   else std::cout << "WARNING: genjets collection is NOT valid" << std::endl;

   // initialize tree
   initTreeStructure();
 
   // setup tree values
   tree_.run		= run;
   tree_.lumi		= lumi;
   tree_.event		= event;
   tree_.ngenjets	= ngenjets;
   for (int ij = 0; ij < ngenjets; ij++){
     tree_.genjet_pt.push_back(genjet_pt[ij]);
     tree_.genjet_e.push_back(genjet_e[ij]);
     tree_.genjet_eta.push_back(genjet_eta[ij]);
     tree_.genjet_phi.push_back(genjet_phi[ij]);
   }

   // fill trees
   tree->Fill();
 
}


//---------------------------------------------------------------------------------------------------
// ------------ method called once each job just before starting event loop  ------------
void 
DisplacedJetsAnalyzer::beginJob()
{

  if (verbose_) std::cout << "Starting job" << std::endl;

  // set up output tree
  tree = fs->make<TTree>("tree","tree");

  // set tree branches
  tree->Branch("run",		&tree_.run,		"run/I");
  tree->Branch("lumi",		&tree_.lumi,		"lumi/I");
  tree->Branch("event",		&tree_.event,		"event/L");
  tree->Branch("ngenjets",	&tree_.ngenjets,	"ngenjets/I");
  tree->Branch("genjet_pt",     &tree_.genjet_pt);
  tree->Branch("genjet_e",      &tree_.genjet_e);
  tree->Branch("genjet_eta",	&tree_.genjet_eta);
  tree->Branch("genjet_phi",	&tree_.genjet_phi);

}


//---------------------------------------------------------------------------------------------------
// initialize trees
void DisplacedJetsAnalyzer::initTreeStructure()
{

  tree_.run		= -500.;
  tree_.lumi		= -500.;
  tree_.event		= 0.;
  tree_.ngenjets	= -500.;

}

//---------------------------------------------------------------------------------------------------
// ------------ method called once each job just after ending the event loop  ------------
void 
DisplacedJetsAnalyzer::endJob() 
{
}

//---------------------------------------------------------------------------------------------------
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedJetsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//---------------------------------------------------------------------------------------------------
//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedJetsAnalyzer);
