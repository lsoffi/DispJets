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

  int				sample;
  long unsigned int		event;
  int				run;
  int				lumi;
  int				ngenjets;
  std::vector<float>            genjet_pt;
  std::vector<float>            genjet_e;
  std::vector<float>            genjet_eta;
  std::vector<float>            genjet_phi;
  std::vector<int>		genjet_ndaug;
  std::vector<int>		genjet_nconst;
  std::vector<float>		genjet_vx;
  std::vector<float>		genjet_vy;
  std::vector<float>		genjet_vz;
  int				ngenpart;
  std::vector<int>		genpar_id;
  std::vector<float>		genpar_pt;
  std::vector<float>		genpar_eta;
  std::vector<float>		genpar_phi;
  std::vector<float>		genpar_vx;
  std::vector<float>		genpar_vy;
  std::vector<float>		genpar_vz;
  std::vector<float>		genpar_Lxy;
  std::vector<float>		genpar_Lxyz;
  std::vector<int>		mom_id;
  std::vector<int>		mom_stat;
  std::vector<float>		mom_e;
  std::vector<float>		mom_m;
  std::vector<float>		mom_pt;
  std::vector<float>		mom_eta;
  std::vector<float>		mom_phi;
  std::vector<float>    	mom_vx;
  std::vector<float>		mom_vy;
  std::vector<float>		mom_vz;
  std::vector<float>		mom_beta;
  std::vector<float>		mom_gama;
  std::vector<float>		mom_Lxy;
  std::vector<float>		mom_Lz;
  std::vector<float>		mom_Lxyz;
  std::vector<float>		mom_ctau;

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
     
      bool	verbose_;
      int	sampleID_;

      // tokens
      edm::EDGetTokenT<std::vector<reco::GenJet> >		genjetToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle> >		genparticleToken_;

      // setup tree;
      TTree* tree;
      tree_struc_ tree_;
     
};

// constants, enums and typedefs
// static data member definitions

const auto pTsort = [](const auto& obj1, const auto& obj2){ return obj1->pt() > obj2->pt(); };

//---------------------------------------------------------------------------------------------------
// ------------ constructor  ------------
DisplacedJetsAnalyzer::DisplacedJetsAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
   verbose_		= iConfig.getUntrackedParameter<bool>("verbose");
   sampleID_		= iConfig.getUntrackedParameter<int>("sampleID",0);
   genjetToken_		= consumes<std::vector<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genjets"));
   genparticleToken_    = consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"));

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
void DisplacedJetsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;

  // --- pick up handles 
  edm::Handle<std::vector<reco::GenJet> > genjets;
  iEvent.getByToken(genjetToken_, genjets);

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByToken(genparticleToken_, genparticles);

  // --- general event info 
  unsigned long int event = iEvent.id().event();   
  int run   = iEvent.id().run();
  int lumi  = iEvent.luminosityBlock();

  // --- genjet info
  int ngenjets = 0;
  std::vector<float>		genjet_pt;
  std::vector<float>		genjet_e;
  std::vector<float>		genjet_eta;
  std::vector<float>		genjet_phi;
  std::vector<int>		genjet_ndaug;
  std::vector<int>		genjet_nconst;
  std::vector<float>		genjet_vx;
  std::vector<float>		genjet_vy;
  std::vector<float>		genjet_vz;

  if (genjets.isValid()){ // make sure have genjet collection
    ngenjets = genjets->size();

    for (const auto & genjet_iter : *genjets){ // loop over genjets
 
      // standard information
      genjet_pt.push_back(genjet_iter.pt()); 
      genjet_e.push_back(genjet_iter.energy()); 
      genjet_eta.push_back(genjet_iter.eta()); 
      genjet_phi.push_back(genjet_iter.phi());
      genjet_ndaug.push_back(genjet_iter.numberOfDaughters());
      genjet_vx.push_back(genjet_iter.vertex().x());
      genjet_vy.push_back(genjet_iter.vertex().y());
      genjet_vz.push_back(genjet_iter.vertex().z());
   
      // constituent information
      std::vector<const reco::GenParticle*> genjet_const = genjet_iter.getGenConstituents();
      genjet_nconst.push_back(genjet_const.size()); 

 
    }// end loop over genjets 
  }// end if genjets.isValid
  else std::cout << "WARNING: genjets collection is NOT valid" << std::endl;

  // --- genparticles info
  int ngenpart = 0;
  int nmothers = 0;
  std::vector<int>	genpar_id;
  std::vector<int>	genpar_stat;
  std::vector<float>	genpar_e;
  std::vector<float>	genpar_pt;
  std::vector<float>	genpar_eta;
  std::vector<float>	genpar_phi;
  std::vector<float>	genpar_vx;
  std::vector<float>	genpar_vy;
  std::vector<float>	genpar_vz;
  std::vector<float>	genpar_Lxy;
  std::vector<float>	genpar_Lxyz;
  std::vector<int>	mom_id;
  std::vector<int>	mom_stat;
  std::vector<float>	mom_e;
  std::vector<float>	mom_m;
  std::vector<float>	mom_pt;
  std::vector<float>	mom_eta;
  std::vector<float>	mom_phi;
  std::vector<float>    mom_vx;
  std::vector<float>	mom_vy;
  std::vector<float>	mom_vz;
  std::vector<float>	mom_beta;
  std::vector<float>	mom_gama;
  std::vector<float>	mom_Lxy;
  std::vector<float>	mom_Lz;
  std::vector<float>	mom_Lxyz;
  std::vector<float>	mom_ctau;

  if (genparticles.isValid()){ // make sure have genparticles collection

    for (const auto & genpar_iter : *genparticles){ // loop over genparticles
      //if (genpar_iter.status() != 1) continue; // scip any non-final state particles
      ngenpart++;  // number stable particles

      // vertex position 
      float vx = genpar_iter.vx();
      float vy = genpar_iter.vy();
      float vz = genpar_iter.vz();
      float gLxy  = std::sqrt(vx*vx + vy*vy);         // xy vertex flight
      float gLxyz = std::sqrt(vx*vx + vy*vy + vz*vz); // xyz vertex flight

      genpar_stat.push_back(genpar_iter.status());
      genpar_id.push_back(genpar_iter.pdgId());
      genpar_pt.push_back(genpar_iter.pt());
      genpar_eta.push_back(genpar_iter.eta());
      genpar_phi.push_back(genpar_iter.phi());
      genpar_vx.push_back(vx);
      genpar_vy.push_back(vy);
      genpar_vz.push_back(vz);
      genpar_Lxy.push_back(gLxy);
      genpar_Lxyz.push_back(gLxyz);

      float mx = 0;
      float my = 0;
      float mz = 0;
      float dx = 0;
      float dy = 0;
      float dz = 0;
      float mBeta = 0; 
      float mGama = 0; 
      float mLxy  = 0; 
      float mLxyz = 0; 
      float mcTau = 0; 

      int n = 0;

      //std::cout << "pdg ID: " << genpar_iter.pdgId() << " w/ status = " << genpar_iter.status() << "=> vx,vy,vz: " << vx << "," << vy << "," << vz << std::endl;
      //if (genpar_iter.mother() != NULL){
      //  for (const auto & mom : *genpar_iter.mother()){
      //    std::cout << "--- mom ID: " << mom.pdgId() << " w/ status = " << mom.status() << "=> vx,vy,vz: " << mom.vx() << "," << mom.vy() << "," << mom.vz() << std::endl; 
      //  }
      //} 

      if (genpar_iter.mother() != NULL){ 
        std::vector< const reco::Candidate * > mother;
        mother.push_back(genpar_iter.mother());
        std::sort(mother.begin(),mother.end(),pTsort);

        for (const auto & mom : *genpar_iter.mother()){
          nmothers++; // number of mothers
          if (verbose_) std::cout << "Number of mothers: " << nmothers << std::endl;
          if (n != 0) continue;
          n++;
          mx = mom.vx();
          my = mom.vy();
          mz = mom.vz();
          dx = vx - mx;
          dy = vy - my;
          dz = vz - mz;
          mBeta = mom.p()/mom.energy();    // mom beta
          mGama = mom.energy()/mom.mass(); // mom gamma
          mLxy  = std::sqrt(dx*dx + dy*dy); 
          mLxyz = std::sqrt(dx*dx + dy*dy + dz*dz); 
          mcTau = std::sqrt(dx*dx + dy*dy + dz*dz) / (mBeta * mGama); 

          mom_id.push_back(mom.pdgId());
          mom_stat.push_back(mom.status());
          mom_m.push_back(mom.mass());
          mom_e.push_back(mom.energy());
          mom_pt.push_back(mom.pt());
          mom_eta.push_back(mom.eta());
          mom_phi.push_back(mom.phi());
          mom_vx.push_back(mx);
          mom_vy.push_back(my);
          mom_vz.push_back(mz);
          mom_beta.push_back(mBeta);
          mom_gama.push_back(mGama);
          mom_Lz.push_back(dz);
          mom_Lxy.push_back(mLxy);
          mom_Lxyz.push_back(mLxyz);
          mom_ctau.push_back(mcTau);
        }
      } 

    }// end loop over genparticles
  }// end if genparticles.isValid
  else std::cout << "WARNING: genparticles collection is NOT valid" << std::endl;


  // --- setup tree values
  initTreeStructure();
  tree_.sample		= sampleID_;
  tree_.run		= run;
  tree_.lumi		= lumi;
  tree_.event		= event;
  tree_.ngenjets	= ngenjets;
  tree_.ngenpart        = ngenpart;

  for (int ij = 0; ij < ngenjets; ij++){
    tree_.genjet_pt.push_back(genjet_pt[ij]);
    tree_.genjet_e.push_back(genjet_e[ij]);
    tree_.genjet_eta.push_back(genjet_eta[ij]);
    tree_.genjet_phi.push_back(genjet_phi[ij]);
    tree_.genjet_ndaug.push_back(genjet_ndaug[ij]);
    tree_.genjet_nconst.push_back(genjet_nconst[ij]);
    tree_.genjet_vx.push_back(genjet_vx[ij]);
    tree_.genjet_vy.push_back(genjet_vy[ij]);
    tree_.genjet_vz.push_back(genjet_vz[ij]);
  }

  if (verbose_) std::cout << "mom_id   " << mom_id.size()    << std::endl; 
  if (verbose_) std::cout << "mom_stat " << mom_stat.size()  << std::endl; 
  if (verbose_) std::cout << "mom_e    " << mom_e.size()     << std::endl; 
  if (verbose_) std::cout << "mom_m    " << mom_m.size()     << std::endl;  
  if (verbose_) std::cout << "mom_pt   " << mom_pt.size()    << std::endl;  
  if (verbose_) std::cout << "mom_eta  " << mom_eta.size()   << std::endl;  
  if (verbose_) std::cout << "mom_phi  " << mom_phi.size()   << std::endl; 
  if (verbose_) std::cout << "mom_vx   " << mom_vx.size()    << std::endl;  
  if (verbose_) std::cout << "mom_vy   " << mom_vy.size()    << std::endl; 
  if (verbose_) std::cout << "mom_vz   " << mom_vz.size()    << std::endl; 
  if (verbose_) std::cout << "mom_beta " << mom_beta.size()  << std::endl;  
  if (verbose_) std::cout << "mom_gama " << mom_gama.size()  << std::endl;  
  if (verbose_) std::cout << "mom_Lxy  " << mom_Lxy.size()   << std::endl;  
  if (verbose_) std::cout << "mom_Lz   " << mom_Lz.size()    << std::endl; 
  if (verbose_) std::cout << "mom_Lxyz " << mom_Lxyz.size()  << std::endl;  
  if (verbose_) std::cout << "mom_ctau " << mom_ctau.size()  << std::endl;  

  for (int ip = 0; ip < ngenpart; ip++){
    tree_.genpar_id.push_back(genpar_id[ip]);
    tree_.genpar_pt.push_back(genpar_pt[ip]);
    tree_.genpar_eta.push_back(genpar_eta[ip]);
    tree_.genpar_phi.push_back(genpar_phi[ip]);
    tree_.genpar_vx.push_back(genpar_vx[ip]);
    tree_.genpar_vy.push_back(genpar_vy[ip]);
    tree_.genpar_vz.push_back(genpar_vz[ip]);
    tree_.genpar_Lxy.push_back(genpar_Lxy[ip]);
    tree_.genpar_Lxyz.push_back(genpar_Lxyz[ip]);
    tree_.mom_id.push_back(mom_id[ip]);
    tree_.mom_stat.push_back(mom_stat[ip]);
    tree_.mom_e.push_back(mom_e[ip]);
    tree_.mom_m.push_back(mom_m[ip]);
    tree_.mom_pt.push_back(mom_pt[ip]);
    tree_.mom_eta.push_back(mom_eta[ip]);
    tree_.mom_phi.push_back(mom_phi[ip]);
    tree_.mom_vx.push_back(mom_vx[ip]);
    tree_.mom_vy.push_back(mom_vy[ip]);
    tree_.mom_vz.push_back(mom_vz[ip]);
    tree_.mom_beta.push_back(mom_beta[ip]);
    tree_.mom_gama.push_back(mom_gama[ip]);
    tree_.mom_Lxy.push_back(mom_Lxy[ip]);
    tree_.mom_Lz.push_back(mom_Lz[ip]);
    tree_.mom_Lxyz.push_back(mom_Lxyz[ip]);
    tree_.mom_ctau.push_back(mom_ctau[ip]);
  }

  // --- fill trees
  tree->Fill();
 
}


//---------------------------------------------------------------------------------------------------
// ------------ method called once each job just before starting event loop  ------------
void DisplacedJetsAnalyzer::beginJob()
{

  if (verbose_) std::cout << "Starting job" << std::endl;

  // --- set up output tree
  tree = fs->make<TTree>("tree","tree");

  // --- set tree branches
  tree->Branch("sample",		&tree_.sample,		"sample/I");
  tree->Branch("run",			&tree_.run,		"run/I");
  tree->Branch("lumi",			&tree_.lumi,		"lumi/I");
  tree->Branch("event",			&tree_.event,		"event/L");

  // gen jet stuff
  tree->Branch("ngenjets",		&tree_.ngenjets,	"ngenjets/I");
  tree->Branch("genjet_pt",     	&tree_.genjet_pt);
  tree->Branch("genjet_e",      	&tree_.genjet_e);
  tree->Branch("genjet_eta",		&tree_.genjet_eta);
  tree->Branch("genjet_phi",		&tree_.genjet_phi);
  tree->Branch("genjet_ndaug",		&tree_.genjet_ndaug);
  tree->Branch("genjet_nconst",		&tree_.genjet_nconst);
  tree->Branch("genjet_vx",		&tree_.genjet_vx);
  tree->Branch("genjet_vy",		&tree_.genjet_vy);
  tree->Branch("genjet_vz",		&tree_.genjet_vz);

  // gen particle stuff
  tree->Branch("ngenpart",		&tree_.ngenpart,	"ngenpart/I");
  tree->Branch("genpar_id",		&tree_.genpar_id);
  tree->Branch("genpar_pt",		&tree_.genpar_pt);
  tree->Branch("genpar_eta",		&tree_.genpar_eta);
  tree->Branch("genpar_phi",		&tree_.genpar_phi);
  tree->Branch("genpar_vx",		&tree_.genpar_vx);
  tree->Branch("genpar_vy",		&tree_.genpar_vy);
  tree->Branch("genpar_vz",		&tree_.genpar_vz);
  tree->Branch("genpar_Lxy",		&tree_.genpar_Lxy);
  tree->Branch("genpar_Lxyz",		&tree_.genpar_Lxyz);

  // gen particle mom stuff
  tree->Branch("mom_id",		&tree_.mom_id);
  tree->Branch("mom_stat",		&tree_.mom_stat);
  tree->Branch("mom_e",			&tree_.mom_e);
  tree->Branch("mom_m",			&tree_.mom_m);
  tree->Branch("mom_pt",		&tree_.mom_pt);
  tree->Branch("mom_eta",		&tree_.mom_eta);
  tree->Branch("mom_phi",		&tree_.mom_phi);
  tree->Branch("mom_vx",		&tree_.mom_vx);
  tree->Branch("mom_vy",		&tree_.mom_vy);
  tree->Branch("mom_vz",		&tree_.mom_vz);
  tree->Branch("mom_beta",		&tree_.mom_beta);
  tree->Branch("mom_gama",		&tree_.mom_gama);
  tree->Branch("mom_Lxy",		&tree_.mom_Lxy);
  tree->Branch("mom_Lz",		&tree_.mom_Lz);
  tree->Branch("mom_Lxyz",		&tree_.mom_Lxyz);
  tree->Branch("mom_ctau",		&tree_.mom_ctau);
}




//---------------------------------------------------------------------------------------------------
// ------------ initialize trees ------------
void DisplacedJetsAnalyzer::initTreeStructure()
{
  
  tree_.sample		= -500.;
  tree_.run		= -500.;
  tree_.lumi		= -500.;
  tree_.event		= 0.;
  tree_.ngenjets	= -500.;
  tree_.ngenpart	= -500.;

}

//---------------------------------------------------------------------------------------------------
// ------------ method called once each job just after ending the event loop  ------------
void DisplacedJetsAnalyzer::endJob() 
{
}

//---------------------------------------------------------------------------------------------------
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DisplacedJetsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//---------------------------------------------------------------------------------------------------
//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedJetsAnalyzer);
