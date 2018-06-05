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
#include "DataFormats/Math/interface/deltaR.h"
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
  std::vector<int>		genpar_qnum;
  std::vector<int>		genpar_match_q1;
  std::vector<int>		genpar_match_q2;
  std::vector<int>		genpar_match_q3;
  std::vector<int>		genpar_match_q4;
  std::vector<float>		genjet_i;
  std::vector<int>		genjet_match;
  std::vector<float>		genjet0_const_st;
  std::vector<float>		genjet0_const_id;
  std::vector<float>		genjet0_const_pt;
  std::vector<float>		genjet0_const_vx;
  std::vector<float>		genjet0_const_vy;
  std::vector<float>		genjet0_const_vz;
  std::vector<float>		genjet1_const_st;
  std::vector<float>		genjet1_const_id;
  std::vector<float>		genjet1_const_pt;
  std::vector<float>		genjet1_const_vx;
  std::vector<float>		genjet1_const_vy;
  std::vector<float>		genjet1_const_vz;
  std::vector<float>		genjet2_const_st;
  std::vector<float>		genjet2_const_id;
  std::vector<float>		genjet2_const_pt;
  std::vector<float>		genjet2_const_vx;
  std::vector<float>		genjet2_const_vy;
  std::vector<float>		genjet2_const_vz;
  std::vector<float>		genjet3_const_st;
  std::vector<float>		genjet3_const_id;
  std::vector<float>		genjet3_const_pt;
  std::vector<float>		genjet3_const_vx;
  std::vector<float>		genjet3_const_vy;
  std::vector<float>		genjet3_const_vz;
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
  std::vector<float>		mom_p;
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

  // --- genjet constituent info
  std::vector<int>              genjet_i;
  std::vector<int>		genjet_match;
  std::vector<float>            genjet0_const_st;
  std::vector<float>            genjet0_const_id;
  std::vector<float>		genjet0_const_pt;
  std::vector<float>		genjet0_const_vx;
  std::vector<float>		genjet0_const_vy;
  std::vector<float>		genjet0_const_vz;
  std::vector<float>		genjet1_const_st;
  std::vector<float>		genjet1_const_id;
  std::vector<float>		genjet1_const_pt;
  std::vector<float>		genjet1_const_vx;
  std::vector<float>		genjet1_const_vy;
  std::vector<float>		genjet1_const_vz;
  std::vector<float>		genjet2_const_st;
  std::vector<float>		genjet2_const_id;
  std::vector<float>		genjet2_const_pt;
  std::vector<float>		genjet2_const_vx;
  std::vector<float>		genjet2_const_vy;
  std::vector<float>		genjet2_const_vz;
  std::vector<float>		genjet3_const_st;
  std::vector<float>		genjet3_const_id;
  std::vector<float>		genjet3_const_pt;
  std::vector<float>		genjet3_const_vx;
  std::vector<float>		genjet3_const_vy;
  std::vector<float>		genjet3_const_vz;

  if (genjets.isValid()){ // make sure have genjet collection
    ngenjets = genjets->size();
 
    int jetiter = 0;
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
      genjet_i.push_back(jetiter); // note: probably not necessary
   
      // constituent information
      std::vector<const reco::GenParticle*> genjet_const = genjet_iter.getGenConstituents();
      genjet_nconst.push_back(genjet_const.size()); 

      int match_LL = 0; // does a constituent come from the hard interaction

      // manually get jet constituents 
      if (genparticles.isValid()){ // make sure have genparticles collection
        for (const auto & genpar_iter : *genparticles){ // loop over genparticles

          // remove non-final state particles
          //if (genpar_iter.status()!=1) continue;         

          // check if deltaR(particle,jet) < 0.4 
          float dR = deltaR(genjet_iter.eta(),genjet_iter.phi(),genpar_iter.eta(),genpar_iter.phi()); 
          if (dR > 0.4) continue; 

          // check if there is a LL particle matching to the jet
          if (genpar_iter.status()==23) match_LL = 1;

          // save genparticle constituents associated with each jet
          if (jetiter==0){
            genjet0_const_st.push_back(genpar_iter.status());
            genjet0_const_id.push_back(genpar_iter.pdgId());
            genjet0_const_pt.push_back(genpar_iter.pt());
            genjet0_const_vx.push_back(genpar_iter.vertex().x());
            genjet0_const_vy.push_back(genpar_iter.vertex().y());
            genjet0_const_vz.push_back(genpar_iter.vertex().z());
          }
          if (jetiter==1){
            genjet1_const_st.push_back(genpar_iter.status());
            genjet1_const_id.push_back(genpar_iter.pdgId());
            genjet1_const_pt.push_back(genpar_iter.pt());
            genjet1_const_vx.push_back(genpar_iter.vertex().x());
            genjet1_const_vy.push_back(genpar_iter.vertex().y());
            genjet1_const_vz.push_back(genpar_iter.vertex().z());
          }
          if (jetiter==2){
            genjet2_const_st.push_back(genpar_iter.status());
            genjet2_const_id.push_back(genpar_iter.pdgId());
            genjet2_const_pt.push_back(genpar_iter.pt());
            genjet2_const_vx.push_back(genpar_iter.vertex().x());
            genjet2_const_vy.push_back(genpar_iter.vertex().y());
            genjet2_const_vz.push_back(genpar_iter.vertex().z());
          }
          if (jetiter==3){
            genjet3_const_st.push_back(genpar_iter.status());
            genjet3_const_id.push_back(genpar_iter.pdgId());
            genjet3_const_pt.push_back(genpar_iter.pt());
            genjet3_const_vx.push_back(genpar_iter.vertex().x());
            genjet3_const_vy.push_back(genpar_iter.vertex().y());
            genjet3_const_vz.push_back(genpar_iter.vertex().z());
          }

        }// end loop over genparticles
      }// end if genparticles.isValid

      genjet_match.push_back(match_LL);
      jetiter++; // count genjet
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
  std::vector<float>	mom_p;
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
  std::vector<int>	quarknum;
  std::vector<int>	match_q1;
  std::vector<int>	match_q2;
  std::vector<int>	match_q3;
  std::vector<int>	match_q4;

  if (genparticles.isValid()){ // make sure have genparticles collection

    float q1_eta = -10000;
    float q2_eta = -10000;
    float q3_eta = -10000;
    float q4_eta = -10000;
    float q1_phi = -10000;
    float q2_phi = -10000;
    float q3_phi = -10000;
    float q4_phi = -10000;
    int   momID  = -10000;
    int   gmomID = -10000;
    bool  findmom = false;


    // store the 4 quarks that are produced in the hard interaction
    int interestingjet = 0;
    for (const auto & genpar_iter : *genparticles){ // loop over genparticles
      // check that the mom or gmom is from the resonance
      if (genpar_iter.mother() != NULL){
        momID   = genpar_iter.mother()->pdgId();
        if (genpar_iter.mother()->mother() != NULL) gmomID  = genpar_iter.mother()->mother()->pdgId();
      }
      findmom = (momID == 35 || momID == 36 || gmomID == 35 || gmomID == 36);
      // only final quarks of hard interaction w/ found mom are interesting 
      if (genpar_iter.status()!=23 || !findmom){ 
         quarknum.push_back(0); 
         continue;
      }
      // save interesting quarks
      interestingjet++;
      quarknum.push_back(interestingjet);
      if (interestingjet==1){ q1_eta = genpar_iter.eta(); q1_phi = genpar_iter.phi(); } 
      if (interestingjet==2){ q2_eta = genpar_iter.eta(); q2_phi = genpar_iter.phi(); } 
      if (interestingjet==3){ q3_eta = genpar_iter.eta(); q3_phi = genpar_iter.phi(); } 
      if (interestingjet==4){ q4_eta = genpar_iter.eta(); q4_phi = genpar_iter.phi(); } 
    }// end loop over genparticles

    std::cout << "here" << std::endl;
    for (const auto & genpar_iter : *genparticles){ // loop over genparticles

      // status: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
      // * 21-29: particles of the hardest process (23 outgoing)
      // * 1: final state particles
      if (genpar_iter.status() < 21 && genpar_iter.status() > 29 && genpar_iter.status() != 1) continue;
      //if (genpar_iter.status() != 23) continue; 
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

      // matching to hard scatter quarks 
      float dR1 = deltaR(q1_eta,q1_phi,genpar_iter.eta(),genpar_iter.phi()); 
      float dR2 = deltaR(q2_eta,q2_phi,genpar_iter.eta(),genpar_iter.phi()); 
      float dR3 = deltaR(q3_eta,q3_phi,genpar_iter.eta(),genpar_iter.phi()); 
      float dR4 = deltaR(q4_eta,q4_phi,genpar_iter.eta(),genpar_iter.phi()); 
      if (dR1 < 0.4) match_q1.push_back(1);
      else           match_q1.push_back(0);
      if (dR2 < 0.4) match_q2.push_back(1);
      else           match_q2.push_back(0);
      if (dR3 < 0.4) match_q3.push_back(1);
      else           match_q3.push_back(0);
      if (dR4 < 0.4) match_q4.push_back(1);
      else           match_q4.push_back(0);

      int n = 0;
      float mx = -10000;
      float my = -10000;
      float mz = -10000;
      float dx = -10000;
      float dy = -10000;
      float dz = -10000;
      float mBeta = -10000; 
      float mGama = -10000; 
      float mLxy  = -10000; 
      float mLxyz = -10000; 
      float mcTau = -10000; 

      // initialize mom info (even for gen particles w/o a mom) 
      int	tmp_mom_id  	= -10000;
      int	tmp_mom_stat 	= -10000;
      float	tmp_mom_e 	= -10000;
      float	tmp_mom_m 	= -10000;
      float	tmp_mom_p 	= -10000;
      float	tmp_mom_pt 	= -10000;
      float	tmp_mom_eta 	= -10000;
      float	tmp_mom_phi 	= -10000;
      float	tmp_mom_vx 	= -10000;
      float	tmp_mom_vy 	= -10000;
      float	tmp_mom_vz 	= -10000;
      float	tmp_mom_beta 	= -10000;
      float	tmp_mom_gama 	= -10000;
      float	tmp_mom_Lxy 	= -10000;
      float	tmp_mom_Lz 	= -10000;
      float	tmp_mom_Lxyz 	= -10000;
      float	tmp_mom_ctau 	= -10000;

      //std::cout << "particle: " << genpar_iter.pdgId() << std::endl;
      //if (genpar_iter.mother() != NULL){
      //  std::vector< const reco::Candidate * > mother;
      //  mother.push_back(genpar_iter.mother());
      //  for (const auto & mom : *genpar_iter.mother()){
      //    std::cout << " - mom:  "  << mom.pdgId()     << std::endl;
      //    std::vector< const reco::Candidate * > gmother;
      //    gmother.push_back(genpar_iter.mother()->mother());
      //    for (const auto & gmom : *genpar_iter.mother()->mother()){
      //      std::cout << " -- gmom: " << gmom.pdgId() << " (status: " << gmom.status() << ")" << std::endl;
      //    } 
      //  }      
      //}

      if (genpar_iter.mother() != NULL){ 
        std::vector< const reco::Candidate * > mother;
        mother.push_back(genpar_iter.mother());
        nmothers = mother.size();
        if (verbose_) std::cout << "Number of mothers: " << nmothers << std::endl;
        std::sort(mother.begin(),mother.end(),pTsort);

        for (const auto & mom : *genpar_iter.mother()){

          if (abs(mom.pdgId()) != 35 && abs(mom.pdgId()) != 36) continue;
          //if (n != 0) continue; // only keep highest pt mother
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

          tmp_mom_id  	= mom.pdgId();
          tmp_mom_stat 	= mom.status();
          tmp_mom_e 	= mom.mass();
          tmp_mom_m 	= mom.energy();
          tmp_mom_p 	= mom.p();
          tmp_mom_pt 	= mom.pt();
          tmp_mom_eta 	= mom.eta();
          tmp_mom_phi 	= mom.phi();
          tmp_mom_vx 	= mx;
          tmp_mom_vy 	= my;
          tmp_mom_vz 	= mz;
          tmp_mom_beta 	= mBeta;
          tmp_mom_gama 	= mGama;
          tmp_mom_Lxy 	= dz;
          tmp_mom_Lz 	= mLxy;
          tmp_mom_Lxyz 	= mLxyz;
          tmp_mom_ctau 	= mcTau;
          
        }// end loop over genparticle mothers
      }// end if genparticles.mother.isValid

      // store only one mom info per gen particle 
      mom_id.push_back(tmp_mom_id);
      mom_stat.push_back(tmp_mom_stat);
      mom_e.push_back(tmp_mom_e);
      mom_m.push_back(tmp_mom_m);
      mom_p.push_back(tmp_mom_p);
      mom_pt.push_back(tmp_mom_pt);
      mom_eta.push_back(tmp_mom_eta);
      mom_phi.push_back(tmp_mom_phi);
      mom_vx.push_back(tmp_mom_vx);
      mom_vy.push_back(tmp_mom_vy);
      mom_vz.push_back(tmp_mom_vz);
      mom_beta.push_back(tmp_mom_beta);
      mom_gama.push_back(tmp_mom_gama);
      mom_Lxy.push_back(tmp_mom_Lxy);
      mom_Lz.push_back(tmp_mom_Lz);
      mom_Lxyz.push_back(tmp_mom_Lxyz);
      mom_ctau.push_back(tmp_mom_ctau);

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
    tree_.genjet_i.push_back(genjet_i[ij]);
    tree_.genjet_match.push_back(genjet_match[ij]);
  }
  for (unsigned int i = 0; i < genjet0_const_id.size(); i++){
    tree_.genjet0_const_st.push_back(genjet0_const_st[i]);
    tree_.genjet0_const_id.push_back(genjet0_const_id[i]);
    tree_.genjet0_const_pt.push_back(genjet0_const_pt[i]);
    tree_.genjet0_const_vx.push_back(genjet0_const_vx[i]);
    tree_.genjet0_const_vy.push_back(genjet0_const_vy[i]);
    tree_.genjet0_const_vz.push_back(genjet0_const_vz[i]);
  } 
  for (unsigned int i = 0; i < genjet1_const_id.size(); i++){
    tree_.genjet1_const_st.push_back(genjet1_const_st[i]);
    tree_.genjet1_const_id.push_back(genjet1_const_id[i]);
    tree_.genjet1_const_pt.push_back(genjet1_const_pt[i]);
    tree_.genjet1_const_vx.push_back(genjet1_const_vx[i]);
    tree_.genjet1_const_vy.push_back(genjet1_const_vy[i]);
    tree_.genjet1_const_vz.push_back(genjet1_const_vz[i]);
  } 
  for (unsigned int i = 0; i < genjet2_const_id.size(); i++){
    tree_.genjet2_const_st.push_back(genjet2_const_st[i]);
    tree_.genjet2_const_id.push_back(genjet2_const_id[i]);
    tree_.genjet2_const_pt.push_back(genjet2_const_pt[i]);
    tree_.genjet2_const_vx.push_back(genjet2_const_vx[i]);
    tree_.genjet2_const_vy.push_back(genjet2_const_vy[i]);
    tree_.genjet2_const_vz.push_back(genjet2_const_vz[i]);
  } 
  for (unsigned int i = 0; i < genjet3_const_id.size(); i++){
    tree_.genjet3_const_st.push_back(genjet3_const_st[i]);
    tree_.genjet3_const_id.push_back(genjet3_const_id[i]);
    tree_.genjet3_const_pt.push_back(genjet3_const_pt[i]);
    tree_.genjet3_const_vx.push_back(genjet3_const_vx[i]);
    tree_.genjet3_const_vy.push_back(genjet3_const_vy[i]);
    tree_.genjet3_const_vz.push_back(genjet3_const_vz[i]);
  } 

  if (verbose_) std::cout << "mom_id   " << mom_id.size()    << std::endl; 
  if (verbose_) std::cout << "mom_stat " << mom_stat.size()  << std::endl; 
  if (verbose_) std::cout << "mom_e    " << mom_e.size()     << std::endl;
  if (verbose_) std::cout << "mom_p    " << mom_p.size()     << std::endl; 
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
    tree_.genpar_qnum.push_back(quarknum[ip]);
    tree_.genpar_match_q1.push_back(match_q1[ip]);
    tree_.genpar_match_q2.push_back(match_q2[ip]);
    tree_.genpar_match_q3.push_back(match_q3[ip]);
    tree_.genpar_match_q4.push_back(match_q4[ip]);
    tree_.mom_id.push_back(mom_id[ip]);
    tree_.mom_stat.push_back(mom_stat[ip]);
    tree_.mom_e.push_back(mom_e[ip]);
    tree_.mom_p.push_back(mom_p[ip]);
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
  tree->Branch("sample",		&tree_.sample,			"sample/I");
  tree->Branch("run",			&tree_.run,			"run/I");
  tree->Branch("lumi",			&tree_.lumi,			"lumi/I");
  tree->Branch("event",			&tree_.event,			"event/L");

  // gen jet stuff
  tree->Branch("ngenjets",		&tree_.ngenjets,		"ngenjets/I");
  tree->Branch("genjet_pt",     	&tree_.genjet_pt);
  tree->Branch("genjet_e",      	&tree_.genjet_e);
  tree->Branch("genjet_eta",		&tree_.genjet_eta);
  tree->Branch("genjet_phi",		&tree_.genjet_phi);
  tree->Branch("genjet_ndaug",		&tree_.genjet_ndaug);
  tree->Branch("genjet_nconst",		&tree_.genjet_nconst);
  tree->Branch("genjet_vx",		&tree_.genjet_vx);
  tree->Branch("genjet_vy",		&tree_.genjet_vy);
  tree->Branch("genjet_vz",		&tree_.genjet_vz);
  tree->Branch("genjet_i",		&tree_.genjet_i);
  tree->Branch("genjet_match",		&tree_.genjet_match);

  // gen jet-particle matching stuff
  tree->Branch("genjet0_const_st",	&tree_.genjet0_const_st);
  tree->Branch("genjet0_const_id",	&tree_.genjet0_const_id);
  tree->Branch("genjet0_const_pt",	&tree_.genjet0_const_pt);
  tree->Branch("genjet0_const_vx",	&tree_.genjet0_const_vx);
  tree->Branch("genjet0_const_vy",	&tree_.genjet0_const_vy);
  tree->Branch("genjet0_const_vz",	&tree_.genjet0_const_vz);
  tree->Branch("genjet1_const_st",	&tree_.genjet1_const_st);
  tree->Branch("genjet1_const_id",	&tree_.genjet1_const_id);
  tree->Branch("genjet1_const_pt",	&tree_.genjet1_const_pt);
  tree->Branch("genjet1_const_vx",	&tree_.genjet1_const_vx);
  tree->Branch("genjet1_const_vy",	&tree_.genjet1_const_vy);
  tree->Branch("genjet1_const_vz",	&tree_.genjet1_const_vz);
  tree->Branch("genjet2_const_st",	&tree_.genjet2_const_st);
  tree->Branch("genjet2_const_id",	&tree_.genjet2_const_id);
  tree->Branch("genjet2_const_pt",	&tree_.genjet2_const_pt);
  tree->Branch("genjet2_const_vx",	&tree_.genjet2_const_vx);
  tree->Branch("genjet2_const_vy",	&tree_.genjet2_const_vy);
  tree->Branch("genjet2_const_vz",	&tree_.genjet2_const_vz);
  tree->Branch("genjet3_const_st",	&tree_.genjet3_const_st);
  tree->Branch("genjet3_const_id",	&tree_.genjet3_const_id);
  tree->Branch("genjet3_const_pt",	&tree_.genjet3_const_pt);
  tree->Branch("genjet3_const_vx",	&tree_.genjet3_const_vx);
  tree->Branch("genjet3_const_vy",	&tree_.genjet3_const_vy);
  tree->Branch("genjet3_const_vz",	&tree_.genjet3_const_vz);

  // gen particle stuff
  tree->Branch("ngenpart",		&tree_.ngenpart,		"ngenpart/I");
  tree->Branch("genpar_id",		&tree_.genpar_id);
  tree->Branch("genpar_pt",		&tree_.genpar_pt);
  tree->Branch("genpar_eta",		&tree_.genpar_eta);
  tree->Branch("genpar_phi",		&tree_.genpar_phi);
  tree->Branch("genpar_vx",		&tree_.genpar_vx);
  tree->Branch("genpar_vy",		&tree_.genpar_vy);
  tree->Branch("genpar_vz",		&tree_.genpar_vz);
  tree->Branch("genpar_Lxy",		&tree_.genpar_Lxy);
  tree->Branch("genpar_Lxyz",		&tree_.genpar_Lxyz);
  tree->Branch("genpar_qnum",		&tree_.genpar_qnum);
  tree->Branch("genpar_match_q1",	&tree_.genpar_match_q1);
  tree->Branch("genpar_match_q2",	&tree_.genpar_match_q2);
  tree->Branch("genpar_match_q3",	&tree_.genpar_match_q3);
  tree->Branch("genpar_match_q4",	&tree_.genpar_match_q4);

  // gen particle mom stuff
  tree->Branch("mom_id",		&tree_.mom_id);
  tree->Branch("mom_stat",		&tree_.mom_stat);
  tree->Branch("mom_e",			&tree_.mom_e);
  tree->Branch("mom_m",			&tree_.mom_m);
  tree->Branch("mom_p",			&tree_.mom_p);
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
