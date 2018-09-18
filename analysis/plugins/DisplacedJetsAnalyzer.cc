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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Math/interface/deltaR.h"
// root include files
#include "TTree.h"

//---------------------------------------------------------------------------------------------------
// ttree structure

struct tree_struc_{

  int				sample;
  int				run;
  int				lumi;
  long unsigned int		event;
  float				weight;
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
  std::vector<int>		genpar_match_q0;
  std::vector<int>		genpar_match_q1;
  std::vector<int>		genpar_match_q2;
  std::vector<int>		genpar_match_q3;
  std::vector<int>		genpar_match_q4;
  std::vector<int>		genpar_match_q1lo;
  std::vector<int>		genpar_match_q2lo;
  std::vector<int>		genpar_match_q3lo;
  std::vector<int>		genpar_match_q4lo;
  std::vector<float>		genpar_dv1;
  std::vector<float>		genpar_dv2;
  std::vector<float>		genpar_dv3;
  std::vector<float>		genpar_dv4;
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
  std::vector<int>		genpar_stat;
  std::vector<float>		genpar_pt;
  std::vector<float>		genpar_eta;
  std::vector<float>		genpar_phi;
  std::vector<float>		genpar_theta;
  std::vector<float>		genpar_vx;
  std::vector<float>		genpar_vy;
  std::vector<float>		genpar_vz;
  std::vector<float>		genpar_Lxy;
  std::vector<float>		genpar_Lxyz;
  std::vector<float>		genpar_la;
  std::vector<float>		genpar_lo;
  std::vector<float>		genpar_laline;
  std::vector<float>		genpar_loline;
  std::vector<float>		genpar_beta;
  std::vector<float>		genpar_q;
  int				nmothers;
  std::vector<int>		mom_id;
  std::vector<int>		mom_stat;
  std::vector<float>		mom_e;
  std::vector<float>		mom_m;
  std::vector<float>		mom_p;
  std::vector<float>		mom_pt;
  std::vector<float>		mom_eta;
  std::vector<float>		mom_phi;
  std::vector<float>		mom_theta;
  std::vector<float>    	mom_vx;
  std::vector<float>		mom_vy;
  std::vector<float>		mom_vz;
  std::vector<float>		mom_beta;
  std::vector<float>		mom_gama;
  std::vector<float>		mom_Lxy;
  std::vector<float>		mom_Lz;
  std::vector<float>		mom_Lxyz;
  std::vector<float>		mom_ctau;
  std::vector<int>		mom_dupl;

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
      void clearVectors(); 
      float getXYZ(float r, float Q, float x0, float y0, float z0, float px0, float py0, float pz0, float & x, float &y, float &z);
      float getLtrue(float s, float z);
      float getL(float x, float y, float z);

      // ----------member data ---------------------------
      edm::Service<TFileService> fs;
     
      bool	verbose_;
      int	sampleID_;

      // tokens
      edm::EDGetTokenT<GenEventInfoProduct>			genInfoToken_;
      edm::EDGetTokenT<std::vector<reco::GenJet> >		genjetToken_;
      edm::EDGetTokenT<std::vector<reco::GenParticle> >		genparticleToken_;
      edm::EDGetTokenT<edm::View<PileupSummaryInfo> >		pileupToken_;
      edm::EDGetTokenT<std::vector<SimVertex> >			verticesToken_;

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
   genInfoToken_	= consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generatorInfo"));
   genjetToken_		= consumes<std::vector<reco::GenJet> >(iConfig.getUntrackedParameter<edm::InputTag>("genjets"));
   genparticleToken_    = consumes<std::vector<reco::GenParticle> >(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"));
   pileupToken_		= consumes<edm::View<PileupSummaryInfo> >(iConfig.getUntrackedParameter<edm::InputTag>("pileupInfo"));
   verticesToken_	= consumes<std::vector<SimVertex> >(iConfig.getUntrackedParameter<edm::InputTag>("vertices"));

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
  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(genInfoToken_, genInfo);  
 
  edm::Handle<std::vector<reco::GenJet> > genjets;
  iEvent.getByToken(genjetToken_, genjets);

  edm::Handle<std::vector<reco::GenParticle> > genparticles;
  iEvent.getByToken(genparticleToken_, genparticles);

  edm::Handle<edm::View<PileupSummaryInfo> > pileupInfo;
  iEvent.getByToken(pileupToken_, pileupInfo);

  edm::Handle<std::vector<SimVertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);

  // --- general event info 
  unsigned long int event = iEvent.id().event();   
  int run                 = iEvent.id().run();
  int lumi                = iEvent.luminosityBlock();
  float weight            = genInfo->weights()[0]; 

  //int nvtx                = -1;   // samples not generated with PU
  //for (unsigned int PV = 0; PV < pileupInfo->size(); ++PV){
  //  int BX = pileupInfo->ptrAt(PV)->getBunchCrossing();
  //  if (BX == 0) nvtx = pileupInfo->ptrAt(PV)->getTrueNumInteractions();
  //}

  // --- genjet info
  int ngenjets = 0;
  std::vector<float>	genjet_pt;
  std::vector<float>	genjet_e;
  std::vector<float>	genjet_eta;
  std::vector<float>	genjet_phi;
  std::vector<int>  	genjet_ndaug;
  std::vector<int>  	genjet_nconst;
  std::vector<float>	genjet_vx;
  std::vector<float>	genjet_vy;
  std::vector<float>	genjet_vz;

  // --- genjet constituent info
  std::vector<int>  	genjet_i;
  std::vector<int>  	genjet_match;
  std::vector<float>	genjet0_const_st;
  std::vector<float>	genjet0_const_id;
  std::vector<float>	genjet0_const_pt;
  std::vector<float>	genjet0_const_vx;
  std::vector<float>	genjet0_const_vy;
  std::vector<float>	genjet0_const_vz;
  std::vector<float>	genjet1_const_st;
  std::vector<float>	genjet1_const_id;
  std::vector<float>	genjet1_const_pt;
  std::vector<float>	genjet1_const_vx;
  std::vector<float>	genjet1_const_vy;
  std::vector<float>	genjet1_const_vz;
  std::vector<float>	genjet2_const_st;
  std::vector<float>	genjet2_const_id;
  std::vector<float>	genjet2_const_pt;
  std::vector<float>	genjet2_const_vx;
  std::vector<float>	genjet2_const_vy;
  std::vector<float>	genjet2_const_vz;
  std::vector<float>	genjet3_const_st;
  std::vector<float>	genjet3_const_id;
  std::vector<float>	genjet3_const_pt;
  std::vector<float>	genjet3_const_vx;
  std::vector<float>	genjet3_const_vy;
  std::vector<float>	genjet3_const_vz;

  // --- genparticles info
  int ngenpart = 0;
  std::vector<int>	genpar_id;
  std::vector<int>	genpar_stat;
  std::vector<float>	genpar_e;
  std::vector<float>	genpar_pt;
  std::vector<float>	genpar_eta;
  std::vector<float>	genpar_phi;
  std::vector<float>	genpar_theta;
  std::vector<float>	genpar_vx;
  std::vector<float>	genpar_vy;
  std::vector<float>	genpar_vz;
  std::vector<float>	genpar_Lxy;
  std::vector<float>	genpar_Lxyz;
  std::vector<float>	genpar_la;
  std::vector<float>	genpar_lo;
  std::vector<float>	genpar_laline;
  std::vector<float>	genpar_loline;
  std::vector<float>	genpar_beta;
  std::vector<float>	genpar_q;
  int nmothers = 0;
  std::vector<int>	mom_id;
  std::vector<int>	mom_stat;
  std::vector<float>	mom_e;
  std::vector<float>	mom_m;
  std::vector<float>	mom_p;
  std::vector<float>	mom_pt;
  std::vector<float>	mom_eta;
  std::vector<float>	mom_phi;
  std::vector<float>	mom_theta;
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
  std::vector<int>	match_q0;
  std::vector<int>	match_q1;
  std::vector<int>	match_q2;
  std::vector<int>	match_q3;
  std::vector<int>	match_q4;
  std::vector<int>	mom_dupl;
  std::vector<int>	match_q1lo;
  std::vector<int>	match_q2lo;
  std::vector<int>	match_q3lo;
  std::vector<int>	match_q4lo;
  std::vector<float>	deltav1;
  std::vector<float>	deltav2;
  std::vector<float>	deltav3;
  std::vector<float>	deltav4;

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

  if (genparticles.isValid()){ // make sure have genparticles collection

    // setup variables to store eta and phi of the quarks of interest
    float q1_eta  = -10000;
    float q2_eta  = -10000;
    float q3_eta  = -10000;
    float q4_eta  = -10000;
    float q1_phi  = -10000;
    float q2_phi  = -10000;
    float q3_phi  = -10000;
    float q4_phi  = -10000;
    // setup variables to store vertex of quarks of interest
    float q1_vx   = -10000;
    float q1_vy   = -10000;
    float q1_vz   = -10000;
    float q2_vx   = -10000;
    float q2_vy   = -10000;
    float q2_vz   = -10000;
    float q3_vx   = -10000;
    float q3_vy   = -10000;
    float q3_vz   = -10000;
    float q4_vx   = -10000;
    float q4_vy   = -10000;
    float q4_vz   = -10000;

    int mom35 = 0; 
    int mom36 = 0;
    float q1_mx   = -10000;
    float q1_my   = -10000;
    float q1_mz   = -10000;

    // store the 4 quarks that are produced in the hard interaction
    int interestingjet = 0;
    for (const auto & genpar_iter : *genparticles){ // loop over genparticles

      // only final quarks of hard interaction 
      if ( genpar_iter.status() != 23 ){
         quarknum.push_back(0); 
         continue;
      }

      //if (genpar_iter.mother(0) != NULL){
      //  std::cout << "Genpar: " << genpar_iter.pdgId() << " " << genpar_iter.status() << " " << genpar_iter.mother(0)->pdgId() << " " << genpar_iter.mother(0)->status() << std::endl;
      //  std::cout << "------- v2: " << genpar_iter.vertex().x() << " " << genpar_iter.vertex().y() << " " << genpar_iter.vertex().z() << std::endl;
      //  std::cout << "------- m1: " << genpar_iter.mother(0)->vx() << " " << genpar_iter.mother(0)->vy() << " " << genpar_iter.mother(0)->vz() << std::endl;
      //}
 
      // save interesting quarks
      interestingjet++;
      quarknum.push_back(interestingjet);

      //std::cout << "    ----- particle " << genpar_iter.pdgId() << std::endl;
      //std::cout << "    ----- momref   " << genpar_iter.mother(0)->pdgId() << std::endl;

      float vx = genpar_iter.vertex().x();
      float vy = genpar_iter.vertex().y();
      float vz = genpar_iter.vertex().z();

      if (interestingjet==1){ q1_eta = genpar_iter.eta(); q1_phi = genpar_iter.phi(); q1_vx = vx; q1_vy = vy; q1_vz = vz;} 
      if (interestingjet==2){ q2_eta = genpar_iter.eta(); q2_phi = genpar_iter.phi(); q2_vx = vx; q2_vy = vy; q2_vz = vz;} 
      if (interestingjet==3){ q3_eta = genpar_iter.eta(); q3_phi = genpar_iter.phi(); q3_vx = vx; q3_vy = vy; q3_vz = vz;} 
      if (interestingjet==4){ q4_eta = genpar_iter.eta(); q4_phi = genpar_iter.phi(); q4_vx = vx; q4_vy = vy; q4_vz = vz;}

      // initialize mom info
      float mx            = -10000;
      float my            = -10000;
      float mz            = -10000;
      float dx            = -10000;
      float dy            = -10000;
      float dz            = -10000;
      float mBeta         = -10000; 
      float mGama         = -10000; 
      float mLxy          = -10000; 
      float mLxyz         = -10000; 
      float mcTau         = -10000; 
      int   tmp_mom_id    = -10000;
      int   tmp_mom_stat  = -10000;
      float tmp_mom_e     = -10000; 
      float tmp_mom_m     = -10000;
      float tmp_mom_p     = -10000;
      float tmp_mom_pt    = -10000;
      float tmp_mom_eta   = -10000;
      float tmp_mom_phi   = -10000;
      float tmp_mom_theta = -10000;
      float tmp_mom_vx    = -10000;
      float tmp_mom_vy    = -10000;
      float tmp_mom_vz    = -10000;
      float tmp_mom_beta  = -10000;
      float tmp_mom_gama  = -10000;
      float tmp_mom_Lxy   = -10000;
      float tmp_mom_Lz    = -10000;
      float tmp_mom_Lxyz  = -10000;
      float tmp_mom_ctau  = -10000;
      float tmp_mom_dupl  = 0;
 
      // now get interesting quarks' mother or gmother information
      std::vector< const reco::Candidate * > mother;
      mother.push_back(genpar_iter.mother(0));
      nmothers = mother.size();

      // now look at mothers/grandmothers
      for (const auto & mom : mother ){
        mx = mom->vx();
        my = mom->vy();
        mz = mom->vz();
        dx = vx - mx;
        dy = vy - my;
        dz = vz - mz;
        mBeta = mom->p()/mom->energy();    // mom beta
        mGama = mom->energy()/mom->mass(); // mom gamma
        mLxy  = std::sqrt(dx*dx + dy*dy); 
        mLxyz = std::sqrt(dx*dx + dy*dy + dz*dz); 
        mcTau = std::sqrt(dx*dx + dy*dy + dz*dz) / (mBeta * mGama); 

        tmp_mom_id  	= mom->pdgId();
        tmp_mom_stat 	= mom->status();
        tmp_mom_e 	= mom->energy();
        tmp_mom_m 	= mom->mass();
        tmp_mom_p 	= mom->p();
        tmp_mom_pt 	= mom->pt();
        tmp_mom_eta 	= mom->eta();
        tmp_mom_phi 	= mom->phi();
        tmp_mom_theta   = mom->theta();
        tmp_mom_vx 	= mom->vx();
        tmp_mom_vy 	= mom->vy();
        tmp_mom_vz 	= mom->vz();
        tmp_mom_beta 	= mBeta;
        tmp_mom_gama 	= mGama;
        tmp_mom_Lxy 	= mLxy;
        tmp_mom_Lz 	= dz;
        tmp_mom_Lxyz 	= mLxyz;
        tmp_mom_ctau 	= mcTau;

        // mom ID for original samples
        //if (std::abs(mom->pdgId())==35) mom35++;
        //if (std::abs(mom->pdgId())==36) mom36++;
        //if (std::abs(mom->pdgId())==35 && mom35>1) tmp_mom_dupl = 1;
        //if (std::abs(mom->pdgId())==36 && mom36>1) tmp_mom_dupl = 1;

        // mom ID for new samples 
        if (std::abs(mom->pdgId())==9000006) mom35++;  
        if (std::abs(mom->pdgId())==9000007) mom36++;
        if (std::abs(mom->pdgId())==9000006 && mom35>1) tmp_mom_dupl = 1;
        if (std::abs(mom->pdgId())==9000007 && mom36>1) tmp_mom_dupl = 1;
        if (interestingjet==1){ q1_mx = mx; q1_my = my; q1_mz = mz;} 
        if (interestingjet==2){ q1_mx = mx; q1_my = my; q1_mz = mz;} 
        if (interestingjet==3){ q1_mx = mx; q1_my = my; q1_mz = mz;} 
        if (interestingjet==4){ q1_mx = mx; q1_my = my; q1_mz = mz;}
      }

      //std::cout << " q1: " << q1_mx << " " << q1_mx << " " << q1_mx << std::endl;
      //std::cout << " q2: " << q2_mx << " " << q2_mx << " " << q2_mx << std::endl;
      //std::cout << " q3: " << q3_mx << " " << q3_mx << " " << q3_mx << std::endl;
      //std::cout << " q4: " << q4_mx << " " << q4_mx << " " << q4_mx << std::endl;

      // store only one mom info per gen particle 
      mom_id.push_back(tmp_mom_id);
      mom_stat.push_back(tmp_mom_stat);
      mom_e.push_back(tmp_mom_e);
      mom_m.push_back(tmp_mom_m);
      mom_p.push_back(tmp_mom_p);
      mom_pt.push_back(tmp_mom_pt);
      mom_eta.push_back(tmp_mom_eta);
      mom_phi.push_back(tmp_mom_phi);
      mom_theta.push_back(tmp_mom_theta);
      mom_vx.push_back(tmp_mom_vx);
      mom_vy.push_back(tmp_mom_vy);
      mom_vz.push_back(tmp_mom_vz);
      mom_beta.push_back(tmp_mom_beta);
      mom_gama.push_back(tmp_mom_gama);
      mom_Lxy.push_back(tmp_mom_Lxy);
      mom_Lz.push_back(tmp_mom_Lz);
      mom_Lxyz.push_back(tmp_mom_Lxyz);
      mom_ctau.push_back(tmp_mom_ctau);
      mom_dupl.push_back(tmp_mom_dupl);
 
    }// end loop over genparticles (for only 4 quarks of interest)

    for (const auto & genpar_iter : *genparticles){ // loop over genparticles

      // status: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
      // * 21-29: particles of the hardest process (23 outgoing)
      // * 1: final state particles
      //if (genpar_iter.status() < 21 && genpar_iter.status() > 29 && genpar_iter.status() != 1) continue;
      //if (genpar_iter.status() != 23) continue; 
      ngenpart++;  // number stable particles

      // variables for timing measurements
      float xT  = -10000;
      float yT  = -10000;
      float zT  = -10000;

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
      genpar_theta.push_back(genpar_iter.theta());
      genpar_vx.push_back(vx);
      genpar_vy.push_back(vy);
      genpar_vz.push_back(vz);
      genpar_Lxy.push_back(gLxy);
      genpar_Lxyz.push_back(gLxyz);

      // matching to hard scatter quarks
      float dR1 = 10; float dR2 = 10; float dR3 = 10; float dR4 = 10;
      float dv1 = 0;  float dv2 = 0;  float dv3 = 0;  float dv4 = 0;
      int match0 = 0; int match1 = 0; int match2 = 0; int match3 = 0; int match4 = 0;
      int match1lo = 0; int match2lo = 0; int match3lo = 0; int match4lo = 0; 

      // dR matching
      dR1 = deltaR(q1_eta,q1_phi,genpar_iter.eta(),genpar_iter.phi()); 
      dR2 = deltaR(q2_eta,q2_phi,genpar_iter.eta(),genpar_iter.phi()); 
      dR3 = deltaR(q3_eta,q3_phi,genpar_iter.eta(),genpar_iter.phi()); 
      dR4 = deltaR(q4_eta,q4_phi,genpar_iter.eta(),genpar_iter.phi());
      // vertex matching
      dv1 = std::fabs(q1_vx-vx)+std::fabs(q1_vy-vy)+std::fabs(q1_vz-vz);    
      dv2 = std::fabs(q2_vx-vx)+std::fabs(q2_vy-vy)+std::fabs(q2_vz-vz);    
      dv3 = std::fabs(q3_vx-vx)+std::fabs(q3_vy-vy)+std::fabs(q3_vz-vz);    
      dv4 = std::fabs(q4_vx-vx)+std::fabs(q4_vy-vy)+std::fabs(q4_vz-vz);    
 
      if (q1_eta < 1.5 && dR1 < 0.4 && dv1 < 10) match1 = 1;
      if (q2_eta < 1.5 && dR2 < 0.4 && dv2 < 10) match2 = 1;
      if (q3_eta < 1.5 && dR3 < 0.4 && dv3 < 10) match3 = 1;
      if (q4_eta < 1.5 && dR4 < 0.4 && dv4 < 10) match4 = 1;
      if (genpar_iter.eta() < 1.5 && dR1 > 0.8 && dR2 > 0.8 && dR3 > 0.8 && dR4 > 0.8) match0 = 1;

      // loose matching (no requirement on matching the SV) -- for use in calculating alpha_jet_PV
      if (q1_eta < 1.5 && dR1 < 0.4) match1lo = 1;
      if (q2_eta < 1.5 && dR2 < 0.4) match2lo = 1;
      if (q3_eta < 1.5 && dR3 < 0.4) match3lo = 1;
      if (q4_eta < 1.5 && dR4 < 0.4) match4lo = 1;

      // if particles match more than one "jet"
      if ( match1 + match2 + match3 + match4 > 1){
        float min_12 = std::min(dR1,dR2);
        float min_34 = std::min(dR3,dR4);
        float min_dR = std::min(min_12,min_34);
        if (dR1 == min_dR){ match1 = 1; match2 = 0; match3 = 0; match4 = 0; } // dR1 is smallest
        if (dR2 == min_dR){ match1 = 0; match2 = 1; match3 = 0; match4 = 0; } // dR2 is smallest
        if (dR3 == min_dR){ match1 = 0; match2 = 0; match3 = 1; match4 = 0; } // dR3 is smallest
        if (dR4 == min_dR){ match1 = 0; match2 = 0; match3 = 0; match4 = 1; } // dR4 is smallest
      }

      // save
      match_q0.push_back(match0);
      match_q1.push_back(match1);
      match_q2.push_back(match2);
      match_q3.push_back(match3);
      match_q4.push_back(match4);
  
      match_q1lo.push_back(match1lo);
      match_q2lo.push_back(match2lo);
      match_q3lo.push_back(match3lo);
      match_q4lo.push_back(match4lo);

      deltav1.push_back(dv1);
      deltav2.push_back(dv2);
      deltav3.push_back(dv3);
      deltav4.push_back(dv4);

      // calculating timing delay
      float px0 = genpar_iter.px();
      float py0 = genpar_iter.py();
      float pz0 = genpar_iter.pz();
      float q   = genpar_iter.charge();
      float s = getXYZ(116.10,q,vx,vy,vz,px0,py0,pz0,xT,yT,zT); // find arclength & xyz coord of point on timing layer
      float loline = getL(xT-q1_mx,yT-q1_my,zT-q1_mz);          // distance from PV to point on timing layer
      float lTline = getL(xT-vx,yT-vy,zT-vz);                   // distance from SV to point on timing layer
      float lo = getLtrue(s,zT-q1_mz);                          // distance from PV to point on timing layer
      float lT = getLtrue(s,zT-vz);                             // distance from SV to point on timing layer

      genpar_lo.push_back(lo);
      genpar_la.push_back(lT);
      genpar_loline.push_back(loline);
      genpar_laline.push_back(lTline);
      genpar_beta.push_back(genpar_iter.p()/genpar_iter.energy());
      genpar_q.push_back(q);
 
    }// end loop over genparticles
  }// end if genparticles.isValid
  else std::cout << "WARNING: genparticles collection is NOT valid" << std::endl;

  // --- setup tree values
  initTreeStructure();
  // --- clear all vectors
  clearVectors();

  tree_.sample		= sampleID_;
  tree_.run		= run;
  tree_.lumi		= lumi;
  tree_.event		= event;
  tree_.weight		= weight;
  tree_.ngenjets	= ngenjets;
  tree_.ngenpart        = ngenpart;
  tree_.nmothers	= nmothers;

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

  if (verbose_) std::cout << "mom_id    " << mom_id.size()    << std::endl; 
  if (verbose_) std::cout << "mom_stat  " << mom_stat.size()  << std::endl; 
  if (verbose_) std::cout << "mom_e     " << mom_e.size()     << std::endl;
  if (verbose_) std::cout << "mom_p     " << mom_p.size()     << std::endl; 
  if (verbose_) std::cout << "mom_m     " << mom_m.size()     << std::endl;  
  if (verbose_) std::cout << "mom_pt    " << mom_pt.size()    << std::endl;  
  if (verbose_) std::cout << "mom_eta   " << mom_eta.size()   << std::endl;  
  if (verbose_) std::cout << "mom_phi   " << mom_phi.size()   << std::endl;
  if (verbose_) std::cout << "mom_theta " << mom_theta.size() << std::endl; 
  if (verbose_) std::cout << "mom_vx    " << mom_vx.size()    << std::endl;  
  if (verbose_) std::cout << "mom_vy    " << mom_vy.size()    << std::endl; 
  if (verbose_) std::cout << "mom_vz    " << mom_vz.size()    << std::endl; 
  if (verbose_) std::cout << "mom_beta  " << mom_beta.size()  << std::endl;  
  if (verbose_) std::cout << "mom_gama  " << mom_gama.size()  << std::endl;  
  if (verbose_) std::cout << "mom_Lxy   " << mom_Lxy.size()   << std::endl;  
  if (verbose_) std::cout << "mom_Lz    " << mom_Lz.size()    << std::endl; 
  if (verbose_) std::cout << "mom_Lxyz  " << mom_Lxyz.size()  << std::endl;  
  if (verbose_) std::cout << "mom_ctau  " << mom_ctau.size()  << std::endl;  

  for (int ip = 0; ip < ngenpart; ip++){
    tree_.genpar_id.push_back(genpar_id[ip]);
    tree_.genpar_stat.push_back(genpar_stat[ip]);
    tree_.genpar_pt.push_back(genpar_pt[ip]);
    tree_.genpar_eta.push_back(genpar_eta[ip]);
    tree_.genpar_phi.push_back(genpar_phi[ip]);
    tree_.genpar_theta.push_back(genpar_theta[ip]);
    tree_.genpar_vx.push_back(genpar_vx[ip]);
    tree_.genpar_vy.push_back(genpar_vy[ip]);
    tree_.genpar_vz.push_back(genpar_vz[ip]);
    tree_.genpar_Lxy.push_back(genpar_Lxy[ip]);
    tree_.genpar_Lxyz.push_back(genpar_Lxyz[ip]);
    tree_.genpar_qnum.push_back(quarknum[ip]);
    tree_.genpar_match_q0.push_back(match_q0[ip]);
    tree_.genpar_match_q1.push_back(match_q1[ip]);
    tree_.genpar_match_q2.push_back(match_q2[ip]);
    tree_.genpar_match_q3.push_back(match_q3[ip]);
    tree_.genpar_match_q4.push_back(match_q4[ip]);
    tree_.genpar_dv1.push_back(deltav1[ip]);
    tree_.genpar_dv2.push_back(deltav2[ip]);
    tree_.genpar_dv3.push_back(deltav3[ip]);
    tree_.genpar_dv4.push_back(deltav4[ip]);
    tree_.genpar_la.push_back(genpar_la[ip]);
    tree_.genpar_lo.push_back(genpar_lo[ip]);
    tree_.genpar_laline.push_back(genpar_laline[ip]);
    tree_.genpar_loline.push_back(genpar_loline[ip]);
    tree_.genpar_beta.push_back(genpar_beta[ip]);
    tree_.genpar_q.push_back(genpar_q[ip]);
    tree_.genpar_match_q1lo.push_back(match_q1lo[ip]);
    tree_.genpar_match_q2lo.push_back(match_q2lo[ip]);
    tree_.genpar_match_q3lo.push_back(match_q3lo[ip]);
    tree_.genpar_match_q4lo.push_back(match_q4lo[ip]);
  }

  // save mom info for the 4 quarks of interest
  for (int ip = 0; ip < 4; ip++){
    tree_.mom_id.push_back(mom_id[ip]);
    tree_.mom_stat.push_back(mom_stat[ip]);
    tree_.mom_e.push_back(mom_e[ip]);
    tree_.mom_p.push_back(mom_p[ip]);
    tree_.mom_m.push_back(mom_m[ip]);
    tree_.mom_pt.push_back(mom_pt[ip]);
    tree_.mom_eta.push_back(mom_eta[ip]);
    tree_.mom_phi.push_back(mom_phi[ip]);
    tree_.mom_theta.push_back(mom_theta[ip]);
    tree_.mom_vx.push_back(mom_vx[ip]);
    tree_.mom_vy.push_back(mom_vy[ip]);
    tree_.mom_vz.push_back(mom_vz[ip]);
    tree_.mom_beta.push_back(mom_beta[ip]);
    tree_.mom_gama.push_back(mom_gama[ip]);
    tree_.mom_Lxy.push_back(mom_Lxy[ip]);
    tree_.mom_Lz.push_back(mom_Lz[ip]);
    tree_.mom_Lxyz.push_back(mom_Lxyz[ip]);
    tree_.mom_ctau.push_back(mom_ctau[ip]);
    tree_.mom_dupl.push_back(mom_dupl[ip]);
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
  tree->Branch("weight",		&tree_.weight,			"weight/F");

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
  tree->Branch("genpar_stat",		&tree_.genpar_stat);
  tree->Branch("genpar_pt",		&tree_.genpar_pt);
  tree->Branch("genpar_eta",		&tree_.genpar_eta);
  tree->Branch("genpar_phi",		&tree_.genpar_phi);
  tree->Branch("genpar_theta",		&tree_.genpar_theta);
  tree->Branch("genpar_vx",		&tree_.genpar_vx);
  tree->Branch("genpar_vy",		&tree_.genpar_vy);
  tree->Branch("genpar_vz",		&tree_.genpar_vz);
  tree->Branch("genpar_Lxy",		&tree_.genpar_Lxy);
  tree->Branch("genpar_Lxyz",		&tree_.genpar_Lxyz);
  tree->Branch("genpar_qnum",		&tree_.genpar_qnum);
  tree->Branch("genpar_match_q0",	&tree_.genpar_match_q0);
  tree->Branch("genpar_match_q1",	&tree_.genpar_match_q1);
  tree->Branch("genpar_match_q2",	&tree_.genpar_match_q2);
  tree->Branch("genpar_match_q3",	&tree_.genpar_match_q3);
  tree->Branch("genpar_match_q4",	&tree_.genpar_match_q4);
  tree->Branch("genpar_dv1",		&tree_.genpar_dv1);
  tree->Branch("genpar_dv2",		&tree_.genpar_dv2);
  tree->Branch("genpar_dv3",		&tree_.genpar_dv3);
  tree->Branch("genpar_dv4",		&tree_.genpar_dv4);
  tree->Branch("genpar_lo",		&tree_.genpar_lo);
  tree->Branch("genpar_la",		&tree_.genpar_la);
  tree->Branch("genpar_loline",		&tree_.genpar_loline);
  tree->Branch("genpar_laline",		&tree_.genpar_laline);
  tree->Branch("genpar_beta",		&tree_.genpar_beta);
  tree->Branch("genpar_q",		&tree_.genpar_q);
  tree->Branch("genpar_match_q1lo",	&tree_.genpar_match_q1lo);
  tree->Branch("genpar_match_q2lo",	&tree_.genpar_match_q2lo);
  tree->Branch("genpar_match_q3lo",	&tree_.genpar_match_q3lo);
  tree->Branch("genpar_match_q4lo",	&tree_.genpar_match_q4lo);

  // gen particle mom stuff
  tree->Branch("nmothers",		&tree_.nmothers,		"nmothers/I");
  tree->Branch("mom_id",		&tree_.mom_id);
  tree->Branch("mom_stat",		&tree_.mom_stat);
  tree->Branch("mom_e",			&tree_.mom_e);
  tree->Branch("mom_m",			&tree_.mom_m);
  tree->Branch("mom_p",			&tree_.mom_p);
  tree->Branch("mom_pt",		&tree_.mom_pt);
  tree->Branch("mom_eta",		&tree_.mom_eta);
  tree->Branch("mom_phi",		&tree_.mom_phi);
  tree->Branch("mom_theta",		&tree_.mom_theta);
  tree->Branch("mom_vx",		&tree_.mom_vx);
  tree->Branch("mom_vy",		&tree_.mom_vy);
  tree->Branch("mom_vz",		&tree_.mom_vz);
  tree->Branch("mom_beta",		&tree_.mom_beta);
  tree->Branch("mom_gama",		&tree_.mom_gama);
  tree->Branch("mom_Lxy",		&tree_.mom_Lxy);
  tree->Branch("mom_Lz",		&tree_.mom_Lz);
  tree->Branch("mom_Lxyz",		&tree_.mom_Lxyz);
  tree->Branch("mom_ctau",		&tree_.mom_ctau);
  tree->Branch("mom_dupl",		&tree_.mom_dupl);
}
//---------------------------------------------------------------------------------------------------
float DisplacedJetsAnalyzer::getXYZ(float r, float Q, float x0, float y0, float z0, float px0, float py0, float pz0, float & x, float &y, float &z)
{
  // constants
  float Bfield = 3.8; // [T]
  float a      = -0.002998*Bfield*Q;           // Paul avery's convention 

  // track params
  float p0     = std::sqrt(px0*px0 + py0*py0); // pt
  float rho    = a/p0;                         // inverse radius of curvature [cm]          
  float dxy    = (-x0*py0 + y0*px0)/p0;        //  
  //std::cout << "dxy: " << dxy << " "  << std::sqrt(x0*x0 + y0*y0) << std::endl; 
  float d0     = -dxy;                         // CMS definition of d0
  float lam    = pz0/p0;                       // 

  // useful variables
  float c      = rho/2.0;
  float B      = c*std::sqrt((r*r - d0*d0)/(1 + 2*c*d0));

  // arc length
  float s      = TMath::ASin(B)*2/rho;

  // xyz coord at r
  float sin_rho_s = TMath::Sin(rho*s);
  float cos_rho_s = TMath::Cos(rho*s);
  x = x0 + (px0/a)*sin_rho_s - (py0/a)*(1-cos_rho_s);
  y = y0 + (py0/a)*sin_rho_s + (px0/a)*(1-cos_rho_s);
  z = z0 + lam*s;

  //float d      = std::sqrt( (x-x0)*(x-x0) + (y-y0)*(y-y0) ); // straight line length
  //float R      = 1/rho; // radius of curvature
  //float salt   = R*TMath::ACos((2*R*R - d*d)/(2*R*R)); // sagita
  //if (s < d) std::cout << " NOOOO S smaller than D " << std::endl; 
  //std::cout << "s: " << s << " " << salt << " vs. d: " << d << std::endl;

  return s;
}

float DisplacedJetsAnalyzer::getLtrue(float s, float z)
{
  float L = std::sqrt(s*s + z*z);
  return L;
}

float DisplacedJetsAnalyzer::getL(float x, float y, float z)
{
  float L = std::sqrt(x*x + y*y + z*z);
  return L;
}

//---------------------------------------------------------------------------------------------------
// ------------ initialize trees ------------
void DisplacedJetsAnalyzer::initTreeStructure()
{
  
  tree_.sample		= -500.;
  tree_.run		= -500.;
  tree_.lumi		= -500.;
  tree_.event		= 0.;
  tree_.weight		= 0.;
  tree_.ngenjets	= -500.;
  tree_.ngenpart	= -500.;
  tree_.nmothers	= -500.;
}

// ------------ clear vectors in trees ------
void DisplacedJetsAnalyzer::clearVectors()
{

  tree_.genjet_pt.clear();
  tree_.genjet_e.clear();
  tree_.genjet_eta.clear();
  tree_.genjet_phi.clear();
  tree_.genjet_ndaug.clear();
  tree_.genjet_nconst.clear();
  tree_.genjet_vx.clear();
  tree_.genjet_vy.clear();
  tree_.genjet_vz.clear();
  tree_.genpar_qnum.clear();
  tree_.genpar_match_q0.clear();
  tree_.genpar_match_q1.clear();
  tree_.genpar_match_q2.clear();
  tree_.genpar_match_q3.clear();
  tree_.genpar_match_q4.clear();
  tree_.genpar_match_q1lo.clear();
  tree_.genpar_match_q2lo.clear();
  tree_.genpar_match_q3lo.clear();
  tree_.genpar_match_q4lo.clear();
  tree_.genpar_dv1.clear();
  tree_.genpar_dv2.clear();
  tree_.genpar_dv3.clear();
  tree_.genpar_dv4.clear();
  tree_.genjet_i.clear();
  tree_.genjet_match.clear();
  tree_.genjet0_const_st.clear();
  tree_.genjet0_const_id.clear();
  tree_.genjet0_const_pt.clear();
  tree_.genjet0_const_vx.clear();
  tree_.genjet0_const_vy.clear();
  tree_.genjet0_const_vz.clear();
  tree_.genjet1_const_st.clear();
  tree_.genjet1_const_id.clear();
  tree_.genjet1_const_pt.clear();
  tree_.genjet1_const_vx.clear();
  tree_.genjet1_const_vy.clear();
  tree_.genjet1_const_vz.clear();
  tree_.genjet2_const_st.clear();
  tree_.genjet2_const_id.clear();
  tree_.genjet2_const_pt.clear();
  tree_.genjet2_const_vx.clear();
  tree_.genjet2_const_vy.clear();
  tree_.genjet2_const_vz.clear();
  tree_.genjet3_const_st.clear();
  tree_.genjet3_const_id.clear();
  tree_.genjet3_const_pt.clear();
  tree_.genjet3_const_vx.clear();
  tree_.genjet3_const_vy.clear();
  tree_.genjet3_const_vz.clear();
  tree_.genpar_id.clear();
  tree_.genpar_stat.clear();
  tree_.genpar_pt.clear();
  tree_.genpar_eta.clear();
  tree_.genpar_phi.clear();
  tree_.genpar_theta.clear();
  tree_.genpar_vx.clear();
  tree_.genpar_vy.clear();
  tree_.genpar_vz.clear();
  tree_.genpar_Lxy.clear();
  tree_.genpar_Lxyz.clear();
  tree_.genpar_la.clear();
  tree_.genpar_lo.clear();
  tree_.genpar_laline.clear();
  tree_.genpar_loline.clear();
  tree_.genpar_beta.clear();
  tree_.genpar_q.clear();
  tree_.mom_id.clear();
  tree_.mom_stat.clear();
  tree_.mom_e.clear();
  tree_.mom_m.clear();
  tree_.mom_p.clear();
  tree_.mom_pt.clear();
  tree_.mom_eta.clear();
  tree_.mom_phi.clear();
  tree_.mom_theta.clear();
  tree_.mom_vx.clear();
  tree_.mom_vy.clear();
  tree_.mom_vz.clear();
  tree_.mom_beta.clear();
  tree_.mom_gama.clear();
  tree_.mom_Lxy.clear();
  tree_.mom_Lz.clear();
  tree_.mom_Lxyz.clear();
  tree_.mom_ctau.clear();
  tree_.mom_dupl.clear();

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
