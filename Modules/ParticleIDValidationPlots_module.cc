////////////////////////////////////////////////////////////////////////
// Class:       ParticleIDValidationPlots
// Plugin Type: analyzer (art v2_05_01)
// File:        ParticleIDValidationPlots_module.cc
//
// Generated at Thu Mar 15 14:40:38 2018 by Kirsty Duffy using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "uboone/ParticleID/Algorithms/WellReconstructedTrackFinder.h"

//#include "uboone/MyClasses/BackTrackerTruthMatch.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "TH1F.h"
#include "TH2F.h"

class ParticleIDValidationPlots;


class ParticleIDValidationPlots : public art::EDAnalyzer {
public:
  explicit ParticleIDValidationPlots(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ParticleIDValidationPlots(ParticleIDValidationPlots const &) = delete;
  ParticleIDValidationPlots(ParticleIDValidationPlots &&) = delete;
  ParticleIDValidationPlots & operator = (ParticleIDValidationPlots const &) = delete;
  ParticleIDValidationPlots & operator = (ParticleIDValidationPlots &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  art::ServiceHandle<art::TFileService> tfs;

  std::string fTrackingAlgo;
  std::string fHitAlgo;
  std::string fHitTrackAssns;
  std::string fTruthMatchingAssns;
  std::string fPIDtag;

  TH1F *WR_truemu_neglogl_mu;
  TH1F *WR_truep_neglogl_mu;
  TH1F *WR_truepi_neglogl_mu;
  TH1F *WR_trueK_neglogl_mu;
  TH1F *WR_truemu_neglogl_p;
  TH1F *WR_truep_neglogl_p;
  TH1F *WR_truepi_neglogl_p;
  TH1F *WR_trueK_neglogl_p;
  TH1F *WR_truemu_neglogl_pi;
  TH1F *WR_truep_neglogl_pi;
  TH1F *WR_truepi_neglogl_pi;
  TH1F *WR_trueK_neglogl_pi;
  TH1F *WR_truemu_neglogl_K;
  TH1F *WR_truep_neglogl_K;
  TH1F *WR_truepi_neglogl_K;
  TH1F *WR_trueK_neglogl_K;

  TH1F *WR_truemu_pull_mu;
  TH1F *WR_truep_pull_mu;
  TH1F *WR_truepi_pull_mu;
  TH1F *WR_trueK_pull_mu;
  TH1F *WR_truemu_pull_p;
  TH1F *WR_truep_pull_p;
  TH1F *WR_truepi_pull_p;
  TH1F *WR_trueK_pull_p;
  TH1F *WR_truemu_pull_pi;
  TH1F *WR_truep_pull_pi;
  TH1F *WR_truepi_pull_pi;
  TH1F *WR_trueK_pull_pi;
  TH1F *WR_truemu_pull_K;
  TH1F *WR_truep_pull_K;
  TH1F *WR_truepi_pull_K;
  TH1F *WR_trueK_pull_K;

  TH1F *WR_truemu_pull_MIP;
  TH1F *WR_truep_pull_MIP;
  TH1F *WR_truepi_pull_MIP;
  TH1F *WR_trueK_pull_MIP;

  TH1F *WR_truemu_PIDA;
  TH1F *WR_truep_PIDA;
  TH1F *WR_truepi_PIDA;
  TH1F *WR_trueK_PIDA;

  TH2F *WR_truemu_dEdxtr_len;
  TH2F *WR_truep_dEdxtr_len;
  TH2F *WR_truepi_dEdxtr_len;
  TH2F *WR_trueK_dEdxtr_len;
  TH2F *WR_truemu_dQdxtr_len;
  TH2F *WR_truep_dQdxtr_len;
  TH2F *WR_truepi_dQdxtr_len;
  TH2F *WR_trueK_dQdxtr_len;

  TH1F *All_Assn_truemu_neglogl_mu;
  TH1F *All_Assn_truep_neglogl_mu;
  TH1F *All_Assn_truepi_neglogl_mu;
  TH1F *All_Assn_trueK_neglogl_mu;
  TH1F *All_Assn_truemu_neglogl_p;
  TH1F *All_Assn_truep_neglogl_p;
  TH1F *All_Assn_truepi_neglogl_p;
  TH1F *All_Assn_trueK_neglogl_p;
  TH1F *All_Assn_truemu_neglogl_pi;
  TH1F *All_Assn_truep_neglogl_pi;
  TH1F *All_Assn_truepi_neglogl_pi;
  TH1F *All_Assn_trueK_neglogl_pi;
  TH1F *All_Assn_truemu_neglogl_K;
  TH1F *All_Assn_truep_neglogl_K;
  TH1F *All_Assn_truepi_neglogl_K;
  TH1F *All_Assn_trueK_neglogl_K;

  TH1F *All_Assn_truemu_pull_mu;
  TH1F *All_Assn_truep_pull_mu;
  TH1F *All_Assn_truepi_pull_mu;
  TH1F *All_Assn_trueK_pull_mu;
  TH1F *All_Assn_truemu_pull_p;
  TH1F *All_Assn_truep_pull_p;
  TH1F *All_Assn_truepi_pull_p;
  TH1F *All_Assn_trueK_pull_p;
  TH1F *All_Assn_truemu_pull_pi;
  TH1F *All_Assn_truep_pull_pi;
  TH1F *All_Assn_truepi_pull_pi;
  TH1F *All_Assn_trueK_pull_pi;
  TH1F *All_Assn_truemu_pull_K;
  TH1F *All_Assn_truep_pull_K;
  TH1F *All_Assn_truepi_pull_K;
  TH1F *All_Assn_trueK_pull_K;

  TH1F *All_Assn_truemu_pull_MIP;
  TH1F *All_Assn_truep_pull_MIP;
  TH1F *All_Assn_truepi_pull_MIP;
  TH1F *All_Assn_trueK_pull_MIP;

  TH1F *All_Assn_truemu_PIDA;
  TH1F *All_Assn_truep_PIDA;
  TH1F *All_Assn_truepi_PIDA;
  TH1F *All_Assn_trueK_PIDA;

  TH2F *All_Assn_truemu_dEdxtr_len;
  TH2F *All_Assn_truep_dEdxtr_len;
  TH2F *All_Assn_truepi_dEdxtr_len;
  TH2F *All_Assn_trueK_dEdxtr_len;
  TH2F *All_Assn_truemu_dQdxtr_len;
  TH2F *All_Assn_truep_dQdxtr_len;
  TH2F *All_Assn_truepi_dQdxtr_len;
  TH2F *All_Assn_trueK_dQdxtr_len;

  /*TH1F *All_AfroBT_truemu_neglogl_mu;
  TH1F *All_AfroBT_truep_neglogl_mu;
  TH1F *All_AfroBT_truepi_neglogl_mu;
  TH1F *All_AfroBT_trueK_neglogl_mu;
  TH1F *All_AfroBT_truemu_neglogl_p;
  TH1F *All_AfroBT_truep_neglogl_p;
  TH1F *All_AfroBT_truepi_neglogl_p;
  TH1F *All_AfroBT_trueK_neglogl_p;
  TH1F *All_AfroBT_truemu_neglogl_pi;
  TH1F *All_AfroBT_truep_neglogl_pi;
  TH1F *All_AfroBT_truepi_neglogl_pi;
  TH1F *All_AfroBT_trueK_neglogl_pi;
  TH1F *All_AfroBT_truemu_neglogl_K;
  TH1F *All_AfroBT_truep_neglogl_K;
  TH1F *All_AfroBT_truepi_neglogl_K;
  TH1F *All_AfroBT_trueK_neglogl_K;

  TH1F *All_AfroBT_truemu_pull_mu;
  TH1F *All_AfroBT_truep_pull_mu;
  TH1F *All_AfroBT_truepi_pull_mu;
  TH1F *All_AfroBT_trueK_pull_mu;
  TH1F *All_AfroBT_truemu_pull_p;
  TH1F *All_AfroBT_truep_pull_p;
  TH1F *All_AfroBT_truepi_pull_p;
  TH1F *All_AfroBT_trueK_pull_p;
  TH1F *All_AfroBT_truemu_pull_pi;
  TH1F *All_AfroBT_truep_pull_pi;
  TH1F *All_AfroBT_truepi_pull_pi;
  TH1F *All_AfroBT_trueK_pull_pi;
  TH1F *All_AfroBT_truemu_pull_K;
  TH1F *All_AfroBT_truep_pull_K;
  TH1F *All_AfroBT_truepi_pull_K;
  TH1F *All_AfroBT_trueK_pull_K;

  TH1F *All_AfroBT_truemu_pull_MIP;
  TH1F *All_AfroBT_truep_pull_MIP;
  TH1F *All_AfroBT_truepi_pull_MIP;
  TH1F *All_AfroBT_trueK_pull_MIP;

  TH1F *All_AfroBT_truemu_PIDA;
  TH1F *All_AfroBT_truep_PIDA;
  TH1F *All_AfroBT_truepi_PIDA;
  TH1F *All_AfroBT_trueK_PIDA;

  TH2F *All_AfroBT_truemu_dEdxtr_len;
  TH2F *All_AfroBT_truep_dEdxtr_len;
  TH2F *All_AfroBT_truepi_dEdxtr_len;
  TH2F *All_AfroBT_trueK_dEdxtr_len;
  TH2F *All_AfroBT_truemu_dQdxtr_len;
  TH2F *All_AfroBT_truep_dQdxtr_len;
  TH2F *All_AfroBT_truepi_dQdxtr_len;
  TH2F *All_AfroBT_trueK_dQdxtr_len;*/
};


ParticleIDValidationPlots::ParticleIDValidationPlots(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fTrackingAlgo = p.get<std::string>("TrackingAlgorithm","pandoraNu::McRecoStage2");
  fHitAlgo = p.get<std::string>("HitProducer","pandoraCosmicHitRemoval::McRecoStage2");
  fHitTrackAssns = p.get<std::string>("HitTrackAssnName","pandoraNu::McRecoStage2");
  fTruthMatchingAssns = p.get<std::string>("HitTruthMatchingAssnName","crHitRemovalTruthMatch::McRecoStage2");
  fPIDtag = p.get<std::string>("ParticleIDProducerModule");
  
  // ---- Well reconstructed
  WR_truemu_neglogl_mu = tfs->make<TH1F>("WR_truemu_neglogl_mu","Well reconstructed tracks, true muons;neg2LL_mu;",200,0,200);
  WR_truep_neglogl_mu = tfs->make<TH1F>("WR_truep_neglogl_mu","Well reconstructed tracks, true protons;neg2LL_mu;",200,0,200);
  WR_truepi_neglogl_mu = tfs->make<TH1F>("WR_truepi_neglogl_mu","Well reconstructed tracks, true pions;neg2LL_mu;",200,0,200);
  WR_trueK_neglogl_mu = tfs->make<TH1F>("WR_trueK_neglogl_mu","Well reconstructed tracks, true kaons;neg2LL_mu;",200,0,200);
  WR_truemu_neglogl_p = tfs->make<TH1F>("WR_truemu_neglogl_p","Well reconstructed tracks, true muons;neg2LL_p;",200,0,200);
  WR_truep_neglogl_p = tfs->make<TH1F>("WR_truep_neglogl_p","Well reconstructed tracks, true protons;neg2LL_p;",200,0,200);
  WR_truepi_neglogl_p = tfs->make<TH1F>("WR_truepi_neglogl_p","Well reconstructed tracks, true pions;neg2LL_p;",200,0,200);
  WR_trueK_neglogl_p = tfs->make<TH1F>("WR_trueK_neglogl_p","Well reconstructed tracks, true kaons;neg2LL_p;",200,0,200);
  WR_truemu_neglogl_pi = tfs->make<TH1F>("WR_truemu_neglogl_pi","Well reconstructed tracks, true muons;neg2LL_pi;",200,0,200);
  WR_truep_neglogl_pi = tfs->make<TH1F>("WR_truep_neglogl_pi","Well reconstructed tracks, true protons;neg2LL_pi;",200,0,200);
  WR_truepi_neglogl_pi = tfs->make<TH1F>("WR_truepi_neglogl_pi","Well reconstructed tracks, true pions;neg2LL_pi;",200,0,200);
  WR_trueK_neglogl_pi = tfs->make<TH1F>("WR_trueK_neglogl_pi","Well reconstructed tracks, true kaons;neg2LL_pi;",200,0,200);
  WR_truemu_neglogl_K = tfs->make<TH1F>("WR_truemu_neglogl_K","Well reconstructed tracks, true muons;neg2LL_K;",200,0,200);
  WR_truep_neglogl_K = tfs->make<TH1F>("WR_truep_neglogl_K","Well reconstructed tracks, true protons;neg2LL_K;",200,0,200);
  WR_truepi_neglogl_K = tfs->make<TH1F>("WR_truepi_neglogl_K","Well reconstructed tracks, true pions;neg2LL_K;",200,0,200);
  WR_trueK_neglogl_K = tfs->make<TH1F>("WR_trueK_neglogl_K","Well reconstructed tracks, true kaons;neg2LL_K;",200,0,200);
  
  WR_truemu_pull_mu = tfs->make<TH1F>("WR_truemu_pull_mu","Well reconstructed tracks, true muons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truep_pull_mu = tfs->make<TH1F>("WR_truep_pull_mu","Well reconstructed tracks, true protons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truepi_pull_mu = tfs->make<TH1F>("WR_truepi_pull_mu","Well reconstructed tracks, true pions;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_trueK_pull_mu = tfs->make<TH1F>("WR_trueK_pull_mu","Well reconstructed tracks, true kaons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truemu_pull_p = tfs->make<TH1F>("WR_truemu_pull_p","Well reconstructed tracks, true muons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truep_pull_p = tfs->make<TH1F>("WR_truep_pull_p","Well reconstructed tracks, true protons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truepi_pull_p = tfs->make<TH1F>("WR_truepi_pull_p","Well reconstructed tracks, true pions;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_trueK_pull_p = tfs->make<TH1F>("WR_trueK_pull_p","Well reconstructed tracks, true kaons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truemu_pull_pi = tfs->make<TH1F>("WR_truemu_pull_pi","Well reconstructed tracks, true muons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truep_pull_pi = tfs->make<TH1F>("WR_truep_pull_pi","Well reconstructed tracks, true protons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truepi_pull_pi = tfs->make<TH1F>("WR_truepi_pull_pi","Well reconstructed tracks, true pions;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_trueK_pull_pi = tfs->make<TH1F>("WR_trueK_pull_pi","Well reconstructed tracks, true kaons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truemu_pull_K = tfs->make<TH1F>("WR_truemu_pull_K","Well reconstructed tracks, true muons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truep_pull_K = tfs->make<TH1F>("WR_truep_pull_K","Well reconstructed tracks, true protons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_truepi_pull_K = tfs->make<TH1F>("WR_truepi_pull_K","Well reconstructed tracks, true pions;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  WR_trueK_pull_K = tfs->make<TH1F>("WR_trueK_pull_K","Well reconstructed tracks, true kaons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);

  WR_truemu_pull_MIP = tfs->make<TH1F>("WR_truemu_pull_MIP","Well reconstructed tracks, true muons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  WR_truep_pull_MIP = tfs->make<TH1F>("WR_truep_pull_MIP","Well reconstructed tracks, true protons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  WR_truepi_pull_MIP = tfs->make<TH1F>("WR_truepi_pull_MIP","Well reconstructed tracks, true pions;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  WR_trueK_pull_MIP = tfs->make<TH1F>("WR_trueK_pull_MIP","Well reconstructed tracks, true kaons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);

  WR_truemu_PIDA = tfs->make<TH1F>("WR_truemu_PIDA","Well reconstructed tracks, true muons;PIDA;",300,0,30);
  WR_truep_PIDA = tfs->make<TH1F>("WR_truep_PIDA","Well reconstructed tracks, true protons;PIDA;",300,0,30);
  WR_truepi_PIDA = tfs->make<TH1F>("WR_truepi_PIDA","Well reconstructed tracks, true pions;PIDA;",300,0,30);
  WR_trueK_PIDA = tfs->make<TH1F>("WR_trueK_PIDA","Well reconstructed tracks, true kaons;PIDA;",300,0,30);

  WR_truemu_dEdxtr_len = tfs->make<TH2F>("WR_truemu_dEdxtr_len","Well reconstructed tracks, true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
  WR_truep_dEdxtr_len = tfs->make<TH2F>("WR_truep_dEdxtr_len","Well reconstructed tracks, true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
  WR_truepi_dEdxtr_len = tfs->make<TH2F>("WR_truepi_dEdxtr_len","Well reconstructed tracks, true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
  WR_trueK_dEdxtr_len = tfs->make<TH2F>("WR_trueK_dEdxtr_len","Well reconstructed tracks, true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);
  WR_truemu_dQdxtr_len = tfs->make<TH2F>("WR_truemu_dQdxtr_len","Well reconstructed tracks, true muons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  WR_truep_dQdxtr_len = tfs->make<TH2F>("WR_truep_dQdxtr_len","Well reconstructed tracks, true protons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  WR_truepi_dQdxtr_len = tfs->make<TH2F>("WR_truepi_dQdxtr_len","Well reconstructed tracks, true pions;Track length (cm);dQ/dx",100,0,700,100,0,500);
  WR_trueK_dQdxtr_len = tfs->make<TH2F>("WR_trueK_dQdxtr_len","Well reconstructed tracks, true kaons;Track length (cm);dQ/dx",100,0,700,100,0,500);



  // ---- All tracks: truth matching using associations
  All_Assn_truemu_neglogl_mu = tfs->make<TH1F>("All_Assn_truemu_neglogl_mu","Well reconstructed tracks (associations truth matching), true muons;neg2LL_mu;",200,0,200);
  All_Assn_truep_neglogl_mu = tfs->make<TH1F>("All_Assn_truep_neglogl_mu","Well reconstructed tracks (associations truth matching), true protons;neg2LL_mu;",200,0,200);
  All_Assn_truepi_neglogl_mu = tfs->make<TH1F>("All_Assn_truepi_neglogl_mu","Well reconstructed tracks (associations truth matching), true pions;neg2LL_mu;",200,0,200);
  All_Assn_trueK_neglogl_mu = tfs->make<TH1F>("All_Assn_trueK_neglogl_mu","Well reconstructed tracks (associations truth matching), true kaons;neg2LL_mu;",200,0,200);
  All_Assn_truemu_neglogl_p = tfs->make<TH1F>("All_Assn_truemu_neglogl_p","Well reconstructed tracks (associations truth matching), true muons;neg2LL_p;",200,0,200);
  All_Assn_truep_neglogl_p = tfs->make<TH1F>("All_Assn_truep_neglogl_p","Well reconstructed tracks (associations truth matching), true protons;neg2LL_p;",200,0,200);
  All_Assn_truepi_neglogl_p = tfs->make<TH1F>("All_Assn_truepi_neglogl_p","Well reconstructed tracks (associations truth matching), true pions;neg2LL_p;",200,0,200);
  All_Assn_trueK_neglogl_p = tfs->make<TH1F>("All_Assn_trueK_neglogl_p","Well reconstructed tracks (associations truth matching), true kaons;neg2LL_p;",200,0,200);
  All_Assn_truemu_neglogl_pi = tfs->make<TH1F>("All_Assn_truemu_neglogl_pi","Well reconstructed tracks (associations truth matching), true muons;neg2LL_pi;",200,0,200);
  All_Assn_truep_neglogl_pi = tfs->make<TH1F>("All_Assn_truep_neglogl_pi","Well reconstructed tracks (associations truth matching), true protons;neg2LL_pi;",200,0,200);
  All_Assn_truepi_neglogl_pi = tfs->make<TH1F>("All_Assn_truepi_neglogl_pi","Well reconstructed tracks (associations truth matching), true pions;neg2LL_pi;",200,0,200);
  All_Assn_trueK_neglogl_pi = tfs->make<TH1F>("All_Assn_trueK_neglogl_pi","Well reconstructed tracks (associations truth matching), true kaons;neg2LL_pi;",200,0,200);
  All_Assn_truemu_neglogl_K = tfs->make<TH1F>("All_Assn_truemu_neglogl_K","Well reconstructed tracks (associations truth matching), true muons;neg2LL_K;",200,0,200);
  All_Assn_truep_neglogl_K = tfs->make<TH1F>("All_Assn_truep_neglogl_K","Well reconstructed tracks (associations truth matching), true protons;neg2LL_K;",200,0,200);
  All_Assn_truepi_neglogl_K = tfs->make<TH1F>("All_Assn_truepi_neglogl_K","Well reconstructed tracks (associations truth matching), true pions;neg2LL_K;",200,0,200);
  All_Assn_trueK_neglogl_K = tfs->make<TH1F>("All_Assn_trueK_neglogl_K","Well reconstructed tracks (associations truth matching), true kaons;neg2LL_K;",200,0,200);
  
  All_Assn_truemu_pull_mu = tfs->make<TH1F>("All_Assn_truemu_pull_mu","Well reconstructed tracks (associations truth matching), true muons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truep_pull_mu = tfs->make<TH1F>("All_Assn_truep_pull_mu","Well reconstructed tracks (associations truth matching), true protons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truepi_pull_mu = tfs->make<TH1F>("All_Assn_truepi_pull_mu","Well reconstructed tracks (associations truth matching), true pions;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_trueK_pull_mu = tfs->make<TH1F>("All_Assn_trueK_pull_mu","Well reconstructed tracks (associations truth matching), true kaons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truemu_pull_p = tfs->make<TH1F>("All_Assn_truemu_pull_p","Well reconstructed tracks (associations truth matching), true muons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truep_pull_p = tfs->make<TH1F>("All_Assn_truep_pull_p","Well reconstructed tracks (associations truth matching), true protons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truepi_pull_p = tfs->make<TH1F>("All_Assn_truepi_pull_p","Well reconstructed tracks (associations truth matching), true pions;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_trueK_pull_p = tfs->make<TH1F>("All_Assn_trueK_pull_p","Well reconstructed tracks (associations truth matching), true kaons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truemu_pull_pi = tfs->make<TH1F>("All_Assn_truemu_pull_pi","Well reconstructed tracks (associations truth matching), true muons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truep_pull_pi = tfs->make<TH1F>("All_Assn_truep_pull_pi","Well reconstructed tracks (associations truth matching), true protons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truepi_pull_pi = tfs->make<TH1F>("All_Assn_truepi_pull_pi","Well reconstructed tracks (associations truth matching), true pions;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_trueK_pull_pi = tfs->make<TH1F>("All_Assn_trueK_pull_pi","Well reconstructed tracks (associations truth matching), true kaons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truemu_pull_K = tfs->make<TH1F>("All_Assn_truemu_pull_K","Well reconstructed tracks (associations truth matching), true muons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truep_pull_K = tfs->make<TH1F>("All_Assn_truep_pull_K","Well reconstructed tracks (associations truth matching), true protons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truepi_pull_K = tfs->make<TH1F>("All_Assn_truepi_pull_K","Well reconstructed tracks (associations truth matching), true pions;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_trueK_pull_K = tfs->make<TH1F>("All_Assn_trueK_pull_K","Well reconstructed tracks (associations truth matching), true kaons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);

  All_Assn_truemu_pull_MIP = tfs->make<TH1F>("All_Assn_truemu_pull_MIP","Well reconstructed tracks (associations truth matching), true muons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_Assn_truep_pull_MIP = tfs->make<TH1F>("All_Assn_truep_pull_MIP","Well reconstructed tracks (associations truth matching), true protons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_Assn_truepi_pull_MIP = tfs->make<TH1F>("All_Assn_truepi_pull_MIP","Well reconstructed tracks (associations truth matching), true pions;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_Assn_trueK_pull_MIP = tfs->make<TH1F>("All_Assn_trueK_pull_MIP","Well reconstructed tracks (associations truth matching), true kaons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);

  All_Assn_truemu_PIDA = tfs->make<TH1F>("All_Assn_truemu_PIDA","Well reconstructed tracks (associations truth matching), true muons;PIDA;",300,0,30);
  All_Assn_truep_PIDA = tfs->make<TH1F>("All_Assn_truep_PIDA","Well reconstructed tracks (associations truth matching), true protons;PIDA;",300,0,30);
  All_Assn_truepi_PIDA = tfs->make<TH1F>("All_Assn_truepi_PIDA","Well reconstructed tracks (associations truth matching), true pions;PIDA;",300,0,30);
  All_Assn_trueK_PIDA = tfs->make<TH1F>("All_Assn_trueK_PIDA","Well reconstructed tracks (associations truth matching), true kaons;PIDA;",300,0,30);

  All_Assn_truemu_dEdxtr_len = tfs->make<TH2F>("All_Assn_truemu_dEdxtr_len","Well reconstructed tracks (associations truth matching), true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_Assn_truep_dEdxtr_len = tfs->make<TH2F>("All_Assn_truep_dEdxtr_len","Well reconstructed tracks (associations truth matching), true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_Assn_truepi_dEdxtr_len = tfs->make<TH2F>("All_Assn_truepi_dEdxtr_len","Well reconstructed tracks (associations truth matching), true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_Assn_trueK_dEdxtr_len = tfs->make<TH2F>("All_Assn_trueK_dEdxtr_len","Well reconstructed tracks (associations truth matching), true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_Assn_truemu_dQdxtr_len = tfs->make<TH2F>("All_Assn_truemu_dQdxtr_len","Well reconstructed tracks (associations truth matching), true muons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_Assn_truep_dQdxtr_len = tfs->make<TH2F>("All_Assn_truep_dQdxtr_len","Well reconstructed tracks (associations truth matching), true protons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_Assn_truepi_dQdxtr_len = tfs->make<TH2F>("All_Assn_truepi_dQdxtr_len","Well reconstructed tracks (associations truth matching), true pions;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_Assn_trueK_dQdxtr_len = tfs->make<TH2F>("All_Assn_trueK_dQdxtr_len","Well reconstructed tracks (associations truth matching), true kaons;Track length (cm);dQ/dx",100,0,700,100,0,500);



  // ---- All tracks: truth matching using associations
  /*All_AfroBT_truemu_neglogl_mu = tfs->make<TH1F>("All_AfroBT_truemu_neglogl_mu","Well reconstructed tracks (Afro's truth matching), true muons;neg2LL_mu;",200,0,200);
  All_AfroBT_truep_neglogl_mu = tfs->make<TH1F>("All_AfroBT_truep_neglogl_mu","Well reconstructed tracks (Afro's truth matching), true protons;neg2LL_mu;",200,0,200);
  All_AfroBT_truepi_neglogl_mu = tfs->make<TH1F>("All_AfroBT_truepi_neglogl_mu","Well reconstructed tracks (Afro's truth matching), true pions;neg2LL_mu;",200,0,200);
  All_AfroBT_trueK_neglogl_mu = tfs->make<TH1F>("All_AfroBT_trueK_neglogl_mu","Well reconstructed tracks (Afro's truth matching), true kaons;neg2LL_mu;",200,0,200);
  All_AfroBT_truemu_neglogl_p = tfs->make<TH1F>("All_AfroBT_truemu_neglogl_p","Well reconstructed tracks (Afro's truth matching), true muons;neg2LL_p;",200,0,200);
  All_AfroBT_truep_neglogl_p = tfs->make<TH1F>("All_AfroBT_truep_neglogl_p","Well reconstructed tracks (Afro's truth matching), true protons;neg2LL_p;",200,0,200);
  All_AfroBT_truepi_neglogl_p = tfs->make<TH1F>("All_AfroBT_truepi_neglogl_p","Well reconstructed tracks (Afro's truth matching), true pions;neg2LL_p;",200,0,200);
  All_AfroBT_trueK_neglogl_p = tfs->make<TH1F>("All_AfroBT_trueK_neglogl_p","Well reconstructed tracks (Afro's truth matching), true kaons;neg2LL_p;",200,0,200);
  All_AfroBT_truemu_neglogl_pi = tfs->make<TH1F>("All_AfroBT_truemu_neglogl_pi","Well reconstructed tracks (Afro's truth matching), true muons;neg2LL_pi;",200,0,200);
  All_AfroBT_truep_neglogl_pi = tfs->make<TH1F>("All_AfroBT_truep_neglogl_pi","Well reconstructed tracks (Afro's truth matching), true protons;neg2LL_pi;",200,0,200);
  All_AfroBT_truepi_neglogl_pi = tfs->make<TH1F>("All_AfroBT_truepi_neglogl_pi","Well reconstructed tracks (Afro's truth matching), true pions;neg2LL_pi;",200,0,200);
  All_AfroBT_trueK_neglogl_pi = tfs->make<TH1F>("All_AfroBT_trueK_neglogl_pi","Well reconstructed tracks (Afro's truth matching), true kaons;neg2LL_pi;",200,0,200);
  All_AfroBT_truemu_neglogl_K = tfs->make<TH1F>("All_AfroBT_truemu_neglogl_K","Well reconstructed tracks (Afro's truth matching), true muons;neg2LL_K;",200,0,200);
  All_AfroBT_truep_neglogl_K = tfs->make<TH1F>("All_AfroBT_truep_neglogl_K","Well reconstructed tracks (Afro's truth matching), true protons;neg2LL_K;",200,0,200);
  All_AfroBT_truepi_neglogl_K = tfs->make<TH1F>("All_AfroBT_truepi_neglogl_K","Well reconstructed tracks (Afro's truth matching), true pions;neg2LL_K;",200,0,200);
  All_AfroBT_trueK_neglogl_K = tfs->make<TH1F>("All_AfroBT_trueK_neglogl_K","Well reconstructed tracks (Afro's truth matching), true kaons;neg2LL_K;",200,0,200);
  
  All_AfroBT_truemu_pull_mu = tfs->make<TH1F>("All_AfroBT_truemu_pull_mu","Well reconstructed tracks (Afro's truth matching), true muons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truep_pull_mu = tfs->make<TH1F>("All_AfroBT_truep_pull_mu","Well reconstructed tracks (Afro's truth matching), true protons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truepi_pull_mu = tfs->make<TH1F>("All_AfroBT_truepi_pull_mu","Well reconstructed tracks (Afro's truth matching), true pions;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_trueK_pull_mu = tfs->make<TH1F>("All_AfroBT_trueK_pull_mu","Well reconstructed tracks (Afro's truth matching), true kaons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truemu_pull_p = tfs->make<TH1F>("All_AfroBT_truemu_pull_p","Well reconstructed tracks (Afro's truth matching), true muons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truep_pull_p = tfs->make<TH1F>("All_AfroBT_truep_pull_p","Well reconstructed tracks (Afro's truth matching), true protons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truepi_pull_p = tfs->make<TH1F>("All_AfroBT_truepi_pull_p","Well reconstructed tracks (Afro's truth matching), true pions;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_trueK_pull_p = tfs->make<TH1F>("All_AfroBT_trueK_pull_p","Well reconstructed tracks (Afro's truth matching), true kaons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truemu_pull_pi = tfs->make<TH1F>("All_AfroBT_truemu_pull_pi","Well reconstructed tracks (Afro's truth matching), true muons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truep_pull_pi = tfs->make<TH1F>("All_AfroBT_truep_pull_pi","Well reconstructed tracks (Afro's truth matching), true protons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truepi_pull_pi = tfs->make<TH1F>("All_AfroBT_truepi_pull_pi","Well reconstructed tracks (Afro's truth matching), true pions;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_trueK_pull_pi = tfs->make<TH1F>("All_AfroBT_trueK_pull_pi","Well reconstructed tracks (Afro's truth matching), true kaons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truemu_pull_K = tfs->make<TH1F>("All_AfroBT_truemu_pull_K","Well reconstructed tracks (Afro's truth matching), true muons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truep_pull_K = tfs->make<TH1F>("All_AfroBT_truep_pull_K","Well reconstructed tracks (Afro's truth matching), true protons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_truepi_pull_K = tfs->make<TH1F>("All_AfroBT_truepi_pull_K","Well reconstructed tracks (Afro's truth matching), true pions;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_AfroBT_trueK_pull_K = tfs->make<TH1F>("All_AfroBT_trueK_pull_K","Well reconstructed tracks (Afro's truth matching), true kaons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);

  All_AfroBT_truemu_pull_MIP = tfs->make<TH1F>("All_AfroBT_truemu_pull_MIP","Well reconstructed tracks (Afro's truth matching), true muons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_AfroBT_truep_pull_MIP = tfs->make<TH1F>("All_AfroBT_truep_pull_MIP","Well reconstructed tracks (Afro's truth matching), true protons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_AfroBT_truepi_pull_MIP = tfs->make<TH1F>("All_AfroBT_truepi_pull_MIP","Well reconstructed tracks (Afro's truth matching), true pions;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_AfroBT_trueK_pull_MIP = tfs->make<TH1F>("All_AfroBT_trueK_pull_MIP","Well reconstructed tracks (Afro's truth matching), true kaons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);

  All_AfroBT_truemu_PIDA = tfs->make<TH1F>("All_AfroBT_truemu_PIDA","Well reconstructed tracks (Afro's truth matching), true muons;PIDA;",300,0,30);
  All_AfroBT_truep_PIDA = tfs->make<TH1F>("All_AfroBT_truep_PIDA","Well reconstructed tracks (Afro's truth matching), true protons;PIDA;",300,0,30);
  All_AfroBT_truepi_PIDA = tfs->make<TH1F>("All_AfroBT_truepi_PIDA","Well reconstructed tracks (Afro's truth matching), true pions;PIDA;",300,0,30);
  All_AfroBT_trueK_PIDA = tfs->make<TH1F>("All_AfroBT_trueK_PIDA","Well reconstructed tracks (Afro's truth matching), true kaons;PIDA;",300,0,30);

  All_AfroBT_truemu_dEdxtr_len = tfs->make<TH2F>("All_AfroBT_truemu_dEdxtr_len","Well reconstructed tracks (Afro's truth matching), true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_AfroBT_truep_dEdxtr_len = tfs->make<TH2F>("All_AfroBT_truep_dEdxtr_len","Well reconstructed tracks (Afro's truth matching), true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_AfroBT_truepi_dEdxtr_len = tfs->make<TH2F>("All_AfroBT_truepi_dEdxtr_len","Well reconstructed tracks (Afro's truth matching), true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_AfroBT_trueK_dEdxtr_len = tfs->make<TH2F>("All_AfroBT_trueK_dEdxtr_len","Well reconstructed tracks (Afro's truth matching), true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_AfroBT_truemu_dQdxtr_len = tfs->make<TH2F>("All_AfroBT_truemu_dQdxtr_len","Well reconstructed tracks (Afro's truth matching), true muons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_AfroBT_truep_dQdxtr_len = tfs->make<TH2F>("All_AfroBT_truep_dQdxtr_len","Well reconstructed tracks (Afro's truth matching), true protons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_AfroBT_truepi_dQdxtr_len = tfs->make<TH2F>("All_AfroBT_truepi_dQdxtr_len","Well reconstructed tracks (Afro's truth matching), true pions;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_AfroBT_trueK_dQdxtr_len = tfs->make<TH2F>("All_AfroBT_trueK_dQdxtr_len","Well reconstructed tracks (Afro's truth matching), true kaons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  */
}

void ParticleIDValidationPlots::analyze(art::Event const & e)
{
  // Get handles to needed information
  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackingAlgo, trackHandle);
  std::vector<art::Ptr<recob::Track>> trackCollection;
  art::fill_ptr_vector(trackCollection, trackHandle);

  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitAlgo, hitHandle);
  
  art::FindManyP<recob::Hit> hits_from_tracks(trackHandle, e, fHitTrackAssns);
  art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitHandle,e,fTruthMatchingAssns);

  const auto& mcpHandle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");

  // --------- Loop over tracks in event ---------- //
  for (auto& track : trackCollection){
    bool isWR = false;
    int WR_pdg = 0;
    int Assn_pdg = 0;
    //int AfroBT_pdg = 0;

    // Check if track is well reconstructed, and get true PDG that way

    for (auto const& thisMCP : (*mcpHandle)){
      if (thisMCP.StatusCode() != 1) continue;

      if (isWellReconstructed((*track), thisMCP)){
	isWR = true;
	WR_pdg = thisMCP.PdgCode();
	break;
      }
    } // for (auto const& thisMCP : (*mcpHandle))


    // Get true PDG from associations
    
    std::unordered_map<int,double> trkide;
    double maxe=-1, tote=0;
    simb::MCParticle const* maxp_me = NULL; //pointer for the particle match we will calculate
    
    std::vector<simb::MCParticle const*> particle_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

    std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track->ID());
    
    //loop only over our hits
    for(size_t i_h=0; i_h<hits_from_track.size(); ++i_h){
      
      particle_vec.clear(); match_vec.clear();
      particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);
      //the .key() gives us the index in the original collection
      
      //loop over particles
      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
	trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
	tote += match_vec[i_p]->energy; //calculate total energy deposited
	if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
	  maxe = trkide[ particle_vec[i_p]->TrackId() ];
	  maxp_me = particle_vec[i_p];
	}
      }//end loop over particles per hit
      
    }
    
    Assn_pdg = maxp_me->PdgCode();
    
      /*std::cout << std::endl;
      std::cout << "Final Match (Assns: pandoraCosmicHitRemoval) is pdg = " << maxp_me->PdgCode() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
		<< " trkid=" << maxp_me->TrackId()
		<< " ke=" << maxp_me->E()-maxp_me->Mass()
		<< "\n\tstart (x,y,z)=(" << maxp_me->Vx()
		<< "," << maxp_me->Vy()
		<< "," << maxp_me->Vz()
		<< ")\tend (x,y,z)=(" << maxp_me->EndX()
		<< "," << maxp_me->EndY()
		<< "," << maxp_me->EndZ() << ")" << std::endl;*/

    // Get true PDG from backtracker
    
    /*BackTrackerTruthMatch backtrackertruthmatch;
    backtrackertruthmatch.MatchToMCParticle(hitHandle,e,hits_from_track);
    art::Ptr<simb::MCParticle> MCP_BT = backtrackertruthmatch.ReturnMCParticle();

    AfroBT_pdg = MCP_BT->PdgCode();*/


    // ------------------- Now calculate PID variables and fill hists ------------------- //
    art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, fPIDtag);
    if (!trackPIDAssn.isValid()){
      std::cout << "trackPIDAssn.isValid() == false. Skipping track." << std::endl;
      continue;
    }
    
    double Bragg_fwd_mu = -999;
    double Bragg_fwd_p = -999;
    double Bragg_fwd_pi = -999;
    double Bragg_fwd_K = -999;
    double Bragg_bwd_mu = -999;
    double Bragg_bwd_p = -999;
    double Bragg_bwd_pi = -999;
    double Bragg_bwd_K = -999;
    double PIDAval = -999;
    double dEdxtruncmean = -999;
    double dQdxtruncmean = -999;
    double trklen = -999;
    
    std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track->ID());
    if (trackPID.size() == 0){
      std::cout << "No track-PID association found for trackID " << track->ID() << ". Skipping track." << std::endl;
      continue;
    }
    std::cout << "trackPID.size() = " << trackPID.size() << std::endl;
    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

    // Loop through AlgScoresVec and find the variables we want
    for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){

      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
      
      if (AlgScore.fAlgName == "BraggPeakLLH"){
	if (anab::kVariableType(AlgScore.fVariableType) == anab::kLogL_fwd){
	  if (AlgScore.fAssumedPdg == 13)   Bragg_fwd_mu = AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 2212) Bragg_fwd_p =  AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 211)  Bragg_fwd_pi = AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 321)  Bragg_fwd_K  = AlgScore.fValue;
	}// if fVariableType == anab::kLogL_fwd
	else if (anab::kVariableType(AlgScore.fVariableType) == anab::kLogL_bwd){
	  if (AlgScore.fAssumedPdg == 13)   Bragg_bwd_mu = AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 2212) Bragg_bwd_p =  AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 211)  Bragg_bwd_pi = AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 321)  Bragg_bwd_K  = AlgScore.fValue;
	} // if fVariableType == anab::kLogL_bwd
      } // if fAlName = BraggPeakLLH

      if (AlgScore.fAlgName == "PIDA" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
	PIDAval = AlgScore.fValue;
      }// if AlgName = PIDA && fVariableType == anab::kPIDA

      if (AlgScore.fAlgName == "TruncatedMean"){
	if (anab::kVariableType(AlgScore.fVariableType) == anab::kdQdxtruncmean) dQdxtruncmean = AlgScore.fValue;
	if (anab::kVariableType(AlgScore.fVariableType) == anab::kdEdxtruncmean) dEdxtruncmean = AlgScore.fValue;
	if (anab::kVariableType(AlgScore.fVariableType) == anab::kTrackLength)   trklen = AlgScore.fValue;
      }// if AlgName = TruncatedMean
      
    } // Loop over AlgScoresVec


    

    // ---------------- Now let's fill some histograms! ------------------- //

    // Use best fit (lowest neg2LogL) from forward vs backward for calculations
    double Bragg_mu = (Bragg_fwd_mu < Bragg_bwd_mu ? Bragg_fwd_mu : Bragg_bwd_mu);
    double Bragg_p  = (Bragg_fwd_p  < Bragg_bwd_p  ? Bragg_fwd_p  : Bragg_bwd_p);
    double Bragg_pi = (Bragg_fwd_pi < Bragg_bwd_pi ? Bragg_fwd_pi : Bragg_bwd_pi);
    double Bragg_K  = (Bragg_fwd_K  < Bragg_bwd_K  ? Bragg_fwd_K  : Bragg_bwd_K);
    
    // Calculate likelihood ratio on a scale of 0 to 1
    double Bragg_pull_mu = Bragg_mu/(Bragg_mu+Bragg_p+Bragg_pi+Bragg_K);
    double Bragg_pull_p  = Bragg_p /(Bragg_mu+Bragg_p+Bragg_pi+Bragg_K);
    double Bragg_pull_pi = Bragg_pi/(Bragg_mu+Bragg_p+Bragg_pi+Bragg_K);
    double Bragg_pull_K  = Bragg_K /(Bragg_mu+Bragg_p+Bragg_pi+Bragg_K);
    double Bragg_pull_MIPproton = (Bragg_mu+Bragg_pi)/(Bragg_mu+Bragg_pi+Bragg_p);

    if (isWR){
      // Well-reconstructed truth matching
      if (TMath::Abs(WR_pdg) == 13){ // True muons
	WR_truemu_neglogl_mu->Fill(Bragg_mu);
	WR_truemu_neglogl_p->Fill(Bragg_p);
	WR_truemu_neglogl_pi->Fill(Bragg_pi);
	WR_truemu_neglogl_K->Fill(Bragg_K);
	WR_truemu_pull_mu->Fill(Bragg_pull_mu);
	WR_truemu_pull_p->Fill(Bragg_pull_p);
	WR_truemu_pull_pi->Fill(Bragg_pull_pi);
	WR_truemu_pull_K->Fill(Bragg_pull_K);
	WR_truemu_pull_MIP->Fill(Bragg_pull_MIPproton);
	WR_truemu_PIDA->Fill(PIDAval);
	WR_truemu_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	WR_truemu_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(WR_pdg) == 2212){ // True protons
	WR_truep_neglogl_mu->Fill(Bragg_mu);
	WR_truep_neglogl_p->Fill(Bragg_p);
	WR_truep_neglogl_pi->Fill(Bragg_pi);
	WR_truep_neglogl_K->Fill(Bragg_K);
	WR_truep_pull_mu->Fill(Bragg_pull_mu);
	WR_truep_pull_p->Fill(Bragg_pull_p);
	WR_truep_pull_pi->Fill(Bragg_pull_pi);
	WR_truep_pull_K->Fill(Bragg_pull_K);
	WR_truep_pull_MIP->Fill(Bragg_pull_MIPproton);
	WR_truep_PIDA->Fill(PIDAval);
	WR_truep_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	WR_truep_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(WR_pdg) == 211){ // True pions
	WR_truepi_neglogl_mu->Fill(Bragg_mu);
	WR_truepi_neglogl_p->Fill(Bragg_p);
	WR_truepi_neglogl_pi->Fill(Bragg_pi);
	WR_truepi_neglogl_K->Fill(Bragg_K);
	WR_truepi_pull_mu->Fill(Bragg_pull_mu);
	WR_truepi_pull_p->Fill(Bragg_pull_p);
	WR_truepi_pull_pi->Fill(Bragg_pull_pi);
	WR_truepi_pull_K->Fill(Bragg_pull_K);
	WR_truepi_pull_MIP->Fill(Bragg_pull_MIPproton);
	WR_truepi_PIDA->Fill(PIDAval);
	WR_truepi_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	WR_truepi_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(WR_pdg) == 321){ // True kaons
	WR_trueK_neglogl_mu->Fill(Bragg_mu);
	WR_trueK_neglogl_p->Fill(Bragg_p);
	WR_trueK_neglogl_pi->Fill(Bragg_pi);
	WR_trueK_neglogl_K->Fill(Bragg_K);
	WR_trueK_pull_mu->Fill(Bragg_pull_mu);
	WR_trueK_pull_p->Fill(Bragg_pull_p);
	WR_trueK_pull_pi->Fill(Bragg_pull_pi);
	WR_trueK_pull_K->Fill(Bragg_pull_K);
	WR_trueK_pull_MIP->Fill(Bragg_pull_MIPproton);
	WR_trueK_PIDA->Fill(PIDAval);
	WR_trueK_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	WR_trueK_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      // Backtracker truth matching (with Afro's converter)
      /*if (TMath::Abs(AfroBT_pdg) == 13){ // True muons
	All_AfroBT_truemu_neglogl_mu->Fill(Bragg_mu);
	All_AfroBT_truemu_neglogl_p->Fill(Bragg_p);
	All_AfroBT_truemu_neglogl_pi->Fill(Bragg_pi);
	All_AfroBT_truemu_neglogl_K->Fill(Bragg_K);
	All_AfroBT_truemu_pull_mu->Fill(Bragg_pull_mu);
	All_AfroBT_truemu_pull_p->Fill(Bragg_pull_p);
	All_AfroBT_truemu_pull_pi->Fill(Bragg_pull_pi);
	All_AfroBT_truemu_pull_K->Fill(Bragg_pull_K);
	All_AfroBT_truemu_pull_MIP->Fill(Bragg_pull_MIPproton);
	All_AfroBT_truemu_PIDA->Fill(PIDAval);
	All_AfroBT_truemu_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	All_AfroBT_truemu_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(AfroBT_pdg) == 2212){ // True protons
	All_AfroBT_truep_neglogl_mu->Fill(Bragg_mu);
	All_AfroBT_truep_neglogl_p->Fill(Bragg_p);
	All_AfroBT_truep_neglogl_pi->Fill(Bragg_pi);
	All_AfroBT_truep_neglogl_K->Fill(Bragg_K);
	All_AfroBT_truep_pull_mu->Fill(Bragg_pull_mu);
	All_AfroBT_truep_pull_p->Fill(Bragg_pull_p);
	All_AfroBT_truep_pull_pi->Fill(Bragg_pull_pi);
	All_AfroBT_truep_pull_K->Fill(Bragg_pull_K);
	All_AfroBT_truep_pull_MIP->Fill(Bragg_pull_MIPproton);
	All_AfroBT_truep_PIDA->Fill(PIDAval);
	All_AfroBT_truep_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	All_AfroBT_truep_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(AfroBT_pdg) == 211){ // True pions
	All_AfroBT_truepi_neglogl_mu->Fill(Bragg_mu);
	All_AfroBT_truepi_neglogl_p->Fill(Bragg_p);
	All_AfroBT_truepi_neglogl_pi->Fill(Bragg_pi);
	All_AfroBT_truepi_neglogl_K->Fill(Bragg_K);
	All_AfroBT_truepi_pull_mu->Fill(Bragg_pull_mu);
	All_AfroBT_truepi_pull_p->Fill(Bragg_pull_p);
	All_AfroBT_truepi_pull_pi->Fill(Bragg_pull_pi);
	All_AfroBT_truepi_pull_K->Fill(Bragg_pull_K);
	All_AfroBT_truepi_pull_MIP->Fill(Bragg_pull_MIPproton);
	All_AfroBT_truepi_PIDA->Fill(PIDAval);
	All_AfroBT_truepi_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	All_AfroBT_truepi_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(AfroBT_pdg) == 321){ // True kaons
	All_AfroBT_trueK_neglogl_mu->Fill(Bragg_mu);
	All_AfroBT_trueK_neglogl_p->Fill(Bragg_p);
	All_AfroBT_trueK_neglogl_pi->Fill(Bragg_pi);
	All_AfroBT_trueK_neglogl_K->Fill(Bragg_K);
	All_AfroBT_trueK_pull_mu->Fill(Bragg_pull_mu);
	All_AfroBT_trueK_pull_p->Fill(Bragg_pull_p);
	All_AfroBT_trueK_pull_pi->Fill(Bragg_pull_pi);
	All_AfroBT_trueK_pull_K->Fill(Bragg_pull_K);
	All_AfroBT_trueK_pull_MIP->Fill(Bragg_pull_MIPproton);
	All_AfroBT_trueK_PIDA->Fill(PIDAval);
	All_AfroBT_trueK_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	All_AfroBT_trueK_dQdxtr_len->Fill(trklen,dQdxtruncmean);
	}*/
      // Associations truth matching
      if (TMath::Abs(Assn_pdg) == 13){ // True muons
	All_Assn_truemu_neglogl_mu->Fill(Bragg_mu);
	All_Assn_truemu_neglogl_p->Fill(Bragg_p);
	All_Assn_truemu_neglogl_pi->Fill(Bragg_pi);
	All_Assn_truemu_neglogl_K->Fill(Bragg_K);
	All_Assn_truemu_pull_mu->Fill(Bragg_pull_mu);
	All_Assn_truemu_pull_p->Fill(Bragg_pull_p);
	All_Assn_truemu_pull_pi->Fill(Bragg_pull_pi);
	All_Assn_truemu_pull_K->Fill(Bragg_pull_K);
	All_Assn_truemu_pull_MIP->Fill(Bragg_pull_MIPproton);
	All_Assn_truemu_PIDA->Fill(PIDAval);
	All_Assn_truemu_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	All_Assn_truemu_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(Assn_pdg) == 2212){ // True protons
	All_Assn_truep_neglogl_mu->Fill(Bragg_mu);
	All_Assn_truep_neglogl_p->Fill(Bragg_p);
	All_Assn_truep_neglogl_pi->Fill(Bragg_pi);
	All_Assn_truep_neglogl_K->Fill(Bragg_K);
	All_Assn_truep_pull_mu->Fill(Bragg_pull_mu);
	All_Assn_truep_pull_p->Fill(Bragg_pull_p);
	All_Assn_truep_pull_pi->Fill(Bragg_pull_pi);
	All_Assn_truep_pull_K->Fill(Bragg_pull_K);
	All_Assn_truep_pull_MIP->Fill(Bragg_pull_MIPproton);
	All_Assn_truep_PIDA->Fill(PIDAval);
	All_Assn_truep_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	All_Assn_truep_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(Assn_pdg) == 211){ // True pions
	All_Assn_truepi_neglogl_mu->Fill(Bragg_mu);
	All_Assn_truepi_neglogl_p->Fill(Bragg_p);
	All_Assn_truepi_neglogl_pi->Fill(Bragg_pi);
	All_Assn_truepi_neglogl_K->Fill(Bragg_K);
	All_Assn_truepi_pull_mu->Fill(Bragg_pull_mu);
	All_Assn_truepi_pull_p->Fill(Bragg_pull_p);
	All_Assn_truepi_pull_pi->Fill(Bragg_pull_pi);
	All_Assn_truepi_pull_K->Fill(Bragg_pull_K);
	All_Assn_truepi_pull_MIP->Fill(Bragg_pull_MIPproton);
	All_Assn_truepi_PIDA->Fill(PIDAval);
	All_Assn_truepi_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	All_Assn_truepi_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
      
      else if (TMath::Abs(Assn_pdg) == 321){ // True kaons
	All_Assn_trueK_neglogl_mu->Fill(Bragg_mu);
	All_Assn_trueK_neglogl_p->Fill(Bragg_p);
	All_Assn_trueK_neglogl_pi->Fill(Bragg_pi);
	All_Assn_trueK_neglogl_K->Fill(Bragg_K);
	All_Assn_trueK_pull_mu->Fill(Bragg_pull_mu);
	All_Assn_trueK_pull_p->Fill(Bragg_pull_p);
	All_Assn_trueK_pull_pi->Fill(Bragg_pull_pi);
	All_Assn_trueK_pull_K->Fill(Bragg_pull_K);
	All_Assn_trueK_pull_MIP->Fill(Bragg_pull_MIPproton);
	All_Assn_trueK_PIDA->Fill(PIDAval);
	All_Assn_trueK_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	All_Assn_trueK_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
    } // if isWR
    
  } // Loop over tracks

  
}

DEFINE_ART_MODULE(ParticleIDValidationPlots)
