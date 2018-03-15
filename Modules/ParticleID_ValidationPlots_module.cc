////////////////////////////////////////////////////////////////////////
// Class:       ParticleID_ValidationPlots
// Plugin Type: analyzer (art v2_05_01)
// File:        ParticleID_ValidationPlots_module.cc
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
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "uboone/ParticleID/Algorithms/WellReconstructedTrackFinder.h"

#include "lardataobj/RecoBase/Track.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "TH1F.h"
#include "TH2F.h"

class ParticleID_ValidationPlots;


class ParticleID_ValidationPlots : public art::EDAnalyzer {
public:
  explicit ParticleID_ValidationPlots(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ParticleID_ValidationPlots(ParticleID_ValidationPlots const &) = delete;
  ParticleID_ValidationPlots(ParticleID_ValidationPlots &&) = delete;
  ParticleID_ValidationPlots & operator = (ParticleID_ValidationPlots const &) = delete;
  ParticleID_ValidationPlots & operator = (ParticleID_ValidationPlots &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  art::ServiceHandle<art::TFileService> tfs;

  std::string fTrackingAlgo;
  std::string fPIDtag;

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

  TH1F *All_BT_truemu_pull_mu;
  TH1F *All_BT_truep_pull_mu;
  TH1F *All_BT_truepi_pull_mu;
  TH1F *All_BT_trueK_pull_mu;
  TH1F *All_BT_truemu_pull_p;
  TH1F *All_BT_truep_pull_p;
  TH1F *All_BT_truepi_pull_p;
  TH1F *All_BT_trueK_pull_p;
  TH1F *All_BT_truemu_pull_pi;
  TH1F *All_BT_truep_pull_pi;
  TH1F *All_BT_truepi_pull_pi;
  TH1F *All_BT_trueK_pull_pi;
  TH1F *All_BT_truemu_pull_K;
  TH1F *All_BT_truep_pull_K;
  TH1F *All_BT_truepi_pull_K;
  TH1F *All_BT_trueK_pull_K;

  TH1F *All_BT_truemu_pull_MIP;
  TH1F *All_BT_truep_pull_MIP;
  TH1F *All_BT_truepi_pull_MIP;
  TH1F *All_BT_trueK_pull_MIP;

  TH1F *All_BT_truemu_PIDA;
  TH1F *All_BT_truep_PIDA;
  TH1F *All_BT_truepi_PIDA;
  TH1F *All_BT_trueK_PIDA;

  TH2F *All_BT_truemu_dEdxtr_len;
  TH2F *All_BT_truep_dEdxtr_len;
  TH2F *All_BT_truepi_dEdxtr_len;
  TH2F *All_BT_trueK_dEdxtr_len;
  TH2F *All_BT_truemu_dQdxtr_len;
  TH2F *All_BT_truep_dQdxtr_len;
  TH2F *All_BT_truepi_dQdxtr_len;
  TH2F *All_BT_trueK_dQdxtr_len;
};


ParticleID_ValidationPlots::ParticleID_ValidationPlots(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fTrackingAlgo = p.get<std::string>("TrackingAlgorithm");
  fPIDtag = p.get<std::string>("ParticleIDProducerModule");
  
  // ---- Well reconstructed
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
  All_Assn_truemu_pull_mu = tfs->make<TH1F>("All_Assn_truemu_pull_mu","Well reconstructed tracks, true muons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truep_pull_mu = tfs->make<TH1F>("All_Assn_truep_pull_mu","Well reconstructed tracks, true protons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truepi_pull_mu = tfs->make<TH1F>("All_Assn_truepi_pull_mu","Well reconstructed tracks, true pions;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_trueK_pull_mu = tfs->make<TH1F>("All_Assn_trueK_pull_mu","Well reconstructed tracks, true kaons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truemu_pull_p = tfs->make<TH1F>("All_Assn_truemu_pull_p","Well reconstructed tracks, true muons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truep_pull_p = tfs->make<TH1F>("All_Assn_truep_pull_p","Well reconstructed tracks, true protons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truepi_pull_p = tfs->make<TH1F>("All_Assn_truepi_pull_p","Well reconstructed tracks, true pions;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_trueK_pull_p = tfs->make<TH1F>("All_Assn_trueK_pull_p","Well reconstructed tracks, true kaons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truemu_pull_pi = tfs->make<TH1F>("All_Assn_truemu_pull_pi","Well reconstructed tracks, true muons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truep_pull_pi = tfs->make<TH1F>("All_Assn_truep_pull_pi","Well reconstructed tracks, true protons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truepi_pull_pi = tfs->make<TH1F>("All_Assn_truepi_pull_pi","Well reconstructed tracks, true pions;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_trueK_pull_pi = tfs->make<TH1F>("All_Assn_trueK_pull_pi","Well reconstructed tracks, true kaons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truemu_pull_K = tfs->make<TH1F>("All_Assn_truemu_pull_K","Well reconstructed tracks, true muons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truep_pull_K = tfs->make<TH1F>("All_Assn_truep_pull_K","Well reconstructed tracks, true protons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_truepi_pull_K = tfs->make<TH1F>("All_Assn_truepi_pull_K","Well reconstructed tracks, true pions;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_Assn_trueK_pull_K = tfs->make<TH1F>("All_Assn_trueK_pull_K","Well reconstructed tracks, true kaons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);

  All_Assn_truemu_pull_MIP = tfs->make<TH1F>("All_Assn_truemu_pull_MIP","Well reconstructed tracks, true muons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_Assn_truep_pull_MIP = tfs->make<TH1F>("All_Assn_truep_pull_MIP","Well reconstructed tracks, true protons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_Assn_truepi_pull_MIP = tfs->make<TH1F>("All_Assn_truepi_pull_MIP","Well reconstructed tracks, true pions;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_Assn_trueK_pull_MIP = tfs->make<TH1F>("All_Assn_trueK_pull_MIP","Well reconstructed tracks, true kaons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);

  All_Assn_truemu_PIDA = tfs->make<TH1F>("All_Assn_truemu_PIDA","Well reconstructed tracks, true muons;PIDA;",300,0,30);
  All_Assn_truep_PIDA = tfs->make<TH1F>("All_Assn_truep_PIDA","Well reconstructed tracks, true protons;PIDA;",300,0,30);
  All_Assn_truepi_PIDA = tfs->make<TH1F>("All_Assn_truepi_PIDA","Well reconstructed tracks, true pions;PIDA;",300,0,30);
  All_Assn_trueK_PIDA = tfs->make<TH1F>("All_Assn_trueK_PIDA","Well reconstructed tracks, true kaons;PIDA;",300,0,30);

  All_Assn_truemu_dEdxtr_len = tfs->make<TH2F>("All_Assn_truemu_dEdxtr_len","Well reconstructed tracks, true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_Assn_truep_dEdxtr_len = tfs->make<TH2F>("All_Assn_truep_dEdxtr_len","Well reconstructed tracks, true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_Assn_truepi_dEdxtr_len = tfs->make<TH2F>("All_Assn_truepi_dEdxtr_len","Well reconstructed tracks, true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_Assn_trueK_dEdxtr_len = tfs->make<TH2F>("All_Assn_trueK_dEdxtr_len","Well reconstructed tracks, true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_Assn_truemu_dQdxtr_len = tfs->make<TH2F>("All_Assn_truemu_dQdxtr_len","Well reconstructed tracks, true muons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_Assn_truep_dQdxtr_len = tfs->make<TH2F>("All_Assn_truep_dQdxtr_len","Well reconstructed tracks, true protons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_Assn_truepi_dQdxtr_len = tfs->make<TH2F>("All_Assn_truepi_dQdxtr_len","Well reconstructed tracks, true pions;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_Assn_trueK_dQdxtr_len = tfs->make<TH2F>("All_Assn_trueK_dQdxtr_len","Well reconstructed tracks, true kaons;Track length (cm);dQ/dx",100,0,700,100,0,500);



  // ---- All tracks: truth matching using associations
  All_BT_truemu_pull_mu = tfs->make<TH1F>("All_BT_truemu_pull_mu","Well reconstructed tracks, true muons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truep_pull_mu = tfs->make<TH1F>("All_BT_truep_pull_mu","Well reconstructed tracks, true protons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truepi_pull_mu = tfs->make<TH1F>("All_BT_truepi_pull_mu","Well reconstructed tracks, true pions;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_trueK_pull_mu = tfs->make<TH1F>("All_BT_trueK_pull_mu","Well reconstructed tracks, true kaons;(neg2LL_mu)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truemu_pull_p = tfs->make<TH1F>("All_BT_truemu_pull_p","Well reconstructed tracks, true muons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truep_pull_p = tfs->make<TH1F>("All_BT_truep_pull_p","Well reconstructed tracks, true protons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truepi_pull_p = tfs->make<TH1F>("All_BT_truepi_pull_p","Well reconstructed tracks, true pions;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_trueK_pull_p = tfs->make<TH1F>("All_BT_trueK_pull_p","Well reconstructed tracks, true kaons;(neg2LL_p)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truemu_pull_pi = tfs->make<TH1F>("All_BT_truemu_pull_pi","Well reconstructed tracks, true muons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truep_pull_pi = tfs->make<TH1F>("All_BT_truep_pull_pi","Well reconstructed tracks, true protons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truepi_pull_pi = tfs->make<TH1F>("All_BT_truepi_pull_pi","Well reconstructed tracks, true pions;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_trueK_pull_pi = tfs->make<TH1F>("All_BT_trueK_pull_pi","Well reconstructed tracks, true kaons;(neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truemu_pull_K = tfs->make<TH1F>("All_BT_truemu_pull_K","Well reconstructed tracks, true muons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truep_pull_K = tfs->make<TH1F>("All_BT_truep_pull_K","Well reconstructed tracks, true protons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_truepi_pull_K = tfs->make<TH1F>("All_BT_truepi_pull_K","Well reconstructed tracks, true pions;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);
  All_BT_trueK_pull_K = tfs->make<TH1F>("All_BT_trueK_pull_K","Well reconstructed tracks, true kaons;(neg2LL_K)/(neg2LL_mu+neg2LL_p+neg2LL_pi+neg2LL_K);",50,0,1);

  All_BT_truemu_pull_MIP = tfs->make<TH1F>("All_BT_truemu_pull_MIP","Well reconstructed tracks, true muons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_BT_truep_pull_MIP = tfs->make<TH1F>("All_BT_truep_pull_MIP","Well reconstructed tracks, true protons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_BT_truepi_pull_MIP = tfs->make<TH1F>("All_BT_truepi_pull_MIP","Well reconstructed tracks, true pions;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);
  All_BT_trueK_pull_MIP = tfs->make<TH1F>("All_BT_trueK_pull_MIP","Well reconstructed tracks, true kaons;(neg2LL_mu+neg2LL_pi)/(neg2LL_mu+neg2LL_p+neg2LL_pi);",50,0,1);

  All_BT_truemu_PIDA = tfs->make<TH1F>("All_BT_truemu_PIDA","Well reconstructed tracks, true muons;PIDA;",300,0,30);
  All_BT_truep_PIDA = tfs->make<TH1F>("All_BT_truep_PIDA","Well reconstructed tracks, true protons;PIDA;",300,0,30);
  All_BT_truepi_PIDA = tfs->make<TH1F>("All_BT_truepi_PIDA","Well reconstructed tracks, true pions;PIDA;",300,0,30);
  All_BT_trueK_PIDA = tfs->make<TH1F>("All_BT_trueK_PIDA","Well reconstructed tracks, true kaons;PIDA;",300,0,30);

  All_BT_truemu_dEdxtr_len = tfs->make<TH2F>("All_BT_truemu_dEdxtr_len","Well reconstructed tracks, true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_BT_truep_dEdxtr_len = tfs->make<TH2F>("All_BT_truep_dEdxtr_len","Well reconstructed tracks, true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_BT_truepi_dEdxtr_len = tfs->make<TH2F>("All_BT_truepi_dEdxtr_len","Well reconstructed tracks, true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_BT_trueK_dEdxtr_len = tfs->make<TH2F>("All_BT_trueK_dEdxtr_len","Well reconstructed tracks, true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_BT_truemu_dQdxtr_len = tfs->make<TH2F>("All_BT_truemu_dQdxtr_len","Well reconstructed tracks, true muons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_BT_truep_dQdxtr_len = tfs->make<TH2F>("All_BT_truep_dQdxtr_len","Well reconstructed tracks, true protons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_BT_truepi_dQdxtr_len = tfs->make<TH2F>("All_BT_truepi_dQdxtr_len","Well reconstructed tracks, true pions;Track length (cm);dQ/dx",100,0,700,100,0,500);
  All_BT_trueK_dQdxtr_len = tfs->make<TH2F>("All_BT_trueK_dQdxtr_len","Well reconstructed tracks, true kaons;Track length (cm);dQ/dx",100,0,700,100,0,500);
  
}

void ParticleID_ValidationPlots::analyze(art::Event const & e)
{
  // Get handles to needed information
  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(ftrackingAlgo, trackHandle);
  std::vector<art::Ptr<recob::Track>> trackCollection;
  art::fill_ptr_vector(trackCollection, trackHandle);

  // --------- Loop over tracks in event ---------- //
  for (auto track : trackCollection){

    // Check if track is well reconstructed, and get true PDG that way
    bool isWR = false;
    int WR_pdg = 0;
    
    const auto& mcpHandle = e.getValidHandle<std::vector<simb::MCParticle>>("largeant");

    for (auto const& thisMCP : (*mcpHandle)){
      if (thisMCP.StatusCode() != 1) continue;

      if (isWellReconstructed(track, thisMCP)){
	isWR = true;
	WR_pdg = thisMCP.PdgCode();
	break;
      }
    } // for (auto const& thisMCP : (*mcpHandle))


    // Get true PDG from associations

    // Get true PDG from backtracker


    // ------------------- Now calculate PID variables and fill hists ------------------- //
    art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, fPIDtag);

    double Bragg_fwd_mu;
    double Bragg_fwd_p;
    double Bragg_fwd_pi;
    double Bragg_fwd_K;
    double Bragg_bwd_mu;
    double Bragg_bwd_p;
    double Bragg_bwd_pi;
    double Bragg_bwd_K;
    double PIDAval;
    double dEdxtruncmean;
    double dQdxtruncmean;
    double trklen;
    
    anab::ParticleID trackPID = trackPIDAssn.at(track.ID());
    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.ParticleIDAlgScores;

    // Loop through AlgScoresVec and find the variables we want
    for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){

      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
      
      if (AlgScore.fAlgName == "BraggPeakLLH"){
	if (AlgScore.fVariableType == anab::kLogL_fwd){
	  if (AlgScore.fAssumedPdg == 13)   Bragg_fwd_mu = AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 2212) Bragg_fwd_p =  AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 211)  Bragg_fwd_pi = AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 321)  Bragg_fwd_K  = AlgScore.fValue;
	}// if fVariableType == anab::kLogL_fwd
	else if (AlgScore.fVariableType == anab::kLogL_bwd){
	  if (AlgScore.fAssumedPdg == 13)   Bragg_bwd_mu = AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 2212) Bragg_bwd_p =  AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 211)  Bragg_bwd_pi = AlgScore.fValue;
	  if (AlgScore.fAssumedPdg == 321)  Bragg_bwd_K  = AlgScore.fValue;
	} // if fVariableType == anab::kLogL_bwd
      } // if fAlName = BraggPeakLLH

      if (AlgScore.fAlgName == "PIDA" && AlgScore.fVariableType == anab::kPIDA){
	PIDAval = AlgScore.fValue;
      }// if AlgName = PIDA && fVariableType == anab::kPIDA

      if (AlgScore.fAlgName == "TruncatedMean"){
	if (AlgScore.fVariableType == anab::kdQdxtruncmean) dQdxtruncmean = AlgScore.fValue;
	if (AlgScore.fVariableType == anab::kdEdxtruncmean) dEdxtruncmean = AlgScore.fValue;
	if (AlgScore.fVariableType == anab::kTrackLength)   trklen = AlgScore.fValue;
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
      if (TMath::Abs(WR_pdg) == 13){ // True muons
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
	WR_trueK_pull_mu->Fill(Bragg_pull_mu);
	WR_trueK_pull_p->Fill(Bragg_pull_p);
	WR_trueK_pull_pi->Fill(Bragg_pull_pi);
	WR_trueK_pull_K->Fill(Bragg_pull_K);
	WR_trueK_pull_MIP->Fill(Bragg_pull_MIPproton);
	WR_trueK_PIDA->Fill(PIDAval);
	WR_trueK_dEdxtr_len->Fill(trklen,dEdxtruncmean);
	WR_trueK_dQdxtr_len->Fill(trklen,dQdxtruncmean);
      }
    } // if isWR
    
  } // Loop over tracks

  
}

DEFINE_ART_MODULE(ParticleID_ValidationPlots)
