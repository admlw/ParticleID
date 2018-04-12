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

#include "uboone/ParticleID/Algorithms/fiducialVolume.h"

//#include "uboone/MyClasses/BackTrackerTruthMatch.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

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
    fidvol::fiducialVolume fid;

    std::vector<double> fv;

    std::string fTrackingAlgo;
    std::string fHitAlgo;
    std::string fHitTrackAssns;
    std::string fTruthMatchingAssns;
    std::string fPIDtag;

    TH1F *TrueBragg_truemu_neglogl_mu;
    TH1F *TrueBragg_truep_neglogl_mu;
    TH1F *TrueBragg_truepi_neglogl_mu;
    TH1F *TrueBragg_trueK_neglogl_mu;
    TH1F *TrueBragg_truemu_neglogl_p;
    TH1F *TrueBragg_truep_neglogl_p;
    TH1F *TrueBragg_truepi_neglogl_p;
    TH1F *TrueBragg_trueK_neglogl_p;
    TH1F *TrueBragg_truemu_neglogl_pi;
    TH1F *TrueBragg_truep_neglogl_pi;
    TH1F *TrueBragg_truepi_neglogl_pi;
    TH1F *TrueBragg_trueK_neglogl_pi;
    TH1F *TrueBragg_truemu_neglogl_K;
    TH1F *TrueBragg_truep_neglogl_K;
    TH1F *TrueBragg_truepi_neglogl_K;
    TH1F *TrueBragg_trueK_neglogl_K;
    TH1F *TrueBragg_truemu_neglogl_MIP;
    TH1F *TrueBragg_truep_neglogl_MIP;
    TH1F *TrueBragg_truepi_neglogl_MIP;
    TH1F *TrueBragg_trueK_neglogl_MIP;
    TH1F *TrueBragg_truemu_neglogl_minmuMIP;
    TH1F *TrueBragg_truep_neglogl_minmuMIP;
    TH1F *TrueBragg_truepi_neglogl_minmuMIP;
    TH1F *TrueBragg_trueK_neglogl_minmuMIP;

    TH2F *TrueBragg_truemu_neglogl_mu_vslength;
    TH2F *TrueBragg_truep_neglogl_mu_vslength;
    TH2F *TrueBragg_truepi_neglogl_mu_vslength;
    TH2F *TrueBragg_trueK_neglogl_mu_vslength;
    TH2F *TrueBragg_truemu_neglogl_MIP_vslength;
    TH2F *TrueBragg_truep_neglogl_MIP_vslength;
    TH2F *TrueBragg_truepi_neglogl_MIP_vslength;
    TH2F *TrueBragg_trueK_neglogl_MIP_vslength;
    TH2F *TrueBragg_truemu_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_truep_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_truepi_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_trueK_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_truemu_neglogl_p_vslength;
    TH2F *TrueBragg_truep_neglogl_p_vslength;
    TH2F *TrueBragg_truepi_neglogl_p_vslength;
    TH2F *TrueBragg_trueK_neglogl_p_vslength;
    TH2F *TrueBragg_truemu_neglogl_mu_vsangle;
    TH2F *TrueBragg_truep_neglogl_mu_vsangle;
    TH2F *TrueBragg_truepi_neglogl_mu_vsangle;
    TH2F *TrueBragg_trueK_neglogl_mu_vsangle;
    TH2F *TrueBragg_truemu_neglogl_MIP_vsangle;
    TH2F *TrueBragg_truep_neglogl_MIP_vsangle;
    TH2F *TrueBragg_truepi_neglogl_MIP_vsangle;
    TH2F *TrueBragg_trueK_neglogl_MIP_vsangle;
    TH2F *TrueBragg_truemu_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_truep_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_truepi_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_trueK_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_truemu_neglogl_p_vsangle;
    TH2F *TrueBragg_truep_neglogl_p_vsangle;
    TH2F *TrueBragg_truepi_neglogl_p_vsangle;
    TH2F *TrueBragg_trueK_neglogl_p_vsangle;
    TH2F *TrueBragg_truemu_neglogl_mu_vsnhits;
    TH2F *TrueBragg_truep_neglogl_mu_vsnhits;
    TH2F *TrueBragg_truepi_neglogl_mu_vsnhits;
    TH2F *TrueBragg_trueK_neglogl_mu_vsnhits;
    TH2F *TrueBragg_truemu_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_truep_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_truepi_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_trueK_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_truemu_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_truep_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_truepi_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_trueK_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_truemu_neglogl_p_vsnhits;
    TH2F *TrueBragg_truep_neglogl_p_vsnhits;
    TH2F *TrueBragg_truepi_neglogl_p_vsnhits;
    TH2F *TrueBragg_trueK_neglogl_p_vsnhits;

    TH2F *TrueBragg_truemu_neglogl_muvsp;
    TH2F *TrueBragg_truep_neglogl_muvsp;
    TH2F *TrueBragg_truepi_neglogl_muvsp;
    TH2F *TrueBragg_trueK_neglogl_muvsp;
    TH2F *TrueBragg_truemu_neglogl_muvspi;
    TH2F *TrueBragg_truep_neglogl_muvspi;
    TH2F *TrueBragg_truepi_neglogl_muvspi;
    TH2F *TrueBragg_trueK_neglogl_muvspi;

    TH2F *TrueBragg_truemu_neglogl_MIPvsp;
    TH2F *TrueBragg_truep_neglogl_MIPvsp;
    TH2F *TrueBragg_truemu_neglogl_minmuMIPvsp;
    TH2F *TrueBragg_truep_neglogl_minmuMIPvsp;

    TH1F *TrueBragg_truemu_neglogl_muoverp;
    TH1F *TrueBragg_truep_neglogl_muoverp;
    TH1F *TrueBragg_truepi_neglogl_muoverp;
    TH1F *TrueBragg_trueK_neglogl_muoverp;
    TH1F *TrueBragg_truemu_neglogl_muminusp;
    TH1F *TrueBragg_truep_neglogl_muminusp;
    TH1F *TrueBragg_truepi_neglogl_muminusp;
    TH1F *TrueBragg_trueK_neglogl_muminusp;

    TH1F *TrueBragg_truemu_neglogl_MIPminusp;
    TH1F *TrueBragg_truep_neglogl_MIPminusp;
    TH1F *TrueBragg_truemu_neglogl_minmuMIPminusp;
    TH1F *TrueBragg_truep_neglogl_minmuMIPminusp;

    TH1F *TrueBragg_truemu_smallest_neglogl;
    TH1F *TrueBragg_truep_smallest_neglogl;
    TH1F *TrueBragg_truepi_smallest_neglogl;
    TH1F *TrueBragg_trueK_smallest_neglogl;

    TH1F *TrueBragg_truemu_PIDA;
    TH1F *TrueBragg_truep_PIDA;
    TH1F *TrueBragg_truepi_PIDA;
    TH1F *TrueBragg_trueK_PIDA;

    TH1F *TrueBragg_chargeEndOverStart_directionCorrect;
    TH1F *TrueBragg_chargeEndOverStart_directionIncorrect;
    TH2F *TrueBragg_chargeEndOverStartVersusNHits_directionCorrect;
    TH2F *TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect;

    TH2F *TrueBragg_chargeEndOverStart_sm0_5_dEdxrr;
    TH2F *TrueBragg_chargeEndOverStart_gr2_dEdxrr;
    TH2F *TrueBragg_chargeEndOverStart_0_5to2_dEdxrr;

    TH2F *TrueBragg_truemu_dEdxtr_len;
    TH2F *TrueBragg_truep_dEdxtr_len;
    TH2F *TrueBragg_truepi_dEdxtr_len;
    TH2F *TrueBragg_trueK_dEdxtr_len;

    TH2F *TrueBragg_correctdirection;

    TH1F *All_truemu_neglogl_mu;
    TH1F *All_truep_neglogl_mu;
    TH1F *All_truepi_neglogl_mu;
    TH1F *All_trueK_neglogl_mu;
    TH1F *All_truemu_neglogl_p;
    TH1F *All_truep_neglogl_p;
    TH1F *All_truepi_neglogl_p;
    TH1F *All_trueK_neglogl_p;
    TH1F *All_truemu_neglogl_pi;
    TH1F *All_truep_neglogl_pi;
    TH1F *All_truepi_neglogl_pi;
    TH1F *All_trueK_neglogl_pi;
    TH1F *All_truemu_neglogl_K;
    TH1F *All_truep_neglogl_K;
    TH1F *All_truepi_neglogl_K;
    TH1F *All_trueK_neglogl_K;
    TH1F *All_truemu_neglogl_MIP;
    TH1F *All_truep_neglogl_MIP;
    TH1F *All_truepi_neglogl_MIP;
    TH1F *All_trueK_neglogl_MIP;
    TH1F *All_truemu_neglogl_minmuMIP;
    TH1F *All_truep_neglogl_minmuMIP;
    TH1F *All_truepi_neglogl_minmuMIP;
    TH1F *All_trueK_neglogl_minmuMIP;

    TH2F *All_truemu_neglogl_mu_vslength;
    TH2F *All_truep_neglogl_mu_vslength;
    TH2F *All_truepi_neglogl_mu_vslength;
    TH2F *All_trueK_neglogl_mu_vslength;
    TH2F *All_truemu_neglogl_MIP_vslength;
    TH2F *All_truep_neglogl_MIP_vslength;
    TH2F *All_truepi_neglogl_MIP_vslength;
    TH2F *All_trueK_neglogl_MIP_vslength;
    TH2F *All_truemu_neglogl_minmuMIP_vslength;
    TH2F *All_truep_neglogl_minmuMIP_vslength;
    TH2F *All_truepi_neglogl_minmuMIP_vslength;
    TH2F *All_trueK_neglogl_minmuMIP_vslength;
    TH2F *All_truemu_neglogl_p_vslength;
    TH2F *All_truep_neglogl_p_vslength;
    TH2F *All_truepi_neglogl_p_vslength;
    TH2F *All_trueK_neglogl_p_vslength;
    TH2F *All_truemu_neglogl_mu_vsangle;
    TH2F *All_truep_neglogl_mu_vsangle;
    TH2F *All_truepi_neglogl_mu_vsangle;
    TH2F *All_trueK_neglogl_mu_vsangle;
    TH2F *All_truemu_neglogl_MIP_vsangle;
    TH2F *All_truep_neglogl_MIP_vsangle;
    TH2F *All_truepi_neglogl_MIP_vsangle;
    TH2F *All_trueK_neglogl_MIP_vsangle;
    TH2F *All_truemu_neglogl_minmuMIP_vsangle;
    TH2F *All_truep_neglogl_minmuMIP_vsangle;
    TH2F *All_truepi_neglogl_minmuMIP_vsangle;
    TH2F *All_trueK_neglogl_minmuMIP_vsangle;
    TH2F *All_truemu_neglogl_p_vsangle;
    TH2F *All_truep_neglogl_p_vsangle;
    TH2F *All_truepi_neglogl_p_vsangle;
    TH2F *All_trueK_neglogl_p_vsangle;
    TH2F *All_truemu_neglogl_mu_vsnhits;
    TH2F *All_truep_neglogl_mu_vsnhits;
    TH2F *All_truepi_neglogl_mu_vsnhits;
    TH2F *All_trueK_neglogl_mu_vsnhits;
    TH2F *All_truemu_neglogl_MIP_vsnhits;
    TH2F *All_truep_neglogl_MIP_vsnhits;
    TH2F *All_truepi_neglogl_MIP_vsnhits;
    TH2F *All_trueK_neglogl_MIP_vsnhits;
    TH2F *All_truemu_neglogl_minmuMIP_vsnhits;
    TH2F *All_truep_neglogl_minmuMIP_vsnhits;
    TH2F *All_truepi_neglogl_minmuMIP_vsnhits;
    TH2F *All_trueK_neglogl_minmuMIP_vsnhits;
    TH2F *All_truemu_neglogl_p_vsnhits;
    TH2F *All_truep_neglogl_p_vsnhits;
    TH2F *All_truepi_neglogl_p_vsnhits;
    TH2F *All_trueK_neglogl_p_vsnhits;

    TH2F *All_truemu_neglogl_muvsp;
    TH2F *All_truep_neglogl_muvsp;
    TH2F *All_truepi_neglogl_muvsp;
    TH2F *All_trueK_neglogl_muvsp;
    TH2F *All_truemu_neglogl_muvspi;
    TH2F *All_truep_neglogl_muvspi;
    TH2F *All_truepi_neglogl_muvspi;
    TH2F *All_trueK_neglogl_muvspi;

    TH2F *All_truemu_neglogl_MIPvsp;
    TH2F *All_truep_neglogl_MIPvsp;
    TH2F *All_truemu_neglogl_minmuMIPvsp;
    TH2F *All_truep_neglogl_minmuMIPvsp;

    TH1F *All_truemu_neglogl_muoverp;
    TH1F *All_truep_neglogl_muoverp;
    TH1F *All_truepi_neglogl_muoverp;
    TH1F *All_trueK_neglogl_muoverp;
    TH1F *All_truemu_neglogl_muminusp;
    TH1F *All_truep_neglogl_muminusp;
    TH1F *All_truepi_neglogl_muminusp;
    TH1F *All_trueK_neglogl_muminusp;

    TH1F *All_truemu_neglogl_MIPminusp;
    TH1F *All_truep_neglogl_MIPminusp;
    TH1F *All_truemu_neglogl_minmuMIPminusp;
    TH1F *All_truep_neglogl_minmuMIPminusp;

    TH1F *All_truemu_smallest_neglogl;
    TH1F *All_truep_smallest_neglogl;
    TH1F *All_truepi_smallest_neglogl;
    TH1F *All_trueK_smallest_neglogl;

    TH1F *All_truemu_PIDA;
    TH1F *All_truep_PIDA;
    TH1F *All_truepi_PIDA;
    TH1F *All_trueK_PIDA;

    TH1F *All_chargeEndOverStart_directionCorrect;
    TH1F *All_chargeEndOverStart_directionIncorrect;
    TH2F *All_chargeEndOverStartVersusNHits_directionCorrect;
    TH2F *All_chargeEndOverStartVersusNHits_directionIncorrect;

    TH1F *Contained_chargeEndOverStart_directionCorrect;
    TH1F *Contained_chargeEndOverStart_directionIncorrect;
    TH2F *Contained_chargeEndOverStartVersusNHits_directionCorrect;
    TH2F *Contained_chargeEndOverStartVersusNHits_directionIncorrect;

    TH2F *All_chargeEndOverStart_sm0_5_dEdxrr;
    TH2F *All_chargeEndOverStart_gr2_dEdxrr;
    TH2F *All_chargeEndOverStart_0_5to2_dEdxrr;

    TH2F *All_truemu_dEdxtr_len;
    TH2F *All_truep_dEdxtr_len;
    TH2F *All_truepi_dEdxtr_len;
    TH2F *All_trueK_dEdxtr_len;

    TH2F *All_correctdirection;
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

  fv = fid.setFiducialVolume(fv, p);
  fid.printFiducialVolume(fv);

  // ---- True Bragg peak
  TrueBragg_truemu_neglogl_mu = tfs->make<TH1F>("TrueBragg_truemu_neglogl_mu","Tracks with true p=0 at end, true muons;neg2LL_mu;",200,0,200);
  TrueBragg_truep_neglogl_mu  = tfs->make<TH1F>("TrueBragg_truep_neglogl_mu","Tracks with true p=0 at end, true protons;neg2LL_mu;",200,0,200);
  TrueBragg_truepi_neglogl_mu = tfs->make<TH1F>("TrueBragg_truepi_neglogl_mu","Tracks with true p=0 at end, true pions;neg2LL_mu;",200,0,200);
  TrueBragg_trueK_neglogl_mu  = tfs->make<TH1F>("TrueBragg_trueK_neglogl_mu","Tracks with true p=0 at end, true kaons;neg2LL_mu;",200,0,200);
  TrueBragg_truemu_neglogl_p  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_p","Tracks with true p=0 at end, true muons;neg2LL_p;",200,0,200);
  TrueBragg_truep_neglogl_p   = tfs->make<TH1F>("TrueBragg_truep_neglogl_p","Tracks with true p=0 at end, true protons;neg2LL_p;",200,0,200);
  TrueBragg_truepi_neglogl_p  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_p","Tracks with true p=0 at end, true pions;neg2LL_p;",200,0,200);
  TrueBragg_trueK_neglogl_p   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_p","Tracks with true p=0 at end, true kaons;neg2LL_p;",200,0,200);
  TrueBragg_truemu_neglogl_pi = tfs->make<TH1F>("TrueBragg_truemu_neglogl_pi","Tracks with true p=0 at end, true muons;neg2LL_pi;",200,0,200);
  TrueBragg_truep_neglogl_pi  = tfs->make<TH1F>("TrueBragg_truep_neglogl_pi","Tracks with true p=0 at end, true protons;neg2LL_pi;",200,0,200);
  TrueBragg_truepi_neglogl_pi = tfs->make<TH1F>("TrueBragg_truepi_neglogl_pi","Tracks with true p=0 at end, true pions;neg2LL_pi;",200,0,200);
  TrueBragg_trueK_neglogl_pi  = tfs->make<TH1F>("TrueBragg_trueK_neglogl_pi","Tracks with true p=0 at end, true kaons;neg2LL_pi;",200,0,200);
  TrueBragg_truemu_neglogl_K  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_K","Tracks with true p=0 at end, true muons;neg2LL_K;",200,0,200);
  TrueBragg_truep_neglogl_K   = tfs->make<TH1F>("TrueBragg_truep_neglogl_K","Tracks with true p=0 at end, true protons;neg2LL_K;",200,0,200);
  TrueBragg_truepi_neglogl_K  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_K","Tracks with true p=0 at end, true pions;neg2LL_K;",200,0,200);
  TrueBragg_trueK_neglogl_K   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_K","Tracks with true p=0 at end, true kaons;neg2LL_K;",200,0,200);
  TrueBragg_truemu_neglogl_MIP  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_MIP","Tracks with true p=0 at end, true muons;neg2LL_MIP;",200,0,200);
  TrueBragg_truep_neglogl_MIP   = tfs->make<TH1F>("TrueBragg_truep_neglogl_MIP","Tracks with true p=0 at end, true protons;neg2LL_MIP;",200,0,200);
  TrueBragg_truepi_neglogl_MIP  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_MIP","Tracks with true p=0 at end, true pions;neg2LL_MIP;",200,0,200);
  TrueBragg_trueK_neglogl_MIP   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_MIP","Tracks with true p=0 at end, true kaons;neg2LL_MIP;",200,0,200);
  TrueBragg_truemu_neglogl_minmuMIP  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_minmuMIP","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);",200,0,200);
  TrueBragg_truep_neglogl_minmuMIP   = tfs->make<TH1F>("TrueBragg_truep_neglogl_minmuMIP","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);",200,0,200);
  TrueBragg_truepi_neglogl_minmuMIP  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_minmuMIP","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);",200,0,200);
  TrueBragg_trueK_neglogl_minmuMIP   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_minmuMIP","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);",200,0,200);

  TrueBragg_truemu_neglogl_mu_vslength = tfs->make<TH2F>("TrueBragg_truemu_neglogl_mu_vslength","Tracks with true p=0 at end, true muons;neg2LL_mu;Track length",200,0,200,500,0,500);
  TrueBragg_truep_neglogl_mu_vslength  = tfs->make<TH2F>("TrueBragg_truep_neglogl_mu_vslength","Tracks with true p=0 at end, true protons;neg2LL_mu;Track length",200,0,200,500,0,500);
  TrueBragg_truepi_neglogl_mu_vslength = tfs->make<TH2F>("TrueBragg_truepi_neglogl_mu_vslength","Tracks with true p=0 at end, true pions;neg2LL_mu;Track length",200,0,200,500,0,500);
  TrueBragg_trueK_neglogl_mu_vslength  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_mu_vslength","Tracks with true p=0 at end, true kaons;neg2LL_mu;Track length",200,0,200,500,0,500);
  TrueBragg_truemu_neglogl_MIP_vslength = tfs->make<TH2F>("TrueBragg_truemu_neglogl_MIP_vslength","Tracks with true p=0 at end, true muons;neg2LL_MIP;Track length",200,0,200,500,0,500);
  TrueBragg_truep_neglogl_MIP_vslength  = tfs->make<TH2F>("TrueBragg_truep_neglogl_MIP_vslength","Tracks with true p=0 at end, true protons;neg2LL_MIP;Track length",200,0,200,500,0,500);
  TrueBragg_truepi_neglogl_MIP_vslength = tfs->make<TH2F>("TrueBragg_truepi_neglogl_MIP_vslength","Tracks with true p=0 at end, true pions;neg2LL_MIP;Track length",200,0,200,500,0,500);
  TrueBragg_trueK_neglogl_MIP_vslength  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_MIP_vslength","Tracks with true p=0 at end, true kaons;neg2LL_MIP;Track length",200,0,200,500,0,500);
  TrueBragg_truemu_neglogl_minmuMIP_vslength = tfs->make<TH2F>("TrueBragg_truemu_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);Track length",200,0,200,500,0,500);
  TrueBragg_truep_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("TrueBragg_truep_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);Track length",200,0,200,500,0,500);
  TrueBragg_truepi_neglogl_minmuMIP_vslength = tfs->make<TH2F>("TrueBragg_truepi_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);Track length",200,0,200,500,0,500);
  TrueBragg_trueK_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);Track length",200,0,200,500,0,500);
  TrueBragg_truemu_neglogl_p_vslength = tfs->make<TH2F>("TrueBragg_truemu_neglogl_p_vslength","Tracks with true p=0 at end, true muons;neg2LL_p;Track length",200,0,200,500,0,500);
  TrueBragg_truep_neglogl_p_vslength  = tfs->make<TH2F>("TrueBragg_truep_neglogl_p_vslength","Tracks with true p=0 at end, true protons;neg2LL_p;Track length",200,0,200,500,0,500);
  TrueBragg_truepi_neglogl_p_vslength = tfs->make<TH2F>("TrueBragg_truepi_neglogl_p_vslength","Tracks with true p=0 at end, true pions;neg2LL_p;Track length",200,0,200,500,0,500);
  TrueBragg_trueK_neglogl_p_vslength  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_p_vslength","Tracks with true p=0 at end, true kaons;neg2LL_p;Track length",200,0,200,500,0,500);
  TrueBragg_truemu_neglogl_mu_vsangle = tfs->make<TH2F>("TrueBragg_truemu_neglogl_mu_vsangle","Tracks with true p=0 at end, true muons;neg2LL_mu;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truep_neglogl_mu_vsangle  = tfs->make<TH2F>("TrueBragg_truep_neglogl_mu_vsangle","Tracks with true p=0 at end, true protons;neg2LL_mu;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truepi_neglogl_mu_vsangle = tfs->make<TH2F>("TrueBragg_truepi_neglogl_mu_vsangle","Tracks with true p=0 at end, true pions;neg2LL_mu;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_trueK_neglogl_mu_vsangle  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_mu_vsangle","Tracks with true p=0 at end, true kaons;neg2LL_mu;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truemu_neglogl_MIP_vsangle = tfs->make<TH2F>("TrueBragg_truemu_neglogl_MIP_vsangle","Tracks with true p=0 at end, true muons;neg2LL_MIP;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truep_neglogl_MIP_vsangle  = tfs->make<TH2F>("TrueBragg_truep_neglogl_MIP_vsangle","Tracks with true p=0 at end, true protons;neg2LL_MIP;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truepi_neglogl_MIP_vsangle = tfs->make<TH2F>("TrueBragg_truepi_neglogl_MIP_vsangle","Tracks with true p=0 at end, true pions;neg2LL_MIP;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_trueK_neglogl_MIP_vsangle  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_MIP_vsangle","Tracks with true p=0 at end, true kaons;neg2LL_MIP;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truemu_neglogl_minmuMIP_vsangle = tfs->make<TH2F>("TrueBragg_truemu_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truep_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("TrueBragg_truep_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truepi_neglogl_minmuMIP_vsangle = tfs->make<TH2F>("TrueBragg_truepi_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_trueK_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truemu_neglogl_p_vsangle = tfs->make<TH2F>("TrueBragg_truemu_neglogl_p_vsangle","Tracks with true p=0 at end, true muons;neg2LL_p;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truep_neglogl_p_vsangle  = tfs->make<TH2F>("TrueBragg_truep_neglogl_p_vsangle","Tracks with true p=0 at end, true protons;neg2LL_p;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truepi_neglogl_p_vsangle = tfs->make<TH2F>("TrueBragg_truepi_neglogl_p_vsangle","Tracks with true p=0 at end, true pions;neg2LL_p;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_trueK_neglogl_p_vsangle  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_p_vsangle","Tracks with true p=0 at end, true kaons;neg2LL_p;Track angle",200,0,200,150,0,TMath::Pi());
  TrueBragg_truemu_neglogl_mu_vsnhits = tfs->make<TH2F>("TrueBragg_truemu_neglogl_mu_vsnhits","Tracks with true p=0 at end, true muons;neg2LL_mu;No. hits",200,0,200,500,0,500);
  TrueBragg_truep_neglogl_mu_vsnhits  = tfs->make<TH2F>("TrueBragg_truep_neglogl_mu_vsnhits","Tracks with true p=0 at end, true protons;neg2LL_mu;No. hits",200,0,200,500,0,500);
  TrueBragg_truepi_neglogl_mu_vsnhits = tfs->make<TH2F>("TrueBragg_truepi_neglogl_mu_vsnhits","Tracks with true p=0 at end, true pions;neg2LL_mu;No. hits",200,0,200,500,0,500);
  TrueBragg_trueK_neglogl_mu_vsnhits  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_mu_vsnhits","Tracks with true p=0 at end, true kaons;neg2LL_mu;No. hits",200,0,200,500,0,500);
  TrueBragg_truemu_neglogl_MIP_vsnhits = tfs->make<TH2F>("TrueBragg_truemu_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true muons;neg2LL_MIP;No. hits",200,0,200,500,0,500);
  TrueBragg_truep_neglogl_MIP_vsnhits  = tfs->make<TH2F>("TrueBragg_truep_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true protons;neg2LL_MIP;No. hits",200,0,200,500,0,500);
  TrueBragg_truepi_neglogl_MIP_vsnhits = tfs->make<TH2F>("TrueBragg_truepi_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true pions;neg2LL_MIP;No. hits",200,0,200,500,0,500);
  TrueBragg_trueK_neglogl_MIP_vsnhits  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true kaons;neg2LL_MIP;No. hits",200,0,200,500,0,500);
  TrueBragg_truemu_neglogl_minmuMIP_vsnhits = tfs->make<TH2F>("TrueBragg_truemu_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);No. hits",200,0,200,500,0,500);
  TrueBragg_truep_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("TrueBragg_truep_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);No. hits",200,0,200,500,0,500);
  TrueBragg_truepi_neglogl_minmuMIP_vsnhits = tfs->make<TH2F>("TrueBragg_truepi_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);No. hits",200,0,200,500,0,500);
  TrueBragg_trueK_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);No. hits",200,0,200,500,0,500);
  TrueBragg_truemu_neglogl_p_vsnhits = tfs->make<TH2F>("TrueBragg_truemu_neglogl_p_vsnhits","Tracks with true p=0 at end, true muons;neg2LL_p;No. hits",200,0,200,500,0,500);
  TrueBragg_truep_neglogl_p_vsnhits  = tfs->make<TH2F>("TrueBragg_truep_neglogl_p_vsnhits","Tracks with true p=0 at end, true protons;neg2LL_p;No. hits",200,0,200,500,0,500);
  TrueBragg_truepi_neglogl_p_vsnhits = tfs->make<TH2F>("TrueBragg_truepi_neglogl_p_vsnhits","Tracks with true p=0 at end, true pions;neg2LL_p;No. hits",200,0,200,500,0,500);
  TrueBragg_trueK_neglogl_p_vsnhits  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_p_vsnhits","Tracks with true p=0 at end, true kaons;neg2LL_p;No. hits",200,0,200,500,0,500);

  TrueBragg_truemu_neglogl_muvsp  = tfs->make<TH2F>("TrueBragg_truemu_neglogl_muvsp","Tracks with true p=0 at end, true muons;neg2LL_mu;neg2LL_p",200,0,200,200,0,200);
  TrueBragg_truep_neglogl_muvsp   = tfs->make<TH2F>("TrueBragg_truep_neglogl_muvsp","Tracks with true p=0 at end, true protons;neg2LL_mu;neg2LL_p",200,0,200,200,0,200);
  TrueBragg_truepi_neglogl_muvsp  = tfs->make<TH2F>("TrueBragg_truepi_neglogl_muvsp","Tracks with true p=0 at end, true pions;neg2LL_mu;neg2LL_p",200,0,200,200,0,200);
  TrueBragg_trueK_neglogl_muvsp   = tfs->make<TH2F>("TrueBragg_trueK_neglogl_muvsp","Tracks with true p=0 at end, true kaons;neg2LL_mu;neg2LL_p",200,0,200,200,0,200);
  TrueBragg_truemu_neglogl_muvspi  = tfs->make<TH2F>("TrueBragg_truemu_neglogl_muvspi","Tracks with true p=0 at end, true muons;neg2LL_mu;neg2LL_pi",200,0,200,200,0,200);
  TrueBragg_truep_neglogl_muvspi   = tfs->make<TH2F>("TrueBragg_truep_neglogl_muvspi","Tracks with true p=0 at end, true protons;neg2LL_mu;neg2LL_pi",200,0,200,200,0,200);
  TrueBragg_truepi_neglogl_muvspi  = tfs->make<TH2F>("TrueBragg_truepi_neglogl_muvspi","Tracks with true p=0 at end, true pions;neg2LL_mu;neg2LL_pi",200,0,200,200,0,200);
  TrueBragg_trueK_neglogl_muvspi   = tfs->make<TH2F>("TrueBragg_trueK_neglogl_muvspi","Tracks with true p=0 at end, true kaons;neg2LL_mu;neg2LL_pi",200,0,200,200,0,200);

  TrueBragg_truemu_neglogl_MIPvsp  = tfs->make<TH2F>("TrueBragg_truemu_neglogl_MIPvsp","Tracks with true p=0 at end, true muons;neg2LL_MIP;neg2LL_p",200,0,200,200,0,200);
  TrueBragg_truep_neglogl_MIPvsp   = tfs->make<TH2F>("TrueBragg_truep_neglogl_MIPvsp","Tracks with true p=0 at end, true protons;neg2LL_MIP;neg2LL_p",200,0,200,200,0,200);
  TrueBragg_truemu_neglogl_minmuMIPvsp  = tfs->make<TH2F>("TrueBragg_truemu_neglogl_minmuMIPvsp","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",200,0,200,200,0,200);
  TrueBragg_truep_neglogl_minmuMIPvsp   = tfs->make<TH2F>("TrueBragg_truep_neglogl_minmuMIPvsp","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",200,0,200,200,0,200);

  TrueBragg_truemu_neglogl_muoverp  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_muoverp","Tracks with true p=0 at end, true muons;neg2LL_mu/neg2LL_p;",200,0,20);
  TrueBragg_truep_neglogl_muoverp   = tfs->make<TH1F>("TrueBragg_truep_neglogl_muoverp","Tracks with true p=0 at end, true protons;neg2LL_mu/neg2LL_p;",200,0,20);
  TrueBragg_truepi_neglogl_muoverp  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_muoverp","Tracks with true p=0 at end, true pions;neg2LL_mu/neg2LL_p;",200,0,20);
  TrueBragg_trueK_neglogl_muoverp   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_muoverp","Tracks with true p=0 at end, true kaons;neg2LL_mu/neg2LL_p;",200,0,20);
  TrueBragg_truemu_neglogl_muminusp  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_muminusp","Tracks with true p=0 at end, true muons;neg2LL_mu-neg2LL_p;",200,-100,100);
  TrueBragg_truep_neglogl_muminusp   = tfs->make<TH1F>("TrueBragg_truep_neglogl_muminusp","Tracks with true p=0 at end, true protons;neg2LL_mu-neg2LL_p;",200,-100,100);
  TrueBragg_truepi_neglogl_muminusp  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_muminusp","Tracks with true p=0 at end, true pions;neg2LL_mu-neg2LL_p;",200,-100,100);
  TrueBragg_trueK_neglogl_muminusp   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_muminusp","Tracks with true p=0 at end, true kaons;neg2LL_mu-neg2LL_p;",200,-100,100);

  TrueBragg_truemu_neglogl_MIPminusp  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_MIPminusp","Tracks with true p=0 at end, true muons;neg2LL_MIP-neg2LL_p;",200,-100,100);
  TrueBragg_truep_neglogl_MIPminusp   = tfs->make<TH1F>("TrueBragg_truep_neglogl_MIPminusp","Tracks with true p=0 at end, true protons;neg2LL_MIP-neg2LL_p;",200,-100,100);
  TrueBragg_truemu_neglogl_minmuMIPminusp  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_minmuMIPminusp","Tracks with true p=0 at end, true muons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-100,100);
  TrueBragg_truep_neglogl_minmuMIPminusp   = tfs->make<TH1F>("TrueBragg_truep_neglogl_minmuMIPminusp","Tracks with true p=0 at end, true protons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-100,100);

  const char* particles[5] = {"#mu", "p", "#pi", "K", "MIP"};
  TrueBragg_truemu_smallest_neglogl  = tfs->make<TH1F>("TrueBragg_truemu_smallest_neglogl","Tracks with true p=0 at end, true muons;Particle type with smallest neg2LL;",5,0,5);
  TrueBragg_truep_smallest_neglogl   = tfs->make<TH1F>("TrueBragg_truep_smallest_neglogl","Tracks with true p=0 at end, true protons;Particle type with smallest neg2LL;",5,0,5);
  TrueBragg_truepi_smallest_neglogl  = tfs->make<TH1F>("TrueBragg_truepi_smallest_neglogl","Tracks with true p=0 at end, true pions;Particle type with smallest neg2LL;",5,0,5);
  TrueBragg_trueK_smallest_neglogl   = tfs->make<TH1F>("TrueBragg_trueK_smallest_neglogl","Tracks with true p=0 at end, true kaons;Particle type with smallest neg2LL;",5,0,5);
  for (size_t i=1; i<=5; i++){
    TrueBragg_truemu_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
    TrueBragg_truep_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
    TrueBragg_truepi_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
    TrueBragg_trueK_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
  }

  TrueBragg_truemu_PIDA = tfs->make<TH1F>("TrueBragg_truemu_PIDA","Tracks with true p=0 at end, true muons;PIDA;",300,0,30);
  TrueBragg_truep_PIDA  = tfs->make<TH1F>("TrueBragg_truep_PIDA","Tracks with true p=0 at end, true protons;PIDA;",300,0,30);
  TrueBragg_truepi_PIDA = tfs->make<TH1F>("TrueBragg_truepi_PIDA","Tracks with true p=0 at end, true pions;PIDA;",300,0,30);
  TrueBragg_trueK_PIDA  = tfs->make<TH1F>("TrueBragg_trueK_PIDA","Tracks with true p=0 at end, true kaons;PIDA;",300,0,30);

  TrueBragg_chargeEndOverStart_directionCorrect   = tfs->make<TH1F>("TrueBragg_chargeEndOverStart_directionCorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
  TrueBragg_chargeEndOverStart_directionIncorrect = tfs->make<TH1F>("TrueBragg_chargeEndOverStart_directionIncorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
  ;
  TrueBragg_chargeEndOverStartVersusNHits_directionCorrect   = tfs->make<TH2F>("TrueBragg_chargeEndOverStartVersusNHits_directionCorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};number of hits used in average", 200, 0, 10, 6, 0, 6);
  TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect = tfs->make<TH2F>("TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 200, 0, 10, 6, 0, 6);
  ;

  TrueBragg_chargeEndOverStart_sm0_5_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_sm0_5_dEdxrr","Tracks with true p=0 at end (end/start average charge < 0.5);Residual range (cm); dEdx",150,0,30,200,0,50);
  TrueBragg_chargeEndOverStart_gr2_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_gr2_dEdxrr","Tracks with true p=0 at end (end/start average charge > 2);Residual range (cm); dEdx",150,0,30,200,0,50);
  TrueBragg_chargeEndOverStart_0_5to2_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_0_5to2_dEdxrr","Tracks with true p=0 at end (end/start average charge = 0.5 - 2);Residual range (cm); dEdx",150,0,30,200,0,50);

  TrueBragg_truemu_dEdxtr_len = tfs->make<TH2F>("TrueBragg_truemu_dEdxtr_len","Tracks with true p=0 at end, true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
  TrueBragg_truep_dEdxtr_len  = tfs->make<TH2F>("TrueBragg_truep_dEdxtr_len","Tracks with true p=0 at end, true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
  TrueBragg_truepi_dEdxtr_len = tfs->make<TH2F>("TrueBragg_truepi_dEdxtr_len","Tracks with true p=0 at end, true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
  TrueBragg_trueK_dEdxtr_len  = tfs->make<TH2F>("TrueBragg_trueK_dEdxtr_len","Tracks with true p=0 at end, true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);

  TrueBragg_correctdirection = tfs->make<TH2F>("TrueBragg_correctdirection","Tracks with true p=0 at end;true particle;PID direction correct?",4,0,4,2,0,2);
  for (size_t i=1; i<=4; i++){
    TrueBragg_correctdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
  }
  TrueBragg_correctdirection->GetYaxis()->SetBinLabel(1,"False");
  TrueBragg_correctdirection->GetYaxis()->SetBinLabel(2,"True");

  // ---- All tracks: truth matching using associations
  All_truemu_neglogl_mu = tfs->make<TH1F>("All_truemu_neglogl_mu","All tracks, true muons;neg2LL_mu;",200,0,200);
  All_truep_neglogl_mu  = tfs->make<TH1F>("All_truep_neglogl_mu","All tracks, true protons;neg2LL_mu;",200,0,200);
  All_truepi_neglogl_mu = tfs->make<TH1F>("All_truepi_neglogl_mu","All tracks, true pions;neg2LL_mu;",200,0,200);
  All_trueK_neglogl_mu  = tfs->make<TH1F>("All_trueK_neglogl_mu","All tracks, true kaons;neg2LL_mu;",200,0,200);
  All_truemu_neglogl_p  = tfs->make<TH1F>("All_truemu_neglogl_p","All tracks, true muons;neg2LL_p;",200,0,200);
  All_truep_neglogl_p   = tfs->make<TH1F>("All_truep_neglogl_p","All tracks, true protons;neg2LL_p;",200,0,200);
  All_truepi_neglogl_p  = tfs->make<TH1F>("All_truepi_neglogl_p","All tracks, true pions;neg2LL_p;",200,0,200);
  All_trueK_neglogl_p   = tfs->make<TH1F>("All_trueK_neglogl_p","All tracks, true kaons;neg2LL_p;",200,0,200);
  All_truemu_neglogl_pi = tfs->make<TH1F>("All_truemu_neglogl_pi","All tracks, true muons;neg2LL_pi;",200,0,200);
  All_truep_neglogl_pi  = tfs->make<TH1F>("All_truep_neglogl_pi","All tracks, true protons;neg2LL_pi;",200,0,200);
  All_truepi_neglogl_pi = tfs->make<TH1F>("All_truepi_neglogl_pi","All tracks, true pions;neg2LL_pi;",200,0,200);
  All_trueK_neglogl_pi  = tfs->make<TH1F>("All_trueK_neglogl_pi","All tracks, true kaons;neg2LL_pi;",200,0,200);
  All_truemu_neglogl_K  = tfs->make<TH1F>("All_truemu_neglogl_K","All tracks, true muons;neg2LL_K;",200,0,200);
  All_truep_neglogl_K   = tfs->make<TH1F>("All_truep_neglogl_K","All tracks, true protons;neg2LL_K;",200,0,200);
  All_truepi_neglogl_K  = tfs->make<TH1F>("All_truepi_neglogl_K","All tracks, true pions;neg2LL_K;",200,0,200);
  All_trueK_neglogl_K   = tfs->make<TH1F>("All_trueK_neglogl_K","All tracks, true kaons;neg2LL_K;",200,0,200);
  All_truemu_neglogl_MIP  = tfs->make<TH1F>("All_truemu_neglogl_MIP","All tracks, true muons;neg2LL_MIP;",200,0,200);
  All_truep_neglogl_MIP   = tfs->make<TH1F>("All_truep_neglogl_MIP","All tracks, true protons;neg2LL_MIP;",200,0,200);
  All_truepi_neglogl_MIP  = tfs->make<TH1F>("All_truepi_neglogl_MIP","All tracks, true pions;neg2LL_MIP;",200,0,200);
  All_trueK_neglogl_MIP   = tfs->make<TH1F>("All_trueK_neglogl_MIP","All tracks, true kaons;neg2LL_MIP;",200,0,200);
  All_truemu_neglogl_minmuMIP  = tfs->make<TH1F>("All_truemu_neglogl_minmuMIP","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);",200,0,200);
  All_truep_neglogl_minmuMIP   = tfs->make<TH1F>("All_truep_neglogl_minmuMIP","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);",200,0,200);
  All_truepi_neglogl_minmuMIP  = tfs->make<TH1F>("All_truepi_neglogl_minmuMIP","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);",200,0,200);
  All_trueK_neglogl_minmuMIP   = tfs->make<TH1F>("All_trueK_neglogl_minmuMIP","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);",200,0,200);

  All_truemu_neglogl_mu_vslength = tfs->make<TH2F>("All_truemu_neglogl_mu_vslength","All tracks, true muons;neg2LL_mu;Track length",200,0,200,500,0,500);
  All_truep_neglogl_mu_vslength  = tfs->make<TH2F>("All_truep_neglogl_mu_vslength","All tracks, true protons;neg2LL_mu;Track length",200,0,200,500,0,500);
  All_truepi_neglogl_mu_vslength = tfs->make<TH2F>("All_truepi_neglogl_mu_vslength","All tracks, true pions;neg2LL_mu;Track length",200,0,200,500,0,500);
  All_trueK_neglogl_mu_vslength  = tfs->make<TH2F>("All_trueK_neglogl_mu_vslength","All tracks, true kaons;neg2LL_mu;Track length",200,0,200,500,0,500);
  All_truemu_neglogl_MIP_vslength = tfs->make<TH2F>("All_truemu_neglogl_MIP_vslength","All tracks, true muons;neg2LL_MIP;Track length",200,0,200,500,0,500);
  All_truep_neglogl_MIP_vslength  = tfs->make<TH2F>("All_truep_neglogl_MIP_vslength","All tracks, true protons;neg2LL_MIP;Track length",200,0,200,500,0,500);
  All_truepi_neglogl_MIP_vslength = tfs->make<TH2F>("All_truepi_neglogl_MIP_vslength","All tracks, true pions;neg2LL_MIP;Track length",200,0,200,500,0,500);
  All_trueK_neglogl_MIP_vslength  = tfs->make<TH2F>("All_trueK_neglogl_MIP_vslength","All tracks, true kaons;neg2LL_MIP;Track length",200,0,200,500,0,500);
  All_truemu_neglogl_minmuMIP_vslength = tfs->make<TH2F>("All_truemu_neglogl_minmuMIP_vslength","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);Track length",200,0,200,500,0,500);
  All_truep_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("All_truep_neglogl_minmuMIP_vslength","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);Track length",200,0,200,500,0,500);
  All_truepi_neglogl_minmuMIP_vslength = tfs->make<TH2F>("All_truepi_neglogl_minmuMIP_vslength","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);Track length",200,0,200,500,0,500);
  All_trueK_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("All_trueK_neglogl_minmuMIP_vslength","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);Track length",200,0,200,500,0,500);
  All_truemu_neglogl_p_vslength = tfs->make<TH2F>("All_truemu_neglogl_p_vslength","All tracks, true muons;neg2LL_p;Track length",200,0,200,500,0,500);
  All_truep_neglogl_p_vslength  = tfs->make<TH2F>("All_truep_neglogl_p_vslength","All tracks, true protons;neg2LL_p;Track length",200,0,200,500,0,500);
  All_truepi_neglogl_p_vslength = tfs->make<TH2F>("All_truepi_neglogl_p_vslength","All tracks, true pions;neg2LL_p;Track length",200,0,200,500,0,500);
  All_trueK_neglogl_p_vslength  = tfs->make<TH2F>("All_trueK_neglogl_p_vslength","All tracks, true kaons;neg2LL_p;Track length",200,0,200,500,0,500);
  All_truemu_neglogl_mu_vsangle = tfs->make<TH2F>("All_truemu_neglogl_mu_vsangle","All tracks, true muons;neg2LL_mu;Track angle",200,0,200,150,0,TMath::Pi());
  All_truep_neglogl_mu_vsangle  = tfs->make<TH2F>("All_truep_neglogl_mu_vsangle","All tracks, true protons;neg2LL_mu;Track angle",200,0,200,150,0,TMath::Pi());
  All_truepi_neglogl_mu_vsangle = tfs->make<TH2F>("All_truepi_neglogl_mu_vsangle","All tracks, true pions;neg2LL_mu;Track angle",200,0,200,150,0,TMath::Pi());
  All_trueK_neglogl_mu_vsangle  = tfs->make<TH2F>("All_trueK_neglogl_mu_vsangle","All tracks, true kaons;neg2LL_mu;Track angle",200,0,200,150,0,TMath::Pi());
  All_truemu_neglogl_MIP_vsangle = tfs->make<TH2F>("All_truemu_neglogl_MIP_vsangle","All tracks, true muons;neg2LL_MIP;Track angle",200,0,200,150,0,TMath::Pi());
  All_truep_neglogl_MIP_vsangle  = tfs->make<TH2F>("All_truep_neglogl_MIP_vsangle","All tracks, true protons;neg2LL_MIP;Track angle",200,0,200,150,0,TMath::Pi());
  All_truepi_neglogl_MIP_vsangle = tfs->make<TH2F>("All_truepi_neglogl_MIP_vsangle","All tracks, true pions;neg2LL_MIP;Track angle",200,0,200,150,0,TMath::Pi());
  All_trueK_neglogl_MIP_vsangle  = tfs->make<TH2F>("All_trueK_neglogl_MIP_vsangle","All tracks, true kaons;neg2LL_MIP;Track angle",200,0,200,150,0,TMath::Pi());
  All_truemu_neglogl_minmuMIP_vsangle = tfs->make<TH2F>("All_truemu_neglogl_minmuMIP_vsangle","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);Track angle",200,0,200,150,0,TMath::Pi());
  All_truep_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("All_truep_neglogl_minmuMIP_vsangle","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);Track angle",200,0,200,150,0,TMath::Pi());
  All_truepi_neglogl_minmuMIP_vsangle = tfs->make<TH2F>("All_truepi_neglogl_minmuMIP_vsangle","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);Track angle",200,0,200,150,0,TMath::Pi());
  All_trueK_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("All_trueK_neglogl_minmuMIP_vsangle","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);Track angle",200,0,200,150,0,TMath::Pi());
  All_truemu_neglogl_p_vsangle = tfs->make<TH2F>("All_truemu_neglogl_p_vsangle","All tracks, true muons;neg2LL_p;Track angle",200,0,200,150,0,TMath::Pi());
  All_truep_neglogl_p_vsangle  = tfs->make<TH2F>("All_truep_neglogl_p_vsangle","All tracks, true protons;neg2LL_p;Track angle",200,0,200,150,0,TMath::Pi());
  All_truepi_neglogl_p_vsangle = tfs->make<TH2F>("All_truepi_neglogl_p_vsangle","All tracks, true pions;neg2LL_p;Track angle",200,0,200,150,0,TMath::Pi());
  All_trueK_neglogl_p_vsangle  = tfs->make<TH2F>("All_trueK_neglogl_p_vsangle","All tracks, true kaons;neg2LL_p;Track angle",200,0,200,150,0,TMath::Pi());
  All_truemu_neglogl_mu_vsnhits = tfs->make<TH2F>("All_truemu_neglogl_mu_vsnhits","All tracks, true muons;neg2LL_mu;No. hits",200,0,200,500,0,500);
  All_truep_neglogl_mu_vsnhits  = tfs->make<TH2F>("All_truep_neglogl_mu_vsnhits","All tracks, true protons;neg2LL_mu;No. hits",200,0,200,500,0,500);
  All_truepi_neglogl_mu_vsnhits = tfs->make<TH2F>("All_truepi_neglogl_mu_vsnhits","All tracks, true pions;neg2LL_mu;No. hits",200,0,200,500,0,500);
  All_trueK_neglogl_mu_vsnhits  = tfs->make<TH2F>("All_trueK_neglogl_mu_vsnhits","All tracks, true kaons;neg2LL_mu;No. hits",200,0,200,500,0,500);
  All_truemu_neglogl_MIP_vsnhits = tfs->make<TH2F>("All_truemu_neglogl_MIP_vsnhits","All tracks, true muons;neg2LL_MIP;No. hits",200,0,200,500,0,500);
  All_truep_neglogl_MIP_vsnhits  = tfs->make<TH2F>("All_truep_neglogl_MIP_vsnhits","All tracks, true protons;neg2LL_MIP;No. hits",200,0,200,500,0,500);
  All_truepi_neglogl_MIP_vsnhits = tfs->make<TH2F>("All_truepi_neglogl_MIP_vsnhits","All tracks, true pions;neg2LL_MIP;No. hits",200,0,200,500,0,500);
  All_trueK_neglogl_MIP_vsnhits  = tfs->make<TH2F>("All_trueK_neglogl_MIP_vsnhits","All tracks, true kaons;neg2LL_MIP;No. hits",200,0,200,500,0,500);
  All_truemu_neglogl_minmuMIP_vsnhits = tfs->make<TH2F>("All_truemu_neglogl_minmuMIP_vsnhits","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);No. hits",200,0,200,500,0,500);
  All_truep_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("All_truep_neglogl_minmuMIP_vsnhits","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);No. hits",200,0,200,500,0,500);
  All_truepi_neglogl_minmuMIP_vsnhits = tfs->make<TH2F>("All_truepi_neglogl_minmuMIP_vsnhits","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);No. hits",200,0,200,500,0,500);
  All_trueK_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("All_trueK_neglogl_minmuMIP_vsnhits","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);No. hits",200,0,200,500,0,500);
  All_truemu_neglogl_p_vsnhits = tfs->make<TH2F>("All_truemu_neglogl_p_vsnhits","All tracks, true muons;neg2LL_p;No. hits",200,0,200,500,0,500);
  All_truep_neglogl_p_vsnhits  = tfs->make<TH2F>("All_truep_neglogl_p_vsnhits","All tracks, true protons;neg2LL_p;No. hits",200,0,200,500,0,500);
  All_truepi_neglogl_p_vsnhits = tfs->make<TH2F>("All_truepi_neglogl_p_vsnhits","All tracks, true pions;neg2LL_p;No. hits",200,0,200,500,0,500);
  All_trueK_neglogl_p_vsnhits  = tfs->make<TH2F>("All_trueK_neglogl_p_vsnhits","All tracks, true kaons;neg2LL_p;No. hits",200,0,200,500,0,500);

  All_truemu_neglogl_muvsp  = tfs->make<TH2F>("All_truemu_neglogl_muvsp","All tracks, true muons;neg2LL_mu;neg2LL_p",200,0,200,200,0,200);
  All_truep_neglogl_muvsp   = tfs->make<TH2F>("All_truep_neglogl_muvsp","All tracks, true protons;neg2LL_mu;neg2LL_p",200,0,200,200,0,200);
  All_truepi_neglogl_muvsp  = tfs->make<TH2F>("All_truepi_neglogl_muvsp","All tracks, true pions;neg2LL_mu;neg2LL_p",200,0,200,200,0,200);
  All_trueK_neglogl_muvsp   = tfs->make<TH2F>("All_trueK_neglogl_muvsp","All tracks, true kaons;neg2LL_mu;neg2LL_p",200,0,200,200,0,200);
  All_truemu_neglogl_muvspi  = tfs->make<TH2F>("All_truemu_neglogl_muvspi","All tracks, true muons;neg2LL_mu;neg2LL_pi",200,0,200,200,0,200);
  All_truep_neglogl_muvspi   = tfs->make<TH2F>("All_truep_neglogl_muvspi","All tracks, true protons;neg2LL_mu;neg2LL_pi",200,0,200,200,0,200);
  All_truepi_neglogl_muvspi  = tfs->make<TH2F>("All_truepi_neglogl_muvspi","All tracks, true pions;neg2LL_mu;neg2LL_pi",200,0,200,200,0,200);
  All_trueK_neglogl_muvspi   = tfs->make<TH2F>("All_trueK_neglogl_muvspi","All tracks, true kaons;neg2LL_mu;neg2LL_pi",200,0,200,200,0,200);

  All_truemu_neglogl_MIPvsp  = tfs->make<TH2F>("All_truemu_neglogl_MIPvsp","All tracks, true muons;neg2LL_MIP;neg2LL_p",200,0,200,200,0,200);
  All_truep_neglogl_MIPvsp   = tfs->make<TH2F>("All_truep_neglogl_MIPvsp","All tracks, true protons;neg2LL_MIP;neg2LL_p",200,0,200,200,0,200);
  All_truemu_neglogl_minmuMIPvsp  = tfs->make<TH2F>("All_truemu_neglogl_minmuMIPvsp","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",200,0,200,200,0,200);
  All_truep_neglogl_minmuMIPvsp   = tfs->make<TH2F>("All_truep_neglogl_minmuMIPvsp","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",200,0,200,200,0,200);

  All_truemu_neglogl_muoverp  = tfs->make<TH1F>("All_truemu_neglogl_muoverp","All tracks, true muons;neg2LL_mu/neg2LL_p;",200,0,50);
  All_truep_neglogl_muoverp   = tfs->make<TH1F>("All_truep_neglogl_muoverp","All tracks, true protons;neg2LL_mu/neg2LL_p;",200,0,50);
  All_truepi_neglogl_muoverp  = tfs->make<TH1F>("All_truepi_neglogl_muoverp","All tracks, true pions;neg2LL_mu/neg2LL_p;",200,0,50);
  All_trueK_neglogl_muoverp   = tfs->make<TH1F>("All_trueK_neglogl_muoverp","All tracks, true kaons;neg2LL_mu/neg2LL_p;",200,0,50);
  All_truemu_neglogl_muminusp  = tfs->make<TH1F>("All_truemu_neglogl_muminusp","All tracks, true muons;neg2LL_mu-neg2LL_p;",200,-100,100);
  All_truep_neglogl_muminusp   = tfs->make<TH1F>("All_truep_neglogl_muminusp","All tracks, true protons;neg2LL_mu-neg2LL_p;",200,-100,100);
  All_truepi_neglogl_muminusp  = tfs->make<TH1F>("All_truepi_neglogl_muminusp","All tracks, true pions;neg2LL_mu-neg2LL_p;",200,-100,100);
  All_trueK_neglogl_muminusp   = tfs->make<TH1F>("All_trueK_neglogl_muminusp","All tracks, true kaons;neg2LL_K;",200,-100,100);

  All_truemu_neglogl_MIPminusp  = tfs->make<TH1F>("All_truemu_neglogl_MIPminusp","All tracks, true muons;neg2LL_MIP-neg2LL_p;",200,-100,100);
  All_truep_neglogl_MIPminusp   = tfs->make<TH1F>("All_truep_neglogl_MIPminusp","All tracks, true protons;neg2LL_MIP-neg2LL_p;",200,-100,100);
  All_truemu_neglogl_minmuMIPminusp  = tfs->make<TH1F>("All_truemu_neglogl_minmuMIPminusp","All tracks, true muons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-100,100);
  All_truep_neglogl_minmuMIPminusp   = tfs->make<TH1F>("All_truep_neglogl_minmuMIPminusp","All tracks, true protons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-100,100);

  All_truemu_smallest_neglogl  = tfs->make<TH1F>("All_truemu_smallest_neglogl","All tracks, true muons;Particle type with smallest neg2LL;",5,0,5);
  All_truep_smallest_neglogl   = tfs->make<TH1F>("All_truep_smallest_neglogl","All tracks, true protons;Particle type with smallest neg2LL;",5,0,5);
  All_truepi_smallest_neglogl  = tfs->make<TH1F>("All_truepi_smallest_neglogl","All tracks, true pions;Particle type with smallest neg2LL;",5,0,5);
  All_trueK_smallest_neglogl   = tfs->make<TH1F>("All_trueK_smallest_neglogl","All tracks, true kaons;Particle type with smallest neg2LL;",5,0,5);
  for (size_t i=1; i<=5; i++){
    All_truemu_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
    All_truep_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
    All_truepi_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
    All_trueK_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
  }

  All_truemu_PIDA = tfs->make<TH1F>("All_truemu_PIDA","All tracks, true muons;PIDA;",300,0,30);
  All_truep_PIDA  = tfs->make<TH1F>("All_truep_PIDA","All tracks, true protons;PIDA;",300,0,30);
  All_truepi_PIDA = tfs->make<TH1F>("All_truepi_PIDA","All tracks, true pions;PIDA;",300,0,30);
  All_trueK_PIDA  = tfs->make<TH1F>("All_trueK_PIDA","All tracks, true kaons;PIDA;",300,0,30);

  All_chargeEndOverStart_directionCorrect   = tfs->make<TH1F>("All_chargeEndOverStart_directionCorrect", "All tracks;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
  All_chargeEndOverStart_directionIncorrect = tfs->make<TH1F>("All_chargeEndOverStart_directionIncorrect", "All tracks;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
  ;
  All_chargeEndOverStartVersusNHits_directionCorrect   = tfs->make<TH2F>("All_chargeEndOverStartVersusNHits_directionCorrect", "All tracks;Charge_{End of track}/Charge_{Start of track};number of hits used in average", 200, 0, 10, 6, 0, 6);
  All_chargeEndOverStartVersusNHits_directionIncorrect = tfs->make<TH2F>("All_chargeEndOverStartVersusNHits_directionIncorrect", "All tracks;Charge_{End of track}/Charge_{Start of track};", 200, 0, 10, 6, 0, 6);
  ;

  Contained_chargeEndOverStart_directionCorrect   = tfs->make<TH1F>("Contained_chargeEndOverStart_directionCorrect", "Contained tracks only;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
  Contained_chargeEndOverStart_directionIncorrect = tfs->make<TH1F>("Contained_chargeEndOverStart_directionIncorrect", "Contained tracks only;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
  ;
  Contained_chargeEndOverStartVersusNHits_directionCorrect   = tfs->make<TH2F>("Contained_chargeEndOverStartVersusNHits_directionCorrect", "Contained tracks only;Charge_{End of track}/Charge_{Start of track};number of hits used in average", 200, 0, 10, 6, 0, 6);
  Contained_chargeEndOverStartVersusNHits_directionIncorrect = tfs->make<TH2F>("Contained_chargeEndOverStartVersusNHits_directionIncorrect", "Contained tracks only;Charge_{End of track}/Charge_{Start of track};", 200, 0, 10, 6, 0, 6);
  ;

  All_chargeEndOverStart_sm0_5_dEdxrr = tfs->make<TH2F>("All_chargeEndOverStart_sm0_5_dEdxrr","All tracks (end/start average charge < 0.5);Residual range (cm); dEdx",150,0,30,200,0,50);
  All_chargeEndOverStart_gr2_dEdxrr = tfs->make<TH2F>("All_chargeEndOverStart_gr2_dEdxrr","All tracks (end/start average charge > 2);Residual range (cm); dEdx",150,0,30,200,0,50);
  All_chargeEndOverStart_0_5to2_dEdxrr = tfs->make<TH2F>("All_chargeEndOverStart_0_5to2_dEdxrr","All tracks (end/start average charge = 0.5 - 2);Residual range (cm); dEdx",150,0,30,200,0,50);

  All_truemu_dEdxtr_len = tfs->make<TH2F>("All_truemu_dEdxtr_len","All tracks, true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_truep_dEdxtr_len  = tfs->make<TH2F>("All_truep_dEdxtr_len","All tracks, true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_truepi_dEdxtr_len = tfs->make<TH2F>("All_truepi_dEdxtr_len","All tracks, true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
  All_trueK_dEdxtr_len  = tfs->make<TH2F>("All_trueK_dEdxtr_len","All tracks, true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);


  All_correctdirection = tfs->make<TH2F>("All_correctdirection","All tracks;true particle;PID direction correct?",4,0,4,2,0,2);
  for (size_t i=1; i<=4; i++){
    All_correctdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
  }
  All_correctdirection->GetYaxis()->SetBinLabel(1,"False");
  All_correctdirection->GetYaxis()->SetBinLabel(2,"True");

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
  art::FindManyP<anab::Calorimetry> calo_from_tracks(trackHandle, e, "pandoraNucali");

  // --------- Loop over tracks in event ---------- //
  for (auto& track : trackCollection){
    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = calo_from_tracks.at(track->ID());

    double angle = track->Theta();
    std::cout << "Looking at track with length " << track->Length() << std::endl;

    // Get true PDG from associations
    int True_pdg = 0;

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

    True_pdg = maxp_me->PdgCode();

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

    // Check if true particle should have a Bragg peak (if the true MCParticle has zero momentum at the end of the track)
    bool TrueBragg = false;
    if (maxp_me->EndPx() == 0 && maxp_me->EndPy() == 0 && maxp_me->EndPz() == 0){
      TrueBragg = true;
    }
    //std::cout << "True particle, EndPx = " << maxp_me->EndPx() << ", EndPy = " << maxp_me->EndPy() << ", EndPz = " << maxp_me->EndPz() << ", TrueBragg = " << TrueBragg << std::endl;

    // ------------------- Get Track end hits over track start hits charge -------------- //
    // for time being, only use Y plane calorimetry
    art::Ptr< anab:: Calorimetry > calo;
    for (auto c : caloFromTrack){
      int planenum = c->PlaneID().Plane;
      if (planenum != 2) continue; // Only use calorimetry from collection plane
      calo = c;
    }
      // Check that caloFromTrack is a valid object? If not, skip track
    if (!calo){
      std::cout << "Did not find a valid calorimetry object. Skipping track." << std::endl;
      continue;
    }
    std::vector<double> dEdx = calo->dEdx();
    std::vector<double> resRange = calo->ResidualRange();

    double nhits = resRange.size();

    // find how many hits to use
    int nHits = dEdx.size();
    int hitsToUse;
    if (nHits >= 10)
      hitsToUse = 5;
    else{
      hitsToUse = std::floor((double)nHits/2.0);
    }

    double averagedEdxTrackStart=0;
    double averagedEdxTrackEnd=0;

    // loop dEdx and take average of first n hits and last n hits
    for (int i = 0; i < (int)dEdx.size(); i++){

      if (i < hitsToUse) averagedEdxTrackStart+=dEdx.at(i);
      else if (i > nHits - hitsToUse -1) averagedEdxTrackEnd+=dEdx.at(i);

    }

    averagedEdxTrackStart = averagedEdxTrackStart/hitsToUse;
    averagedEdxTrackEnd   = averagedEdxTrackEnd/hitsToUse;

    double dEdxStartEndRatio = averagedEdxTrackEnd/averagedEdxTrackStart;

    // ------------------- Make plots of dEdx vs residual range for different track start/end dEdx ratios ------------------- //

    if (dEdxStartEndRatio < 0.5){
      for (int i=0; i < (int)resRange.size(); i++){
        All_chargeEndOverStart_sm0_5_dEdxrr->Fill(resRange.at(i),dEdx.at(i));
        if (TrueBragg){
          TrueBragg_chargeEndOverStart_sm0_5_dEdxrr->Fill(resRange.at(i),dEdx.at(i));
        }
      }
    }
    else if (dEdxStartEndRatio > 2.0){
      for (int i=0; i < (int)resRange.size(); i++){
        All_chargeEndOverStart_gr2_dEdxrr->Fill(resRange.at(i),dEdx.at(i));
        if (TrueBragg){
          TrueBragg_chargeEndOverStart_gr2_dEdxrr->Fill(resRange.at(i),dEdx.at(i));
        }
      }

    }
    else {
      for (int i=0; i < (int)resRange.size(); i++){
        All_chargeEndOverStart_0_5to2_dEdxrr->Fill(resRange.at(i),dEdx.at(i));
        if (TrueBragg){
          TrueBragg_chargeEndOverStart_0_5to2_dEdxrr->Fill(resRange.at(i),dEdx.at(i));
        }
      }
    }
    //std::cout << dEdxStartEndRatio << std::endl;

    // ------------------- Check track direction ---------------------------------------- //

    TVector3 RecoTrackDir(track->End().X()-track->Start().X(),
                          track->End().Y()-track->Start().Y(),
                          track->End().Z()-track->Start().Z());
    RecoTrackDir = RecoTrackDir.Unit();

    TVector3 True_start(maxp_me->Vx(), maxp_me->Vy(), maxp_me->Vz());
    TVector3 True_end(maxp_me->EndX(), maxp_me->EndY(), maxp_me->EndZ());
    TVector3 TrueTrackDir = True_end - True_start;
    TrueTrackDir = TrueTrackDir.Unit();

    /*std::cout << "RecoTrackDir = " ;
    RecoTrackDir.Print();
    std::cout << std::endl
              << "TrueTrackDir = " ;
    TrueTrackDir.Print();
    std::cout << std::endl
              << "RecoTrackDir.Dot(TrueTrackDir) = " << RecoTrackDir.Dot(TrueTrackDir) << std::endl;*/

    // If dot product < 0, track direction is wrong
    if (RecoTrackDir.Dot(TrueTrackDir) < 0){
      All_chargeEndOverStart_directionIncorrect->Fill(dEdxStartEndRatio);
      All_chargeEndOverStartVersusNHits_directionIncorrect->Fill(dEdxStartEndRatio, hitsToUse);

      if (TrueBragg){
        TrueBragg_chargeEndOverStart_directionIncorrect->Fill(dEdxStartEndRatio);
        TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect->Fill(dEdxStartEndRatio, hitsToUse);
      }


      if (fid.isInFiducialVolume(True_start, fv) && fid.isInFiducialVolume(True_end, fv)){
        Contained_chargeEndOverStart_directionIncorrect->Fill(dEdxStartEndRatio);
        Contained_chargeEndOverStartVersusNHits_directionIncorrect->Fill(dEdxStartEndRatio, hitsToUse);
      }
    }
    else{ // else, track direction is correct
      All_chargeEndOverStart_directionCorrect->Fill(dEdxStartEndRatio);
      All_chargeEndOverStartVersusNHits_directionCorrect->Fill(dEdxStartEndRatio, hitsToUse);

      if (TrueBragg){
        TrueBragg_chargeEndOverStart_directionCorrect->Fill(dEdxStartEndRatio);
        TrueBragg_chargeEndOverStartVersusNHits_directionCorrect->Fill(dEdxStartEndRatio, hitsToUse);
      }

      if (fid.isInFiducialVolume(True_start, fv) && fid.isInFiducialVolume(True_end, fv)){
        Contained_chargeEndOverStart_directionCorrect->Fill(dEdxStartEndRatio);
        Contained_chargeEndOverStartVersusNHits_directionCorrect->Fill(dEdxStartEndRatio, hitsToUse);
      }
    }



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
    double noBragg_fwd_MIP = -999;
    double PIDAval = -999;
    double dEdxtruncmean = -999;
    double trklen = -999;

    std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track->ID());
    if (trackPID.size() == 0){
      std::cout << "No track-PID association found for trackID " << track->ID() << ". Skipping track." << std::endl;
      continue;
    }
    //std::cout << "trackPID.size() = " << trackPID.size() << std::endl;
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
          if (AlgScore.fAssumedPdg == 0)    noBragg_fwd_MIP = AlgScore.fValue;
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
        if (anab::kVariableType(AlgScore.fVariableType) == anab::kdEdxtruncmean) dEdxtruncmean = AlgScore.fValue;
        if (anab::kVariableType(AlgScore.fVariableType) == anab::kTrackLength)   trklen = AlgScore.fValue;
      }// if AlgName = TruncatedMean

    } // Loop over AlgScoresVec

    // Some couts for debugging
    // std::cout << "From analyzer module:" << std::endl
    //   << "neg2LogL mu fwd = " << Bragg_fwd_mu << std::endl
    //   << "neg2LogL p fwd = " << Bragg_fwd_p << std::endl
    //   << "neg2LogL pi fwd = " << Bragg_fwd_pi << std::endl
    //   << "neg2LogL K fwd = " << Bragg_fwd_K << std::endl
    //   << "neg2LogL mu bwd = " << Bragg_bwd_mu << std::endl
    //   << "neg2LogL p bwd = " << Bragg_bwd_p << std::endl
    //   << "neg2LogL pi bwd = " << Bragg_bwd_pi << std::endl
    //   << "neg2LogL K bwd = " << Bragg_bwd_K << std::endl;


    // ---------------- Now let's fill some histograms! ------------------- //

    // Use best fit (lowest neg2LogL) from forward vs backward for calculations
    double Bragg_mu = (Bragg_fwd_mu < Bragg_bwd_mu ? Bragg_fwd_mu : Bragg_bwd_mu);
    double Bragg_p  = (Bragg_fwd_p  < Bragg_bwd_p  ? Bragg_fwd_p  : Bragg_bwd_p);
    double Bragg_pi = (Bragg_fwd_pi < Bragg_bwd_pi ? Bragg_fwd_pi : Bragg_bwd_pi);
    double Bragg_K  = (Bragg_fwd_K  < Bragg_bwd_K  ? Bragg_fwd_K  : Bragg_bwd_K);

    //double Bragg_smallest = std::min({Bragg_mu, Bragg_p, Bragg_pi, Bragg_K, noBragg_fwd_MIP});
    double Bragg_smallest = std::min({Bragg_mu, Bragg_p, noBragg_fwd_MIP});
    if (Bragg_smallest == Bragg_mu) std::cout << "PID estimate: muon" << std::endl;
    else if (Bragg_smallest == Bragg_p) std::cout << "PID estimate: proton" << std::endl;
    else if (Bragg_smallest == Bragg_pi) std::cout << "PID estimate: pion" << std::endl;
    else if (Bragg_smallest == Bragg_K) std::cout << "PID estimate: kaon" << std::endl;
    else if (Bragg_smallest == noBragg_fwd_MIP) std::cout << "PID estimate: muon (no Bragg peak)" << std::endl;

    bool PID_fwd = false;
    if (Bragg_smallest == Bragg_fwd_mu || Bragg_smallest == Bragg_fwd_p || Bragg_smallest == Bragg_fwd_pi || Bragg_smallest == Bragg_fwd_K || Bragg_smallest == noBragg_fwd_MIP) PID_fwd = true;

    if (TrueBragg){
      // Well-reconstructed truth matching
      if (TMath::Abs(True_pdg) == 13){ // True muons
        TrueBragg_truemu_neglogl_mu->Fill(Bragg_mu);
        TrueBragg_truemu_neglogl_p->Fill(Bragg_p);
        TrueBragg_truemu_neglogl_pi->Fill(Bragg_pi);
        TrueBragg_truemu_neglogl_K->Fill(Bragg_K);
        TrueBragg_truemu_PIDA->Fill(PIDAval);
        TrueBragg_truemu_dEdxtr_len->Fill(trklen,dEdxtruncmean);
        TrueBragg_truemu_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
        TrueBragg_truemu_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
        TrueBragg_truemu_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
        TrueBragg_truemu_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
        TrueBragg_truemu_neglogl_MIPvsp->Fill(noBragg_fwd_MIP,Bragg_p);
        TrueBragg_truemu_neglogl_minmuMIPvsp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),Bragg_p);
        TrueBragg_truemu_neglogl_MIPminusp->Fill(noBragg_fwd_MIP-Bragg_p);
        TrueBragg_truemu_neglogl_minmuMIPminusp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP)-Bragg_p);
        TrueBragg_truemu_neglogl_MIP->Fill(noBragg_fwd_MIP);
        TrueBragg_truemu_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
        TrueBragg_truemu_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
        TrueBragg_truemu_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
        TrueBragg_truemu_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
        TrueBragg_truemu_neglogl_p_vslength->Fill(Bragg_p,trklen);
        TrueBragg_truemu_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
        TrueBragg_truemu_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
        TrueBragg_truemu_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
        TrueBragg_truemu_neglogl_p_vsangle->Fill(Bragg_p,angle);
        TrueBragg_truemu_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
        TrueBragg_truemu_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
        TrueBragg_truemu_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
        TrueBragg_truemu_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

        if (Bragg_smallest == Bragg_mu) TrueBragg_truemu_smallest_neglogl->Fill(0.5);
        else if (Bragg_smallest == Bragg_p) TrueBragg_truemu_smallest_neglogl->Fill(1.5);
        else if (Bragg_smallest == Bragg_pi) TrueBragg_truemu_smallest_neglogl->Fill(2.5);
        else if (Bragg_smallest == Bragg_K) TrueBragg_truemu_smallest_neglogl->Fill(3.5);
        else if (Bragg_smallest == noBragg_fwd_MIP)
        TrueBragg_truemu_smallest_neglogl->Fill(4.5);

        if ((PID_fwd && RecoTrackDir.Dot(TrueTrackDir) > 0) || (!PID_fwd && RecoTrackDir.Dot(TrueTrackDir) < 0)){ // correct track direction
          TrueBragg_correctdirection->Fill(0.5,1.5);
        }
        else{ // right track direction
          TrueBragg_correctdirection->Fill(0.5,0.5);
        }
      }

      else if (TMath::Abs(True_pdg) == 2212){ // True protons
        TrueBragg_truep_neglogl_mu->Fill(Bragg_mu);
        TrueBragg_truep_neglogl_p->Fill(Bragg_p);
        TrueBragg_truep_neglogl_pi->Fill(Bragg_pi);
        TrueBragg_truep_neglogl_K->Fill(Bragg_K);
        TrueBragg_truep_PIDA->Fill(PIDAval);
        TrueBragg_truep_dEdxtr_len->Fill(trklen,dEdxtruncmean);
        TrueBragg_truep_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
        TrueBragg_truep_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
        TrueBragg_truep_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
        TrueBragg_truep_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
        TrueBragg_truep_neglogl_MIPvsp->Fill(noBragg_fwd_MIP,Bragg_p);
        TrueBragg_truep_neglogl_minmuMIPvsp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),Bragg_p);
        TrueBragg_truep_neglogl_MIPminusp->Fill(noBragg_fwd_MIP-Bragg_p);
        TrueBragg_truep_neglogl_minmuMIPminusp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP)-Bragg_p);
        TrueBragg_truep_neglogl_MIP->Fill(noBragg_fwd_MIP);
        TrueBragg_truep_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
        TrueBragg_truep_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
        TrueBragg_truep_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
        TrueBragg_truep_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
        TrueBragg_truep_neglogl_p_vslength->Fill(Bragg_p,trklen);
        TrueBragg_truep_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
        TrueBragg_truep_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
        TrueBragg_truep_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
        TrueBragg_truep_neglogl_p_vsangle->Fill(Bragg_p,angle);
        TrueBragg_truep_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
        TrueBragg_truep_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
        TrueBragg_truep_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
        TrueBragg_truep_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

        if (Bragg_smallest == Bragg_mu) TrueBragg_truep_smallest_neglogl->Fill(0.5);
        else if (Bragg_smallest == Bragg_p) TrueBragg_truep_smallest_neglogl->Fill(1.5);
        else if (Bragg_smallest == Bragg_pi) TrueBragg_truep_smallest_neglogl->Fill(2.5);
        else if (Bragg_smallest == Bragg_K) TrueBragg_truep_smallest_neglogl->Fill(3.5);
        else if (Bragg_smallest == noBragg_fwd_MIP)
        TrueBragg_truep_smallest_neglogl->Fill(4.5);

        if ((PID_fwd && RecoTrackDir.Dot(TrueTrackDir) > 0) || (!PID_fwd && RecoTrackDir.Dot(TrueTrackDir) < 0)){ // correct track direction
          TrueBragg_correctdirection->Fill(1.5,1.5);
        }
        else{ // right track direction
          TrueBragg_correctdirection->Fill(1.5,0.5);
        }
      }

      else if (TMath::Abs(True_pdg) == 211){ // True pions
        TrueBragg_truepi_neglogl_mu->Fill(Bragg_mu);
        TrueBragg_truepi_neglogl_p->Fill(Bragg_p);
        TrueBragg_truepi_neglogl_pi->Fill(Bragg_pi);
        TrueBragg_truepi_neglogl_K->Fill(Bragg_K);
        TrueBragg_truepi_PIDA->Fill(PIDAval);
        TrueBragg_truepi_dEdxtr_len->Fill(trklen,dEdxtruncmean);
        TrueBragg_truepi_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
        TrueBragg_truepi_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
        TrueBragg_truepi_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
        TrueBragg_truepi_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
        TrueBragg_truepi_neglogl_MIP->Fill(noBragg_fwd_MIP);
        TrueBragg_truepi_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
        TrueBragg_truepi_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
        TrueBragg_truepi_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
        TrueBragg_truepi_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
        TrueBragg_truepi_neglogl_p_vslength->Fill(Bragg_p,trklen);
        TrueBragg_truepi_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
        TrueBragg_truepi_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
        TrueBragg_truepi_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
        TrueBragg_truepi_neglogl_p_vsangle->Fill(Bragg_p,angle);
        TrueBragg_truepi_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
        TrueBragg_truepi_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
        TrueBragg_truepi_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
        TrueBragg_truepi_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

        if (Bragg_smallest == Bragg_mu) TrueBragg_truepi_smallest_neglogl->Fill(0.5);
        else if (Bragg_smallest == Bragg_p) TrueBragg_truepi_smallest_neglogl->Fill(1.5);
        else if (Bragg_smallest == Bragg_pi) TrueBragg_truepi_smallest_neglogl->Fill(2.5);
        else if (Bragg_smallest == Bragg_K) TrueBragg_truepi_smallest_neglogl->Fill(3.5);
        else if (Bragg_smallest == noBragg_fwd_MIP)
        TrueBragg_truepi_smallest_neglogl->Fill(4.5);

        if ((PID_fwd && RecoTrackDir.Dot(TrueTrackDir) > 0) || (!PID_fwd && RecoTrackDir.Dot(TrueTrackDir) < 0)){ // correct track direction
          TrueBragg_correctdirection->Fill(2.5,1.5);
        }
        else{ // right track direction
          TrueBragg_correctdirection->Fill(2.5,0.5);
        }
      }

      else if (TMath::Abs(True_pdg) == 321){ // True kaons
        TrueBragg_trueK_neglogl_mu->Fill(Bragg_mu);
        TrueBragg_trueK_neglogl_p->Fill(Bragg_p);
        TrueBragg_trueK_neglogl_pi->Fill(Bragg_pi);
        TrueBragg_trueK_neglogl_K->Fill(Bragg_K);
        TrueBragg_trueK_PIDA->Fill(PIDAval);
        TrueBragg_trueK_dEdxtr_len->Fill(trklen,dEdxtruncmean);
        TrueBragg_trueK_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
        TrueBragg_trueK_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
        TrueBragg_trueK_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
        TrueBragg_trueK_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
        TrueBragg_trueK_neglogl_MIP->Fill(noBragg_fwd_MIP);
        TrueBragg_trueK_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
        TrueBragg_trueK_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
        TrueBragg_trueK_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
        TrueBragg_trueK_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
        TrueBragg_trueK_neglogl_p_vslength->Fill(Bragg_p,trklen);
        TrueBragg_trueK_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
        TrueBragg_trueK_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
        TrueBragg_trueK_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
        TrueBragg_trueK_neglogl_p_vsangle->Fill(Bragg_p,angle);
        TrueBragg_trueK_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
        TrueBragg_trueK_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
        TrueBragg_trueK_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
        TrueBragg_trueK_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

        if (Bragg_smallest == Bragg_mu) TrueBragg_trueK_smallest_neglogl->Fill(0.5);
        else if (Bragg_smallest == Bragg_p) TrueBragg_trueK_smallest_neglogl->Fill(1.5);
        else if (Bragg_smallest == Bragg_pi) TrueBragg_trueK_smallest_neglogl->Fill(2.5);
        else if (Bragg_smallest == Bragg_K) TrueBragg_trueK_smallest_neglogl->Fill(3.5);
        else if (Bragg_smallest == noBragg_fwd_MIP)
        TrueBragg_trueK_smallest_neglogl->Fill(4.5);

        if ((PID_fwd && RecoTrackDir.Dot(TrueTrackDir) > 0) || (!PID_fwd && RecoTrackDir.Dot(TrueTrackDir) < 0)){ // correct track direction
          TrueBragg_correctdirection->Fill(3.5,1.5);
        }
        else{ // right track direction
          TrueBragg_correctdirection->Fill(3.5,0.5);
        }
      }
    } // end if(TrueBragg)

    // All particles
    if (TMath::Abs(True_pdg) == 13){ // True muons
      All_truemu_neglogl_mu->Fill(Bragg_mu);
      All_truemu_neglogl_p->Fill(Bragg_p);
      All_truemu_neglogl_pi->Fill(Bragg_pi);
      All_truemu_neglogl_K->Fill(Bragg_K);
      All_truemu_PIDA->Fill(PIDAval);
      All_truemu_dEdxtr_len->Fill(trklen,dEdxtruncmean);
      All_truemu_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
      All_truemu_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
      All_truemu_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
      All_truemu_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
      All_truemu_neglogl_MIPvsp->Fill(noBragg_fwd_MIP,Bragg_p);
      All_truemu_neglogl_minmuMIPvsp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),Bragg_p);
      All_truemu_neglogl_MIPminusp->Fill(noBragg_fwd_MIP-Bragg_p);
      All_truemu_neglogl_minmuMIPminusp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP)-Bragg_p);
      All_truemu_neglogl_MIP->Fill(noBragg_fwd_MIP);
      All_truemu_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
      All_truemu_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
      All_truemu_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
      All_truemu_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
      All_truemu_neglogl_p_vslength->Fill(Bragg_p,trklen);
      All_truemu_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
      All_truemu_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
      All_truemu_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
      All_truemu_neglogl_p_vsangle->Fill(Bragg_p,angle);
      All_truemu_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
      All_truemu_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
      All_truemu_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
      All_truemu_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

      if (Bragg_smallest == Bragg_mu) All_truemu_smallest_neglogl->Fill(0.5);
      else if (Bragg_smallest == Bragg_p) All_truemu_smallest_neglogl->Fill(1.5);
      else if (Bragg_smallest == Bragg_pi) All_truemu_smallest_neglogl->Fill(2.5);
      else if (Bragg_smallest == Bragg_K) All_truemu_smallest_neglogl->Fill(3.5);
      else if (Bragg_smallest == noBragg_fwd_MIP)
      All_truemu_smallest_neglogl->Fill(4.5);

      if ((PID_fwd && RecoTrackDir.Dot(TrueTrackDir) > 0) || (!PID_fwd && RecoTrackDir.Dot(TrueTrackDir) < 0)){ // correct track direction
        All_correctdirection->Fill(0.5,1.5);
      }
      else{ // right track direction
        All_correctdirection->Fill(0.5,0.5);
      }
    }

    else if (TMath::Abs(True_pdg) == 2212){ // True protons
      All_truep_neglogl_mu->Fill(Bragg_mu);
      All_truep_neglogl_p->Fill(Bragg_p);
      All_truep_neglogl_pi->Fill(Bragg_pi);
      All_truep_neglogl_K->Fill(Bragg_K);
      All_truep_PIDA->Fill(PIDAval);
      All_truep_dEdxtr_len->Fill(trklen,dEdxtruncmean);
      All_truep_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
      All_truep_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
      All_truep_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
      All_truep_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
      All_truep_neglogl_MIPvsp->Fill(noBragg_fwd_MIP,Bragg_p);
      All_truep_neglogl_minmuMIPvsp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),Bragg_p);
      All_truep_neglogl_MIPminusp->Fill(noBragg_fwd_MIP-Bragg_p);
      All_truep_neglogl_minmuMIPminusp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP)-Bragg_p);
      All_truep_neglogl_MIP->Fill(noBragg_fwd_MIP);
      All_truep_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
      All_truep_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
      All_truep_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
      All_truep_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
      All_truep_neglogl_p_vslength->Fill(Bragg_p,trklen);
      All_truep_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
      All_truep_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
      All_truep_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
      All_truep_neglogl_p_vsangle->Fill(Bragg_p,angle);
      All_truep_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
      All_truep_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
      All_truep_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
      All_truep_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

      if (Bragg_smallest == Bragg_mu) All_truep_smallest_neglogl->Fill(0.5);
      else if (Bragg_smallest == Bragg_p) All_truep_smallest_neglogl->Fill(1.5);
      else if (Bragg_smallest == Bragg_pi) All_truep_smallest_neglogl->Fill(2.5);
      else if (Bragg_smallest == Bragg_K) All_truep_smallest_neglogl->Fill(3.5);
      else if (Bragg_smallest == noBragg_fwd_MIP)
      All_truep_smallest_neglogl->Fill(4.5);

      if ((PID_fwd && RecoTrackDir.Dot(TrueTrackDir) > 0) || (!PID_fwd && RecoTrackDir.Dot(TrueTrackDir) < 0)){ // correct track direction
        All_correctdirection->Fill(1.5,1.5);
      }
      else{ // right track direction
        All_correctdirection->Fill(1.5,0.5);
      }
    }

    else if (TMath::Abs(True_pdg) == 211){ // True pions
      All_truepi_neglogl_mu->Fill(Bragg_mu);
      All_truepi_neglogl_p->Fill(Bragg_p);
      All_truepi_neglogl_pi->Fill(Bragg_pi);
      All_truepi_neglogl_K->Fill(Bragg_K);
      All_truepi_PIDA->Fill(PIDAval);
      All_truepi_dEdxtr_len->Fill(trklen,dEdxtruncmean);
      All_truepi_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
      All_truepi_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
      All_truepi_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
      All_truepi_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
      All_truepi_neglogl_MIP->Fill(noBragg_fwd_MIP);
      All_truepi_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
      All_truepi_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
      All_truepi_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
      All_truepi_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
      All_truepi_neglogl_p_vslength->Fill(Bragg_p,trklen);
      All_truepi_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
      All_truepi_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
      All_truepi_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
      All_truepi_neglogl_p_vsangle->Fill(Bragg_p,angle);
      All_truepi_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
      All_truepi_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
      All_truepi_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
      All_truepi_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

      if (Bragg_smallest == Bragg_mu) All_truepi_smallest_neglogl->Fill(0.5);
      else if (Bragg_smallest == Bragg_p) All_truepi_smallest_neglogl->Fill(1.5);
      else if (Bragg_smallest == Bragg_pi) All_truepi_smallest_neglogl->Fill(2.5);
      else if (Bragg_smallest == Bragg_K) All_truepi_smallest_neglogl->Fill(3.5);
      else if (Bragg_smallest == noBragg_fwd_MIP)
      All_truepi_smallest_neglogl->Fill(4.5);

      if ((PID_fwd && RecoTrackDir.Dot(TrueTrackDir) > 0) || (!PID_fwd && RecoTrackDir.Dot(TrueTrackDir) < 0)){ // correct track direction
        All_correctdirection->Fill(2.5,1.5);
      }
      else{ // right track direction
        All_correctdirection->Fill(2.5,0.5);
      }
    }

    else if (TMath::Abs(True_pdg) == 321){ // True kaons
      All_trueK_neglogl_mu->Fill(Bragg_mu);
      All_trueK_neglogl_p->Fill(Bragg_p);
      All_trueK_neglogl_pi->Fill(Bragg_pi);
      All_trueK_neglogl_K->Fill(Bragg_K);
      All_trueK_PIDA->Fill(PIDAval);
      All_trueK_dEdxtr_len->Fill(trklen,dEdxtruncmean);
      All_trueK_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
      All_trueK_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
      All_trueK_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
      All_trueK_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
      All_trueK_neglogl_MIP->Fill(noBragg_fwd_MIP);
      All_trueK_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
      All_trueK_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
      All_trueK_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
      All_trueK_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
      All_trueK_neglogl_p_vslength->Fill(Bragg_p,trklen);
      All_trueK_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
      All_trueK_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
      All_trueK_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
      All_trueK_neglogl_p_vsangle->Fill(Bragg_p,angle);
      All_trueK_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
      All_trueK_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
      All_trueK_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
      All_trueK_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

      if (Bragg_smallest == Bragg_mu) All_trueK_smallest_neglogl->Fill(0.5);
      else if (Bragg_smallest == Bragg_p) All_trueK_smallest_neglogl->Fill(1.5);
      else if (Bragg_smallest == Bragg_pi) All_trueK_smallest_neglogl->Fill(2.5);
      else if (Bragg_smallest == Bragg_K) All_trueK_smallest_neglogl->Fill(3.5);
      else if (Bragg_smallest == noBragg_fwd_MIP)
      All_trueK_smallest_neglogl->Fill(4.5);

      if ((PID_fwd && RecoTrackDir.Dot(TrueTrackDir) > 0) || (!PID_fwd && RecoTrackDir.Dot(TrueTrackDir) < 0)){ // correct track direction
        All_correctdirection->Fill(3.5,1.5);
      }
      else{ // right track direction
        All_correctdirection->Fill(3.5,0.5);
      }
    }

  } // Loop over tracks


}

DEFINE_ART_MODULE(ParticleIDValidationPlots)
