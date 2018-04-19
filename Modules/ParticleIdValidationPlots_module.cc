////////////////////////////////////////////////////////////////////////
// Class:       ParticleIdValidationPlots
// Plugin Type: analyzer (art v2_05_01)
// File:        ParticleIdValidationPlots_module.cc
//
// Generated at Thu Mar 15 14:40:38 2018 by Kirsty Duffy using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 *
 * \class ParticleIdValidationPlots
 *
 * \brief ParticleId analyzer module
 *
 * \author Kirsty Duffy (kduffy@fnal.gov), Adam Lister (alister1@lancaster.ac.uk)
 *
 * \date 2018/04/18
 *
 */

// Art
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

// LArSoft
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

// Algorithms
#include "uboone/ParticleID/Algorithms/FiducialVolume.h"

class ParticleIdValidationPlots;


class ParticleIdValidationPlots : public art::EDAnalyzer {
  public:
    explicit ParticleIdValidationPlots(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ParticleIdValidationPlots(ParticleIdValidationPlots const &) = delete;
    ParticleIdValidationPlots(ParticleIdValidationPlots &&) = delete;
    ParticleIdValidationPlots & operator = (ParticleIdValidationPlots const &) = delete;
    ParticleIdValidationPlots & operator = (ParticleIdValidationPlots &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;

  private:

    art::ServiceHandle<art::TFileService> tfs;
    fidvol::fiducialVolume fid;

    std::vector<double> fv;

    bool fIsData;
    std::string fTrackLabel;
    std::string fHitLabel;
    std::string fHitTrackAssns;
    std::string fCaloTrackAssns;
    std::string fHitTruthAssns;
    std::string fPIDLabel;

    /** Setup root tree  */
    TTree *pidTree;

    int true_PDG = -999;
    double true_start_momentum = -999;
    double true_start_x = -999; 
    double true_start_y = -999; 
    double true_start_z = -999; 
    double true_end_momentum = -999; 
    double true_end_x = -999; 
    double true_end_y = -999; 
    double true_end_z = -999; 
    double track_start_x;
    double track_start_y;
    double track_start_z;
    double track_end_x;      
    double track_end_y;
    double track_end_z;
    double track_neglogl_fwd_mu;
    double track_neglogl_fwd_p;
    double track_neglogl_fwd_pi;
    double track_neglogl_fwd_k;
    double track_neglogl_fwd_mip;
    double track_neglogl_bwd_mu;
    double track_neglogl_bwd_p;
    double track_neglogl_bwd_pi;
    double track_neglogl_bwd_k;
    // double track_neglogl_bwd_mip; // not used right now
    double track_PIDA_mean;
    double track_PIDA_median;
    double track_PIDA_kde;
    double track_length;
    double track_dEdx;
    double track_theta;
    double track_phi;
    double track_nhits;

    /** Histograms for all tracks, i.e. can be used by data */
    TH1F* AllTracks_PIDA;
    TH1F *AllTracks_neglogl_mu;
    TH1F *AllTracks_neglogl_p;
    TH1F *AllTracks_neglogl_pi;
    TH1F *AllTracks_neglogl_K;
    TH1F *AllTracks_neglogl_MIP;
    TH1F *AllTracks_neglogl_minmuMIP;
    TH2F *AllTracks_neglogl_mu_vslength;
    TH2F *AllTracks_neglogl_p_vslength;
    TH2F *AllTracks_neglogl_pi_vslength;
    TH2F *AllTracks_neglogl_K_vslength;
    TH2F *AllTracks_neglogl_MIP_vslength;
    TH2F *AllTracks_neglogl_minmuMIP_vslength;
    TH2F *AllTracks_neglogl_mu_vsangle;
    TH2F *AllTracks_neglogl_p_vsangle;
    TH2F *AllTracks_neglogl_pi_vsangle;
    TH2F *AllTracks_neglogl_K_vsangle;
    TH2F *AllTracks_neglogl_MIP_vsangle;
    TH2F *AllTracks_neglogl_minmuMIP_vsangle;
    TH2F *AllTracks_neglogl_mu_vsnhits;
    TH2F *AllTracks_neglogl_p_vsnhits;
    TH2F *AllTracks_neglogl_pi_vsnhits;
    TH2F *AllTracks_neglogl_K_vsnhits;
    TH2F *AllTracks_neglogl_MIP_vsnhits;
    TH2F *AllTracks_neglogl_minmuMIP_vsnhits;
    TH2F *AllTracks_neglogl_muvsp;
    TH2F *AllTracks_neglogl_muvspi;
    TH2F *AllTracks_neglogl_MIPvsp;
    TH2F *AllTracks_neglogl_minmuMIPvsp;
    TH1F *AllTracks_neglogl_muoverp;
    TH1F *AllTracks_neglogl_muminusp;
    TH1F *AllTracks_neglogl_MIPminusp;
    TH1F *AllTracks_neglogl_minmuMIPminusp;
    TH2F* AllTracks_dEdxtr_len;
    TH1F *AllTracks_smallest_neglogl;

    TH2F *All_chargeEndOverStart_sm0_5_dEdxrr;
    TH2F *All_chargeEndOverStart_gr2_dEdxrr;
    TH2F *All_chargeEndOverStart_0_5to2_dEdxrr;

    /** Histograms for tracks expected to have Bragg peaks */
    TH1F *TrueBragg_truemu_neglogl_mu;
    TH1F *TrueBragg_truep_neglogl_mu;
    TH1F *TrueBragg_truepi_neglogl_mu;
    TH1F *TrueBragg_trueK_neglogl_mu;
    TH1F *TrueBragg_truee_neglogl_mu;
    TH1F *TrueBragg_truemu_neglogl_p;
    TH1F *TrueBragg_truep_neglogl_p;
    TH1F *TrueBragg_truepi_neglogl_p;
    TH1F *TrueBragg_trueK_neglogl_p;
    TH1F *TrueBragg_truee_neglogl_p;
    TH1F *TrueBragg_truemu_neglogl_pi;
    TH1F *TrueBragg_truep_neglogl_pi;
    TH1F *TrueBragg_truepi_neglogl_pi;
    TH1F *TrueBragg_trueK_neglogl_pi;
    TH1F *TrueBragg_truee_neglogl_pi;
    TH1F *TrueBragg_truemu_neglogl_K;
    TH1F *TrueBragg_truep_neglogl_K;
    TH1F *TrueBragg_truepi_neglogl_K;
    TH1F *TrueBragg_trueK_neglogl_K;
    TH1F *TrueBragg_truee_neglogl_K;
    TH1F *TrueBragg_truemu_neglogl_MIP;
    TH1F *TrueBragg_truep_neglogl_MIP;
    TH1F *TrueBragg_truepi_neglogl_MIP;
    TH1F *TrueBragg_trueK_neglogl_MIP;
    TH1F *TrueBragg_truee_neglogl_MIP;
    TH1F *TrueBragg_truemu_neglogl_minmuMIP;
    TH1F *TrueBragg_truep_neglogl_minmuMIP;
    TH1F *TrueBragg_truepi_neglogl_minmuMIP;
    TH1F *TrueBragg_trueK_neglogl_minmuMIP;
    TH1F *TrueBragg_truee_neglogl_minmuMIP;

    TH2F *TrueBragg_truemu_neglogl_mu_vslength;
    TH2F *TrueBragg_truep_neglogl_mu_vslength;
    TH2F *TrueBragg_truepi_neglogl_mu_vslength;
    TH2F *TrueBragg_trueK_neglogl_mu_vslength;
    TH2F *TrueBragg_truee_neglogl_mu_vslength;
    TH2F *TrueBragg_truemu_neglogl_MIP_vslength;
    TH2F *TrueBragg_truep_neglogl_MIP_vslength;
    TH2F *TrueBragg_truepi_neglogl_MIP_vslength;
    TH2F *TrueBragg_trueK_neglogl_MIP_vslength;
    TH2F *TrueBragg_truee_neglogl_MIP_vslength;
    TH2F *TrueBragg_truemu_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_truep_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_truepi_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_trueK_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_truee_neglogl_minmuMIP_vslength;
    TH2F *TrueBragg_truemu_neglogl_p_vslength;
    TH2F *TrueBragg_truep_neglogl_p_vslength;
    TH2F *TrueBragg_truepi_neglogl_p_vslength;
    TH2F *TrueBragg_trueK_neglogl_p_vslength;
    TH2F *TrueBragg_truee_neglogl_p_vslength;
    TH2F *TrueBragg_truemu_neglogl_mu_vsangle;
    TH2F *TrueBragg_truep_neglogl_mu_vsangle;
    TH2F *TrueBragg_truepi_neglogl_mu_vsangle;
    TH2F *TrueBragg_trueK_neglogl_mu_vsangle;
    TH2F *TrueBragg_truee_neglogl_mu_vsangle;
    TH2F *TrueBragg_truemu_neglogl_MIP_vsangle;
    TH2F *TrueBragg_truep_neglogl_MIP_vsangle;
    TH2F *TrueBragg_truepi_neglogl_MIP_vsangle;
    TH2F *TrueBragg_trueK_neglogl_MIP_vsangle;
    TH2F *TrueBragg_truee_neglogl_MIP_vsangle;
    TH2F *TrueBragg_truemu_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_truep_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_truepi_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_trueK_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_truee_neglogl_minmuMIP_vsangle;
    TH2F *TrueBragg_truemu_neglogl_p_vsangle;
    TH2F *TrueBragg_truep_neglogl_p_vsangle;
    TH2F *TrueBragg_truepi_neglogl_p_vsangle;
    TH2F *TrueBragg_trueK_neglogl_p_vsangle;
    TH2F *TrueBragg_truee_neglogl_p_vsangle;
    TH2F *TrueBragg_truemu_neglogl_mu_vsnhits;
    TH2F *TrueBragg_truep_neglogl_mu_vsnhits;
    TH2F *TrueBragg_truepi_neglogl_mu_vsnhits;
    TH2F *TrueBragg_trueK_neglogl_mu_vsnhits;
    TH2F *TrueBragg_truee_neglogl_mu_vsnhits;
    TH2F *TrueBragg_truemu_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_truep_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_truepi_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_trueK_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_truee_neglogl_MIP_vsnhits;
    TH2F *TrueBragg_truemu_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_truep_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_truepi_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_trueK_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_truee_neglogl_minmuMIP_vsnhits;
    TH2F *TrueBragg_truemu_neglogl_p_vsnhits;
    TH2F *TrueBragg_truep_neglogl_p_vsnhits;
    TH2F *TrueBragg_truepi_neglogl_p_vsnhits;
    TH2F *TrueBragg_trueK_neglogl_p_vsnhits;
    TH2F *TrueBragg_truee_neglogl_p_vsnhits;

    TH2F *TrueBragg_truemu_neglogl_muvsp;
    TH2F *TrueBragg_truep_neglogl_muvsp;
    TH2F *TrueBragg_truepi_neglogl_muvsp;
    TH2F *TrueBragg_trueK_neglogl_muvsp;
    TH2F *TrueBragg_truee_neglogl_muvsp;
    TH2F *TrueBragg_truemu_neglogl_muvspi;
    TH2F *TrueBragg_truep_neglogl_muvspi;
    TH2F *TrueBragg_truepi_neglogl_muvspi;
    TH2F *TrueBragg_trueK_neglogl_muvspi;
    TH2F *TrueBragg_truee_neglogl_muvspi;

    TH2F *TrueBragg_truemu_neglogl_MIPvsp;
    TH2F *TrueBragg_truep_neglogl_MIPvsp;
    TH2F *TrueBragg_truepi_neglogl_MIPvsp;
    TH2F *TrueBragg_trueK_neglogl_MIPvsp;
    TH2F *TrueBragg_truee_neglogl_MIPvsp;
    TH2F *TrueBragg_truemu_neglogl_minmuMIPvsp;
    TH2F *TrueBragg_truep_neglogl_minmuMIPvsp;
    TH2F *TrueBragg_truepi_neglogl_minmuMIPvsp;
    TH2F *TrueBragg_trueK_neglogl_minmuMIPvsp;
    TH2F *TrueBragg_truee_neglogl_minmuMIPvsp;

    TH1F *TrueBragg_truemu_neglogl_muoverp;
    TH1F *TrueBragg_truep_neglogl_muoverp;
    TH1F *TrueBragg_truepi_neglogl_muoverp;
    TH1F *TrueBragg_trueK_neglogl_muoverp;
    TH1F *TrueBragg_truee_neglogl_muoverp;
    TH1F *TrueBragg_truemu_neglogl_muminusp;
    TH1F *TrueBragg_truep_neglogl_muminusp;
    TH1F *TrueBragg_truepi_neglogl_muminusp;
    TH1F *TrueBragg_trueK_neglogl_muminusp;
    TH1F *TrueBragg_truee_neglogl_muminusp;

    TH1F *TrueBragg_truemu_neglogl_MIPminusp;
    TH1F *TrueBragg_truep_neglogl_MIPminusp;
    TH1F *TrueBragg_truepi_neglogl_MIPminusp;
    TH1F *TrueBragg_trueK_neglogl_MIPminusp;
    TH1F *TrueBragg_truee_neglogl_MIPminusp;
    TH1F *TrueBragg_truemu_neglogl_minmuMIPminusp;
    TH1F *TrueBragg_truep_neglogl_minmuMIPminusp;
    TH1F *TrueBragg_truepi_neglogl_minmuMIPminusp;
    TH1F *TrueBragg_trueK_neglogl_minmuMIPminusp;
    TH1F *TrueBragg_truee_neglogl_minmuMIPminusp;

    TH1F *TrueBragg_truemu_smallest_neglogl;
    TH1F *TrueBragg_truep_smallest_neglogl;
    TH1F *TrueBragg_truepi_smallest_neglogl;
    TH1F *TrueBragg_trueK_smallest_neglogl;
    TH1F *TrueBragg_truee_smallest_neglogl;

    TH1F *TrueBragg_truemu_PIDA;
    TH1F *TrueBragg_truep_PIDA;
    TH1F *TrueBragg_truepi_PIDA;
    TH1F *TrueBragg_trueK_PIDA;
    TH1F *TrueBragg_truee_PIDA;

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
    TH2F *TrueBragg_truee_dEdxtr_len;

    TH2F *TrueBragg_correctdirection;
    TH2F *TrueBragg_incorrectdirection;
    TH1F *TrueBragg_PIDdir;

    TH2F *TrueBragg_correctdirection_PIDdir;
    TH2F *TrueBragg_incorrectdirection_PIDdir;

    /** All tracks which are matched to an MCParticle */
    TH1F *All_truemu_neglogl_mu;
    TH1F *All_truep_neglogl_mu;
    TH1F *All_truepi_neglogl_mu;
    TH1F *All_trueK_neglogl_mu;
    TH1F *All_truee_neglogl_mu;
    TH1F *All_truemu_neglogl_p;
    TH1F *All_truep_neglogl_p;
    TH1F *All_truepi_neglogl_p;
    TH1F *All_trueK_neglogl_p;
    TH1F *All_truee_neglogl_p;
    TH1F *All_truemu_neglogl_pi;
    TH1F *All_truep_neglogl_pi;
    TH1F *All_truepi_neglogl_pi;
    TH1F *All_trueK_neglogl_pi;
    TH1F *All_truee_neglogl_pi;
    TH1F *All_truemu_neglogl_K;
    TH1F *All_truep_neglogl_K;
    TH1F *All_truepi_neglogl_K;
    TH1F *All_trueK_neglogl_K;
    TH1F *All_truee_neglogl_K;
    TH1F *All_truemu_neglogl_MIP;
    TH1F *All_truep_neglogl_MIP;
    TH1F *All_truepi_neglogl_MIP;
    TH1F *All_trueK_neglogl_MIP;
    TH1F *All_truee_neglogl_MIP;
    TH1F *All_truemu_neglogl_minmuMIP;
    TH1F *All_truep_neglogl_minmuMIP;
    TH1F *All_truepi_neglogl_minmuMIP;
    TH1F *All_trueK_neglogl_minmuMIP;
    TH1F *All_truee_neglogl_minmuMIP;

    TH2F *All_truemu_neglogl_mu_vslength;
    TH2F *All_truep_neglogl_mu_vslength;
    TH2F *All_truepi_neglogl_mu_vslength;
    TH2F *All_trueK_neglogl_mu_vslength;
    TH2F *All_truee_neglogl_mu_vslength;
    TH2F *All_truemu_neglogl_MIP_vslength;
    TH2F *All_truep_neglogl_MIP_vslength;
    TH2F *All_truepi_neglogl_MIP_vslength;
    TH2F *All_trueK_neglogl_MIP_vslength;
    TH2F *All_truee_neglogl_MIP_vslength;
    TH2F *All_truemu_neglogl_minmuMIP_vslength;
    TH2F *All_truep_neglogl_minmuMIP_vslength;
    TH2F *All_truepi_neglogl_minmuMIP_vslength;
    TH2F *All_trueK_neglogl_minmuMIP_vslength;
    TH2F *All_truee_neglogl_minmuMIP_vslength;
    TH2F *All_truemu_neglogl_p_vslength;
    TH2F *All_truep_neglogl_p_vslength;
    TH2F *All_truepi_neglogl_p_vslength;
    TH2F *All_trueK_neglogl_p_vslength;
    TH2F *All_truee_neglogl_p_vslength;
    TH2F *All_truemu_neglogl_mu_vsangle;
    TH2F *All_truep_neglogl_mu_vsangle;
    TH2F *All_truepi_neglogl_mu_vsangle;
    TH2F *All_trueK_neglogl_mu_vsangle;
    TH2F *All_truee_neglogl_mu_vsangle;
    TH2F *All_truemu_neglogl_MIP_vsangle;
    TH2F *All_truep_neglogl_MIP_vsangle;
    TH2F *All_truepi_neglogl_MIP_vsangle;
    TH2F *All_trueK_neglogl_MIP_vsangle;
    TH2F *All_truee_neglogl_MIP_vsangle;
    TH2F *All_truemu_neglogl_minmuMIP_vsangle;
    TH2F *All_truep_neglogl_minmuMIP_vsangle;
    TH2F *All_truepi_neglogl_minmuMIP_vsangle;
    TH2F *All_trueK_neglogl_minmuMIP_vsangle;
    TH2F *All_truee_neglogl_minmuMIP_vsangle;
    TH2F *All_truemu_neglogl_p_vsangle;
    TH2F *All_truep_neglogl_p_vsangle;
    TH2F *All_truepi_neglogl_p_vsangle;
    TH2F *All_trueK_neglogl_p_vsangle;
    TH2F *All_truee_neglogl_p_vsangle;
    TH2F *All_truemu_neglogl_mu_vsnhits;
    TH2F *All_truep_neglogl_mu_vsnhits;
    TH2F *All_truepi_neglogl_mu_vsnhits;
    TH2F *All_trueK_neglogl_mu_vsnhits;
    TH2F *All_truee_neglogl_mu_vsnhits;
    TH2F *All_truemu_neglogl_MIP_vsnhits;
    TH2F *All_truep_neglogl_MIP_vsnhits;
    TH2F *All_truepi_neglogl_MIP_vsnhits;
    TH2F *All_trueK_neglogl_MIP_vsnhits;
    TH2F *All_truee_neglogl_MIP_vsnhits;
    TH2F *All_truemu_neglogl_minmuMIP_vsnhits;
    TH2F *All_truep_neglogl_minmuMIP_vsnhits;
    TH2F *All_truepi_neglogl_minmuMIP_vsnhits;
    TH2F *All_trueK_neglogl_minmuMIP_vsnhits;
    TH2F *All_truee_neglogl_minmuMIP_vsnhits;
    TH2F *All_truemu_neglogl_p_vsnhits;
    TH2F *All_truep_neglogl_p_vsnhits;
    TH2F *All_truepi_neglogl_p_vsnhits;
    TH2F *All_trueK_neglogl_p_vsnhits;
    TH2F *All_truee_neglogl_p_vsnhits;

    TH2F *All_truemu_neglogl_muvsp;
    TH2F *All_truep_neglogl_muvsp;
    TH2F *All_truepi_neglogl_muvsp;
    TH2F *All_trueK_neglogl_muvsp;
    TH2F *All_truee_neglogl_muvsp;
    TH2F *All_truemu_neglogl_muvspi;
    TH2F *All_truep_neglogl_muvspi;
    TH2F *All_truepi_neglogl_muvspi;
    TH2F *All_trueK_neglogl_muvspi;
    TH2F *All_truee_neglogl_muvspi;

    TH2F *All_truemu_neglogl_MIPvsp;
    TH2F *All_truep_neglogl_MIPvsp;
    TH2F *All_truepi_neglogl_MIPvsp;
    TH2F *All_trueK_neglogl_MIPvsp;
    TH2F *All_truee_neglogl_MIPvsp;
    TH2F *All_truemu_neglogl_minmuMIPvsp;
    TH2F *All_truep_neglogl_minmuMIPvsp;
    TH2F *All_truepi_neglogl_minmuMIPvsp;
    TH2F *All_trueK_neglogl_minmuMIPvsp;
    TH2F *All_truee_neglogl_minmuMIPvsp;

    TH1F *All_truemu_neglogl_muoverp;
    TH1F *All_truep_neglogl_muoverp;
    TH1F *All_truepi_neglogl_muoverp;
    TH1F *All_trueK_neglogl_muoverp;
    TH1F *All_truee_neglogl_muoverp;
    TH1F *All_truemu_neglogl_muminusp;
    TH1F *All_truep_neglogl_muminusp;
    TH1F *All_truepi_neglogl_muminusp;
    TH1F *All_trueK_neglogl_muminusp;
    TH1F *All_truee_neglogl_muminusp;

    TH1F *All_truemu_neglogl_MIPminusp;
    TH1F *All_truep_neglogl_MIPminusp;
    TH1F *All_truepi_neglogl_MIPminusp;
    TH1F *All_trueK_neglogl_MIPminusp;
    TH1F *All_truee_neglogl_MIPminusp;
    TH1F *All_truemu_neglogl_minmuMIPminusp;
    TH1F *All_truep_neglogl_minmuMIPminusp;
    TH1F *All_truepi_neglogl_minmuMIPminusp;
    TH1F *All_trueK_neglogl_minmuMIPminusp;
    TH1F *All_truee_neglogl_minmuMIPminusp;

    TH1F *All_truemu_smallest_neglogl;
    TH1F *All_truep_smallest_neglogl;
    TH1F *All_truepi_smallest_neglogl;
    TH1F *All_trueK_smallest_neglogl;
    TH1F *All_truee_smallest_neglogl;

    TH1F *All_truemu_PIDA;
    TH1F *All_truep_PIDA;
    TH1F *All_truepi_PIDA;
    TH1F *All_trueK_PIDA;
    TH1F *All_truee_PIDA;

    TH1F *All_chargeEndOverStart_directionCorrect;
    TH1F *All_chargeEndOverStart_directionIncorrect;
    TH2F *All_chargeEndOverStartVersusNHits_directionCorrect;
    TH2F *All_chargeEndOverStartVersusNHits_directionIncorrect;

    TH1F *Contained_chargeEndOverStart_directionCorrect;
    TH1F *Contained_chargeEndOverStart_directionIncorrect;
    TH2F *Contained_chargeEndOverStartVersusNHits_directionCorrect;
    TH2F *Contained_chargeEndOverStartVersusNHits_directionIncorrect;

    TH2F *All_truemu_dEdxtr_len;
    TH2F *All_truep_dEdxtr_len;
    TH2F *All_truepi_dEdxtr_len;
    TH2F *All_trueK_dEdxtr_len;
    TH2F *All_truee_dEdxtr_len;

    TH2F *All_correctdirection;
    TH2F *All_incorrectdirection;

};


ParticleIdValidationPlots::ParticleIdValidationPlots(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fhicl::ParameterSet const p_fv     = p.get<fhicl::ParameterSet>("FiducialVolume");
  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");

  fIsData = p.get<bool>("IsDataPlotsOnly", "false");
  fTrackLabel = p_labels.get<std::string>("TrackLabel","pandoraNu::McRecoStage2");
  fHitLabel = p_labels.get<std::string>("HitLabel","pandoraCosmicHitRemoval::McRecoStage2");
  fHitTrackAssns = p_labels.get<std::string>("HitTrackAssn","pandoraNu::McRecoStage2");
  fCaloTrackAssns = p_labels.get<std::string>("CaloTrackAssn", "pandoraNucali::McRecoStage2");
  fHitTruthAssns = p_labels.get<std::string>("HitTruthAssn","crHitRemovalTruthMatch::McRecoStage2");
  fPIDLabel = p_labels.get<std::string>("ParticleIdLabel");

  fv = fid.setFiducialVolume(fv, p_fv);
  fid.printFiducialVolume(fv);


  AllTracks_PIDA                   = tfs->make<TH1F>("AllTracks_PIDA"                   , ";PIDA;"                , 400                    , 0    , 50);
  AllTracks_neglogl_mu             = tfs->make<TH1F>("AllTracks_neglogl_mu"             , ";neg2LL_mu;"           , 400                    , 0    , 400);
  AllTracks_neglogl_p              = tfs->make<TH1F>("AllTracks_neglogl_p"              , ";neg2LL_p;"            , 400                    , 0    , 400);
  AllTracks_neglogl_pi             = tfs->make<TH1F>("AllTracks_neglogl_pi"             , ";neg2LL_pi;"           , 400                    , 0    , 400);
  AllTracks_neglogl_K              = tfs->make<TH1F>("AllTracks_neglogl_K"              , ";neg2LL_K;"            , 400                    , 0    , 400);
  AllTracks_neglogl_MIP            = tfs->make<TH1F>("AllTracks_neglogl_MIP"            , ";neg2LL_MIP;"          , 400                    , 0    , 400);
  AllTracks_neglogl_minmuMIP       = tfs->make<TH1F>("AllTracks_neglogl_minmuMIP"       , ";neg2LL_minmuMIP;"     , 400                    , 0    , 400);
  AllTracks_neglogl_muoverp        = tfs->make<TH1F>("AllTracks_neglogl_muoverp"        , ";neg2LL_mu/neg2LL_p;"  , 400                    , 0    , 20);
  AllTracks_neglogl_muminusp       = tfs->make<TH1F>("AllTracks_neglogl_muminusp"       , ";neg2LL_mu-neg2LL_p;"  , 400                    , -200 , 200);
  AllTracks_neglogl_MIPminusp      = tfs->make<TH1F>("AllTracks_neglogl_MIPminusp"      , ";neg2LL_MIP-neg2LL_p;" , 400                    , -200 , 200);
  AllTracks_neglogl_minmuMIPminusp = tfs->make<TH1F>("AllTracks_neglogl_minmuMIPminusp" , ";min(neg2LL_mu         , neg2LL_MIP)-neg2LL_p;" , 400  , -200  , 200);

  AllTracks_neglogl_mu_vslength       = tfs->make<TH2F>("AllTracks_neglogl_mu_vslength"       , ";neg2LL_mu;Track length (cm)"       , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_p_vslength        = tfs->make<TH2F>("AllTracks_neglogl_p_vslength"        , ";neg2LL_p;Track length (cm)"        , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_pi_vslength       = tfs->make<TH2F>("AllTracks_neglogl_pi_vslength"       , ";neg2LL_pi;Track length (cm)"       , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_K_vslength        = tfs->make<TH2F>("AllTracks_neglogl_K_vslength"        , ";neg2LL_K;Track length (cm)"        , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_MIP_vslength      = tfs->make<TH2F>("AllTracks_neglogl_MIP_vslength"      , ";neg2LL_MIP;Track length (cm)"      , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_minmuMIP_vslength = tfs->make<TH2F>("AllTracks_neglogl_minmuMIP_vslength" , ";neg2LL_minmuMIP;Track length (cm)" , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_mu_vsangle        = tfs->make<TH2F>("AllTracks_neglogl_mu_vsangle"        , ";neg2LL_mu;Track angle"             , 400                   , 0   , 400 , 150 , 0   , TMath::Pi());
  AllTracks_neglogl_p_vsangle         = tfs->make<TH2F>("AllTracks_neglogl_p_vsangle"         , ";neg2LL_p;Track angle"              , 400                   , 0   , 400 , 150 , 0   , TMath::Pi());
  AllTracks_neglogl_pi_vsangle        = tfs->make<TH2F>("AllTracks_neglogl_pi_vsangle"        , ";neg2LL_pi;Track angle"             , 400                   , 0   , 400 , 150 , 0   , TMath::Pi());
  AllTracks_neglogl_K_vsangle         = tfs->make<TH2F>("AllTracks_neglogl_K_vsangle"         , ";neg2LL_K;Track angle"              , 400                   , 0   , 400 , 150 , 0   , TMath::Pi());
  AllTracks_neglogl_MIP_vsangle       = tfs->make<TH2F>("AllTracks_neglogl_MIP_vsangle"       , ";neg2LL_MIP;Track angle"            , 400                   , 0   , 400 , 150 , 0   , TMath::Pi());
  AllTracks_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("AllTracks_neglogl_minmuMIP_vsangle"  , ";neg2LL_minmuMIP;Track angle"       , 400                   , 0   , 400 , 150 , 0   , TMath::Pi());
  AllTracks_neglogl_mu_vsnhits        = tfs->make<TH2F>("AllTracks_neglogl_mu_vsnhits"        , ";neg2LL_mu;N. hits"                 , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_p_vsnhits         = tfs->make<TH2F>("AllTracks_neglogl_p_vsnhits"         , ";neg2LL_p;N. hits"                  , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_pi_vsnhits        = tfs->make<TH2F>("AllTracks_neglogl_pi_vsnhits"        , ";neg2LL_pi;N. hits"                 , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_K_vsnhits         = tfs->make<TH2F>("AllTracks_neglogl_K_vsnhits"         , ";neg2LL_K;N. hits"                  , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_MIP_vsnhits       = tfs->make<TH2F>("AllTracks_neglogl_MIP_vsnhits"       , ";neg2LL_MIP;N. hits"                , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("AllTracks_neglogl_minmuMIP_vsnhits"  , ";neg2LL_minmuMIP;N. hits"           , 400                   , 0   , 400 , 500 , 0   , 500);
  AllTracks_neglogl_muvsp             = tfs->make<TH2F>("AllTracks_neglogl_muvsp"             , ";neg2LL_mu;neg2LL_p"                , 400                   , 0   , 400 , 400 , 0   , 400);
  AllTracks_neglogl_muvspi            = tfs->make<TH2F>("AllTracks_neglogl_muvspi"            , ";neg2LL_mu;neg2LL_pi"               , 400                   , 0   , 400 , 400 , 0   , 400);
  AllTracks_neglogl_MIPvsp            = tfs->make<TH2F>("AllTracks_neglogl_MIPvsp"            , ";neg2LL_MIP;neg2LL_p"               , 400                   , 0   , 400 , 400 , 0   , 400);
  AllTracks_neglogl_minmuMIPvsp       = tfs->make<TH2F>("AllTracks_neglogl_minmuMIPvsp"       , ";min(neg2LL_mu, neg2LL_MIP);neg2LL_p" , 400 , 0   , 400 , 400 , 0             , 400);
  AllTracks_dEdxtr_len                = tfs->make<TH2F>("AllTracks_dEdxtr_len"                , ";Track length (cm);dE/dx"           , 100                   , 0   , 700 , 100 , 0   , 50);

  const char* particles[5] = {"#mu", "p", "#pi", "K", "MIP"};
  AllTracks_smallest_neglogl  = tfs->make<TH1F>("AllTracks_smallest_neglogl","All tracks;Particle type with smallest neg2LL;",5,0,5);
  for (size_t i=1; i<=5; i++){
    AllTracks_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
  }

  All_chargeEndOverStart_sm0_5_dEdxrr = tfs->make<TH2F>("All_chargeEndOverStart_sm0_5_dEdxrr","All tracks (end/start average charge < 0.5);Residual range (cm); dEdx",150,0,30,400,0,50);
  All_chargeEndOverStart_gr2_dEdxrr = tfs->make<TH2F>("All_chargeEndOverStart_gr2_dEdxrr","All tracks (end/start average charge > 2);Residual range (cm); dEdx",150,0,30,400,0,50);
  All_chargeEndOverStart_0_5to2_dEdxrr = tfs->make<TH2F>("All_chargeEndOverStart_0_5to2_dEdxrr","All tracks (end/start average charge = 0.5 - 2);Residual range (cm); dEdx",150,0,30,400,0,50);

  if (!fIsData){
    // ---- True Bragg peak
    TrueBragg_truemu_neglogl_mu = tfs->make<TH1F>("TrueBragg_truemu_neglogl_mu","Tracks with true p=0 at end, true muons;neg2LL_mu;",400,0,400);
    TrueBragg_truep_neglogl_mu  = tfs->make<TH1F>("TrueBragg_truep_neglogl_mu","Tracks with true p=0 at end, true protons;neg2LL_mu;",400,0,400);
    TrueBragg_truepi_neglogl_mu = tfs->make<TH1F>("TrueBragg_truepi_neglogl_mu","Tracks with true p=0 at end, true pions;neg2LL_mu;",400,0,400);
    TrueBragg_trueK_neglogl_mu  = tfs->make<TH1F>("TrueBragg_trueK_neglogl_mu","Tracks with true p=0 at end, true kaons;neg2LL_mu;",400,0,400);
    TrueBragg_truee_neglogl_mu  = tfs->make<TH1F>("TrueBragg_truee_neglogl_mu","Tracks with true p=0 at end, true electrons;neg2LL_mu;",400,0,400);
    TrueBragg_truemu_neglogl_p  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_p","Tracks with true p=0 at end, true muons;neg2LL_p;",400,0,400);
    TrueBragg_truep_neglogl_p   = tfs->make<TH1F>("TrueBragg_truep_neglogl_p","Tracks with true p=0 at end, true protons;neg2LL_p;",400,0,400);
    TrueBragg_truepi_neglogl_p  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_p","Tracks with true p=0 at end, true pions;neg2LL_p;",400,0,400);
    TrueBragg_trueK_neglogl_p   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_p","Tracks with true p=0 at end, true kaons;neg2LL_p;",400,0,400);
    TrueBragg_truee_neglogl_p   = tfs->make<TH1F>("TrueBragg_truee_neglogl_p","Tracks with true p=0 at end, true electrons;neg2LL_p;",400,0,400);
    TrueBragg_truemu_neglogl_pi = tfs->make<TH1F>("TrueBragg_truemu_neglogl_pi","Tracks with true p=0 at end, true muons;neg2LL_pi;",400,0,400);
    TrueBragg_truep_neglogl_pi  = tfs->make<TH1F>("TrueBragg_truep_neglogl_pi","Tracks with true p=0 at end, true protons;neg2LL_pi;",400,0,400);
    TrueBragg_truepi_neglogl_pi = tfs->make<TH1F>("TrueBragg_truepi_neglogl_pi","Tracks with true p=0 at end, true pions;neg2LL_pi;",400,0,400);
    TrueBragg_trueK_neglogl_pi  = tfs->make<TH1F>("TrueBragg_trueK_neglogl_pi","Tracks with true p=0 at end, true kaons;neg2LL_pi;",400,0,400);
    TrueBragg_truee_neglogl_pi  = tfs->make<TH1F>("TrueBragg_truee_neglogl_pi","Tracks with true p=0 at end, true electrons;neg2LL_pi;",400,0,400);
    TrueBragg_truemu_neglogl_K  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_K","Tracks with true p=0 at end, true muons;neg2LL_K;",400,0,400);
    TrueBragg_truep_neglogl_K   = tfs->make<TH1F>("TrueBragg_truep_neglogl_K","Tracks with true p=0 at end, true protons;neg2LL_K;",400,0,400);
    TrueBragg_truepi_neglogl_K  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_K","Tracks with true p=0 at end, true pions;neg2LL_K;",400,0,400);
    TrueBragg_trueK_neglogl_K   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_K","Tracks with true p=0 at end, true kaons;neg2LL_K;",400,0,400);
    TrueBragg_truee_neglogl_K   = tfs->make<TH1F>("TrueBragg_truee_neglogl_K","Tracks with true p=0 at end, true electrons;neg2LL_K;",400,0,400);
    TrueBragg_truemu_neglogl_MIP  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_MIP","Tracks with true p=0 at end, true muons;neg2LL_MIP;",400,0,400);
    TrueBragg_truep_neglogl_MIP   = tfs->make<TH1F>("TrueBragg_truep_neglogl_MIP","Tracks with true p=0 at end, true protons;neg2LL_MIP;",400,0,400);
    TrueBragg_truepi_neglogl_MIP  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_MIP","Tracks with true p=0 at end, true pions;neg2LL_MIP;",400,0,400);
    TrueBragg_trueK_neglogl_MIP   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_MIP","Tracks with true p=0 at end, true kaons;neg2LL_MIP;",400,0,400);
    TrueBragg_truee_neglogl_MIP   = tfs->make<TH1F>("TrueBragg_truee_neglogl_MIP","Tracks with true p=0 at end, true electrons;neg2LL_MIP;",400,0,400);
    TrueBragg_truemu_neglogl_minmuMIP  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_minmuMIP","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);",400,0,400);
    TrueBragg_truep_neglogl_minmuMIP   = tfs->make<TH1F>("TrueBragg_truep_neglogl_minmuMIP","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);",400,0,400);
    TrueBragg_truepi_neglogl_minmuMIP  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_minmuMIP","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);",400,0,400);
    TrueBragg_trueK_neglogl_minmuMIP   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_minmuMIP","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);",400,0,400);
    TrueBragg_truee_neglogl_minmuMIP   = tfs->make<TH1F>("TrueBragg_truee_neglogl_minmuMIP","Tracks with true p=0 at end, true electrons;min(neg2LL_mu, neg2LL_MIP);",400,0,400);

    TrueBragg_truemu_neglogl_mu_vslength = tfs->make<TH2F>("TrueBragg_truemu_neglogl_mu_vslength","Tracks with true p=0 at end, true muons;neg2LL_mu;Track length",400,0,400,500,0,500);
    TrueBragg_truep_neglogl_mu_vslength  = tfs->make<TH2F>("TrueBragg_truep_neglogl_mu_vslength","Tracks with true p=0 at end, true protons;neg2LL_mu;Track length",400,0,400,500,0,500);
    TrueBragg_truepi_neglogl_mu_vslength = tfs->make<TH2F>("TrueBragg_truepi_neglogl_mu_vslength","Tracks with true p=0 at end, true pions;neg2LL_mu;Track length",400,0,400,500,0,500);
    TrueBragg_trueK_neglogl_mu_vslength  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_mu_vslength","Tracks with true p=0 at end, true kaons;neg2LL_mu;Track length",400,0,400,500,0,500);
    TrueBragg_truee_neglogl_mu_vslength  = tfs->make<TH2F>("TrueBragg_truee_neglogl_mu_vslength","Tracks with true p=0 at end, true electrons;neg2LL_mu;Track length",400,0,400,500,0,500);
    TrueBragg_truemu_neglogl_MIP_vslength = tfs->make<TH2F>("TrueBragg_truemu_neglogl_MIP_vslength","Tracks with true p=0 at end, true muons;neg2LL_MIP;Track length",400,0,400,500,0,500);
    TrueBragg_truep_neglogl_MIP_vslength  = tfs->make<TH2F>("TrueBragg_truep_neglogl_MIP_vslength","Tracks with true p=0 at end, true protons;neg2LL_MIP;Track length",400,0,400,500,0,500);
    TrueBragg_truepi_neglogl_MIP_vslength = tfs->make<TH2F>("TrueBragg_truepi_neglogl_MIP_vslength","Tracks with true p=0 at end, true pions;neg2LL_MIP;Track length",400,0,400,500,0,500);
    TrueBragg_trueK_neglogl_MIP_vslength  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_MIP_vslength","Tracks with true p=0 at end, true kaons;neg2LL_MIP;Track length",400,0,400,500,0,500);
    TrueBragg_truee_neglogl_MIP_vslength  = tfs->make<TH2F>("TrueBragg_truee_neglogl_MIP_vslength","Tracks with true p=0 at end, true electrons;neg2LL_MIP;Track length",400,0,400,500,0,500);
    TrueBragg_truemu_neglogl_minmuMIP_vslength = tfs->make<TH2F>("TrueBragg_truemu_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    TrueBragg_truep_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("TrueBragg_truep_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    TrueBragg_truepi_neglogl_minmuMIP_vslength = tfs->make<TH2F>("TrueBragg_truepi_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    TrueBragg_trueK_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    TrueBragg_truee_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("TrueBragg_truee_neglogl_minmuMIP_vslength","Tracks with true p=0 at end, true electrons;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    TrueBragg_truemu_neglogl_p_vslength = tfs->make<TH2F>("TrueBragg_truemu_neglogl_p_vslength","Tracks with true p=0 at end, true muons;neg2LL_p;Track length",400,0,400,500,0,500);
    TrueBragg_truep_neglogl_p_vslength  = tfs->make<TH2F>("TrueBragg_truep_neglogl_p_vslength","Tracks with true p=0 at end, true protons;neg2LL_p;Track length",400,0,400,500,0,500);
    TrueBragg_truepi_neglogl_p_vslength = tfs->make<TH2F>("TrueBragg_truepi_neglogl_p_vslength","Tracks with true p=0 at end, true pions;neg2LL_p;Track length",400,0,400,500,0,500);
    TrueBragg_trueK_neglogl_p_vslength  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_p_vslength","Tracks with true p=0 at end, true kaons;neg2LL_p;Track length",400,0,400,500,0,500);
    TrueBragg_truee_neglogl_p_vslength  = tfs->make<TH2F>("TrueBragg_truee_neglogl_p_vslength","Tracks with true p=0 at end, true electrons;neg2LL_p;Track length",400,0,400,500,0,500);
    TrueBragg_truemu_neglogl_mu_vsangle = tfs->make<TH2F>("TrueBragg_truemu_neglogl_mu_vsangle","Tracks with true p=0 at end, true muons;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truep_neglogl_mu_vsangle  = tfs->make<TH2F>("TrueBragg_truep_neglogl_mu_vsangle","Tracks with true p=0 at end, true protons;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truepi_neglogl_mu_vsangle = tfs->make<TH2F>("TrueBragg_truepi_neglogl_mu_vsangle","Tracks with true p=0 at end, true pions;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_trueK_neglogl_mu_vsangle  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_mu_vsangle","Tracks with true p=0 at end, true kaons;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truee_neglogl_mu_vsangle  = tfs->make<TH2F>("TrueBragg_truee_neglogl_mu_vsangle","Tracks with true p=0 at end, true electrons;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truemu_neglogl_MIP_vsangle = tfs->make<TH2F>("TrueBragg_truemu_neglogl_MIP_vsangle","Tracks with true p=0 at end, true muons;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truep_neglogl_MIP_vsangle  = tfs->make<TH2F>("TrueBragg_truep_neglogl_MIP_vsangle","Tracks with true p=0 at end, true protons;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truepi_neglogl_MIP_vsangle = tfs->make<TH2F>("TrueBragg_truepi_neglogl_MIP_vsangle","Tracks with true p=0 at end, true pions;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_trueK_neglogl_MIP_vsangle  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_MIP_vsangle","Tracks with true p=0 at end, true kaons;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truee_neglogl_MIP_vsangle  = tfs->make<TH2F>("TrueBragg_truee_neglogl_MIP_vsangle","Tracks with true p=0 at end, true electrons;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truemu_neglogl_minmuMIP_vsangle = tfs->make<TH2F>("TrueBragg_truemu_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truep_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("TrueBragg_truep_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truepi_neglogl_minmuMIP_vsangle = tfs->make<TH2F>("TrueBragg_truepi_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_trueK_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truee_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("TrueBragg_truee_neglogl_minmuMIP_vsangle","Tracks with true p=0 at end, true electrons;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truemu_neglogl_p_vsangle = tfs->make<TH2F>("TrueBragg_truemu_neglogl_p_vsangle","Tracks with true p=0 at end, true muons;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truep_neglogl_p_vsangle  = tfs->make<TH2F>("TrueBragg_truep_neglogl_p_vsangle","Tracks with true p=0 at end, true protons;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truepi_neglogl_p_vsangle = tfs->make<TH2F>("TrueBragg_truepi_neglogl_p_vsangle","Tracks with true p=0 at end, true pions;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_trueK_neglogl_p_vsangle  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_p_vsangle","Tracks with true p=0 at end, true kaons;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truee_neglogl_p_vsangle  = tfs->make<TH2F>("TrueBragg_truee_neglogl_p_vsangle","Tracks with true p=0 at end, true electrons;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    TrueBragg_truemu_neglogl_mu_vsnhits = tfs->make<TH2F>("TrueBragg_truemu_neglogl_mu_vsnhits","Tracks with true p=0 at end, true muons;neg2LL_mu;No. hits",400,0,400,500,0,500);
    TrueBragg_truep_neglogl_mu_vsnhits  = tfs->make<TH2F>("TrueBragg_truep_neglogl_mu_vsnhits","Tracks with true p=0 at end, true protons;neg2LL_mu;No. hits",400,0,400,500,0,500);
    TrueBragg_truepi_neglogl_mu_vsnhits = tfs->make<TH2F>("TrueBragg_truepi_neglogl_mu_vsnhits","Tracks with true p=0 at end, true pions;neg2LL_mu;No. hits",400,0,400,500,0,500);
    TrueBragg_trueK_neglogl_mu_vsnhits  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_mu_vsnhits","Tracks with true p=0 at end, true kaons;neg2LL_mu;No. hits",400,0,400,500,0,500);
    TrueBragg_truee_neglogl_mu_vsnhits  = tfs->make<TH2F>("TrueBragg_truee_neglogl_mu_vsnhits","Tracks with true p=0 at end, true electrons;neg2LL_mu;No. hits",400,0,400,500,0,500);
    TrueBragg_truemu_neglogl_MIP_vsnhits = tfs->make<TH2F>("TrueBragg_truemu_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true muons;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    TrueBragg_truep_neglogl_MIP_vsnhits  = tfs->make<TH2F>("TrueBragg_truep_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true protons;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    TrueBragg_truepi_neglogl_MIP_vsnhits = tfs->make<TH2F>("TrueBragg_truepi_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true pions;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    TrueBragg_trueK_neglogl_MIP_vsnhits  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true kaons;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    TrueBragg_truee_neglogl_MIP_vsnhits  = tfs->make<TH2F>("TrueBragg_truee_neglogl_MIP_vsnhits","Tracks with true p=0 at end, true electrons;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    TrueBragg_truemu_neglogl_minmuMIP_vsnhits = tfs->make<TH2F>("TrueBragg_truemu_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    TrueBragg_truep_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("TrueBragg_truep_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    TrueBragg_truepi_neglogl_minmuMIP_vsnhits = tfs->make<TH2F>("TrueBragg_truepi_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    TrueBragg_trueK_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    TrueBragg_truee_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("TrueBragg_truee_neglogl_minmuMIP_vsnhits","Tracks with true p=0 at end, true electrons;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    TrueBragg_truemu_neglogl_p_vsnhits = tfs->make<TH2F>("TrueBragg_truemu_neglogl_p_vsnhits","Tracks with true p=0 at end, true muons;neg2LL_p;No. hits",400,0,400,500,0,500);
    TrueBragg_truep_neglogl_p_vsnhits  = tfs->make<TH2F>("TrueBragg_truep_neglogl_p_vsnhits","Tracks with true p=0 at end, true protons;neg2LL_p;No. hits",400,0,400,500,0,500);
    TrueBragg_truepi_neglogl_p_vsnhits = tfs->make<TH2F>("TrueBragg_truepi_neglogl_p_vsnhits","Tracks with true p=0 at end, true pions;neg2LL_p;No. hits",400,0,400,500,0,500);
    TrueBragg_trueK_neglogl_p_vsnhits  = tfs->make<TH2F>("TrueBragg_trueK_neglogl_p_vsnhits","Tracks with true p=0 at end, true kaons;neg2LL_p;No. hits",400,0,400,500,0,500);
    TrueBragg_truee_neglogl_p_vsnhits  = tfs->make<TH2F>("TrueBragg_truee_neglogl_p_vsnhits","Tracks with true p=0 at end, true electrons;neg2LL_p;No. hits",400,0,400,500,0,500);

    TrueBragg_truemu_neglogl_muvsp  = tfs->make<TH2F>("TrueBragg_truemu_neglogl_muvsp","Tracks with true p=0 at end, true muons;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truep_neglogl_muvsp   = tfs->make<TH2F>("TrueBragg_truep_neglogl_muvsp","Tracks with true p=0 at end, true protons;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truepi_neglogl_muvsp  = tfs->make<TH2F>("TrueBragg_truepi_neglogl_muvsp","Tracks with true p=0 at end, true pions;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_trueK_neglogl_muvsp   = tfs->make<TH2F>("TrueBragg_trueK_neglogl_muvsp","Tracks with true p=0 at end, true kaons;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truee_neglogl_muvsp   = tfs->make<TH2F>("TrueBragg_truee_neglogl_muvsp","Tracks with true p=0 at end, true electrons;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truemu_neglogl_muvspi  = tfs->make<TH2F>("TrueBragg_truemu_neglogl_muvspi","Tracks with true p=0 at end, true muons;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);
    TrueBragg_truep_neglogl_muvspi   = tfs->make<TH2F>("TrueBragg_truep_neglogl_muvspi","Tracks with true p=0 at end, true protons;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);
    TrueBragg_truepi_neglogl_muvspi  = tfs->make<TH2F>("TrueBragg_truepi_neglogl_muvspi","Tracks with true p=0 at end, true pions;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);
    TrueBragg_trueK_neglogl_muvspi   = tfs->make<TH2F>("TrueBragg_trueK_neglogl_muvspi","Tracks with true p=0 at end, true kaons;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);
    TrueBragg_truee_neglogl_muvspi   = tfs->make<TH2F>("TrueBragg_truee_neglogl_muvspi","Tracks with true p=0 at end, true electrons;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);

    TrueBragg_truemu_neglogl_MIPvsp  = tfs->make<TH2F>("TrueBragg_truemu_neglogl_MIPvsp","Tracks with true p=0 at end, true muons;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truep_neglogl_MIPvsp   = tfs->make<TH2F>("TrueBragg_truep_neglogl_MIPvsp","Tracks with true p=0 at end, true protons;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truepi_neglogl_MIPvsp   = tfs->make<TH2F>("TrueBragg_truepi_neglogl_MIPvsp","Tracks with true p=0 at end, true pions;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_trueK_neglogl_MIPvsp   = tfs->make<TH2F>("TrueBragg_trueK_neglogl_MIPvsp","Tracks with true p=0 at end, true kaons;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truee_neglogl_MIPvsp   = tfs->make<TH2F>("TrueBragg_truee_neglogl_MIPvsp","Tracks with true p=0 at end, true electrons;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truemu_neglogl_minmuMIPvsp  = tfs->make<TH2F>("TrueBragg_truemu_neglogl_minmuMIPvsp","Tracks with true p=0 at end, true muons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truep_neglogl_minmuMIPvsp   = tfs->make<TH2F>("TrueBragg_truep_neglogl_minmuMIPvsp","Tracks with true p=0 at end, true protons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truepi_neglogl_minmuMIPvsp   = tfs->make<TH2F>("TrueBragg_truepi_neglogl_minmuMIPvsp","Tracks with true p=0 at end, true pions;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);
    TrueBragg_trueK_neglogl_minmuMIPvsp   = tfs->make<TH2F>("TrueBragg_trueK_neglogl_minmuMIPvsp","Tracks with true p=0 at end, true kaons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);
    TrueBragg_truee_neglogl_minmuMIPvsp   = tfs->make<TH2F>("TrueBragg_truee_neglogl_minmuMIPvsp","Tracks with true p=0 at end, true electrons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);

    TrueBragg_truemu_neglogl_muoverp  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_muoverp","Tracks with true p=0 at end, true muons;neg2LL_mu/neg2LL_p;",400,0,20);
    TrueBragg_truep_neglogl_muoverp   = tfs->make<TH1F>("TrueBragg_truep_neglogl_muoverp","Tracks with true p=0 at end, true protons;neg2LL_mu/neg2LL_p;",400,0,20);
    TrueBragg_truepi_neglogl_muoverp  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_muoverp","Tracks with true p=0 at end, true pions;neg2LL_mu/neg2LL_p;",400,0,20);
    TrueBragg_trueK_neglogl_muoverp   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_muoverp","Tracks with true p=0 at end, true kaons;neg2LL_mu/neg2LL_p;",400,0,20);
    TrueBragg_truee_neglogl_muoverp   = tfs->make<TH1F>("TrueBragg_truee_neglogl_muoverp","Tracks with true p=0 at end, true electrons;neg2LL_mu/neg2LL_p;",400,0,20);
    TrueBragg_truemu_neglogl_muminusp  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_muminusp","Tracks with true p=0 at end, true muons;neg2LL_mu-neg2LL_p;",400,-200,200);
    TrueBragg_truep_neglogl_muminusp   = tfs->make<TH1F>("TrueBragg_truep_neglogl_muminusp","Tracks with true p=0 at end, true protons;neg2LL_mu-neg2LL_p;",400,-200,200);
    TrueBragg_truepi_neglogl_muminusp  = tfs->make<TH1F>("TrueBragg_truepi_neglogl_muminusp","Tracks with true p=0 at end, true pions;neg2LL_mu-neg2LL_p;",400,-200,200);
    TrueBragg_trueK_neglogl_muminusp   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_muminusp","Tracks with true p=0 at end, true kaons;neg2LL_mu-neg2LL_p;",400,-200,200);
    TrueBragg_truee_neglogl_muminusp   = tfs->make<TH1F>("TrueBragg_truee_neglogl_muminusp","Tracks with true p=0 at end, true electrons;neg2LL_mu-neg2LL_p;",400,-200,200);

    TrueBragg_truemu_neglogl_MIPminusp  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_MIPminusp","Tracks with true p=0 at end, true muons;neg2LL_MIP-neg2LL_p;",400,-200,200);
    TrueBragg_truep_neglogl_MIPminusp   = tfs->make<TH1F>("TrueBragg_truep_neglogl_MIPminusp","Tracks with true p=0 at end, true protons;neg2LL_MIP-neg2LL_p;",400,-200,200);
    TrueBragg_truepi_neglogl_MIPminusp   = tfs->make<TH1F>("TrueBragg_truepi_neglogl_MIPminusp","Tracks with true p=0 at end, true pions;neg2LL_MIP-neg2LL_p;",400,-200,200);
    TrueBragg_trueK_neglogl_MIPminusp   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_MIPminusp","Tracks with true p=0 at end, true kaons;neg2LL_MIP-neg2LL_p;",400,-200,200);
    TrueBragg_truee_neglogl_MIPminusp   = tfs->make<TH1F>("TrueBragg_truee_neglogl_MIPminusp","Tracks with true p=0 at end, true electrons;neg2LL_MIP-neg2LL_p;",400,-200,200);
    TrueBragg_truemu_neglogl_minmuMIPminusp  = tfs->make<TH1F>("TrueBragg_truemu_neglogl_minmuMIPminusp","Tracks with true p=0 at end, true muons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",400,-200,200);
    TrueBragg_truep_neglogl_minmuMIPminusp   = tfs->make<TH1F>("TrueBragg_truep_neglogl_minmuMIPminusp","Tracks with true p=0 at end, true protons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",400,-200,200);
    TrueBragg_truepi_neglogl_minmuMIPminusp   = tfs->make<TH1F>("TrueBragg_truepi_neglogl_minmuMIPminusp","Tracks with true p=0 at end, true pions;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",400,-200,200);
    TrueBragg_trueK_neglogl_minmuMIPminusp   = tfs->make<TH1F>("TrueBragg_trueK_neglogl_minmuMIPminusp","Tracks with true p=0 at end, true kaons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",400,-200,200);
    TrueBragg_truee_neglogl_minmuMIPminusp   = tfs->make<TH1F>("TrueBragg_truee_neglogl_minmuMIPminusp","Tracks with true p=0 at end, true electrons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",400,-200,200);

    TrueBragg_truemu_smallest_neglogl  = tfs->make<TH1F>("TrueBragg_truemu_smallest_neglogl","Tracks with true p=0 at end, true muons;Particle type with smallest neg2LL;",5,0,5);
    TrueBragg_truep_smallest_neglogl   = tfs->make<TH1F>("TrueBragg_truep_smallest_neglogl","Tracks with true p=0 at end, true protons;Particle type with smallest neg2LL;",5,0,5);
    TrueBragg_truepi_smallest_neglogl  = tfs->make<TH1F>("TrueBragg_truepi_smallest_neglogl","Tracks with true p=0 at end, true pions;Particle type with smallest neg2LL;",5,0,5);
    TrueBragg_trueK_smallest_neglogl   = tfs->make<TH1F>("TrueBragg_trueK_smallest_neglogl","Tracks with true p=0 at end, true kaons;Particle type with smallest neg2LL;",5,0,5);
    TrueBragg_truee_smallest_neglogl   = tfs->make<TH1F>("TrueBragg_truee_smallest_neglogl","Tracks with true p=0 at end, true electrons;Particle type with smallest neg2LL;",5,0,5);
    for (size_t i=1; i<=5; i++){
      TrueBragg_truemu_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
      TrueBragg_truep_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
      TrueBragg_truepi_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
      TrueBragg_trueK_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
      TrueBragg_truee_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
    }

    TrueBragg_truemu_PIDA = tfs->make<TH1F>("TrueBragg_truemu_PIDA","Tracks with true p=0 at end, true muons;PIDA;",300,0,30);
    TrueBragg_truep_PIDA  = tfs->make<TH1F>("TrueBragg_truep_PIDA","Tracks with true p=0 at end, true protons;PIDA;",300,0,30);
    TrueBragg_truepi_PIDA = tfs->make<TH1F>("TrueBragg_truepi_PIDA","Tracks with true p=0 at end, true pions;PIDA;",300,0,30);
    TrueBragg_trueK_PIDA  = tfs->make<TH1F>("TrueBragg_trueK_PIDA","Tracks with true p=0 at end, true kaons;PIDA;",300,0,30);
    TrueBragg_truee_PIDA  = tfs->make<TH1F>("TrueBragg_truee_PIDA","Tracks with true p=0 at end, true electrons;PIDA;",300,0,30);

    TrueBragg_chargeEndOverStart_directionCorrect   = tfs->make<TH1F>("TrueBragg_chargeEndOverStart_directionCorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
    TrueBragg_chargeEndOverStart_directionIncorrect = tfs->make<TH1F>("TrueBragg_chargeEndOverStart_directionIncorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 100, 0, 10);
    ;
    TrueBragg_chargeEndOverStartVersusNHits_directionCorrect   = tfs->make<TH2F>("TrueBragg_chargeEndOverStartVersusNHits_directionCorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};number of hits used in average", 400, 0, 10, 6, 0, 6);
    TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect = tfs->make<TH2F>("TrueBragg_chargeEndOverStartVersusNHits_directionIncorrect", "Tracks with true p=0 at end;Charge_{End of track}/Charge_{Start of track};", 400, 0, 10, 6, 0, 6);
    ;

    TrueBragg_chargeEndOverStart_sm0_5_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_sm0_5_dEdxrr","Tracks with true p=0 at end (end/start average charge < 0.5);Residual range (cm); dEdx",150,0,30,400,0,50);
    TrueBragg_chargeEndOverStart_gr2_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_gr2_dEdxrr","Tracks with true p=0 at end (end/start average charge > 2);Residual range (cm); dEdx",150,0,30,400,0,50);
    TrueBragg_chargeEndOverStart_0_5to2_dEdxrr = tfs->make<TH2F>("TrueBragg_chargeEndOverStart_0_5to2_dEdxrr","Tracks with true p=0 at end (end/start average charge = 0.5 - 2);Residual range (cm); dEdx",150,0,30,400,0,50);

    TrueBragg_truemu_dEdxtr_len = tfs->make<TH2F>("TrueBragg_truemu_dEdxtr_len","Tracks with true p=0 at end, true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
    TrueBragg_truep_dEdxtr_len  = tfs->make<TH2F>("TrueBragg_truep_dEdxtr_len","Tracks with true p=0 at end, true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
    TrueBragg_truepi_dEdxtr_len = tfs->make<TH2F>("TrueBragg_truepi_dEdxtr_len","Tracks with true p=0 at end, true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
    TrueBragg_trueK_dEdxtr_len  = tfs->make<TH2F>("TrueBragg_trueK_dEdxtr_len","Tracks with true p=0 at end, true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);
    TrueBragg_truee_dEdxtr_len  = tfs->make<TH2F>("TrueBragg_truee_dEdxtr_len","Tracks with true p=0 at end, true electrons;Track length (cm);dE/dx",100,0,700,100,0,50);

    TrueBragg_correctdirection = tfs->make<TH2F>("TrueBragg_correctdirection","Tracks with true p=0 at end, reconstructed correct direction;true particle;PID direction correct?",4,0,4,2,0,2);
    TrueBragg_incorrectdirection = tfs->make<TH2F>("TrueBragg_incorrectdirection","Tracks with true p=0 at end, reconstructed incorrect direction;true particle;PID direction correct?",4,0,4,2,0,2);
    for (size_t i=1; i<=4; i++){
      TrueBragg_correctdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
      TrueBragg_incorrectdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
    }
    TrueBragg_correctdirection->GetYaxis()->SetBinLabel(1,"False");
    TrueBragg_correctdirection->GetYaxis()->SetBinLabel(2,"True");
    TrueBragg_incorrectdirection->GetYaxis()->SetBinLabel(1,"False");
    TrueBragg_incorrectdirection->GetYaxis()->SetBinLabel(2,"True");

    // ---- All tracks: truth matching using associations
    All_truemu_neglogl_mu = tfs->make<TH1F>("All_truemu_neglogl_mu","All tracks, true muons;neg2LL_mu;",400,0,400);
    All_truep_neglogl_mu  = tfs->make<TH1F>("All_truep_neglogl_mu","All tracks, true protons;neg2LL_mu;",400,0,400);
    All_truepi_neglogl_mu = tfs->make<TH1F>("All_truepi_neglogl_mu","All tracks, true pions;neg2LL_mu;",400,0,400);
    All_trueK_neglogl_mu  = tfs->make<TH1F>("All_trueK_neglogl_mu","All tracks, true kaons;neg2LL_mu;",400,0,400);
    All_truee_neglogl_mu  = tfs->make<TH1F>("All_truee_neglogl_mu","All tracks, true electrons;neg2LL_mu;",400,0,400);
    All_truemu_neglogl_p  = tfs->make<TH1F>("All_truemu_neglogl_p","All tracks, true muons;neg2LL_p;",400,0,400);
    All_truep_neglogl_p   = tfs->make<TH1F>("All_truep_neglogl_p","All tracks, true protons;neg2LL_p;",400,0,400);
    All_truepi_neglogl_p  = tfs->make<TH1F>("All_truepi_neglogl_p","All tracks, true pions;neg2LL_p;",400,0,400);
    All_trueK_neglogl_p   = tfs->make<TH1F>("All_trueK_neglogl_p","All tracks, true kaons;neg2LL_p;",400,0,400);
    All_truee_neglogl_p   = tfs->make<TH1F>("All_truee_neglogl_p","All tracks, true electrons;neg2LL_p;",400,0,400);
    All_truemu_neglogl_pi = tfs->make<TH1F>("All_truemu_neglogl_pi","All tracks, true muons;neg2LL_pi;",400,0,400);
    All_truep_neglogl_pi  = tfs->make<TH1F>("All_truep_neglogl_pi","All tracks, true protons;neg2LL_pi;",400,0,400);
    All_truepi_neglogl_pi = tfs->make<TH1F>("All_truepi_neglogl_pi","All tracks, true pions;neg2LL_pi;",400,0,400);
    All_trueK_neglogl_pi  = tfs->make<TH1F>("All_trueK_neglogl_pi","All tracks, true kaons;neg2LL_pi;",400,0,400);
    All_truee_neglogl_pi  = tfs->make<TH1F>("All_truee_neglogl_pi","All tracks, true electrons;neg2LL_pi;",400,0,400);
    All_truemu_neglogl_K  = tfs->make<TH1F>("All_truemu_neglogl_K","All tracks, true muons;neg2LL_K;",400,0,400);
    All_truep_neglogl_K   = tfs->make<TH1F>("All_truep_neglogl_K","All tracks, true protons;neg2LL_K;",400,0,400);
    All_truepi_neglogl_K  = tfs->make<TH1F>("All_truepi_neglogl_K","All tracks, true pions;neg2LL_K;",400,0,400);
    All_trueK_neglogl_K   = tfs->make<TH1F>("All_trueK_neglogl_K","All tracks, true kaons;neg2LL_K;",400,0,400);
    All_truee_neglogl_K   = tfs->make<TH1F>("All_truee_neglogl_K","All tracks, true electrons;neg2LL_K;",400,0,400);
    All_truemu_neglogl_MIP  = tfs->make<TH1F>("All_truemu_neglogl_MIP","All tracks, true muons;neg2LL_MIP;",400,0,400);
    All_truep_neglogl_MIP   = tfs->make<TH1F>("All_truep_neglogl_MIP","All tracks, true protons;neg2LL_MIP;",400,0,400);
    All_truepi_neglogl_MIP  = tfs->make<TH1F>("All_truepi_neglogl_MIP","All tracks, true pions;neg2LL_MIP;",400,0,400);
    All_trueK_neglogl_MIP   = tfs->make<TH1F>("All_trueK_neglogl_MIP","All tracks, true kaons;neg2LL_MIP;",400,0,400);
    All_truee_neglogl_MIP   = tfs->make<TH1F>("All_truee_neglogl_MIP","All tracks, true electrons;neg2LL_MIP;",400,0,400);
    All_truemu_neglogl_minmuMIP  = tfs->make<TH1F>("All_truemu_neglogl_minmuMIP","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);",400,0,400);
    All_truep_neglogl_minmuMIP   = tfs->make<TH1F>("All_truep_neglogl_minmuMIP","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);",400,0,400);
    All_truepi_neglogl_minmuMIP  = tfs->make<TH1F>("All_truepi_neglogl_minmuMIP","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);",400,0,400);
    All_trueK_neglogl_minmuMIP   = tfs->make<TH1F>("All_trueK_neglogl_minmuMIP","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);",400,0,400);
    All_truee_neglogl_minmuMIP   = tfs->make<TH1F>("All_truee_neglogl_minmuMIP","All tracks, true electrons;min(neg2LL_mu, neg2LL_MIP);",400,0,400);

    All_truemu_neglogl_mu_vslength = tfs->make<TH2F>("All_truemu_neglogl_mu_vslength","All tracks, true muons;neg2LL_mu;Track length",400,0,400,500,0,500);
    All_truep_neglogl_mu_vslength  = tfs->make<TH2F>("All_truep_neglogl_mu_vslength","All tracks, true protons;neg2LL_mu;Track length",400,0,400,500,0,500);
    All_truepi_neglogl_mu_vslength = tfs->make<TH2F>("All_truepi_neglogl_mu_vslength","All tracks, true pions;neg2LL_mu;Track length",400,0,400,500,0,500);
    All_trueK_neglogl_mu_vslength  = tfs->make<TH2F>("All_trueK_neglogl_mu_vslength","All tracks, true kaons;neg2LL_mu;Track length",400,0,400,500,0,500);
    All_truee_neglogl_mu_vslength  = tfs->make<TH2F>("All_truee_neglogl_mu_vslength","All tracks, true electrons;neg2LL_mu;Track length",400,0,400,500,0,500);
    All_truemu_neglogl_MIP_vslength = tfs->make<TH2F>("All_truemu_neglogl_MIP_vslength","All tracks, true muons;neg2LL_MIP;Track length",400,0,400,500,0,500);
    All_truep_neglogl_MIP_vslength  = tfs->make<TH2F>("All_truep_neglogl_MIP_vslength","All tracks, true protons;neg2LL_MIP;Track length",400,0,400,500,0,500);
    All_truepi_neglogl_MIP_vslength = tfs->make<TH2F>("All_truepi_neglogl_MIP_vslength","All tracks, true pions;neg2LL_MIP;Track length",400,0,400,500,0,500);
    All_trueK_neglogl_MIP_vslength  = tfs->make<TH2F>("All_trueK_neglogl_MIP_vslength","All tracks, true kaons;neg2LL_MIP;Track length",400,0,400,500,0,500);
    All_truee_neglogl_MIP_vslength  = tfs->make<TH2F>("All_truee_neglogl_MIP_vslength","All tracks, true electrons;neg2LL_MIP;Track length",400,0,400,500,0,500);
    All_truemu_neglogl_minmuMIP_vslength = tfs->make<TH2F>("All_truemu_neglogl_minmuMIP_vslength","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    All_truep_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("All_truep_neglogl_minmuMIP_vslength","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    All_truepi_neglogl_minmuMIP_vslength = tfs->make<TH2F>("All_truepi_neglogl_minmuMIP_vslength","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    All_trueK_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("All_trueK_neglogl_minmuMIP_vslength","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    All_truee_neglogl_minmuMIP_vslength  = tfs->make<TH2F>("All_truee_neglogl_minmuMIP_vslength","All tracks, true electrons;min(neg2LL_mu, neg2LL_MIP);Track length",400,0,400,500,0,500);
    All_truemu_neglogl_p_vslength = tfs->make<TH2F>("All_truemu_neglogl_p_vslength","All tracks, true muons;neg2LL_p;Track length",400,0,400,500,0,500);
    All_truep_neglogl_p_vslength  = tfs->make<TH2F>("All_truep_neglogl_p_vslength","All tracks, true protons;neg2LL_p;Track length",400,0,400,500,0,500);
    All_truepi_neglogl_p_vslength = tfs->make<TH2F>("All_truepi_neglogl_p_vslength","All tracks, true pions;neg2LL_p;Track length",400,0,400,500,0,500);
    All_trueK_neglogl_p_vslength  = tfs->make<TH2F>("All_trueK_neglogl_p_vslength","All tracks, true kaons;neg2LL_p;Track length",400,0,400,500,0,500);
    All_truee_neglogl_p_vslength  = tfs->make<TH2F>("All_truee_neglogl_p_vslength","All tracks, true electrons;neg2LL_p;Track length",400,0,400,500,0,500);
    All_truemu_neglogl_mu_vsangle = tfs->make<TH2F>("All_truemu_neglogl_mu_vsangle","All tracks, true muons;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    All_truep_neglogl_mu_vsangle  = tfs->make<TH2F>("All_truep_neglogl_mu_vsangle","All tracks, true protons;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    All_truepi_neglogl_mu_vsangle = tfs->make<TH2F>("All_truepi_neglogl_mu_vsangle","All tracks, true pions;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    All_trueK_neglogl_mu_vsangle  = tfs->make<TH2F>("All_trueK_neglogl_mu_vsangle","All tracks, true kaons;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    All_truee_neglogl_mu_vsangle  = tfs->make<TH2F>("All_truee_neglogl_mu_vsangle","All tracks, true electrons;neg2LL_mu;Track angle",400,0,400,150,0,TMath::Pi());
    All_truemu_neglogl_MIP_vsangle = tfs->make<TH2F>("All_truemu_neglogl_MIP_vsangle","All tracks, true muons;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    All_truep_neglogl_MIP_vsangle  = tfs->make<TH2F>("All_truep_neglogl_MIP_vsangle","All tracks, true protons;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    All_truepi_neglogl_MIP_vsangle = tfs->make<TH2F>("All_truepi_neglogl_MIP_vsangle","All tracks, true pions;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    All_trueK_neglogl_MIP_vsangle  = tfs->make<TH2F>("All_trueK_neglogl_MIP_vsangle","All tracks, true kaons;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    All_truee_neglogl_MIP_vsangle  = tfs->make<TH2F>("All_truee_neglogl_MIP_vsangle","All tracks, true electrons;neg2LL_MIP;Track angle",400,0,400,150,0,TMath::Pi());
    All_truemu_neglogl_minmuMIP_vsangle = tfs->make<TH2F>("All_truemu_neglogl_minmuMIP_vsangle","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    All_truep_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("All_truep_neglogl_minmuMIP_vsangle","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    All_truepi_neglogl_minmuMIP_vsangle = tfs->make<TH2F>("All_truepi_neglogl_minmuMIP_vsangle","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    All_trueK_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("All_trueK_neglogl_minmuMIP_vsangle","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    All_truee_neglogl_minmuMIP_vsangle  = tfs->make<TH2F>("All_truee_neglogl_minmuMIP_vsangle","All tracks, true electrons;min(neg2LL_mu, neg2LL_MIP);Track angle",400,0,400,150,0,TMath::Pi());
    All_truemu_neglogl_p_vsangle = tfs->make<TH2F>("All_truemu_neglogl_p_vsangle","All tracks, true muons;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    All_truep_neglogl_p_vsangle  = tfs->make<TH2F>("All_truep_neglogl_p_vsangle","All tracks, true protons;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    All_truepi_neglogl_p_vsangle = tfs->make<TH2F>("All_truepi_neglogl_p_vsangle","All tracks, true pions;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    All_trueK_neglogl_p_vsangle  = tfs->make<TH2F>("All_trueK_neglogl_p_vsangle","All tracks, true kaons;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    All_truee_neglogl_p_vsangle  = tfs->make<TH2F>("All_truee_neglogl_p_vsangle","All tracks, true electrons;neg2LL_p;Track angle",400,0,400,150,0,TMath::Pi());
    All_truemu_neglogl_mu_vsnhits = tfs->make<TH2F>("All_truemu_neglogl_mu_vsnhits","All tracks, true muons;neg2LL_mu;No. hits",400,0,400,500,0,500);
    All_truep_neglogl_mu_vsnhits  = tfs->make<TH2F>("All_truep_neglogl_mu_vsnhits","All tracks, true protons;neg2LL_mu;No. hits",400,0,400,500,0,500);
    All_truepi_neglogl_mu_vsnhits = tfs->make<TH2F>("All_truepi_neglogl_mu_vsnhits","All tracks, true pions;neg2LL_mu;No. hits",400,0,400,500,0,500);
    All_trueK_neglogl_mu_vsnhits  = tfs->make<TH2F>("All_trueK_neglogl_mu_vsnhits","All tracks, true kaons;neg2LL_mu;No. hits",400,0,400,500,0,500);
    All_truee_neglogl_mu_vsnhits  = tfs->make<TH2F>("All_truee_neglogl_mu_vsnhits","All tracks, true electrons;neg2LL_mu;No. hits",400,0,400,500,0,500);
    All_truemu_neglogl_MIP_vsnhits = tfs->make<TH2F>("All_truemu_neglogl_MIP_vsnhits","All tracks, true muons;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    All_truep_neglogl_MIP_vsnhits  = tfs->make<TH2F>("All_truep_neglogl_MIP_vsnhits","All tracks, true protons;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    All_truepi_neglogl_MIP_vsnhits = tfs->make<TH2F>("All_truepi_neglogl_MIP_vsnhits","All tracks, true pions;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    All_trueK_neglogl_MIP_vsnhits  = tfs->make<TH2F>("All_trueK_neglogl_MIP_vsnhits","All tracks, true kaons;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    All_truee_neglogl_MIP_vsnhits  = tfs->make<TH2F>("All_truee_neglogl_MIP_vsnhits","All tracks, true electrons;neg2LL_MIP;No. hits",400,0,400,500,0,500);
    All_truemu_neglogl_minmuMIP_vsnhits = tfs->make<TH2F>("All_truemu_neglogl_minmuMIP_vsnhits","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    All_truep_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("All_truep_neglogl_minmuMIP_vsnhits","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    All_truepi_neglogl_minmuMIP_vsnhits = tfs->make<TH2F>("All_truepi_neglogl_minmuMIP_vsnhits","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    All_trueK_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("All_trueK_neglogl_minmuMIP_vsnhits","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    All_truee_neglogl_minmuMIP_vsnhits  = tfs->make<TH2F>("All_truee_neglogl_minmuMIP_vsnhits","All tracks, true electrons;min(neg2LL_mu, neg2LL_MIP);No. hits",400,0,400,500,0,500);
    All_truemu_neglogl_p_vsnhits = tfs->make<TH2F>("All_truemu_neglogl_p_vsnhits","All tracks, true muons;neg2LL_p;No. hits",400,0,400,500,0,500);
    All_truep_neglogl_p_vsnhits  = tfs->make<TH2F>("All_truep_neglogl_p_vsnhits","All tracks, true protons;neg2LL_p;No. hits",400,0,400,500,0,500);
    All_truepi_neglogl_p_vsnhits = tfs->make<TH2F>("All_truepi_neglogl_p_vsnhits","All tracks, true pions;neg2LL_p;No. hits",400,0,400,500,0,500);
    All_trueK_neglogl_p_vsnhits  = tfs->make<TH2F>("All_trueK_neglogl_p_vsnhits","All tracks, true kaons;neg2LL_p;No. hits",400,0,400,500,0,500);
    All_truee_neglogl_p_vsnhits  = tfs->make<TH2F>("All_truee_neglogl_p_vsnhits","All tracks, true electrons;neg2LL_p;No. hits",400,0,400,500,0,500);

    All_truemu_neglogl_muvsp  = tfs->make<TH2F>("All_truemu_neglogl_muvsp","All tracks, true muons;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    All_truep_neglogl_muvsp   = tfs->make<TH2F>("All_truep_neglogl_muvsp","All tracks, true protons;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    All_truepi_neglogl_muvsp  = tfs->make<TH2F>("All_truepi_neglogl_muvsp","All tracks, true pions;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    All_trueK_neglogl_muvsp   = tfs->make<TH2F>("All_trueK_neglogl_muvsp","All tracks, true kaons;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    All_truee_neglogl_muvsp   = tfs->make<TH2F>("All_truee_neglogl_muvsp","All tracks, true electrons;neg2LL_mu;neg2LL_p",400,0,400,400,0,400);
    All_truemu_neglogl_muvspi  = tfs->make<TH2F>("All_truemu_neglogl_muvspi","All tracks, true muons;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);
    All_truep_neglogl_muvspi   = tfs->make<TH2F>("All_truep_neglogl_muvspi","All tracks, true protons;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);
    All_truepi_neglogl_muvspi  = tfs->make<TH2F>("All_truepi_neglogl_muvspi","All tracks, true pions;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);
    All_trueK_neglogl_muvspi   = tfs->make<TH2F>("All_trueK_neglogl_muvspi","All tracks, true kaons;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);
    All_truee_neglogl_muvspi   = tfs->make<TH2F>("All_truee_neglogl_muvspi","All tracks, true electrons;neg2LL_mu;neg2LL_pi",400,0,400,400,0,400);

    All_truemu_neglogl_MIPvsp  = tfs->make<TH2F>("All_truemu_neglogl_MIPvsp","All tracks, true muons;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    All_truep_neglogl_MIPvsp   = tfs->make<TH2F>("All_truep_neglogl_MIPvsp","All tracks, true protons;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    All_truepi_neglogl_MIPvsp   = tfs->make<TH2F>("All_truepi_neglogl_MIPvsp","All tracks, true pions;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    All_trueK_neglogl_MIPvsp   = tfs->make<TH2F>("All_trueK_neglogl_MIPvsp","All tracks, true kaons;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    All_truee_neglogl_MIPvsp   = tfs->make<TH2F>("All_truee_neglogl_MIPvsp","All tracks, true electrons;neg2LL_MIP;neg2LL_p",400,0,400,400,0,400);
    All_truemu_neglogl_minmuMIPvsp  = tfs->make<TH2F>("All_truemu_neglogl_minmuMIPvsp","All tracks, true muons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);
    All_truep_neglogl_minmuMIPvsp   = tfs->make<TH2F>("All_truep_neglogl_minmuMIPvsp","All tracks, true protons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);
    All_truepi_neglogl_minmuMIPvsp   = tfs->make<TH2F>("All_truepi_neglogl_minmuMIPvsp","All tracks, true pions;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);
    All_trueK_neglogl_minmuMIPvsp   = tfs->make<TH2F>("All_trueK_neglogl_minmuMIPvsp","All tracks, true kaons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);
    All_truee_neglogl_minmuMIPvsp   = tfs->make<TH2F>("All_truee_neglogl_minmuMIPvsp","All tracks, true electrons;min(neg2LL_mu, neg2LL_MIP);neg2LL_p",400,0,400,400,0,400);

    All_truemu_neglogl_muoverp  = tfs->make<TH1F>("All_truemu_neglogl_muoverp","All tracks, true muons;neg2LL_mu/neg2LL_p;",200,0,50);
    All_truep_neglogl_muoverp   = tfs->make<TH1F>("All_truep_neglogl_muoverp","All tracks, true protons;neg2LL_mu/neg2LL_p;",200,0,50);
    All_truepi_neglogl_muoverp  = tfs->make<TH1F>("All_truepi_neglogl_muoverp","All tracks, true pions;neg2LL_mu/neg2LL_p;",200,0,50);
    All_trueK_neglogl_muoverp   = tfs->make<TH1F>("All_trueK_neglogl_muoverp","All tracks, true kaons;neg2LL_mu/neg2LL_p;",200,0,50);
    All_truee_neglogl_muoverp   = tfs->make<TH1F>("All_truee_neglogl_muoverp","All tracks, true electrons;neg2LL_mu/neg2LL_p;",200,0,50);
    All_truemu_neglogl_muminusp  = tfs->make<TH1F>("All_truemu_neglogl_muminusp","All tracks, true muons;neg2LL_mu-neg2LL_p;",200,-200,200);
    All_truep_neglogl_muminusp   = tfs->make<TH1F>("All_truep_neglogl_muminusp","All tracks, true protons;neg2LL_mu-neg2LL_p;",200,-200,200);
    All_truepi_neglogl_muminusp  = tfs->make<TH1F>("All_truepi_neglogl_muminusp","All tracks, true pions;neg2LL_mu-neg2LL_p;",200,-200,200);
    All_trueK_neglogl_muminusp   = tfs->make<TH1F>("All_trueK_neglogl_muminusp","All tracks, true kaons;neg2LL_K;",200,-200,200);
    All_truee_neglogl_muminusp   = tfs->make<TH1F>("All_truee_neglogl_muminusp","All tracks, true electrons;neg2LL_K;",200,-200,200);

    All_truemu_neglogl_MIPminusp  = tfs->make<TH1F>("All_truemu_neglogl_MIPminusp","All tracks, true muons;neg2LL_MIP-neg2LL_p;",200,-200,200);
    All_truep_neglogl_MIPminusp   = tfs->make<TH1F>("All_truep_neglogl_MIPminusp","All tracks, true protons;neg2LL_MIP-neg2LL_p;",200,-200,200);
    All_truepi_neglogl_MIPminusp   = tfs->make<TH1F>("All_truepi_neglogl_MIPminusp","All tracks, true pions;neg2LL_MIP-neg2LL_p;",200,-200,200);
    All_trueK_neglogl_MIPminusp   = tfs->make<TH1F>("All_trueK_neglogl_MIPminusp","All tracks, true kaons;neg2LL_MIP-neg2LL_p;",200,-200,200);
    All_truee_neglogl_MIPminusp   = tfs->make<TH1F>("All_truee_neglogl_MIPminusp","All tracks, true electrons;neg2LL_MIP-neg2LL_p;",200,-200,200);
    All_truemu_neglogl_minmuMIPminusp  = tfs->make<TH1F>("All_truemu_neglogl_minmuMIPminusp","All tracks, true muons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-200,200);
    All_truep_neglogl_minmuMIPminusp   = tfs->make<TH1F>("All_truep_neglogl_minmuMIPminusp","All tracks, true protons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-200,200);
    All_truepi_neglogl_minmuMIPminusp   = tfs->make<TH1F>("All_truepi_neglogl_minmuMIPminusp","All tracks, true pions;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-200,200);
    All_trueK_neglogl_minmuMIPminusp   = tfs->make<TH1F>("All_trueK_neglogl_minmuMIPminusp","All tracks, true kaons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-200,200);
    All_truee_neglogl_minmuMIPminusp   = tfs->make<TH1F>("All_truee_neglogl_minmuMIPminusp","All tracks, true electrons;min(neg2LL_MIP, neg2LL_mu)-neg2LL_p;",200,-200,200);

    All_truemu_smallest_neglogl  = tfs->make<TH1F>("All_truemu_smallest_neglogl","All tracks, true muons;Particle type with smallest neg2LL;",5,0,5);
    All_truep_smallest_neglogl   = tfs->make<TH1F>("All_truep_smallest_neglogl","All tracks, true protons;Particle type with smallest neg2LL;",5,0,5);
    All_truepi_smallest_neglogl  = tfs->make<TH1F>("All_truepi_smallest_neglogl","All tracks, true pions;Particle type with smallest neg2LL;",5,0,5);
    All_trueK_smallest_neglogl   = tfs->make<TH1F>("All_trueK_smallest_neglogl","All tracks, true kaons;Particle type with smallest neg2LL;",5,0,5);
    All_truee_smallest_neglogl   = tfs->make<TH1F>("All_truee_smallest_neglogl","All tracks, true electrons;Particle type with smallest neg2LL;",5,0,5);

    for (size_t i=1; i<=5; i++){
      All_truemu_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
      All_truep_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
      All_truepi_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
      All_trueK_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
      All_truee_smallest_neglogl->GetXaxis()->SetBinLabel(i,particles[i-1]);
    }

    All_truemu_PIDA = tfs->make<TH1F>("All_truemu_PIDA","All tracks, true muons;PIDA;",300,0,30);
    All_truep_PIDA  = tfs->make<TH1F>("All_truep_PIDA","All tracks, true protons;PIDA;",300,0,30);
    All_truepi_PIDA = tfs->make<TH1F>("All_truepi_PIDA","All tracks, true pions;PIDA;",300,0,30);
    All_trueK_PIDA  = tfs->make<TH1F>("All_trueK_PIDA","All tracks, true kaons;PIDA;",300,0,30);
    All_truee_PIDA  = tfs->make<TH1F>("All_truee_PIDA","All tracks, true electrons;PIDA;",300,0,30);

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

    All_truemu_dEdxtr_len = tfs->make<TH2F>("All_truemu_dEdxtr_len","All tracks, true muons;Track length (cm);dE/dx",100,0,700,100,0,50);
    All_truep_dEdxtr_len  = tfs->make<TH2F>("All_truep_dEdxtr_len","All tracks, true protons;Track length (cm);dE/dx",100,0,700,100,0,50);
    All_truepi_dEdxtr_len = tfs->make<TH2F>("All_truepi_dEdxtr_len","All tracks, true pions;Track length (cm);dE/dx",100,0,700,100,0,50);
    All_trueK_dEdxtr_len  = tfs->make<TH2F>("All_trueK_dEdxtr_len","All tracks, true kaons;Track length (cm);dE/dx",100,0,700,100,0,50);
    All_truee_dEdxtr_len  = tfs->make<TH2F>("All_truee_dEdxtr_len","All tracks, true electrons;Track length (cm);dE/dx",100,0,700,100,0,50);

    All_correctdirection = tfs->make<TH2F>("All_correctdirection","All tracks, reconstructed correct direction;true particle;PID direction correct?",4,0,4,2,0,2);
    All_incorrectdirection = tfs->make<TH2F>("All_incorrectdirection","All tracks, reconstructed incorrect direction;true particle;PID direction correct?",4,0,4,2,0,2);
    for (size_t i=1; i<=4; i++){
      All_correctdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
      All_incorrectdirection->GetXaxis()->SetBinLabel(i,particles[i-1]);
    }
    All_correctdirection->GetYaxis()->SetBinLabel(1,"False");
    All_correctdirection->GetYaxis()->SetBinLabel(2,"True");
    All_incorrectdirection->GetYaxis()->SetBinLabel(1,"False");
    All_incorrectdirection->GetYaxis()->SetBinLabel(2,"True");
  }
}

void ParticleIdValidationPlots::beginJob(){

  pidTree = tfs->make<TTree>("pidTree" , "pidTree");

  pidTree->Branch( "true_PDG"            , &true_PDG            ) ;
  pidTree->Branch( "true_start_momentum" , &true_start_momentum ) ;
  pidTree->Branch( "true_start_x"        , &true_start_x        ) ;
  pidTree->Branch( "true_start_y"        , &true_start_y        ) ;
  pidTree->Branch( "true_start_z"        , &true_start_z        ) ;
  pidTree->Branch( "true_end_momentum"   , &true_end_momentum   ) ;
  pidTree->Branch( "true_end_x"          , &true_end_x          ) ;
  pidTree->Branch( "true_end_y"          , &true_end_y          ) ;
  pidTree->Branch( "true_end_z"          , &true_end_z          ) ;

  pidTree->Branch( "track_start_x"           , &track_start_x       ) ;
  pidTree->Branch( "track_start_y"           , &track_start_y       ) ;
  pidTree->Branch( "track_start_z"           , &track_start_z       ) ;
  pidTree->Branch( "track_end_x"             , &track_end_x         ) ;
  pidTree->Branch( "track_end_y"             , &track_end_y         ) ;
  pidTree->Branch( "track_end_z"             , &track_end_z         ) ;
  pidTree->Branch( "track_neglogl_mu"        , &track_neglogl_fwd_mu    ) ;
  pidTree->Branch( "track_neglogl_p"         , &track_neglogl_fwd_p     ) ;
  pidTree->Branch( "track_neglogl_pi"        , &track_neglogl_fwd_pi    ) ;
  pidTree->Branch( "track_neglogl_k"         , &track_neglogl_fwd_k     ) ;
  pidTree->Branch( "track_neglogl_mip"       , &track_neglogl_fwd_mip   ) ;
  pidTree->Branch( "track_neglogl_mu"        , &track_neglogl_bwd_mu    ) ;
  pidTree->Branch( "track_neglogl_p"         , &track_neglogl_bwd_p     ) ;
  pidTree->Branch( "track_neglogl_pi"        , &track_neglogl_bwd_pi    ) ;
  pidTree->Branch( "track_neglogl_k"         , &track_neglogl_bwd_k     ) ;
  //pidTree->Branch( " track_neglogl_mip   " , track_neglogl_bwd_mip   ) ;
  pidTree->Branch( "track_PIDA_mean"              , &track_PIDA_mean          ) ;
  pidTree->Branch( "track_PIDA_median"              , &track_PIDA_median          ) ;
  pidTree->Branch( "track_PIDA_kde"              , &track_PIDA_kde          ) ;
  pidTree->Branch( "track_length"            , &track_length        ) ;
  pidTree->Branch( "track_dEdx"              , &track_dEdx          ) ;
  pidTree->Branch( "track_theta"             , &track_theta         ) ;
  pidTree->Branch( "track_phi"               , &track_phi           ) ;
  pidTree->Branch( "track_nhits"             , &track_nhits         ) ;


}

void ParticleIdValidationPlots::analyze(art::Event const & e)
{
  fIsData = e.isRealData();

  std::cout << "[PARTICLEIDVALID] Using IsData: " << fIsData << " settings!" << std::endl;
  std::cout << "[PARTICLEIDVALID] Using track/calorimetry associations " << fCaloTrackAssns << std::endl;

  // Get handles to needed information
  art::Handle<std::vector<recob::Track>> trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector<art::Ptr<recob::Track>> trackCollection;
  art::fill_ptr_vector(trackCollection, trackHandle);

  art::Handle<std::vector<recob::Hit>> hitHandle;
  e.getByLabel(fHitLabel, hitHandle);

  art::FindManyP<recob::Hit> hits_from_tracks(trackHandle, e, fHitTrackAssns);
  art::FindManyP<anab::Calorimetry> calo_from_tracks(trackHandle, e, fCaloTrackAssns);

  // --------- Loop over tracks in event ---------- //
  for (auto& track : trackCollection){
    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = calo_from_tracks.at(track->ID());

    double angle = track->Theta();
    std::cout << "[PARTICLEIDVALID] Looking at track with length " << track->Length() << std::endl;

    bool TrueBragg = false;
    int True_pdg = 0;
    simb::MCParticle const* maxp_me = NULL; //pointer for the particle match we will calculate

    track_length = track->Length();
    track_theta = track->Theta();
    track_phi = track->Phi();
    track_start_x = track->Start().X();
    track_start_y = track->Start().Y();
    track_start_z = track->Start().Z();
    track_end_x = track->End().X();
    track_end_y = track->End().Y();
    track_end_z = track->End().Z();

    if (!fIsData){
      // Get true PDG from associations
      std::unordered_map<int,double> trkide;
      double maxe=-1, tote=0;

      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track->ID());

      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitHandle,e,fHitTruthAssns);

      //loop only over our hits
      for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
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
      true_PDG = maxp_me->PdgCode();
      true_start_momentum = maxp_me->P();
      true_start_x = maxp_me->Vx();
      true_start_y = maxp_me->Vy();
      true_start_z = maxp_me->Vz();
      true_end_momentum = std::sqrt( std::pow(maxp_me->EndMomentum().X(),2) + std::pow(maxp_me->EndMomentum().Y(),2) + std::pow(maxp_me->EndMomentum().Z(),2));
      true_end_x = maxp_me->EndX();
      true_end_y = maxp_me->EndY();
      true_end_z = maxp_me->EndZ();

      // Check if true particle should have a Bragg peak (if the true MCParticle has zero momentum at the end of the track)
      TrueBragg = false;
      if (maxp_me->EndPx() == 0 && maxp_me->EndPy() == 0 && maxp_me->EndPz() == 0){
        TrueBragg = true;
      }

    } // end if(!fIsData)

    // ------------------- Get Track end hits over track start hits charge -------------- //

    std::cout << "[PARTICLEIDVALID] Getting calorimetry information" << std::endl;

    // for time being, only use Y plane calorimetry
    art::Ptr< anab:: Calorimetry > calo;
    for (auto c : caloFromTrack){
      int planenum = c->PlaneID().Plane;
      if (planenum != 2) continue; // Only use calorimetry from collection plane
      calo = c;
    }
    // Check that caloFromTrack is a valid object? If not, skip track
    if (!calo){
      std::cout << "[PARTICLEIDVALID] Did not find a valid calorimetry object. Skipping track." << std::endl;
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

    std::cout << "[PARTICLEIDVALID] Making plots! " << std::endl;

    /**
     * Note that here we don't need to check if fIsData because TrueBragg is
     * false by construction if fIsData = true
     */

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
    TVector3 RecoTrackDir, TrueTrackDir;
    if (!fIsData){
      RecoTrackDir = {track->End().X()-track->Start().X(),
        track->End().Y()-track->Start().Y(),
        track->End().Z()-track->Start().Z()};
      RecoTrackDir = RecoTrackDir.Unit();

      TVector3 True_start(maxp_me->Vx(), maxp_me->Vy(), maxp_me->Vz());
      TVector3 True_end(maxp_me->EndX(), maxp_me->EndY(), maxp_me->EndZ());
      TrueTrackDir = True_end - True_start;
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
    }

    // ------------------- Now calculate PID variables and fill hists ------------------- //
    art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, e, fPIDLabel);
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
    double PIDAval_mean = -999;
    double PIDAval_median = -999;
    double PIDAval_kde = -999;
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

      if (AlgScore.fAlgName == "PIDA_mean" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
        PIDAval_mean = AlgScore.fValue;
      }// if AlgName = PIDA && fVariableType == anab::kPIDA
      if (AlgScore.fAlgName == "PIDA_medeian" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
        PIDAval_median = AlgScore.fValue;
      }// if AlgName = PIDA && fVariableType == anab::kPIDA
      if (AlgScore.fAlgName == "PIDA_kde" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
        PIDAval_kde = AlgScore.fValue;
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

    track_neglogl_fwd_mu = Bragg_fwd_mu;
    track_neglogl_bwd_mu = Bragg_bwd_mu;
    track_neglogl_fwd_p = Bragg_fwd_p;
    track_neglogl_bwd_p = Bragg_bwd_p;
    track_neglogl_fwd_pi = Bragg_fwd_pi;
    track_neglogl_bwd_pi = Bragg_bwd_pi;
    track_neglogl_fwd_k = Bragg_fwd_K;
    track_neglogl_bwd_k = Bragg_bwd_K;
    track_neglogl_fwd_mip = noBragg_fwd_MIP;
    track_dEdx = dEdxtruncmean; 
    track_PIDA_mean = PIDAval_mean;
    track_PIDA_median = PIDAval_median;
    track_PIDA_kde = PIDAval_kde;
    track_nhits = nhits;

    //double Bragg_smallest = std::min({Bragg_mu, Bragg_p, Bragg_pi, Bragg_K, noBragg_fwd_MIP});
    double Bragg_smallest = std::min({Bragg_mu, Bragg_p, noBragg_fwd_MIP});
    if (Bragg_smallest == Bragg_mu) std::cout << "PID estimate: muon" << std::endl;
    else if (Bragg_smallest == Bragg_p) std::cout << "PID estimate: proton" << std::endl;
    else if (Bragg_smallest == Bragg_pi) std::cout << "PID estimate: pion" << std::endl;
    else if (Bragg_smallest == Bragg_K) std::cout << "PID estimate: kaon" << std::endl;
    else if (Bragg_smallest == noBragg_fwd_MIP) std::cout << "PID estimate: muon (no Bragg peak)" << std::endl;

    bool PID_fwd = false;
    if (Bragg_smallest == Bragg_fwd_mu || Bragg_smallest == Bragg_fwd_p || Bragg_smallest == Bragg_fwd_pi || Bragg_smallest == Bragg_fwd_K || Bragg_smallest == noBragg_fwd_MIP) PID_fwd = true;
    std::cout << "Filling alltracks" << std::endl;
    AllTracks_neglogl_mu->Fill(Bragg_mu);
    AllTracks_neglogl_p->Fill(Bragg_p);
    AllTracks_neglogl_pi->Fill(Bragg_pi);
    AllTracks_neglogl_K->Fill(Bragg_K);
    AllTracks_neglogl_MIP->Fill(noBragg_fwd_MIP);
    AllTracks_neglogl_minmuMIP->Fill(std::min(Bragg_mu, noBragg_fwd_MIP));
    AllTracks_neglogl_mu_vslength->Fill(Bragg_mu, trklen);
    AllTracks_neglogl_p_vslength->Fill(Bragg_p, trklen);
    AllTracks_neglogl_pi_vslength->Fill(Bragg_pi, trklen);
    AllTracks_neglogl_K_vslength->Fill(Bragg_K, trklen);
    AllTracks_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP, trklen);
    AllTracks_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu, noBragg_fwd_MIP), trklen);
    AllTracks_neglogl_mu_vsangle->Fill(Bragg_mu, angle);
    AllTracks_neglogl_p_vsangle->Fill(Bragg_p, angle);
    AllTracks_neglogl_pi_vsangle->Fill(Bragg_pi, angle);
    AllTracks_neglogl_K_vsangle->Fill(Bragg_K, angle);
    AllTracks_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP, angle);
    AllTracks_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu, noBragg_fwd_MIP), angle);
    AllTracks_neglogl_mu_vsnhits->Fill(Bragg_mu, nhits);
    AllTracks_neglogl_p_vsnhits->Fill(Bragg_p, nhits);
    AllTracks_neglogl_pi_vsnhits->Fill(Bragg_pi, nhits);
    AllTracks_neglogl_K_vsnhits->Fill(Bragg_K, nhits);
    AllTracks_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP, nhits);
    AllTracks_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu, noBragg_fwd_MIP), nhits);
    AllTracks_neglogl_muvsp->Fill(Bragg_mu, Bragg_p);
    AllTracks_neglogl_muvspi->Fill(Bragg_mu, Bragg_pi);
    AllTracks_neglogl_MIPvsp->Fill(noBragg_fwd_MIP, Bragg_p);
    AllTracks_neglogl_minmuMIPvsp->Fill(std::min(Bragg_mu, noBragg_fwd_MIP), Bragg_p);
    AllTracks_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
    AllTracks_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
    AllTracks_neglogl_MIPminusp->Fill(noBragg_fwd_MIP-Bragg_p);
    AllTracks_neglogl_minmuMIPminusp->Fill(std::min(Bragg_mu, noBragg_fwd_MIP)-Bragg_p);
    AllTracks_dEdxtr_len->Fill(trklen, dEdxtruncmean);

    if (Bragg_smallest == Bragg_mu) AllTracks_smallest_neglogl->Fill(0.5);
    else if (Bragg_smallest == Bragg_p) AllTracks_smallest_neglogl->Fill(1.5);
    else if (Bragg_smallest == Bragg_pi) AllTracks_smallest_neglogl->Fill(2.5);
    else if (Bragg_smallest == Bragg_K) AllTracks_smallest_neglogl->Fill(3.5);
    else if (Bragg_smallest == noBragg_fwd_MIP)
      AllTracks_smallest_neglogl->Fill(4.5);

    if (!fIsData){
      if (TrueBragg){
        // Well-reconstructed truth matching
        if (TMath::Abs(True_pdg) == 13){ // True muons
          TrueBragg_truemu_neglogl_mu->Fill(Bragg_mu);
          TrueBragg_truemu_neglogl_p->Fill(Bragg_p);
          TrueBragg_truemu_neglogl_pi->Fill(Bragg_pi);
          TrueBragg_truemu_neglogl_K->Fill(Bragg_K);
          TrueBragg_truemu_PIDA->Fill(PIDAval_mean);
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

          if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
            if (PID_fwd) TrueBragg_correctdirection->Fill(0.5,1.5);
            else TrueBragg_correctdirection->Fill(0.5,0.5);
          }
          else{ // reco dir is wrong
            if (!PID_fwd) TrueBragg_incorrectdirection->Fill(0.5,1.5);
            else TrueBragg_incorrectdirection->Fill(0.5,0.5);
          }
        }

        else if (TMath::Abs(True_pdg) == 2212){ // True protons
          TrueBragg_truep_neglogl_mu->Fill(Bragg_mu);
          TrueBragg_truep_neglogl_p->Fill(Bragg_p);
          TrueBragg_truep_neglogl_pi->Fill(Bragg_pi);
          TrueBragg_truep_neglogl_K->Fill(Bragg_K);
          TrueBragg_truep_PIDA->Fill(PIDAval_mean);
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

          if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
            if (PID_fwd) TrueBragg_correctdirection->Fill(1.5,1.5);
            else TrueBragg_correctdirection->Fill(1.5,0.5);
          }
          else{ // reco dir is wrong
            if (!PID_fwd) TrueBragg_incorrectdirection->Fill(1.5,1.5);
            else TrueBragg_incorrectdirection->Fill(1.5,0.5);
          }
        }

        else if (TMath::Abs(True_pdg) == 211){ // True pions
          TrueBragg_truepi_neglogl_mu->Fill(Bragg_mu);
          TrueBragg_truepi_neglogl_p->Fill(Bragg_p);
          TrueBragg_truepi_neglogl_pi->Fill(Bragg_pi);
          TrueBragg_truepi_neglogl_K->Fill(Bragg_K);
          TrueBragg_truepi_PIDA->Fill(PIDAval_mean);
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

          if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
            if (PID_fwd) TrueBragg_correctdirection->Fill(2.5,1.5);
            else TrueBragg_correctdirection->Fill(2.5,0.5);
          }
          else{ // reco dir is wrong
            if (!PID_fwd) TrueBragg_incorrectdirection->Fill(2.5,1.5);
            else TrueBragg_incorrectdirection->Fill(2.5,0.5);
          }
        }

        else if (TMath::Abs(True_pdg) == 321){ // True kaons
          TrueBragg_trueK_neglogl_mu->Fill(Bragg_mu);
          TrueBragg_trueK_neglogl_p->Fill(Bragg_p);
          TrueBragg_trueK_neglogl_pi->Fill(Bragg_pi);
          TrueBragg_trueK_neglogl_K->Fill(Bragg_K);
          TrueBragg_trueK_PIDA->Fill(PIDAval_mean);
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

          if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
            if (PID_fwd) TrueBragg_correctdirection->Fill(3.5,1.5);
            else TrueBragg_correctdirection->Fill(3.5,0.5);
          }
          else{ // reco dir is wrong
            if (!PID_fwd) TrueBragg_incorrectdirection->Fill(3.5,1.5);
            else TrueBragg_incorrectdirection->Fill(3.5,0.5);
          }
        }


        else if (TMath::Abs(True_pdg) == 11){ // True electrons
          TrueBragg_truee_neglogl_mu->Fill(Bragg_mu);
          TrueBragg_truee_neglogl_p->Fill(Bragg_p);
          TrueBragg_truee_neglogl_pi->Fill(Bragg_pi);
          TrueBragg_truee_neglogl_K->Fill(Bragg_K);
          TrueBragg_truee_PIDA->Fill(PIDAval_mean);
          TrueBragg_truee_dEdxtr_len->Fill(trklen,dEdxtruncmean);
          TrueBragg_truee_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
          TrueBragg_truee_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
          TrueBragg_truee_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
          TrueBragg_truee_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
          TrueBragg_truee_neglogl_MIPvsp->Fill(noBragg_fwd_MIP,Bragg_p);
          TrueBragg_truee_neglogl_minmuMIPvsp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),Bragg_p);
          TrueBragg_truee_neglogl_MIPminusp->Fill(noBragg_fwd_MIP-Bragg_p);
          TrueBragg_truee_neglogl_minmuMIPminusp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP)-Bragg_p);
          TrueBragg_truee_neglogl_MIP->Fill(noBragg_fwd_MIP);
          TrueBragg_truee_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
          TrueBragg_truee_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
          TrueBragg_truee_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
          TrueBragg_truee_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
          TrueBragg_truee_neglogl_p_vslength->Fill(Bragg_p,trklen);
          TrueBragg_truee_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
          TrueBragg_truee_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
          TrueBragg_truee_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
          TrueBragg_truee_neglogl_p_vsangle->Fill(Bragg_p,angle);
          TrueBragg_truee_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
          TrueBragg_truee_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
          TrueBragg_truee_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
          TrueBragg_truee_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

          if (Bragg_smallest == Bragg_mu) TrueBragg_truee_smallest_neglogl->Fill(0.5);
          else if (Bragg_smallest == Bragg_p) TrueBragg_truee_smallest_neglogl->Fill(1.5);
          else if (Bragg_smallest == Bragg_pi) TrueBragg_truee_smallest_neglogl->Fill(2.5);
          else if (Bragg_smallest == Bragg_K) TrueBragg_truee_smallest_neglogl->Fill(3.5);
          else if (Bragg_smallest == noBragg_fwd_MIP)
            TrueBragg_truee_smallest_neglogl->Fill(4.5);
        }
      } // end if(TrueBragg)

      // All particles
      if (TMath::Abs(True_pdg) == 13){ // True muons
        All_truemu_neglogl_mu->Fill(Bragg_mu);
        All_truemu_neglogl_p->Fill(Bragg_p);
        All_truemu_neglogl_pi->Fill(Bragg_pi);
        All_truemu_neglogl_K->Fill(Bragg_K);
        All_truemu_PIDA->Fill(PIDAval_mean);
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

        if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
          if (PID_fwd) All_correctdirection->Fill(0.5,1.5);
          else All_correctdirection->Fill(0.5,0.5);
        }
        else{ // reco dir is wrong
          if (!PID_fwd) All_incorrectdirection->Fill(0.5,1.5);
          else All_incorrectdirection->Fill(0.5,0.5);
        }
      }

      else if (TMath::Abs(True_pdg) == 2212){ // True protons
        All_truep_neglogl_mu->Fill(Bragg_mu);
        All_truep_neglogl_p->Fill(Bragg_p);
        All_truep_neglogl_pi->Fill(Bragg_pi);
        All_truep_neglogl_K->Fill(Bragg_K);
        All_truep_PIDA->Fill(PIDAval_mean);
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

        if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
          if (PID_fwd) All_correctdirection->Fill(1.5,1.5);
          else All_correctdirection->Fill(1.5,0.5);
        }
        else{ // reco dir is wrong
          if (!PID_fwd) All_incorrectdirection->Fill(1.5,1.5);
          else All_incorrectdirection->Fill(1.5,0.5);
        }
      }

      else if (TMath::Abs(True_pdg) == 211){ // True pions
        All_truepi_neglogl_mu->Fill(Bragg_mu);
        All_truepi_neglogl_p->Fill(Bragg_p);
        All_truepi_neglogl_pi->Fill(Bragg_pi);
        All_truepi_neglogl_K->Fill(Bragg_K);
        All_truepi_PIDA->Fill(PIDAval_mean);
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

        if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
          if (PID_fwd) All_correctdirection->Fill(20.5,1.5);
          else All_correctdirection->Fill(2.5,0.5);
        }
        else{ // reco dir is wrong
          if (!PID_fwd) All_incorrectdirection->Fill(2.5,1.5);
          else All_incorrectdirection->Fill(2.5,0.5);
        }
      }

      else if (TMath::Abs(True_pdg) == 321){ // True kaons
        All_trueK_neglogl_mu->Fill(Bragg_mu);
        All_trueK_neglogl_p->Fill(Bragg_p);
        All_trueK_neglogl_pi->Fill(Bragg_pi);
        All_trueK_neglogl_K->Fill(Bragg_K);
        All_trueK_PIDA->Fill(PIDAval_mean);
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

        if (RecoTrackDir.Dot(TrueTrackDir) > 0){ // reco dir is right
          if (PID_fwd) All_correctdirection->Fill(3.5,1.5);
          else All_correctdirection->Fill(3.5,0.5);
        }
        else{ // reco dir is wrong
          if (!PID_fwd) All_incorrectdirection->Fill(3.5,1.5);
          else All_incorrectdirection->Fill(3.5,0.5);
        }
      }

      else if (TMath::Abs(True_pdg) == 11){ // True electrons
        All_truee_neglogl_mu->Fill(Bragg_mu);
        All_truee_neglogl_p->Fill(Bragg_p);
        All_truee_neglogl_pi->Fill(Bragg_pi);
        All_truee_neglogl_K->Fill(Bragg_K);
        All_truee_PIDA->Fill(PIDAval_mean);
        All_truee_dEdxtr_len->Fill(trklen,dEdxtruncmean);
        All_truee_neglogl_muvsp->Fill(Bragg_mu,Bragg_p);
        All_truee_neglogl_muvspi->Fill(Bragg_mu,Bragg_pi);
        All_truee_neglogl_muoverp->Fill(Bragg_mu/Bragg_p);
        All_truee_neglogl_muminusp->Fill(Bragg_mu-Bragg_p);
        All_truee_neglogl_MIPvsp->Fill(noBragg_fwd_MIP,Bragg_p);
        All_truee_neglogl_minmuMIPvsp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),Bragg_p);
        All_truee_neglogl_MIPminusp->Fill(noBragg_fwd_MIP-Bragg_p);
        All_truee_neglogl_minmuMIPminusp->Fill(std::min(Bragg_mu,noBragg_fwd_MIP)-Bragg_p);
        All_truee_neglogl_MIP->Fill(noBragg_fwd_MIP);
        All_truee_neglogl_minmuMIP->Fill(std::min(Bragg_mu,noBragg_fwd_MIP));
        All_truee_neglogl_mu_vslength->Fill(Bragg_mu,trklen);
        All_truee_neglogl_MIP_vslength->Fill(noBragg_fwd_MIP,trklen);
        All_truee_neglogl_minmuMIP_vslength->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),trklen);
        All_truee_neglogl_p_vslength->Fill(Bragg_p,trklen);
        All_truee_neglogl_mu_vsangle->Fill(Bragg_mu,angle);
        All_truee_neglogl_MIP_vsangle->Fill(noBragg_fwd_MIP,angle);
        All_truee_neglogl_minmuMIP_vsangle->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),angle);
        All_truee_neglogl_p_vsangle->Fill(Bragg_p,angle);
        All_truee_neglogl_mu_vsnhits->Fill(Bragg_mu,nhits);
        All_truee_neglogl_MIP_vsnhits->Fill(noBragg_fwd_MIP,nhits);
        All_truee_neglogl_minmuMIP_vsnhits->Fill(std::min(Bragg_mu,noBragg_fwd_MIP),nhits);
        All_truee_neglogl_p_vsnhits->Fill(Bragg_p,nhits);

        if (Bragg_smallest == Bragg_mu) All_truee_smallest_neglogl->Fill(0.5);
        else if (Bragg_smallest == Bragg_p) All_truee_smallest_neglogl->Fill(1.5);
        else if (Bragg_smallest == Bragg_pi) All_truee_smallest_neglogl->Fill(2.5);
        else if (Bragg_smallest == Bragg_K) All_truee_smallest_neglogl->Fill(3.5);
        else if (Bragg_smallest == noBragg_fwd_MIP)
          All_truee_smallest_neglogl->Fill(4.5);
      }
    }// end !fIsData

    std::cout << "[PARTICLEIDVALID] Filling tree. " << std::endl;
    pidTree->Fill();

  } // Loop over tracks


}

DEFINE_ART_MODULE(ParticleIdValidationPlots)
