////////////////////////////////////////////////////////////////////////
// Class:       ParticleId
// Plugin Type: producer (art v2_05_01)
// File:        ParticleId_module.cc
//
// Generated at Wed Jan 31 11:25:52 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 *
 * \class ParticleId
 *
 * \brief ParticleId producer module
 *
 * \author Kirsty Duffy (kduffy@fnal.gov), Adam Lister (alister1@lancaster.ac.uk)
 *
 * \date 2018/04/18
 *
 */

// base includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"

// data products
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larana/TruncatedMean/Algorithm/TruncMean.h"

// local includes
#include "uboone/ParticleID/Algorithms/GetDaughterTracksShowers.h"
#include "uboone/ParticleID/Algorithms/FiducialVolume.h"
#include "uboone/ParticleID/Algorithms/PIDA.h"
#include "uboone/ParticleID/Algorithms/Bragg_negLogL_Estimator.h"
#include "uboone/ParticleID/Algorithms/LandauGaussian.h"

// root includes
#include "TVector3.h"

// cpp includes
#include <memory>

namespace UBPID{
  class ParticleId;
}

class UBPID::ParticleId : public art::EDProducer {
  public:
    explicit ParticleId(fhicl::ParameterSet const & p);

    ParticleId();
    ParticleId(ParticleId const &) = delete;
    ParticleId(ParticleId &&) = delete;
    ParticleId & operator = (ParticleId const &) = delete;
    ParticleId & operator = (ParticleId &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    std::vector<double> fv;

  private:

    // fcl
    std::string fTrackLabel;
    std::string fCaloLabel;
    double fCutDistance;
    double fCutFraction;

    // fidvol related
    fidvol::fiducialVolume fid;
    
    // for PIDA
    particleid::PIDA pida;

    // for likelihood-based PID
    particleid::Bragg_negLogL_Estimator braggcalc;

    // For truncated mean
    TruncMean trm;

    //other
    bool isData;
};


UBPID::ParticleId::ParticleId(fhicl::ParameterSet const & p)
{

  fhicl::ParameterSet const p_fv     = p.get<fhicl::ParameterSet>("FiducialVolume");
  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");
  fhicl::ParameterSet const p_bragg  = p.get<fhicl::ParameterSet>("BraggAlgo");

  // fcl parameters
  fTrackLabel = p_labels.get< std::string > ("TrackLabel");
  fCaloLabel = p_labels.get< std::string > ("CalorimetryLabel");
  fCutDistance  = p.get< double > ("DaughterFinderCutDistance");
  fCutFraction  = p.get< double > ("DaughterFinderCutFraction");

  fv = fid.setFiducialVolume(fv, p_fv);
  fid.printFiducialVolume(fv);
  braggcalc.configure(p_bragg);
  braggcalc.printConfiguration();

  // this module produces a anab::ParticleID object and
  // an association to the track which produced it
  produces< std::vector<anab::ParticleID> >();
  produces< art::Assns< recob::Track, anab::ParticleID> >();

}

void UBPID::ParticleId::produce(art::Event & e)
{

  bool isData = e.isRealData();

  if (!isData) std::cout << "[ParticleID] Running simulated data." << std::endl;
  else std::cout << "[ParticleID] Running physics data." << std::endl;

  // produce collection of particleID objects
  std::unique_ptr< std::vector<anab::ParticleID> > particleIDCollection( new std::vector<anab::ParticleID> );
  std::unique_ptr< art::Assns <recob::Track, anab::ParticleID> > trackParticleIdAssn( new art::Assns<recob::Track, anab::ParticleID> );

  //
  // get handles to needed information
  //

  // tracks...
  art::Handle < std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector< art::Ptr<recob::Track> > trackCollection;
  art::fill_ptr_vector(trackCollection, trackHandle);

  // calorimetry object...
  art::FindManyP<anab::Calorimetry> caloFromTracks(trackHandle, e, fCaloLabel);

  for (auto& track : trackCollection){

    // Skip tracks/events with no valid calorimetry object associated to it
    if (!caloFromTracks.isValid()){
      std::cout << "[ParticleID] Calorimetry<->Track associations are not valid for this track. Skipping." << std::endl;
      continue;
    }

    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = caloFromTracks.at(track->ID());

    // for time being, only use Y plane calorimetry
    art::Ptr< anab:: Calorimetry > calo;
    int planenum = -1;
    for (auto c : caloFromTrack){
      planenum = c->PlaneID().Plane;
      if (planenum != 2) continue; // Only use calorimetry from collection plane
      calo = c;
    }
    
    // Check that caloFromTrack is a valid object
    if (!calo){
      std::cout << "[ParticleID] Calorimetry on plane " << planenum << " is unavailable. Skipping." << std::endl;
      continue;
    }

    std::vector<double> dEdx = calo->dEdx();
    std::vector<double> resRange = calo->ResidualRange();

    // int nDaughters = GetNDaughterTracks((*trackHandle), track->ID(), fCutDistance, fCutFraction);
    // std::cout << "[ParticleID]  Found track with " << nDaughters << " reconstructed daughters." << std::endl;

    // Vairables for ParticleID Class
    std::vector<anab::sParticleIDAlgScores> AlgScoresVec;
    anab::sParticleIDAlgScores Bragg_fwd_mu;
    anab::sParticleIDAlgScores Bragg_fwd_p;
    anab::sParticleIDAlgScores Bragg_fwd_pi;
    anab::sParticleIDAlgScores Bragg_fwd_k;
    anab::sParticleIDAlgScores Bragg_bwd_mu;
    anab::sParticleIDAlgScores Bragg_bwd_p;
    anab::sParticleIDAlgScores Bragg_bwd_pi;
    anab::sParticleIDAlgScores Bragg_bwd_k;
    anab::sParticleIDAlgScores noBragg_fwd_MIP;
    anab::sParticleIDAlgScores PIDAval_mean;
    anab::sParticleIDAlgScores PIDAval_median;
    anab::sParticleIDAlgScores PIDAval_kde;
    anab::sParticleIDAlgScores dEdxtruncmean;
    anab::sParticleIDAlgScores trklen;

    // bool   isContained = false;

    // Evaluate PID only for fully-contained particles
    // TVector3 trackStart = track->Vertex();
    // TVector3 trackEnd = track->End();

    // if (fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){
    //   isContained = true;
    // }

    // if (isContained){

    // Check if particle has reconstructed "daughters" - if it does, there may be no Bragg peak
    // and PID might not be accurate
    // if (nDaughters == 0){

    // std::cout << "[ParticleID]  >> Track is fully contained and has no daughters " << std::endl;

    /**
     * Algorithm 1: PIDA
     * This makes use of Bruce home-brewed PIDA calculation, which can be 
     * calculated via three methods:
     * (1) mean (original implementation from B. Baller)
     * (2) median (T. Yang & V. Meddage)
     * (3) kernel density estimator (A. Lister)
     */

    // mean
    PIDAval_mean.fAlgName = "PIDA_mean";
    PIDAval_mean.fVariableType = anab::kPIDA;
    PIDAval_mean.fValue = pida.getPida(dEdx, resRange, "mean");
    AlgScoresVec.push_back(PIDAval_mean);

    // median
    PIDAval_median.fAlgName = "PIDA_median";
    PIDAval_median.fVariableType = anab::kPIDA;
    PIDAval_median.fValue = pida.getPida(dEdx, resRange, "median");
    AlgScoresVec.push_back(PIDAval_median);

    // median
    PIDAval_kde.fAlgName = "PIDA_kde";
    PIDAval_kde.fVariableType = anab::kPIDA;
    PIDAval_kde.fValue = pida.getPida(dEdx, resRange, "kde");
    AlgScoresVec.push_back(PIDAval_kde);

    /**
     * Algorithm 2: BraggPeakLLH
     * Uses B. Ballers theory, along with landau-gaussian distributions with
     * widths measured from data and simulation to estimate the likelihood for 
     * each hit in a track to have come from each particle species.
     */
    Bragg_fwd_mu.fAlgName = "BraggPeakLLH";
    Bragg_fwd_p.fAlgName  = "BraggPeakLLH";
    Bragg_fwd_pi.fAlgName = "BraggPeakLLH";
    Bragg_fwd_k.fAlgName  = "BraggPeakLLH";
    Bragg_bwd_mu.fAlgName = "BraggPeakLLH";
    Bragg_bwd_p.fAlgName  = "BraggPeakLLH";
    Bragg_bwd_pi.fAlgName = "BraggPeakLLH";
    Bragg_bwd_k.fAlgName  = "BraggPeakLLH";
    Bragg_fwd_mu.fVariableType = anab::kLogL_fwd;
    Bragg_fwd_p.fVariableType  = anab::kLogL_fwd;
    Bragg_fwd_pi.fVariableType = anab::kLogL_fwd;
    Bragg_fwd_k.fVariableType  = anab::kLogL_fwd;
    Bragg_bwd_mu.fVariableType = anab::kLogL_bwd;
    Bragg_bwd_p.fVariableType  = anab::kLogL_bwd;
    Bragg_bwd_pi.fVariableType = anab::kLogL_bwd;
    Bragg_bwd_k.fVariableType  = anab::kLogL_bwd;
    Bragg_fwd_mu.fAssumedPdg = 13;
    Bragg_fwd_p.fAssumedPdg = 2212;
    Bragg_fwd_pi.fAssumedPdg = 211;
    Bragg_fwd_k.fAssumedPdg = 321;
    Bragg_bwd_mu.fAssumedPdg = 13;
    Bragg_bwd_p.fAssumedPdg = 2212;
    Bragg_bwd_pi.fAssumedPdg = 211;
    Bragg_bwd_k.fAssumedPdg = 321;
    Bragg_fwd_mu.fValue = braggcalc.getNegLogL(dEdx, resRange, Bragg_fwd_mu.fAssumedPdg, true);
    Bragg_fwd_p.fValue  = braggcalc.getNegLogL(dEdx, resRange, Bragg_fwd_p.fAssumedPdg,  true);
    Bragg_fwd_pi.fValue = braggcalc.getNegLogL(dEdx, resRange, Bragg_fwd_pi.fAssumedPdg, true);
    Bragg_fwd_k.fValue  = braggcalc.getNegLogL(dEdx, resRange, Bragg_fwd_k.fAssumedPdg,  true);
    Bragg_bwd_mu.fValue = braggcalc.getNegLogL(dEdx, resRange, Bragg_bwd_mu.fAssumedPdg, false);
    Bragg_bwd_p.fValue  = braggcalc.getNegLogL(dEdx, resRange, Bragg_bwd_p.fAssumedPdg,  false);
    Bragg_bwd_pi.fValue = braggcalc.getNegLogL(dEdx, resRange, Bragg_bwd_pi.fAssumedPdg, false);
    Bragg_bwd_k.fValue  = braggcalc.getNegLogL(dEdx, resRange, Bragg_bwd_k.fAssumedPdg,  false);

    // Special case: MIP-like probability
    // fit to the flat MIP region of dEdx with residual range > 15 cm
    noBragg_fwd_MIP.fAlgName = "BraggPeakLLH";
    noBragg_fwd_MIP.fVariableType = anab::kLogL_fwd;
    noBragg_fwd_MIP.fAssumedPdg = 0;
    noBragg_fwd_MIP.fValue = braggcalc.getNegLogL(dEdx, resRange, noBragg_fwd_MIP.fAssumedPdg, true);

    AlgScoresVec.push_back(Bragg_fwd_mu);
    AlgScoresVec.push_back(Bragg_fwd_p);
    AlgScoresVec.push_back(Bragg_fwd_pi);
    AlgScoresVec.push_back(Bragg_fwd_k);
    AlgScoresVec.push_back(Bragg_bwd_mu);
    AlgScoresVec.push_back(Bragg_bwd_p);
    AlgScoresVec.push_back(Bragg_bwd_pi);
    AlgScoresVec.push_back(Bragg_bwd_k);
    AlgScoresVec.push_back(noBragg_fwd_MIP);

    /**
     * Algorithm 3: Truncated mean dE/dx versus track length
     * Makes use of the "Truncated Mean" algorithm developed by D. Caratelli
     * to plot the truncated mean dE/dx  of a track
     * versus its length for separation.
     */
    dEdxtruncmean.fAlgName = "TruncatedMean";
    dEdxtruncmean.fVariableType = anab::kdEdxtruncmean;
    trklen.fAlgName = "TruncatedMean";
    trklen.fVariableType = anab::kTrackLength;

    size_t nmin = 1;
    size_t nmax = 1;
    const size_t currentiteration = 0;
    const size_t lmin = 1;
    const float convergencelimit = 0.1;
    const float nsigma = 1.0;
    dEdxtruncmean.fValue = (double)trm.CalcIterativeTruncMean(dEdx, nmin, nmax, currentiteration, lmin, convergencelimit, nsigma);
    trklen.fValue = track->Length();

    AlgScoresVec.push_back(dEdxtruncmean);
    AlgScoresVec.push_back(trklen);
      /*  }

          } // end if(isContained)
          else{
      // If particle is *not* contained, assume it is a muon
      // Set pdg=13. No other PID variables will be set because those are only filled for contained particles
      pdg = 13;
      }*/


      // -------------------------------------------------------------------------- //
      // Finally, fill product with the variables that we calculated above and make assns

      std::cout << "[ParticleID] >> Making particleIDCollection... " << std::endl;
      anab::ParticleID PID_object(AlgScoresVec);
    particleIDCollection->push_back(PID_object);

    std::cout << "[ParticleID] Making assn... " << std::endl;
    util::CreateAssn(*this, e, *particleIDCollection, track, *trackParticleIdAssn);

  }

  e.put(std::move(particleIDCollection));
  e.put(std::move(trackParticleIdAssn));

}

DEFINE_ART_MODULE(UBPID::ParticleId)
