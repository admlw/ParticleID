////////////////////////////////////////////////////////////////////////
// Class:       ParticleId
// Plugin Type: producer (art v2_05_01)
// File:        ParticleId_module.cc
//
// Generated at Wed Jan 31 11:25:52 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

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

// Truncated mean calculator
//#include "larana/TruncatedMean/Algorithm/TruncMean.h"
// UPDATE WITH NEW CHECKOUT!
#include "uboone/ParticleID/TruncatedMeanCopy/Algorithm/TruncMean.h"


// local includes
#include "uboone/ParticleID/Algorithms/GetDaughterTracksShowers.h"
#include "uboone/ParticleID/Algorithms/fiducialVolume.h"
#include "uboone/ParticleID/Algorithms/PIDA.h"
#include "uboone/ParticleID/Algorithms/Bragg_negLogL_Estimator.h"

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
  
  // Selected optional functions.
  void beginJob() override;
  void endJob() override;
  
  std::vector<double> fv;
  
private:
  
  // fcl
  std::string fTrackingAlgo; 
  std::string fCaloLabel;
  std::string fPidaType;
  double fCutDistance;
  double fCutFraction;
  double fBraggWidthMu;
  double fBraggWidthP;
  double fBraggWidthPi;
  double fBraggWidthK;
  
  // fidvol related
  fidvol::fiducialVolume fid;
  particleid::PIDA pida;
  
  // for likelihood-based PID
  particleid::Bragg_negLogL_Estimator braggcalc;
  
  // For truncated mean
  //TruncMean trm;
  
  //other
  bool isData;
};


UBPID::ParticleId::ParticleId(fhicl::ParameterSet const & p)
{
  
  // fcl parameters
  fTrackingAlgo = p.get< std::string > ("TrackingAlgorithm");
  fCaloLabel = p.get< std::string > ("CalorimetryModule");
  fPidaType = p.get< std::string > ("PIDACalcType");
  fCutDistance  = p.get< double > ("DaughterFinderCutDistance");
  fCutFraction  = p.get< double > ("DaughterFinderCutFraction");
  fBraggWidthMu = p.get< double > ("dEdxWidthMu",0.1);
  fBraggWidthP  = p.get< double > ("dEdxWidthP",0.2);
  fBraggWidthPi = p.get< double > ("dEdxWidthPi",0.1);
  fBraggWidthK  = p.get< double > ("dEdxWidthK",0.1);
  
  fv = fid.setFiducialVolume(fv, p);
  fid.printFiducialVolume(fv);
  
  braggcalc.setWidthMu(fBraggWidthMu);
  braggcalc.setWidthP(fBraggWidthP);
  braggcalc.setWidthPi(fBraggWidthPi);
  braggcalc.setWidthK(fBraggWidthK);
  
  // this module produces a anab::ParticleID object and
  // an association to the track which produced it
  produces< std::vector<anab::ParticleID> >();
  produces< art::Assns< recob::Track, anab::ParticleID> >();
  
}

void UBPID::ParticleId::beginJob()
{
  // Implementation of optional member function here.
}

void UBPID::ParticleId::produce(art::Event & e)
{
  
  bool isData = e.isRealData();
  
  if (!isData) std::cout << "[ParticleID]  Running Simulated Data" << std::endl;
  
  // produce collection of particleID objects
  std::unique_ptr< std::vector<anab::ParticleID> > particleIDCollection( new std::vector<anab::ParticleID> );
  std::unique_ptr< art::Assns <recob::Track, anab::ParticleID> > trackParticleIdAssn( new art::Assns<recob::Track, anab::ParticleID> );
  
  //
  // get handles to needed information
  //

  // tracks...
  art::Handle < std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackingAlgo, trackHandle);
  std::vector< art::Ptr<recob::Track> > trackCollection;
  art::fill_ptr_vector(trackCollection, trackHandle);
  
  // calorimetry object...
  art::FindManyP<anab::Calorimetry> caloFromTracks(trackHandle, e, fCaloLabel);
  
  for (auto& track : trackCollection){
						      
    // Skip tracks/events wtih no valid calorimetry object associated to it
    if (!caloFromTracks.isValid()){
      std::cout << "Did not find valid calorimetry object for this event. Skipping track..." << std::endl;
      continue;
    }
												
    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = caloFromTracks.at(track->ID());
		   
    // for time being, only use Y plane calorimetry
    art::Ptr< anab:: Calorimetry > calo;
    for (auto c : caloFromTrack){
      int planenum = c->PlaneID().Plane;
      if (planenum != 2) continue; // Only use calorimetry from collection plane
      calo = c;
    }
    // Check that caloFromTrack is a valid object? If not, do what? Skip track? Return nothing?
    if (!calo){
      std::cout << "Did not find a valid calorimetry object. Skipping track." << std::endl;
      continue;
    }
    
    std::vector<double> dEdx = calo->dEdx();
    std::vector<double> dQdx = calo->dQdx();
    std::vector<double> resRange = calo->ResidualRange();

    int nDaughters = GetNDaughterTracks((*trackHandle), track->ID(), fCutDistance, fCutFraction);

    std::cout << "[ParticleID]  Found track with " << nDaughters << " reconstructed daughters." << std::endl;

    // Vairables for ParticleID Class
    std::vector<anab::sParticleIDAlgScores> AlgScoresVec;
    anab::sParticleIDAlgScores Bragg_fwd_mu;
    anab::sParticleIDAlgScores Bragg_fwd_p;
    anab::sParticleIDAlgScores Bragg_fwd_pi;
    anab::sParticleIDAlgScores Bragg_fwd_K;
    anab::sParticleIDAlgScores Bragg_bwd_mu;
    anab::sParticleIDAlgScores Bragg_bwd_p;
    anab::sParticleIDAlgScores Bragg_bwd_pi;
    anab::sParticleIDAlgScores Bragg_bwd_K;
    anab::sParticleIDAlgScores PIDAval;
    //anab::sParticleIDAlgScores dEdxtruncmean;
    //anab::sParticleIDAlgScores dQdxtruncmean;
    anab::sParticleIDAlgScores trklen;
    // bool   isContained = false;

    
    // -------------------------------------------------------------------------- //
    // Start calculating PID variables

    // Evaluate PID only for fully-contained particles
    TVector3 trackStart = track->Vertex();
    TVector3 trackEnd = track->End();

    //if (fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){
    //  isContained = true;
    //}
    
    //if (isContained){

      // Check if particle has reconstructed "daughters" - if it does, there may be no Bragg peak
      // and PID might not be accurate
      //if (nDaughters == 0){

    //std::cout << "[ParticleID]  >> Track is fully contained and has no daughters " << std::endl;

	// ------ Algorithm 1:
	// ------ PIDA ------ //
	PIDAval.fAlgName = "PIDA";
	PIDAval.fVariableType = anab::kPIDA;
	PIDAval.fValue = pida.getPida(dEdx, resRange, fPidaType);

	AlgScoresVec.push_back(PIDAval);
	
        //std::cout << "[ParticleID] >> PIDA value: " << PIDAval.fValue << std::endl;

	// ------ Algorithm 2:
	// ------ Likelihood compared to Bragg peak theoretical prediction ------ //
	Bragg_fwd_mu.fAlgName = "BraggPeakLLH";
	Bragg_fwd_p.fAlgName  = "BraggPeakLLH";
	Bragg_fwd_pi.fAlgName = "BraggPeakLLH";
	Bragg_fwd_K.fAlgName  = "BraggPeakLLH";
	Bragg_bwd_mu.fAlgName = "BraggPeakLLH";
	Bragg_bwd_p.fAlgName  = "BraggPeakLLH";
	Bragg_bwd_pi.fAlgName = "BraggPeakLLH";
	Bragg_bwd_K.fAlgName  = "BraggPeakLLH";
	Bragg_fwd_mu.fVariableType = anab::kLogL_fwd;
	Bragg_fwd_p.fVariableType  = anab::kLogL_fwd;
	Bragg_fwd_pi.fVariableType = anab::kLogL_fwd;
	Bragg_fwd_K.fVariableType  = anab::kLogL_fwd;
	Bragg_bwd_mu.fVariableType = anab::kLogL_bwd;
	Bragg_bwd_p.fVariableType  = anab::kLogL_bwd;
	Bragg_bwd_pi.fVariableType = anab::kLogL_bwd;
	Bragg_bwd_K.fVariableType  = anab::kLogL_bwd;
	Bragg_fwd_mu.fAssumedPdg = 13;
	Bragg_fwd_p.fAssumedPdg = 2212;
	Bragg_fwd_pi.fAssumedPdg = 211;
	Bragg_fwd_K.fAssumedPdg = 321;
	Bragg_bwd_mu.fAssumedPdg = 13;
	Bragg_bwd_p.fAssumedPdg = 2212;
	Bragg_bwd_pi.fAssumedPdg = 211;
	Bragg_bwd_K.fAssumedPdg = 321;
	Bragg_fwd_mu.fValue = braggcalc.getNegLogL(dEdx, resRange, Bragg_fwd_mu.fAssumedPdg, true);
	Bragg_fwd_p.fValue  = braggcalc.getNegLogL(dEdx, resRange, Bragg_fwd_p.fAssumedPdg,  true);
	Bragg_fwd_pi.fValue = braggcalc.getNegLogL(dEdx, resRange, Bragg_fwd_pi.fAssumedPdg, true);
	Bragg_fwd_K.fValue  = braggcalc.getNegLogL(dEdx, resRange, Bragg_fwd_K.fAssumedPdg,  true);
	Bragg_bwd_mu.fValue = braggcalc.getNegLogL(dEdx, resRange, Bragg_bwd_mu.fAssumedPdg, false);
	Bragg_bwd_p.fValue  = braggcalc.getNegLogL(dEdx, resRange, Bragg_bwd_p.fAssumedPdg,  false);
	Bragg_bwd_pi.fValue = braggcalc.getNegLogL(dEdx, resRange, Bragg_bwd_pi.fAssumedPdg, false);
	Bragg_bwd_K.fValue  = braggcalc.getNegLogL(dEdx, resRange, Bragg_bwd_K.fAssumedPdg,  false);

	AlgScoresVec.push_back(Bragg_fwd_mu);
	AlgScoresVec.push_back(Bragg_fwd_p);
	AlgScoresVec.push_back(Bragg_fwd_pi);
	AlgScoresVec.push_back(Bragg_fwd_K);
	AlgScoresVec.push_back(Bragg_bwd_mu);
	AlgScoresVec.push_back(Bragg_bwd_p);
	AlgScoresVec.push_back(Bragg_bwd_pi);
	AlgScoresVec.push_back(Bragg_bwd_K);

	// Use best fit (lowest neg2LogL) from forward vs backward for calculations
	double Bragg_mu = (Bragg_fwd_mu.fValue < Bragg_bwd_mu.fValue ? Bragg_fwd_mu.fValue : Bragg_bwd_mu.fValue);
	double Bragg_p  = (Bragg_fwd_p.fValue  < Bragg_bwd_p.fValue  ? Bragg_fwd_p.fValue  : Bragg_bwd_p.fValue);
	double Bragg_pi = (Bragg_fwd_pi.fValue < Bragg_bwd_pi.fValue ? Bragg_fwd_pi.fValue : Bragg_bwd_pi.fValue);
	double Bragg_K  = (Bragg_fwd_K.fValue  < Bragg_bwd_K.fValue  ? Bragg_fwd_K.fValue  : Bragg_bwd_K.fValue);

	// Calculate likelihood ratio on a scale of 0 to 1
        double Bragg_pull_mu = Bragg_mu/(Bragg_mu+Bragg_p+Bragg_pi+Bragg_K);
	double Bragg_pull_p  = Bragg_p /(Bragg_mu+Bragg_p+Bragg_pi+Bragg_K);
	double Bragg_pull_pi = Bragg_pi/(Bragg_mu+Bragg_p+Bragg_pi+Bragg_K);
	double Bragg_pull_K  = Bragg_K /(Bragg_mu+Bragg_p+Bragg_pi+Bragg_K);
	double Bragg_pull_MIPproton = (Bragg_mu+Bragg_pi)/(Bragg_mu+Bragg_pi+Bragg_p);

	// Return PDG code based on most likely particle
	std::cout << "neglogl fwd mu = " << Bragg_fwd_mu.fValue << std::endl;
	std::cout << "neglogl fwd p  = " << Bragg_fwd_p.fValue  << std::endl;
	std::cout << "neglogl fwd pi = " << Bragg_fwd_pi.fValue << std::endl;
	std::cout << "neglogl fwd K  = " << Bragg_fwd_K.fValue  << std::endl;
	std::cout << "neglogl fwd mu = " << Bragg_bwd_mu.fValue << std::endl;
	std::cout << "neglogl bwd p  = " << Bragg_bwd_p.fValue  << std::endl;
	std::cout << "neglogl bwd pi = " << Bragg_bwd_pi.fValue << std::endl;
	std::cout << "neglogl bwd K  = " << Bragg_bwd_K.fValue  << std::endl;
								    
	std::cout << "Bragg_pull_mu = " << Bragg_pull_mu << std::endl;
	std::cout << "Bragg_pull_p  = " << Bragg_pull_p  << std::endl;
	std::cout << "Bragg_pull_pi = " << Bragg_pull_pi << std::endl;
	std::cout << "Bragg_pull_K  = " << Bragg_pull_K  << std::endl << std::endl;
							    
        std::cout << "Bragg_pull_MIP/proton = " << Bragg_pull_MIPproton << std::endl;
	
	// ------ Algorithm 3:
	// ------ Truncated mean dE/dx vs track length ------ //
	/*dQdxtruncmean.fAlgName = "TruncatedMean";
	dQdxtruncmean.fVariableType = anab::kdQdxtruncmean;
        dEdxtruncmean.fAlgName = "TruncatedMean";
        dEdxtruncmean.fVariableType = anab::kdEdxtruncmean;
        trklen.fAlgName = "TruncatedMean";
        trklen.fVariableType = anab::kTrackLength;

	// CalcIterativeTruncMean(std::vector<float> v, const size_t& nmin, const size_t& nmax, const size_t& currentiteration, const size_t lmin, const float& convergencelimit, const float& nsigma)
	// Setting nmin=1, nmax=1, currentiteration=0 does one iteration only
	// Setting lmin=1 and convergencelimit=0.1 don't really affect it if only doing one iteration
	// Setting nsigma=1 truncates at +/- 1 sigma (as done by Marco and Libo)
        size_t nmin = 1;
        size_t nmax = 1;
        const size_t currentiteration = 0;
        const size_t lmin = 1;
        const float convergencelimit = 0.1;
        const float nsigma = 1.0;
									   dQdxtruncmean.fValue = (double)trm.CalcIterativeTruncMean(dQdx, nmin, nmax, currentiteration, lmin, convergencelimit, nsigma);
        dEdxtruncmean.fValue = (double)trm.CalcIterativeTruncMean(dEdx, nmin, nmax, currentiteration, lmin, convergencelimit, nsigma);
        trklen.fValue = track->Length();

        AlgScoresVec.push_back(dQdxtruncmean);
        AlgScoresVec.push_back(dEdxtruncmean);
        AlgScoresVec.push_back(trklen);
	*/
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


        std::cout << "[ParticleID]  >> Making assn... " << std::endl;
        util::CreateAssn(*this, e, *particleIDCollection, track, *trackParticleIdAssn);

  }
      
  e.put(std::move(particleIDCollection));
  e.put(std::move(trackParticleIdAssn));
  
}


void UBPID::ParticleId::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(UBPID::ParticleId)
