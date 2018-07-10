////////////////////////////////////////////////////////////////////////
// Class:       UBPID::CalibrateSimulatedData
// Plugin Type: producer (art v2_05_01)
// File:        UBPID::CalibrateSimulatedData_module.cc
//
// Generated at Thu Jun  7 10:24:02 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 *
 * \class ParticleId
 *
 * \brief Calibration producer module
 *
 * \author Kirsty Duffy (kduffy@fnal.gov), Adam Lister (a.lister1@lancaster.ac.uk)
 *
 * \date 2018/04/18
 *
 */

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

// local includes
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"

#include "TRandom3.h"

#include <memory>

namespace UBPID{
  class CalibrateSimulatedData;
}

class UBPID::CalibrateSimulatedData : public art::EDProducer {
  public:
    explicit CalibrateSimulatedData(fhicl::ParameterSet const & p);

    // Plugins should not be copied or assigned.
    CalibrateSimulatedData(UBPID::CalibrateSimulatedData const &) = delete;
    CalibrateSimulatedData(UBPID::CalibrateSimulatedData &&) = delete;
    CalibrateSimulatedData & operator = (UBPID::CalibrateSimulatedData const &) = delete;
    CalibrateSimulatedData & operator = (UBPID::CalibrateSimulatedData &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

  private:

    // Declare member data here.
    std::string fCaloLabel;
    std::string fTrackLabel;
    bool fIsSimSmear;
    bool fIsSetSeed;
    int fRandomSeed;

    std::vector<double> fSimGausSmearWidth;

};


UBPID::CalibrateSimulatedData::CalibrateSimulatedData(fhicl::ParameterSet const & p)
  // :
  // Initialize member data here.
{

  fhicl::ParameterSet p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");

  fCaloLabel = p_labels.get< std::string > ("CalorimetryLabel");
  fTrackLabel = p_labels.get< std::string > ("TrackLabel");
  fIsSimSmear = p.get< bool > ("IsSimulationSmearing");
  fSimGausSmearWidth = p.get< std::vector<double> > ("SimulationGausSmearWidth");
  fIsSetSeed = p.get< bool > ("IsSetSeed");
  fRandomSeed = p.get< int > ("RandomSeed");

  produces< std::vector<anab::Calorimetry> >();
  produces< art::Assns< recob::Track, anab::Calorimetry> >();

  std::cout << "[CalibrateSimulatedData] >> Track Label: " << fTrackLabel << std::endl;
  std::cout << "[CalibrateSimulatedData] >> Calo Label: " << fCaloLabel << std::endl;
  std::cout << "[CalibrateSimulatedData] The following is only useful if running simulated data" << std::endl;
  std::cout << "[CalibrateSimulatedData] >> Do simulation smearing? " << fIsSimSmear << std::endl;
  std::cout << "[CalibrateSimulatedData] >> Smearing simulated data by [" 
    << fSimGausSmearWidth.at(0) << ", " << fSimGausSmearWidth.at(1) 
    << ", " << fSimGausSmearWidth.at(2) << "]" << std::endl;

}

void UBPID::CalibrateSimulatedData::produce(art::Event & e)
{

  bool isData = e.isRealData();

  // get seed
  int run = e.run();
  int sub_run = e.subRun();
  int event = e.event();

  // setup random 
  TRandom3 r;
  
  if (fIsSetSeed){
    r.SetSeed(fRandomSeed);
  }
  else{
    std::string randomSeedString =  Form("%i%i%i", event, sub_run, run);
    r.SetSeed(std::atoi(randomSeedString.c_str()));
  }
  
  std::cout << "[CalibrateSimulatedData] Run " << run << " SubRun " << sub_run << " Event " << event << std::endl;
  std::cout << "[CalibrateSimulatedData] The random number seed is " << r.GetSeed() << std::endl;

  std::unique_ptr< std::vector<anab::Calorimetry> > calorimetryCollection( new std::vector<anab::Calorimetry> );
  std::unique_ptr< art::Assns <recob::Track, anab::Calorimetry> > trackCalorimetryAssn( new art::Assns<recob::Track, anab::Calorimetry> );

  // tracks...
  art::Handle < std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector< art::Ptr<recob::Track> > trackCollection;
  art::fill_ptr_vector(trackCollection, trackHandle);

  // calorimetry object...
  art::FindManyP<anab::Calorimetry> caloFromTracks(trackHandle, e, fCaloLabel);

  for (auto& track : trackCollection){

    if (!caloFromTracks.isValid()){
      std::cout << "[CalibrateSimulatedData] Calorimetry<->Track associations are not valid for this track. Skipping." << std::endl;
      continue;
    }

    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = caloFromTracks.at(track->ID());

    art::Ptr< anab:: Calorimetry > calo;
    int planenum;
    for (auto c : caloFromTrack){
      planenum = c->PlaneID().Plane;
      calo = c;

      if (!calo || planenum < 0 || planenum > 2){
        std::cout << "[CalibrateSimulatedData] Calorimetry on plane " << planenum << " is unavailable. Skipping." << std::endl;
        continue;
      }

      std::vector<double> dEdx = calo->dEdx();

      /**
       * If this is simulated data, then we want to smear the dE/dx values
       * by a gaussian with a width defined in fcl
       */

      if (isData == false && fIsSimSmear == true){

        for (size_t i = 0; i < dEdx.size(); i++){

          double simulationSmear = r.Gaus(1., fSimGausSmearWidth.at(planenum)); 
          dEdx.at(i) = dEdx.at(i) * simulationSmear;

        }

      }

      /** build calorimetry object */
      anab::Calorimetry caloObj(
          calo->KineticEnergy(),
          dEdx,
          calo->dQdx(),
          calo->ResidualRange(),
          calo->DeadWireResRC(),
          calo->Range(),
          calo->TrkPitchVec(),
          calo->PlaneID()
          );

      calorimetryCollection->push_back(caloObj);

      util::CreateAssn(*this, e, *calorimetryCollection, track, *trackCalorimetryAssn);


    }


  }

  e.put(std::move(calorimetryCollection));
  e.put(std::move(trackCalorimetryAssn));
}

DEFINE_ART_MODULE(UBPID::CalibrateSimulatedData)
