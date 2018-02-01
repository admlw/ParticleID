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

// data products
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/Track.h"

// UBXSec includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"

// local includes
#include "uboone/ParticleID/Algorithms/GetDaughterTracksShowers.h"
#include "uboone/ParticleID/Algorithms/fiducialVolume.h"

// cpp includes
#include <memory>

class ParticleId;


class ParticleId : public art::EDProducer {
  public:
    explicit ParticleId(fhicl::ParameterSet const & p);

    ParticleId(ParticleId const &) = delete;
    ParticleId(ParticleId &&) = delete;
    ParticleId & operator = (ParticleId const &) = delete;
    ParticleId & operator = (ParticleId &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;
    
    double xl, xh, yl, yh, zl, zh;
  
  private:

    // fcl
    std::string fTrackingAlgo; 
    double fCutDistance;
    double fCutFraction;

    // fidvol related
    fidvol::fiducialVolume fid;

    // ubxsec related
    ubana::SelectionResult selRes;

    //other
    bool isData;
};


ParticleId::ParticleId(fhicl::ParameterSet const & p)
{


  
  // fcl parameters
  fTrackingAlgo = p.get< std::string > ("TrackingAlgorithm");
  fCutDistance  = p.get< double > ("DaughterFinderCutDistance");
  fCutFraction  = p.get< double > ("DaughterFinderCutFraction");

  fid.setFiducialVolume(xl, xh, yl, yh, zl, zh, p);
  fid.printFiducialVolume(xl, xh, yl, yh, zl, zh);

  // this module produces a anab::ParticleID object and
  // an association to the track which produced it
  produces< std::vector<anab::ParticleID> >();

}

void ParticleId::beginJob()
{
  // Implementation of optional member function here.
}

void ParticleId::produce(art::Event & e)
{

  bool isData = e.isRealData();

  selRes.SetSelectionStatus(false);
  selRes.SetFailureReason("NoUBXSec");

  // produce collection of particleID objects
  std::unique_ptr< std::vector<anab::ParticleID> > particleIDCollection( new std::vector<anab::ParticleID> );

  //
  // get handles to needed information
  //
  
  // selection result...
  art::Handle< std::vector< ubana::SelectionResult> > selectionHandle;
  e.getByLabel("UBXSec", selectionHandle);

  if (!isData) std::cout << ">> Running Simulated Data" << std::endl;

  if (!selectionHandle.isValid()){

    std::cout << " >> Selection product not found. " << std::endl;
    mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found" 
      << std::endl;

    throw std::exception();
  }

  art::FindManyP< ubana::TPCObject > selectedTpcObjects(selectionHandle, e, "UBXSec");

  // and tracks...
  //art::Handle < std::vector<recob::Track> > trackHandle;
  //e.getByLabel(fTrackingAlgo, trackHandle);

  // def vars for loop
  art::Ptr< ubana::TPCObject > selectedTpcObject;
  
  if (selectedTpcObjects.at(0).size() == 1){

    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const std::vector< recob::Track >& selectedTracks = selectedTpcObject->GetTracks();

    for (size_t i = 0; i < selectedTracks.size(); i++){

      const recob::Track& track = selectedTracks.at(i);

      int nDaughters = GetNDaughterTracks(selectedTracks, track.ID(), fCutDistance, fCutFraction);

      std::cout << nDaughters << std::endl;

    }
  
  }

}

void ParticleId::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ParticleId)
