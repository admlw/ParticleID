////////////////////////////////////////////////////////////////////////
// Class:       ParticleIdAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        ParticleIdAnalyzer_module.cc
//
// Generated at Fri Feb  2 10:56:38 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"

// data products
#include "lardataobj/RecoBase/Track.h"

// UBXSec includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"


class ParticleIdAnalyzer;


class ParticleIdAnalyzer : public art::EDAnalyzer {
  public:
    explicit ParticleIdAnalyzer(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ParticleIdAnalyzer(ParticleIdAnalyzer const &) = delete;
    ParticleIdAnalyzer(ParticleIdAnalyzer &&) = delete;
    ParticleIdAnalyzer & operator = (ParticleIdAnalyzer const &) = delete;
    ParticleIdAnalyzer & operator = (ParticleIdAnalyzer &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

  private:

   std::string fTrackingAlgo;

   // ubxsec related
    ubana::SelectionResult selRes;

};


ParticleIdAnalyzer::ParticleIdAnalyzer(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fTrackingAlgo = p.get< std::string > ("TrackingAlgorithm");

}

void ParticleIdAnalyzer::beginJob()
{
  // Implementation of optional member function here.
}

void ParticleIdAnalyzer::analyze(art::Event const & e)
{

  bool isData = e.isRealData();

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
      TVector3 trackStart = track.Vertex();
      TVector3 trackEnd = track.End();

      std::cout << trackStart.X() << " " << trackEnd.X() << std::endl;


    }
  
  }

}

void ParticleIdAnalyzer::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ParticleIdAnalyzer)
