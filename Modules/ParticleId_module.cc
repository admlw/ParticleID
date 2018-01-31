////////////////////////////////////////////////////////////////////////
// Class:       ParticleId
// Plugin Type: producer (art v2_05_01)
// File:        ParticleId_module.cc
//
// Generated at Wed Jan 31 11:25:52 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

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

private:

  std::string fTrackingAlgo; 

};


ParticleId::ParticleId(fhicl::ParameterSet const & p)
{

  // fcl parameters
  fTrackingAlgo = p.get< std::string > ("TrackingAlgorithm");

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

  // get handles to needed information
  art::Handle< std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackingAlgo, trackHandle);

  // produce collection of particleID objects
  std::unique_ptr< std::vector<anab::ParticleID> > particleIDCollection( new std::vector<anab::ParticleID> );

  for (auto const& track : (*trackHandle)){

    track.Length();

  }
  
}

void ParticleId::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(ParticleId)
