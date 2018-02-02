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
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// local includes
#include "uboone/ParticleID/Algorithms/GetDaughterTracksShowers.h"
#include "uboone/ParticleID/Algorithms/fiducialVolume.h"
#include "uboone/ParticleID/Algorithms/PIDA.h"

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
    double fCutDistance;
    double fCutFraction;

    // fidvol related
    fidvol::fiducialVolume fid;
    particleid::PIDA pida;

    //other
    bool isData;
};


UBPID::ParticleId::ParticleId(fhicl::ParameterSet const & p)
{

  // fcl parameters
  fTrackingAlgo = p.get< std::string > ("TrackingAlgorithm");
  fCutDistance  = p.get< double > ("DaughterFinderCutDistance");
  fCutFraction  = p.get< double > ("DaughterFinderCutFraction");

  fv = fid.setFiducialVolume(fv, p);
  fid.printFiducialVolume(fv);

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

  std::vector< anab::ParticleID > pidVector;

  for (auto& track : trackCollection){

    TVector3 trackStart = track->Vertex();
    TVector3 trackEnd = track->End();

    int nDaughters = GetNDaughterTracks((*trackHandle), track->ID(), fCutDistance, fCutFraction);

    std::cout << "[ParticleID]  Found track with " << nDaughters << " reconstructed daughters." << std::endl;

    int pdg = 0;
    int ndf = 0;
    double minchi2 = 0.0;
    double deltachi2 = 0.0;
    double chi2proton = 0.0;
    double chi2kaon = 0.0;
    double chi2pion = 0.0;
    double chi2muon = 0.0;
    double missingE = 0.0;
    double missingEAvg = 0.0;
    double pidaVal = -1.0;

    // if track is fully contained and is not reinteracting then fill PIDA
    if (nDaughters == 0 && fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){

      std::cout << "[ParticleID]  >> Track is fully contained and has no daughters " << std::endl;

      pidaVal = pida.getPida();

    }

    geo::PlaneID planeid(0,0,0);

    particleIDCollection->push_back(anab::ParticleID(
          pdg,
          ndf,
          minchi2,
          deltachi2,
          chi2proton,
          chi2kaon,
          chi2pion,
          chi2muon,
          missingE,
          missingEAvg,
          pidaVal,
          planeid));


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
