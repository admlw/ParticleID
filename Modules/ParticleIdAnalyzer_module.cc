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
#include "art/Framework/Services/Optional/TFileService.h"

// data products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// UBXSec includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"

// local includes
#include "uboone/ParticleID/Algorithms/GetDaughterTracksShowers.h"
#include "uboone/ParticleID/Algorithms/fiducialVolume.h"
#include "uboone/ParticleID/Algorithms/PIDA.h"

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"

#include <iostream>

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

    std::vector<double> fv;

  private:

    art::ServiceHandle< art::TFileService > tfs; 

    std::string fTrackingAlgo;
    std::string fCaloLabel;

    double fCutDistance;
    double fCutFraction;

    // fidvol related
    fidvol::fiducialVolume fid;
    particleid::PIDA pida;

    // ubxsec related
    ubana::SelectionResult selRes;

    // cc1munp related
    int cc1munpRun;
    int cc1munpSubRun;
    int cc1munpEvent;
    int cc1munpMuonID;
    int cc1munpProtonID;
    TTree* ttree;

    int nProtons = 0;

    TH1D* hMuonPreCutMean;
    TH1D* hMuonPreCutMedian;
    TH1D* hMuonPreCutKde;
    TH1D* hProtonPreCutMean;
    TH1D* hProtonPreCutMedian;
    TH1D* hProtonPreCutKde;

    TH1D* hPreCutMean;
    TH1D* hPreCutMedian;
    TH1D* hPreCutKde;

    TH1D* hMuonPostCutMean;
    TH1D* hMuonPostCutMedian;
    TH1D* hMuonPostCutKde;
    TH1D* hProtonPostCutMean;
    TH1D* hProtonPostCutMedian;
    TH1D* hProtonPostCutKde;

    TH1D* hPostCutMean;
    TH1D* hPostCutMedian;
    TH1D* hPostCutKde;

    TH2D* hProtonStartYZ;
    TH2D* hMuonStartYZ;

    std::vector<TH1D*> protonPIDAVals;

};


ParticleIdAnalyzer::ParticleIdAnalyzer(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fTrackingAlgo = p.get< std::string > ("TrackingAlgorithm");
  fCaloLabel = p.get< std::string > ("CalorimetryModule");
  fCutDistance = p.get< double > ("DaughterFinderCutDistance");
  fCutFraction = p.get< double > ("DaughterFinderCutFraction");

  fv = fid.setFiducialVolume(fv, p);
  fid.printFiducialVolume(fv);

}

void ParticleIdAnalyzer::beginJob()
{

  std::cout << "Opening file..." << std::endl;
  TFile *file = new TFile("/uboone/app/users/alister1/particleID/ubcode_v06_26_01_10/srcs/uboonecode/uboone/ParticleID/cc1unp_output.root", "READ");
  std::cout << "...Done." << std::endl;

  if (file->IsZombie())
    std::cout << "File is zombie!" << std::endl;

  ttree = (TTree*)file->Get("cc1unpselana/fMC_TrunMean");

  ttree->SetBranchAddress("fRun", &cc1munpRun);
  ttree->SetBranchAddress("fSubRun", &cc1munpSubRun);
  ttree->SetBranchAddress("fEvent", &cc1munpEvent);
  ttree->SetBranchAddress("trkmuoncandid", &cc1munpMuonID);
  ttree->SetBranchAddress("trkprotoncandid", &cc1munpProtonID);

  hMuonPreCutMean = tfs->make<TH1D>("hMuonPreCutMean", ";;", 100, 0, 30);
  hMuonPreCutMedian = tfs->make<TH1D>("hMuonPreCutMedian", ";;", 100, 0, 30);
  hMuonPreCutKde = tfs->make<TH1D>("hMuonPreCutKde", ";;", 100, 0, 30);
  hProtonPreCutMean = tfs->make<TH1D>("hProtonPreCutMean", ";;", 100, 0, 30);
  hProtonPreCutMedian = tfs->make<TH1D>("hProtonPreCutMedian", ";;", 100, 0, 30);
  hProtonPreCutKde = tfs->make<TH1D>("hProtonPreCutKde", ";;", 100, 0, 30);

  hPreCutMean = tfs->make<TH1D>("hPreCutMean", ";;", 100, 0, 30);
  hPreCutMedian = tfs->make<TH1D>("hPreCutMedian", ";;", 100, 0, 30);
  hPreCutKde = tfs->make<TH1D>("hPreCutKde", ";;", 100, 0, 30);

  hMuonPostCutMean = tfs->make<TH1D>("hMuonPostCutMean", ";;", 100, 0, 30);
  hMuonPostCutMedian = tfs->make<TH1D>("hMuonPostCutMedian", ";;", 100, 0, 30);
  hMuonPostCutKde = tfs->make<TH1D>("hMuonPostCutKde", ";;", 100, 0, 30);
  hProtonPostCutMean = tfs->make<TH1D>("hProtonPostCutMean", ";;", 100, 0, 30);
  hProtonPostCutMedian = tfs->make<TH1D>("hProtonPostCutMedian", ";;", 100, 0, 30);
  hProtonPostCutKde = tfs->make<TH1D>("hProtonPostCutKde", ";;", 100, 0, 30);

  hPostCutMean = tfs->make<TH1D>("hPostCutMean", ";;", 100, 0, 30);
  hPostCutMedian = tfs->make<TH1D>("hPostCutMedian", ";;", 100, 0, 30);
  hPostCutKde = tfs->make<TH1D>("hPostCutKde", ";;", 100, 0, 30);

  hProtonStartYZ = tfs->make<TH2D>("hProtonStartYZ", ";;", 50, 0, 1036, 50, -116.5, 116.5);
  hMuonStartYZ = tfs->make<TH2D>("hMuonStartYZ", ";;", 50, 0, 1036, 50, -116.5, 116.5);

  for (int i = 0; i < 50; i ++){

    TString th1name = Form("protonPIDAVals_%i", i);
    protonPIDAVals.push_back(tfs->make<TH1D>(th1name, ";PIDA vals;", 1000, 0, 50));

  }

  std::cout << "Histograms created." << std::endl;
}

void ParticleIdAnalyzer::analyze(art::Event const & e)
{

  //bool isData = e.isRealData();
  bool isSelected = false;

  int fRun = e.run();
  int fSubRun = e.subRun();
  int fEvent = e.event();

  std::cout << "----- " << fRun << "." << fSubRun << "." << fEvent << std::endl;

  std::vector<int> muonIds;
  std::vector<int> protonIds;

  for (int i = 0; i < ttree->GetEntries(); i++){

    ttree->GetEntry(i);


    if (cc1munpRun == fRun && cc1munpSubRun == fSubRun && cc1munpEvent == fEvent){

      std::cout << "Found an event!" << std::endl;
      isSelected = true;
      muonIds.push_back(cc1munpMuonID);
      protonIds.push_back(cc1munpProtonID);
      continue;
    }

  }

  if (isSelected == false) return;
  //
  // get handles to needed information
  //

  // selection result...
  art::Handle< std::vector< recob::Track > > trackHandle;
  e.getByLabel(fTrackingAlgo, trackHandle);
  std::vector< art::Ptr< recob::Track > > trackPtrs;
  art::fill_ptr_vector(trackPtrs, trackHandle);

  art::FindManyP< anab::Calorimetry > caloFromTracks(trackHandle, e, "pandoraNucalo");

  for (size_t j = 0; j < trackPtrs.size(); j++){

    art::Ptr< recob::Track > track     = trackPtrs.at(j);
    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = caloFromTracks.at(track->ID());
    art::Ptr< anab::Calorimetry > calo = caloFromTrack.at(0);
    std::cout << "USING PLANE: " << calo->PlaneID() << std::endl;

    std::vector< double > dEdx         = calo->dEdx();
    std::vector< double > resRange     = calo->ResidualRange();

    int trackID          = track->ID();
    int nDaughters       = GetNDaughterTracks((*trackHandle), trackID, fCutDistance, fCutFraction);

    TVector3 trackStart  = track->Vertex();
    TVector3 trackEnd    = track->End();

    double pidaValMean   = pida.getPida(dEdx, resRange, "mean");
    double pidaValMedian = pida.getPida(dEdx, resRange, "median");
    double pidaValKde    = pida.getPida(dEdx, resRange, "kde");

    for (size_t i = 0; i < muonIds.size(); i++){

      if (trackID == muonIds.at(i)){

        std::cout << ">> Found Candidate Muon!" << std::endl;

        hMuonPreCutMean->Fill(pidaValMean);
        hMuonPreCutMedian->Fill(pidaValMedian);
        hMuonPreCutKde->Fill(pidaValKde);

        hPreCutMean->Fill(pidaValMean);
        hPreCutMedian->Fill(pidaValMean);
        hPreCutKde->Fill(pidaValMean);

        hMuonStartYZ->Fill(trackStart.Z(), trackStart.Y());

        if (nDaughters == 0 && fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){

          hMuonPostCutMean->Fill(pidaValMean);
          hMuonPostCutMedian->Fill(pidaValMedian);
          hMuonPostCutKde->Fill(pidaValKde);

          hPostCutMean->Fill(pidaValMean);
          hPostCutMedian->Fill(pidaValMean);
          hPostCutKde->Fill(pidaValMean);

        }

      }

    }

    for (size_t i = 0; i < protonIds.size(); i++){

      if (trackID == protonIds.at(i)){

        std::cout << ">> Found Candidate Proton!" << std::endl;

        hProtonPreCutMean->Fill(pidaValMean);
        hProtonPreCutMedian->Fill(pidaValMedian);
        hProtonPreCutKde->Fill(pidaValKde);

        hPreCutMean->Fill(pidaValMean);
        hPreCutMedian->Fill(pidaValMean);
        hPreCutKde->Fill(pidaValMean);

        hProtonStartYZ->Fill(trackStart.Z(), trackStart.Y());

        if (nDaughters == 0 && fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){

          hProtonPostCutMean->Fill(pidaValMean);
          hProtonPostCutMedian->Fill(pidaValMedian);
          hProtonPostCutKde->Fill(pidaValKde);

          hPostCutMean->Fill(pidaValMean);
          hPostCutMedian->Fill(pidaValMean);
          hPostCutKde->Fill(pidaValMean);

          if (nProtons < 50){
            TString th1name = Form("run%i_event%i_trkid%i", fRun, fEvent, trackID);
            protonPIDAVals.at(nProtons)->SetName(th1name);
            for (size_t j = 0; j < dEdx.size(); j++){

              std::cout << "dEdx "  << j << ": " << dEdx.at(j) << " ResRg: " << resRange.at(j) <<  " particleId: " << dEdx.at(j)*std::pow(resRange.at(j), 0.42) << std::endl;
              protonPIDAVals.at(nProtons)->Fill(dEdx.at(j)*std::pow(resRange.at(j),0.42));
            }

            nProtons++;
          }
        }

      }

    }


    /*
       art::Handle< std::vector< ubana::SelectionResult> > selectionHandle;
       e.getByLabel("UBXSec", selectionHandle);

       if (!isData) std::cout << ">> Running Simulated Data" << std::endl;

       if (!selectionHandle.isValid()){

       std::cout << " >> Selection product not found. " << std::endl;
       mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found" 
       << std::endl;
       o
       throw std::exception();
       }

       art::FindManyP< ubana::TPCObject > selectedTpcObjects(selectionHandle, e, "UBXSec");
       art::FindManyP< anab::ParticleID > particleIdFromTracks(trackHandle, e, "pid");

    // def vars for loop
    art::Ptr< ubana::TPCObject > selectedTpcObject;

    if (selectedTpcObjects.at(0).size() == 1){

    selectedTpcObject = selectedTpcObjects.at(0).at(0);

    const std::vector< recob::Track >& selectedTracks = selectedTpcObject->GetTracks();

    for (size_t i = 0; i < selectedTracks.size(); i++){

    const recob::Track track = selectedTracks.at(i);
    art::Ptr< recob::Track > trackP = lar::util::make_ptr(track);

    std::vector< art::Ptr<anab::ParticleID> > pid = particleIdFromTracks.at(track.key());

    std::cout << pid.at(0)->PIDA() << std::endl;

    }

    }
    */
  }
}

void ParticleIdAnalyzer::endJob()
{
  /*
     hMuonPreCutMean->Write();
     hMuonPreCutMedian->Write();
     hMuonPreCutKde->Write();
     hProtonPreCutMean->Write();
     hProtonPreCutMedian->Write();
     hProtonPreCutKde->Write();
     hPreCutMean->Write();
     hPreCutMedian->Write();
     hPreCutKde->Write();
     hMuonPostCutMean->Write();
     hMuonPostCutMedian->Write();
     hMuonPostCutKde->Write();
     hProtonPostCutMean->Write();
     hProtonPostCutMedian->Write();
     hProtonPostCutKde->Write();
     hPostCutMean->Write();
     hPostCutMedian->Write();
     hPostCutKde->Write();
     */
}

DEFINE_ART_MODULE(ParticleIdAnalyzer)
