/*************************************************************
 * 
 * This is a quick macro for counting how many daughter tracks a
 * reconstructed track has, using the number of tracks with start
 * or end points within a given distance of the end point of the
 * first track.
 *
 * Don't forget to set up gallery first:
 * setup gallery v1_03_08 -q e10:nu:prof
 *
 * root -b
 * root [0] .L Test_countdaughtertracks.C++
 * root [1] MakePlots("path/to/reco2/file.root or path/to/list/of/input/files.txt")
 *
 * Kirsty Duffy (kduffy@fnal.gov), Fermilab, Jan 28 2018
 * 
 *************************************************************/


//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// Include GetNDaughterTracks/Showers functions
#include "../Algorithms/GetDaughterTracksShowers.h"
#include "../Algorithms/GetDaughterTracksShowers.cxx"

// Function to get root file directly *or* get root file names from a text file
#include <boost/algorithm/string/predicate.hpp>
std::vector<std::string> GetFileList(std::string input_file)
{
  std::vector<std::string> filenames;
  if(boost::algorithm::ends_with(input_file,".root")){
    filenames.emplace_back(input_file);
    std::cout << "Added file " << filenames.back() << " to be processed." << std::endl;
  }
  else{
    std::ifstream input(input_file);
    for(std::string line; std::getline(input,line);){
      filenames.emplace_back(line);
      std::cout << "Added file " << filenames.back() << " to be processed." << std::endl;
    }
  }
  return filenames;
}




// ---------------------- Main Function ------------------------ //

void CountDaughterTracks(gallery::Event *ev, double cutdist, double cutfrac, std::vector<std::string> trktags, TH2F &h_reco, TH2F &h_eff, TH2F &h_effnorm)
{
  
    // Loop through track producers
    for (std::string trktag_str : trktags){
      //std::cout << "--- " << trktag_str << " ---" << std::endl;
      art::InputTag const trktag(trktag_str);
      
      // Get tracks
      auto const& trk_handle = ev->getValidHandle<std::vector<recob::Track>>(trktag);
      //art::Handle<std::vector<recob::Track>> trk_handle;
      //ev.getByLabel(trktag,trk_handle);
      const auto& trackVec(*trk_handle);
      //std::cout << trk_handle->size() << " tracks in event" << std::endl;

      
      // Loop through tracks
      // For each one, look at other tracks in the event and try to work out if it
      // has "daughters"
      for (auto const& track : trackVec){
	
	// Does it have daughters?
	int ndaughters = GetNDaughterTracks(trackVec,track.ID(),cutdist,cutfrac);
	if (trktag_str == "pandoraNu" || trktag_str == "pandoraCosmic"){
	  // Include showers for pandoraNu and pandoraCosmic
	  auto const& shwr_handle = ev->getValidHandle<std::vector<recob::Shower>>(trktag);
	  ndaughters += GetNDaughterShowers(trackVec, track.ID(), (*shwr_handle), cutdist, cutfrac);
	}
	//std::cout << "  Track " << track.ID() << ": " << ndaughters << " daughters within 10 " << std::endl;
	h_reco.Fill(cutdist,cutfrac,ndaughters);

	
	// If using pandora, check PFPs to see if this looks right
	if (trktag_str == "pandoraNu" || trktag_str == "pandoraCosmic"){
	  art::FindManyP<recob::PFParticle> PFP_from_tracks(trk_handle, (*ev), trktag);
	  auto pfp = PFP_from_tracks.at(track.ID());
	  //if (pfp.size() > 1) std::cout << "[WARNING] pfp size = " << pfp.size() << std::endl;
	  //std::cout << "           PFP daughters: " << pfp.at(0)->NumDaughters() << std::endl;

	  //if (pfp.at(0)->NumDaughters() != ndaughters) std::cout << "PROBLEM: PFP daughters = " << pfp.at(0)->NumDaughters() << ", track/shower daughters = " << ndaughters << std::endl;

	  // Add to bin content in efficiency histogram
	  double efficiency;
	  //if (ndaughters == 0 && pfp.at(0)->NumDaughters() == 0){ efficiency = 0.;}
	  //else if (pfp.at(0)->NumDaughters == 0){ efficiency = ndaughters;}
	  //else { efficiency = double(ndaughters)/double(pfp.at(0)->NumDaughters());}
	  efficiency = ndaughters - pfp.at(0)->NumDaughters();

	  //if (efficiency != 0.) std::cout << "cutdist = " << cutdist << ", cutfrac = " << cutfrac << ", NDaughters = " << ndaughters << ", nPFPdaughters = " << pfp.at(0)->NumDaughters() << ", efficiency = " << efficiency << std::endl;
	  
	  h_eff.Fill(cutdist,cutfrac,efficiency);
	  h_effnorm.Fill(cutdist,cutfrac); // To normalise efficiency histogram per bin

	} // end if (pandoraNu || pandoraCosmic)
	
      } // end loop over tracks

    } // end loop over track producers
  
}


// ---------------------- This function makes some plots ------------------------ //

void MakePlots(std::string input_files)
{
  // Tracking algorithms
  std::vector<std::string> trktags;
  trktags.push_back("pandoraNu");
  trktags.push_back("pandoraCosmic");
  //trktags.push_back("pandoraNuKalmanTrack");
  //trktags.push_back("pandoraNuKalmanShower");
  //trktags.push_back("pandoraCosmicKalmanTrack");
  trktags.push_back("pmtrack");
  trktags.push_back("pandoraNuPMA");
  trktags.push_back("pandoraNuKHit");
  trktags.push_back("pandoraCosmicKHit");

  // Output file
  TFile *fout = new TFile("out_DaughterTracks.root","recreate");

  
  gStyle->SetOptStat(0);

  // Format files list
  std::vector<std::string> filenames = GetFileList(input_files);

  // Make histograms
  const int n_producers = trktags.size();
  TH2F *hists_reco[n_producers];
  TH2F *hists_effnorm[n_producers];
  TH2F *hists_eff[n_producers];
  for (int i_producer=0; i_producer<n_producers; i_producer++){
    hists_reco[i_producer] = new TH2F(std::string("hreco_"+trktags.at(i_producer)).c_str(),";Cut distance (cm);Fraction of track length (for maximum cut)",20,0,20,20,0,1);
    hists_effnorm[i_producer]  = (TH2F*)hists_reco[i_producer]->Clone(std::string("heffnorm_"+trktags.at(i_producer)).c_str());
    hists_eff[i_producer]  = (TH2F*)hists_reco[i_producer]->Clone(std::string("heff_"+trktags.at(i_producer)).c_str());
  }
  
  // Loop through events
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){

    std::cout << "Processing " 
	      << "Run " << ev.eventAuxiliary().run() << ", "
	      << "Event " << ev.eventAuxiliary().event() << std::endl;

    // Loop through track producers and make plots
    for (int i_producer=0; i_producer<n_producers; i_producer++){
      
      std::vector<std::string> trktags_thisproducer;
      trktags_thisproducer.push_back(trktags.at(i_producer));
      
      for (int bin_cutdist=1; bin_cutdist<hists_reco[i_producer]->GetXaxis()->GetNbins()+1; bin_cutdist++){
	for (int bin_cutfrac=1; bin_cutfrac<hists_reco[i_producer]->GetYaxis()->GetNbins()+1; bin_cutfrac++){
	  
	  double cutdist = hists_reco[i_producer]->GetXaxis()->GetBinCenter(bin_cutdist);
	  double cutfrac = hists_reco[i_producer]->GetYaxis()->GetBinCenter(bin_cutfrac);
	  
	  // Fill histograms
	  CountDaughterTracks(&ev, cutdist, cutfrac, trktags_thisproducer, *hists_reco[i_producer], *hists_eff[i_producer], *hists_effnorm[i_producer]);
	  
	}// bin_cutfrac
      }// bin_cutdist
      
    }// i_producer

  } // loop through gallery events
      
  // Loop through track producers again and write histograms to file and to pdfs
  fout->cd();
  TCanvas *c1 = new TCanvas();
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  
  // Do all the reco ones first so they have the same colour scheme
  for (int i_producer=0; i_producer<n_producers; i_producer++){
    hists_reco[i_producer]->Draw("colz");
    c1->Print(std::string(std::string(hists_reco[i_producer]->GetName())+".pdf").c_str());
    hists_reco[i_producer]->Write();
  }

    // Now change colour scheme for efficiency histograms
    const Int_t NRGBs = 3;
    const Int_t NCont = 30;
    Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 1.00, 1.00 };
    Double_t green[NRGBs] = { 0.00, 1.00, 0.00 };
    Double_t blue[NRGBs]  = { 1.00, 1.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs,stops,red,green,blue,NCont);
    gStyle->SetNumberContours(NCont+1);
    
  // Now do efficiency hists (pandoraNu and pandoraCosmic only)
  for (int i_producer=0; i_producer<n_producers; i_producer++){
    if (!(trktags.at(i_producer) == "pandoraNu" || trktags.at(i_producer) == "pandoraCosmic")){ continue; }

    // Normalise in each bin by dividing by h_effnorm (was filled the same number of times but with
    // weight 1 -- the bin content will just be the number of tracks)
    hists_eff[i_producer]->Divide(hists_effnorm[i_producer]);
    
    hists_eff[i_producer]->SetTitle("No. reco daughters - no. PFP daughters (per track, averaged over all tracks)");
    hists_eff[i_producer]->GetZaxis()->SetRangeUser(-1,1.);

    hists_eff[i_producer]->Draw("colz");
    c1->Print(std::string(std::string(hists_eff[i_producer]->GetName())+".pdf").c_str());
    hists_eff[i_producer]->Write();
    hists_effnorm[i_producer]->Write();

    std::cout << "hists_reco[" << trktags.at(i_producer) << "] NEntries = " << hists_reco[i_producer]->GetEntries() << std::endl;
    std::cout << "hists_effnorm[" << trktags.at(i_producer) << "] NEntries = " << hists_effnorm[i_producer]->GetEntries() << std::endl;
    std::cout << "hists_eff[" << trktags.at(i_producer) << "] NEntries = " << hists_eff[i_producer]->GetEntries() << std::endl;
  }
}
