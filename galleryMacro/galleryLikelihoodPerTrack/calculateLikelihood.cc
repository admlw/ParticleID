/*****************************************
* calculate likelihood for each track compared to predicted likelihood
****************************************/

// c++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <exception>

// root includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TString.h"
#include "TRandom3.h"

// art includes
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTrajectory.h"
#include "/uboone/app/users/alister1/particleID/ubcode_v06_26_01_10/srcs/uboonecode/uboone/ParticleID/Algorithms/Theory_dEdx_resrange.cxx"
//#include "/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/ParticleID/Algorithms/Bragg_negLogL_Estimator.cxx"
#include "Bragg_negLogL_Estimator_new.cxx"

//const simb::MCParticle* getMatchedParticle(art::FindManyP<recob::Hit> trackHitAssn, art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit, art::InputTag hitMatcherTag, size_t it){

/* get MCParticle from a given track... */
/*
std::vector< art::Ptr<recob::Hit> > trkHitPtrs = trackHitAssn.at(it);
std::unordered_map<int, double> trkide;
double maxe=-1, tote=0;
const simb::MCParticle* mcp;

std::vector<simb::MCParticle const*> particleVec;
std::vector<anab::BackTrackerHitMatchingData const*> matchVec;

// loop hits
for (size_t ih = 0; ih < trkHitPtrs.size(); ih++){

particleVec.clear();
matchVec.clear();
particlesPerHit.get(trkHitPtrs[ih].key(), particleVec, matchVec);

// loop over particles
for (size_t ip = 0; ip < particleVec.size(); ip++){

trkide[particleVec[ip]->TrackId()]+=matchVec[ip]->energy;
tote += matchVec[ip]->energy;
if (trkide[particleVec[ip]->TrackId()] > maxe){
maxe = trkide[ particleVec[ip]->TrackId()];
mcp = particleVec[ip];

}
}
}

if (maxe < 0.0*tote)
mcp = 0;

return mcp;
}
*/

bool isInTPC(double x, double y, double z){

  double DET_X_HIGH = 256.0;
  double DET_X_LOW = 0.0;
  double DET_Y_HIGH = 116.5;
  double DET_Y_LOW = -116.5;
  double DET_Z_HIGH = 1040;
  double DET_Z_LOW = 0.0;

  if (x > DET_X_LOW && x < DET_X_HIGH &&
    y > DET_Y_LOW && y < DET_Y_HIGH &&
    z > DET_Z_LOW && z < DET_Z_HIGH)
    return true;
    else return false;

  }
  /*
  void initialiseLandau(TF1* landau, double meanVal, std::string ptype){

  double width;
  if (ptype == "muon")
  width = 0.099;
  else if (ptype == "proton")
  width = 0.2113;
  else
  throw std::runtime_error("NO NO NO, Why are you trying to use a particle type that isn't a proton or a muon? Are you stupid, Adam?");

  landau->SetParameters(1, meanVal, width);
  double integral = landau->Integral(0,width*100);
  double mean = landau->Mean(0, width*100);
  landau->SetParameters(1./integral, meanVal - (mean - meanVal), width);

}
*/
double likelihood(double obs, double trues, std::string ptype){

  TF1* landau = new TF1("landau", "landau", 0, 500);
  //  initialiseLandau(landau, trues, ptype);

  double width;
  if (ptype == "muon")
  width = 0.099;
  else if (ptype == "proton")
  width = 0.211;
  else
  throw std::runtime_error("NO NO NO, Why are you trying to use a particle type that isn't a proton or a muon? Are you stupid, Adam?");

  landau->SetParameters(1, trues, width);
  double integral = landau->Integral(0,500);
  double mean = landau->Mean(trues-100, trues+100);

  //std::cout << "[CALCLIKELIHOOD] MEAN: " << mean << " INTEGRAL: " << integral << std::endl;

  landau->SetParameters(1./integral, trues - (mean - trues), width);

  //std::cout << "[CALCLIKELIHOOD] VAL: " << obs << " LIKELIHOOD " << landau->Eval(obs) << std::endl;

  double landauVal = landau->Eval(obs);
  double likelihoodval;

  if (landauVal == 0)
  likelihoodval = -2*std::log(1e-50);
  else{
    likelihoodval = -2*std::log(landau->Eval(obs));
  }
  //std::cout << "[CALCLIKELIHOOD] LIKELIHOOD VAL: " << likelihoodval << std::endl;

  return likelihoodval;

}

bool isWellReconstructed(recob::Track track, simb::MCParticle mcp){

  double tsy = track.Start().Y();
  double tey = track.End().Y();
  double tsz = track.Start().Z();
  double tez = track.End().Z();

  double msy = mcp.Vy();
  double mey = mcp.EndY();
  double msz = mcp.Vz();
  double mez = mcp.EndZ();

  double twoDStartRes = std::sqrt(std::pow(msy-tsy,2)+std::pow(msz-tsz,2));
  double twoDStartResFlip = std::sqrt(std::pow(msy-tey,2)+std::pow(msz-tez,2));
  double twoDEndRes = std::sqrt(std::pow(mey-tey,2)+std::pow(mez-tez,2));
  double twoDEndResFlip = std::sqrt(std::pow(mey-tsy,2)+std::pow(mez-tsz,2));

  if ((twoDStartRes < 2.0 && twoDEndRes < 2.0) ||
  (twoDStartResFlip < 2.0 && twoDEndResFlip < 2.0)){
    std::cout << "[CALCLIKELIHOOD] Found a match!" << std::endl;
    std::cout << "[CALCLIKELIHOOD] Start Res Fwd : " << twoDStartRes << std::endl;
    std::cout << "[CALCLIKELIHOOD] End Res Fwd   : " << twoDEndRes << std::endl;
    std::cout << "[CALCLIKELIHOOD] Start Res Bwd : " << twoDStartResFlip << std::endl;
    std::cout << "[CALCLIKELIHOOD] End Res Bwd   : " << twoDEndResFlip << std::endl;
    return true;
  }
  else {
    /*    std::cout << "[CALCLIKELIHOOD] Start Res Fwd : " << twoDStartRes << std::endl;
    std::cout << "[CALCLIKELIHOOD] End Res Fwd   : " << twoDEndRes << std::endl;
    std::cout << "[CALCLIKELIHOOD] Start Res Bwd : " << twoDStartResFlip << std::endl;
    std::cout << "[CALCLIKELIHOOD] End Res Bwd   : " << twoDEndResFlip << std::endl;
    */    return false;
  }
}

int main(int argv, char** argc){

  std::vector<std::string> filenames;

  // extract input filename from input textfile
  std::string file_name;
  std::ifstream input_file(argc[1]);
  while (getline(input_file,file_name))
  filenames.push_back(file_name);

  art::InputTag trackTag { "pandoraNu::McRecoStage2" };
  art::InputTag hitTag { "gaushit" };
  art::InputTag hitMatcherTag { "gaushitTruthMatch" };
  art::InputTag caloTag { "pandoraNucali" }; // calibrated

  TFile *fOutput = new TFile("outputfile.root", "RECREATE");
  TFile *fBadLikelihood = new TFile("badLikelihood.root", "RECREATE");
  TFile *fGoodLikelihood = new TFile("goodLikelihood.root", "RECREATE");
  TFile *fBadLikelihoodMcpDists = new TFile("badLikelihoodMcpDists.root", "RECREATE");
  //TFile *fGoodLikelihoodMcpDists = new TFile("goodLikelihoodMcpDists.root", "RECREATE");

  // introduce randomness
  //TH1D* hTrueProtonLikelihoodProton = new TH1D("hTrueProtonLikelihoodProton", ";Likelihood(Proton);", 1000, -100, 100);
  //TH1D* hTrueProtonLikelihoodMuon = new TH1D("hTrueProtonLikelihoodMuon", ";Likelihood(Muon);", 1000, -100, 100);
  //TH2D* hTrueProtonLikelihoodProtonLikelihoodMuon = new TH2D("hTrueProtonLikelihoodProtonLikelihoodMuon", ";Likelihood Muon; Likelihood Proton", 1000, -100, 100, 1000, -100, 100);
  //TH1D* hTrueMuonLikelihoodMuon = new TH1D("hTrueMuonLikelihoodMuon", ";Likelihood(Muon);", 1000, -10, 10);
  //TH1D* hTrueMuonLikelihoodProton = new TH1D("hTrueMuonLikelihoodProton", ";Likelihood(Proton);", 1000, -10, 10);
  //TH2D* hTrueMuonLikelihoodProtonLikelihoodMuon = new TH2D("hTrueMuonLikelihoodProtonLikelihoodMuon", ";Likelihood Muon; Likelihood Proton", 1000, -100, 100, 1000, -100, 100);
  TH1D* hTrueMuonMuonLikelihoodRatio = new TH1D("hTrueMuonMuonLikelihoodRatio", ";Score;", 20, 0, 1);
  TH1D* hTrueProtonMuonLikelihoodRatio = new TH1D("hTrueProtonMuonLikelihoodRatio", ";Score;", 20, 0, 1);
  TH2D* hTrueMuonMuonLikelihoodRatioRecoLength = new TH2D("hTrueMuonMuonLikelihoodRatioRecoLength", ";Score;Track Length", 20, 0, 1, 700, 0, 700);
  TH2D* hTrueProtonMuonLikelihoodRatioRecoLength = new TH2D("hTrueProtonMuonLikelihoodRatioRecoLength", ";Score;Track Length", 20, 0, 1, 700, 0, 700);
  TH1D* hTrueMuonMcpThetaXZ = new TH1D("hTrueMuonMcpThetaXZ", ";Theta;", 100, -3.15, 3.15);
  TH1D* hTrueMuonMcpThetaYZ = new TH1D("hTrueMuonMcpThetaYZ", ";Theta;", 200, -3.15, 3.15);
  TH2D* hTrueMuonMcpThetaXZThetaYZ = new TH2D("hTrueMuonMcpThetaXZThetaYZ", ";#theta_{XZ};#theta_{YZ}",100, -3.15, 3.15, 200, -3.15, 3.15);
  TH2D* hTrueMuonMcpYZStart = new TH2D("hTrueMuonMcpYZStart", ";Z (cm);Y (cm)", 200, 0, 1040, 100, -116.5, 116.5);
  TH2D* hTrueMuonMcpYZEnd = new TH2D("hTrueMuonMcpYZEnd", ";Z (cm);Y (cm)", 200, 0, 1040, 100, -116.5, 116.5);
  TH1D* hTrueProtonMcpThetaXZ = new TH1D("hTrueProtonMcpThetaXZ", ";Theta;", 100, 0, 3.15);
  TH1D* hTrueProtonMcpThetaYZ = new TH1D("hTrueProtonMcpThetaYZ", ";Theta;", 200, -3.15, 3.15);
  TH2D* hTrueProtonMcpThetaXZThetaYZ = new TH2D("hTrueProtonMcpThetaXZThetaYZ", ";#theta_{XZ};#theta_{YZ}",100, 0, 3.15, 200, -3.15, 3.15);
  TH2D* hTrueProtonMcpYZStart = new TH2D("hTrueProtonMcpYZStart", ";Z (cm);Y (cm)", 200, 0, 1040, 100, -116.5, 116.5);
  TH2D* hTrueProtonMcpYZEnd = new TH2D("hTrueProtonMcpYZEnd", ";Z (cm);Y (cm)", 200, 0, 1040, 100, -116.5, 116.5);

  // new Histograms
  TH1D* h_trueProton_protonLikelihoodZeroToOne = new TH1D("h_trueProton_protonLikelihoodZeroToOne", "h_trueProton_protonLikelihoodZeroToOne", 50, 0, 1);
  TH1D* h_trueProton_muonLikelihoodZeroToOne = new TH1D("h_trueProton_muonLikelihoodZeroToOne", "h_trueProton_muonLikelihoodZeroToOne", 50, 0, 1);
  TH1D* h_trueMuon_protonLikelihoodZeroToOne = new TH1D("h_trueMuon_protonLikelihoodZeroToOne", "h_trueMuon_protonLikelihoodZeroToOne", 50, 0, 1);
  TH1D* h_trueMuon_muonLikelihoodZeroToOne = new TH1D("h_trueMuon_muonLikelihoodZeroToOne", "h_trueMuon_muonLikelihoodZeroToOne", 50, 0, 1);
  TH1D* h_trueProton_MipConsistencyPull = new TH1D("h_trueProton_MipConsistencyPull", "h_trueProton_MipConsistencyPull", 50, 0, 1);
  TH1D* h_trueMuon_MipConsistencyPull = new TH1D("h_trueMuon_MipConsistencyPull", "h_trueMuon_MipConsistencyPull", 50, 0, 1);

  particleid::Theory_dEdx_resrange blah;
  TGraph* protonTheory = (TGraph*)blah.g_ThdEdxRR_Proton;
  TGraph* muonTheory   = (TGraph*)blah.g_ThdEdxRR_Muon;
  fBadLikelihood->cd();
  protonTheory->Write();
  muonTheory->Write();

  fGoodLikelihood->cd();
  protonTheory->Write();
  muonTheory->Write();

  fOutput->cd();
  protonTheory->Write();
  muonTheory->Write();

  int n_correct = 0;
  int n_total = 0;


  // begin event loop
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()){

    std::cout << "------ Processing "
    << "Run " << ev.eventAuxiliary().run() << ", "
    << "Event " << ev.eventAuxiliary().event() << ", "
    << "Time " << ev.eventAuxiliary().time().timeHigh() << "------" << std::endl;

    const auto& trackHandle = ev.getValidHandle< std::vector<recob::Track> >(trackTag);
    const auto& trackVec(*trackHandle);
    const auto& mcpHandle = ev.getValidHandle< std::vector< simb::MCParticle > >("largeant");
    //std::vector< art::Ptr<simb::MCParticle> > mcpCollection;
    //art::fill_ptr_vector(mcpCollection, mcpHandle);
    const auto& hitHandle = ev.getValidHandle< std::vector<recob::Hit> >(hitTag);

    art::FindManyP<recob::Hit> trackHitAssn(trackHandle, ev, trackTag);
    art::FindMany<anab::Calorimetry> trackCaloAssn(trackHandle, ev, caloTag);

    for (size_t it = 0; it < trackVec.size(); it++){

      // get tracks and mcparticles
      art::FindMany<simb::MCParticle, anab::BackTrackerHitMatchingData> particlesPerHit(hitHandle, ev, hitMatcherTag);
      recob::Track track = trackVec.at(it);
      //const simb::MCParticle* mcp = getMatchedParticle(trackHitAssn, particlesPerHit, hitMatcherTag, it)

      std::cout << "Looking at track with length " << track.Length() << std::endl;

      simb::MCParticle mcp;
      bool isDefined = false;
      std::cout << "[CALCLIKELIHOOD] Looping MCParticles..." << std::endl;
      for(auto const& thisMcp : (*mcpHandle)){

        if (thisMcp.StatusCode() != 1 || !isInTPC(thisMcp.EndX(), thisMcp.EndY(), thisMcp.EndZ()) || !isInTPC(thisMcp.Vx(), thisMcp.Vy(), thisMcp.Vz())) continue;

        if (isWellReconstructed(track, thisMcp)){
          isDefined = true;
          mcp = thisMcp;
          break;
        }

      }


      if (!isDefined) continue;
      std::cout << "[CALCLIKELIHOOD] Found good match" << std::endl;

      //if (!isWellReconstructed(track, mcp)) continue;

      //if (!isInTPC(mcp.EndX(), mcp.EndY(), mcp.EndZ())) continue;

      std::vector< const anab::Calorimetry* > calos;
      trackCaloAssn.get(it, calos);

      int pdgCode = std::abs(mcp.PdgCode());

      std::cout << "[CALCLIKELIHOOD] Number of Calorimetry objects: " << calos.size() << std::endl;
      for (unsigned int i = 0; i < calos.size(); i++){

        const anab::Calorimetry* calo = calos.at(i);

        if (calo->PlaneID().Plane != 2) continue;

        const std::vector<double>& dedx = calo->dEdx();
        const std::vector<double>& resrg = calo->ResidualRange();

        // Determine which way the track is going - take 10 hits (or 1/3 track)
        // from start and 10 hits (or 1/3 track) from end and see which has
        // higher average dEdx (ignore first and last hits)
        bool ForwardsGoing = false;
        bool BackwardsGoing = false;
        int onethirdtrack = std::floor((resrg.size()-2.0)/3.0);
        int n_startend = std::min(onethirdtrack, 10);
        double mean_smallrr = 0;
        double mean_highrr = 0;
        for (int i_hit=1; i_hit < n_startend-1; i_hit++){
          mean_smallrr += dedx.at(i_hit);
          mean_highrr  += dedx.at(resrg.size()-i_hit-1);
        }
        mean_smallrr/= n_startend;
        mean_highrr/= n_startend;

        std::cout << "mean_smallrr = " << mean_smallrr << ", mean_highrr = " << mean_highrr << std::endl;
        if (mean_smallrr/mean_highrr < 0.5){
          std::cout << "Track is backwards going" << std::endl;
          BackwardsGoing = true;
        }
        else if (mean_smallrr/mean_highrr > 2.0){
          std::cout << "Track is forwards going" << std::endl;
          ForwardsGoing = true;
        }
        else{
          std::cout << "Cannot tell track direction reliably" << std::endl;
        }

        // double muonWidth = 0.1;
        // double protonWidth = 0.2;

        particleid::Bragg_negLogL_Estimator algorithm;
        // algorithm.setWidthMu(muonWidth);
        // algorithm.setWidthP(protonWidth);
        algorithm.setResRangeMaxShift(5.0);
        algorithm.setResRangeMinShift(-5.0);

        double rr_shift_mufwd=0;
        double rr_shift_pfwd=0;
        double rr_shift_pifwd=0;
        double rr_shift_mubwd=0;
        double rr_shift_pbwd=0;
        double rr_shift_pibwd=0;

        /*std::cout << "[CALCLIKELIHOOD] PDG: " << pdgCode << std::endl;
        double totalProtonLogLikelihoodBwd = algorithm.getNegLogL(dedx, resrg, 2212, 0, &rr_shift_pbwd);
        std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Proton Bwd: " << totalProtonLogLikelihoodBwd << std::endl;
        std::cout << "                    rr shift = " << rr_shift_pbwd << std::endl;
        double totalMuonLogLikelihoodBwd   = algorithm.getNegLogL(dedx, resrg, 13, 0, &rr_shift_mubwd);
        std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Muon Bwd:   " << totalMuonLogLikelihoodBwd << std::endl;
        std::cout << "                    rr shift = " << rr_shift_mubwd << std::endl;
        double totalPionLogLikelihoodBwd   = algorithm.getNegLogL(dedx, resrg, 211, 0, &rr_shift_pibwd);
        std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Pion Bwd:   " << totalPionLogLikelihoodBwd << std::endl;
        std::cout << "                    rr shift = " << rr_shift_pibwd << std::endl;
        double totalProtonLogLikelihoodFwd = algorithm.getNegLogL(dedx, resrg, 2212, 1, &rr_shift_pfwd);
        std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Proton Fwd: " << totalProtonLogLikelihoodFwd << std::endl;
        std::cout << "                    rr shift = " << rr_shift_pfwd << std::endl;
        double totalMuonLogLikelihoodFwd   = algorithm.getNegLogL(dedx, resrg, 13, 1, &rr_shift_mufwd);
        std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Muon Fwd:   " << totalMuonLogLikelihoodFwd << std::endl;
        std::cout << "                    rr shift = " << rr_shift_mufwd << std::endl;
        double totalPionLogLikelihoodFwd   = algorithm.getNegLogL(dedx, resrg, 211, 1, &rr_shift_pifwd);
        std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Pion Fwd:   " << totalPionLogLikelihoodFwd << std::endl;
        std::cout << "                    rr shift = " << rr_shift_pifwd << std::endl;

        double Bragg_mu = std::min(totalMuonLogLikelihoodFwd, totalMuonLogLikelihoodBwd);
        double Bragg_p  = std::min(totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd);
        double Bragg_pi = std::min(totalPionLogLikelihoodFwd, totalPionLogLikelihoodBwd);*/

        double Bragg_mu = 0;
        double Bragg_p = 0;
        double Bragg_pi = 0;
        if (ForwardsGoing){
          Bragg_p = algorithm.getNegLogL(dedx, resrg, 2212, 1, &rr_shift_pfwd);
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Proton Fwd: " << Bragg_p << std::endl;
          std::cout << "                    rr shift = " << rr_shift_pfwd << std::endl;
          Bragg_mu   = algorithm.getNegLogL(dedx, resrg, 13, 1, &rr_shift_mufwd);
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Muon Fwd:   " << Bragg_mu << std::endl;
          std::cout << "                    rr shift = " << rr_shift_mufwd << std::endl;
          Bragg_pi   = algorithm.getNegLogL(dedx, resrg, 211, 1, &rr_shift_pifwd);
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Pion Fwd:   " << Bragg_pi << std::endl;
          std::cout << "                    rr shift = " << rr_shift_pifwd << std::endl;
        }
        if (BackwardsGoing){
          std::cout << "[CALCLIKELIHOOD] PDG: " << pdgCode << std::endl;
          Bragg_p = algorithm.getNegLogL(dedx, resrg, 2212, 0, &rr_shift_pbwd);
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Proton Bwd: " << Bragg_p << std::endl;
          std::cout << "                    rr shift = " << rr_shift_pbwd << std::endl;
          Bragg_mu   = algorithm.getNegLogL(dedx, resrg, 13, 0, &rr_shift_mubwd);
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Muon Bwd:   " << Bragg_mu << std::endl;
          std::cout << "                    rr shift = " << rr_shift_mubwd << std::endl;
          Bragg_pi   = algorithm.getNegLogL(dedx, resrg, 211, 0, &rr_shift_pibwd);
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Pion Bwd:   " << Bragg_pi << std::endl;
          std::cout << "                    rr shift = " << rr_shift_pibwd << std::endl;
        }

        if (!ForwardsGoing && !BackwardsGoing){
          // No Bragg peak - what to do?
        }

        double mipConsistencyPull = (Bragg_mu+Bragg_pi)/(Bragg_p+Bragg_mu+Bragg_pi);
        /*
        std::cout << "true PDG: " << pdgCode << std::endl;
        std::cout << "proton likelihood Fwd: " << totalProtonLogLikelihoodFwd << std::endl;
        std::cout << "proton likelihood Bwd: " << totalProtonLogLikelihoodBwd << std::endl;
        std::cout << "muon likelihood Fwd: " << totalMuonLogLikelihoodFwd << std::endl;
        std::cout << "muon likelihood Bwd: " << totalMuonLogLikelihoodBwd << std::endl;
        */
        bool correctPDG = false;
        if (pdgCode == 2212){
          if (std::min(totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd) >= std::min(totalMuonLogLikelihoodFwd, totalMuonLogLikelihoodBwd)){
            correctPDG = false;
          }
          else{
            correctPDG = true;
          }
        }
        else if (pdgCode == 13 || pdgCode == -13 || pdgCode == 211 || pdgCode == -211)
        {
          if (std::min(totalMuonLogLikelihoodFwd, totalMuonLogLikelihoodBwd) >= std::min(totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd)){
            correctPDG = false;
          }
          else{
            correctPDG = true;
          }
        }
        if (correctPDG) n_correct++;
        n_total++;

        if (pdgCode == 2212){

          h_trueProton_MipConsistencyPull->Fill(mipConsistencyPull);

          double likelihoodProtonUnderMuonAssmp   = totalMuonLogLikelihoodFwd;
          double likelihoodProtonUnderProtonAssmp = totalProtonLogLikelihoodFwd;
          double likelihoodProtonUnderMuonBwdAssmp = totalMuonLogLikelihoodBwd;
          double likelihoodProtonUnderProtonBwdAssmp = totalProtonLogLikelihoodBwd;

          double likelihoodMuMin = std::min(likelihoodProtonUnderMuonAssmp, likelihoodProtonUnderMuonBwdAssmp);
          double likelihoodPrMin = std::min(likelihoodProtonUnderProtonAssmp, likelihoodProtonUnderProtonBwdAssmp);

          //hTrueProtonLikelihoodProton->Fill(likelihoodPrMin);
          //hTrueProtonLikelihoodMuon->Fill(likelihoodMuMin);
          //hTrueProtonLikelihoodProtonLikelihoodMuon->Fill(likelihoodMuMin, likelihoodPrMin);
          hTrueProtonMuonLikelihoodRatio->Fill((likelihoodPrMin/(likelihoodMuMin+likelihoodPrMin)));
          hTrueProtonMuonLikelihoodRatioRecoLength->Fill((likelihoodPrMin/(likelihoodMuMin+likelihoodPrMin)), track.Length());

          h_trueProton_protonLikelihoodZeroToOne->Fill(likelihoodPrMin/(likelihoodPrMin+likelihoodMuMin));
          h_trueProton_muonLikelihoodZeroToOne->Fill(likelihoodMuMin/(likelihoodPrMin+likelihoodMuMin));

          if (std::min(likelihoodProtonUnderProtonAssmp, likelihoodProtonUnderProtonBwdAssmp) >= std::min(likelihoodProtonUnderMuonAssmp, likelihoodProtonUnderMuonBwdAssmp)){

            fBadLikelihood->cd();
            TH2D* h = new TH2D(Form("h_ProtonMistakenForMuon_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

            for (size_t i = 0; i < dedx.size(); i++){

              if (likelihoodProtonUnderMuonAssmp < likelihoodProtonUnderMuonBwdAssmp)
              h->Fill(resrg.at(i), dedx.at(i));
              else
              h->Fill(resrg.at(resrg.size()-i-1), dedx.at(i));

            }
            h->Write();

            fBadLikelihoodMcpDists->cd();
            double mcpThetaXZ = std::atan2(mcp.Px(), mcp.Pz());
            double mcpThetaYZ = std::atan2(mcp.Py(), mcp.Pz());
            hTrueProtonMcpThetaXZ->Fill(mcpThetaXZ);
            hTrueProtonMcpThetaYZ->Fill(mcpThetaYZ);
            hTrueProtonMcpThetaXZThetaYZ->Fill(mcpThetaXZ, mcpThetaYZ);
            hTrueProtonMcpYZStart->Fill(mcp.Vz(), mcp.Vy());
            hTrueProtonMcpYZEnd->Fill(mcp.EndZ(), mcp.EndY());

          }
          else{

            fGoodLikelihood->cd();
            TH2D* h = new TH2D(Form("h_Proton_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

            for (size_t i = 0; i < dedx.size(); i++){

              if (likelihoodProtonUnderProtonAssmp < likelihoodProtonUnderProtonBwdAssmp)
              h->Fill(resrg.at(i), dedx.at(i));
              else
              h->Fill(resrg.at(resrg.size() - i -1), dedx.at(i));

            }
            h->Write();

          }


        }

        if (pdgCode == 13 || pdgCode == 211){
          h_trueMuon_MipConsistencyPull->Fill(mipConsistencyPull);

          double likelihoodMuonUnderMuonAssmp   = totalMuonLogLikelihoodFwd;
          double likelihoodMuonUnderProtonAssmp = totalProtonLogLikelihoodFwd;
          double likelihoodMuonUnderMuonBwdAssmp = totalMuonLogLikelihoodFwd;
          double likelihoodMuonUnderProtonBwdAssmp = totalMuonLogLikelihoodBwd;

          double likelihoodMuMin = std::min(likelihoodMuonUnderMuonAssmp, likelihoodMuonUnderMuonBwdAssmp);
          double likelihoodPrMin = std::min(likelihoodMuonUnderProtonAssmp, likelihoodMuonUnderProtonBwdAssmp);

          //hTrueMuonLikelihoodMuon->Fill(likelihoodMuMin);
          //hTrueMuonLikelihoodProton->Fill(likelihoodPrMin);
          //hTrueMuonLikelihoodProtonLikelihoodMuon->Fill(likelihoodMuMin, likelihoodPrMin);
          hTrueMuonMuonLikelihoodRatio->Fill(likelihoodMuMin/(likelihoodMuMin+likelihoodPrMin));
          hTrueMuonMuonLikelihoodRatioRecoLength->Fill((likelihoodMuMin/(likelihoodMuMin+likelihoodPrMin)), track.Length());

          h_trueMuon_protonLikelihoodZeroToOne->Fill(likelihoodPrMin/(likelihoodPrMin+likelihoodMuMin));
          h_trueMuon_muonLikelihoodZeroToOne->Fill(likelihoodMuMin/(likelihoodPrMin+likelihoodMuMin));

          std::cout << "[CALCLIKELIHOOD] True Muon." << std::endl;
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Muon:   " << likelihoodMuonUnderMuonAssmp << std::endl;
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of Proton: " << likelihoodMuonUnderProtonAssmp << std::endl;
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of MuonBwd:   " << likelihoodMuonUnderMuonBwdAssmp << std::endl;
          std::cout << "[CALCLIKELIHOOD] >> Likelihood Under assumption of ProtonBwd: " << likelihoodMuonUnderProtonBwdAssmp << std::endl;

          if (std::min(likelihoodMuonUnderProtonAssmp, likelihoodMuonUnderProtonBwdAssmp) <= std::min(likelihoodMuonUnderMuonAssmp, likelihoodMuonUnderMuonBwdAssmp)){

            fBadLikelihood->cd();
            TH2D* h = new TH2D(Form("h_MuonMistakenForProton_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

            for (size_t i = 0; i < dedx.size(); i++){

              if (likelihoodMuonUnderProtonAssmp < likelihoodMuonUnderProtonBwdAssmp)
              h->Fill(resrg.at(i), dedx.at(i));
              else
              h->Fill(resrg.at(resrg.size() - i -1), dedx.at(i));

            }
            h->Write();

            fBadLikelihoodMcpDists->cd();
            double mcpThetaXZ = std::atan2(mcp.Px(), mcp.Pz());
            double mcpThetaYZ = std::atan2(mcp.Py(), mcp.Pz());
            hTrueMuonMcpThetaXZ->Fill(mcpThetaXZ);
            hTrueMuonMcpThetaYZ->Fill(mcpThetaYZ);
            hTrueMuonMcpThetaXZThetaYZ->Fill(mcpThetaXZ, mcpThetaYZ);
            hTrueMuonMcpYZStart->Fill(mcp.Vz(), mcp.Vy());
            hTrueMuonMcpYZEnd->Fill(mcp.EndZ(), mcp.EndY());


          }
          else{

            fGoodLikelihood->cd();
            TH2D* h = new TH2D(Form("h_Muon_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

            for (size_t i = 0; i < dedx.size(); i++){

              if (likelihoodMuonUnderMuonAssmp < likelihoodMuonUnderMuonBwdAssmp)
              h->Fill(resrg.at(i), dedx.at(i));
              else
              h->Fill(resrg.at(resrg.size() - i -1), dedx.at(i));

            }
            h->Write();
          }

        }

        // Save plots for all tracks to root file (for debugging, for now)
        TH2D* h = new TH2D(Form("h_run%i_event%i_track%i", ev.eventAuxiliary().run(), ev.eventAuxiliary().event(), (int)it), ";Residual Range (cm); dE/dx (MeV/cm)", 200, 0, 50, 200, 0, 50);

        for (size_t i = 0; i < dedx.size(); i++){

          // skip first and last hits in plots
          if (i == 0) continue;
          if (i == dedx.size()-1) continue;

          double min_neg2logl = std::min({totalMuonLogLikelihoodFwd, totalMuonLogLikelihoodBwd, totalProtonLogLikelihoodFwd, totalProtonLogLikelihoodBwd, totalPionLogLikelihoodFwd, totalPionLogLikelihoodBwd});

          double rr_shift;
          if (min_neg2logl == totalMuonLogLikelihoodFwd) rr_shift = rr_shift_mufwd;
          if (min_neg2logl == totalProtonLogLikelihoodFwd) rr_shift = rr_shift_pfwd;
          if (min_neg2logl == totalPionLogLikelihoodFwd) rr_shift = rr_shift_pifwd;
          if (min_neg2logl == totalMuonLogLikelihoodBwd) rr_shift = rr_shift_mubwd;
          if (min_neg2logl == totalProtonLogLikelihoodBwd) rr_shift = rr_shift_pbwd;
          if (min_neg2logl == totalPionLogLikelihoodBwd) rr_shift = rr_shift_pibwd;

          if ((min_neg2logl == totalMuonLogLikelihoodFwd) || (min_neg2logl == totalPionLogLikelihoodFwd) || (min_neg2logl == totalProtonLogLikelihoodFwd))
            h->Fill(resrg.at(i)+rr_shift, dedx.at(i));
          else
            h->Fill(resrg.at(resrg.size() - i -1)+rr_shift, dedx.at(i));

        }
        fOutput->cd();
        h->Write();

      }

    }
  }

  fOutput->cd();
  hTrueMuonMuonLikelihoodRatio->Write();
  hTrueMuonMuonLikelihoodRatioRecoLength->Write();
  //hTrueMuonLikelihoodProton->Write();
  //hTrueMuonLikelihoodMuon->Write();
  //hTrueMuonLikelihoodProtonLikelihoodMuon->Write();
  //hTrueProtonLikelihoodMuon->Write();
  //hTrueProtonLikelihoodProton->Write();
  //hTrueProtonLikelihoodProtonLikelihoodMuon->Write();
  hTrueProtonMuonLikelihoodRatio->Write();
  hTrueProtonMuonLikelihoodRatioRecoLength->Write();
  h_trueMuon_muonLikelihoodZeroToOne->Write();
  h_trueMuon_protonLikelihoodZeroToOne->Write();
  h_trueProton_muonLikelihoodZeroToOne->Write();
  h_trueProton_protonLikelihoodZeroToOne->Write();
  h_trueMuon_MipConsistencyPull->Write();
  h_trueProton_MipConsistencyPull->Write();

  fBadLikelihoodMcpDists->cd();
  hTrueProtonMcpThetaXZ->Write();
  hTrueProtonMcpThetaYZ->Write();
  hTrueProtonMcpThetaXZThetaYZ->Write();
  hTrueProtonMcpYZStart->Write();
  hTrueProtonMcpYZEnd->Write();
  hTrueMuonMcpThetaXZ->Write();
  hTrueMuonMcpThetaYZ->Write();
  hTrueMuonMcpThetaXZThetaYZ->Write();
  hTrueMuonMcpYZStart->Write();
  hTrueMuonMcpYZEnd->Write();

  std::cout << "Number of tracks correctly identified: " << n_correct << std::endl;
  std::cout << "Number of total tracks: " << n_total << std::endl;

  return 0;

}
