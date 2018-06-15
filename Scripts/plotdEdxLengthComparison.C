/**
 * Dear future Adam and Kirsty:
 * I'm sorry this is a mess. I'm sorry that it uses a copy/pasted version of the truncated mean code.
 * I wanted a quick and dirty test of how well the truncated mean/length PID worked.
 *
 * summary looks to be that the separation is *really* excellent but the data/mc agreement is pretty
 * bad, so I don't trust it.
 *
 * Anyway to run, use:
 * root -l bnbcos.root onbeam.root offbeam.root plotdEdxLengthComparison.C
 *
 * You need to update the scale factors in the code but I've put them right at the top so they're
 * easy to find
 */ 

#include "truncatedMean.C"

/**
 * code to convert dQdx to dE/dx. Not actually used in this script
 * except to check that I did the dE/dx to dQ/dx conversion correctly
 */ 
double convertdQdxTodEdx(double dQdx){
  double rho      = 1.383;                    // LAr density in g/cm^3
  double Wion     = 1000./4.237e7;        // 23.6 eV = 1e, Wion in MeV/e
  double E_field  = 0.273;                           // Electric Field in the drift region in KV/cm
  double Beta     = 0.212 / (rho * E_field);
  double Alpha    = 0.930;
  double dEdx = (exp(Beta * Wion * dQdx ) - Alpha) / Beta;

  return dEdx;
}

/**
 * convert dE/dx (MeV/cm) to dQ/dx (e-/cm). This is done because we don't save dQ/dx because it
 * is a silly variable that nobody should ever use, but this PID makes use of the <dQ/dx>_{tr}
 */ 
double convertdEdxTodQdx(double dEdx){
  double rho      = 1.383;                    // LAr density in g/cm^3
  double Wion     = 1000./4.237e7;        // 23.6 eV = 1e, Wion in MeV/e
  double E_field  = 0.273;                           // Electric Field in the drift region in KV/cm
  double Beta     = 0.212 / (rho * E_field);
  double Alpha    = 0.930;

  double dQdx = (std::log((dEdx*Beta)+Alpha))/(Beta*Wion);

  return dQdx;
}

/**
 * helper function which converts a vector of dE/dx values to a vector of
 * dQ/dx values
 */
std::vector<double> getdQdxVector(std::vector<double>* dEdxVals){

  std::vector<double> dQdxVals;

  for (size_t i = 0; i < dEdxVals->size(); i++){

    dQdxVals.push_back(convertdEdxTodQdx(dEdxVals->at(i)));

  }

  return dQdxVals;

}

/**
 * helper function to convert ADC/cm to e-/cm. This isn't needed because
 * the above function already converts to e-/cm.
 */ 
double convertADCToE(double dQdx){

  return dQdx; //* 196.98;

}

/**
 * main
 */ 
void plotdEdxLengthComparison(){

  // length cuts are used to look at a slice of the
  // 2D distribution in dQ/dx space.
  double track_length_high_cut = 300;
  double track_length_low_cut = 0;
  double off_beam_scaling = 0.77;
  double mc_scaling = 1.4;
  double proton_dqdx_scaling = 1.1;
  double proton_normalisation_scaling = 0.5;
  int highval = 300000;
  
  TTree *tree_mc = (TTree*)_file0->Get("pidvalid/pidTree");

  /** BNB+COS MC plots and variables */
  double track_length_mc;
  std::vector<double> *track_dEdx_mc = 0;
  std::vector<double> *track_dEdx_mc_perhit_y = 0;
  int true_PDG;

  tree_mc->SetBranchAddress("track_length", &track_length_mc);
  tree_mc->SetBranchAddress("track_dEdx"  , &track_dEdx_mc);
  tree_mc->SetBranchAddress("track_dEdx_perhit_y", &track_dEdx_mc_perhit_y);
  tree_mc->SetBranchAddress("true_PDG"    , &true_PDG);

  TTree *tree_data_onbeam = (TTree*)_file1->Get("pidvalid/pidTree");

  TH2D* h_dEdx_length_proton = new TH2D("h_dEdx_length_proton", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 300);
  TH2D* h_dEdx_length_muon   = new TH2D("h_dEdx_length_muon", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 300);
  TH2D* h_dEdx_length_pion   = new TH2D("h_dEdx_length_pion", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 300);
  TH2D* h_dEdx_length_kaon   = new TH2D("h_dEdx_length_kaon", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 300);
  TH2D* h_dEdx_length_other  = new TH2D("h_dEdx_length_other", ";<dQ/dx>_{tr} (e-/cm);track length (cm)", 50, -1, highval, 50, 0, 300);
  TH2D* h_dEdx_length_total  = new TH2D("h_dEdx_length_total", ";<dQ/dx>_{tr} (e-/cm);track length(cm)", 50, -1, highval, 50, 0, 300);

  TH1D* h_muID_truemu = new TH1D("h_muID_truemu", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_muID_truep = new TH1D("h_muID_truep", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_muID_truepi = new TH1D("h_muID_truepi", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_muID_truek = new TH1D("h_muID_truek", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_muID_trueother = new TH1D("h_muID_trueother", ";<dQdx>_{tr};", 50, 0, highval);

  TH1D* h_pID_truemu = new TH1D("h_pID_truemu", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_pID_truep = new TH1D("h_pID_truep", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_pID_truepi = new TH1D("h_pID_truepi", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_pID_truek = new TH1D("h_pID_truek", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_pID_trueother = new TH1D("h_pID_trueother", ";<dQdx>_{tr};", 50, 0, highval);

  /** On-beam plots and variables */
  double track_length_data_onbeam;
  std::vector<double> *track_dEdx_data_onbeam = 0;
  std::vector<double> *track_dEdx_data_onbeam_perhit_y = 0;

  tree_data_onbeam->SetBranchAddress("track_length", &track_length_data_onbeam);
  tree_data_onbeam->SetBranchAddress("track_dEdx"  , &track_dEdx_data_onbeam);
  tree_data_onbeam->SetBranchAddress("track_dEdx_perhit_y", &track_dEdx_data_onbeam_perhit_y);

  TTree *tree_data_offbeam = (TTree*)_file2->Get("pidvalid/pidTree");

  TH2D* h_dEdx_length_data_onbeam = new TH2D("h_dEdx_length_data_onbeam", ";<dQdx>_{tr} (e-/cm)", 50, -1, highval, 50, 0, 300);
  TH1D* h_muID_data_onbeam = new TH1D("h_muID_data_onbeam", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_pID_data_onbeam = new TH1D("h_pID_data_onbeam", ";<dQdx>_{tr};", 50, 0, highval);

  /** off-beam plots and variables */
  double track_length_data_offbeam;
  std::vector<double> *track_dEdx_data_offbeam = 0;
  std::vector<double> *track_dEdx_data_offbeam_perhit_y = 0;

  tree_data_offbeam->SetBranchAddress("track_length", &track_length_data_offbeam);
  tree_data_offbeam->SetBranchAddress("track_dEdx"  , &track_dEdx_data_offbeam);
  tree_data_offbeam->SetBranchAddress("track_dEdx_perhit_y", &track_dEdx_data_offbeam_perhit_y);

  TH2D* h_dEdx_length_data_offbeam = new TH2D("h_dEdx_length_data_offbeam", ";<dQdx>_{tr} (e-/cm)", 50, -1, highval, 50, 0, 300);
  TH1D* h_muID_data_offbeam = new TH1D("h_muID_data_offbeam", ";<dQdx>_{tr};", 50, 0, highval);
  TH1D* h_pID_data_offbeam = new TH1D("h_pID_data_offbeam", ";<dQdx>_{tr};", 50, 0, highval);

  /** loop bnb+cos mc */
  for (int i = 0; i < tree_mc->GetEntries(); i++){

    std::fstream infile("dQdxSeparator2.txt");
    tree_mc->GetEntry(i);

    /** convert dEdx values to dQdx values */
    std::vector<double> track_dQdx_mc_perhit_y = getdQdxVector(track_dEdx_mc_perhit_y);

    double truncatedMeandQdx = (double)CalcIterativeTruncMean(track_dQdx_mc_perhit_y, 1, 1, 0, 1, 0.1, 1.0, std::numeric_limits<double>::max());
    double dQdxInE = convertADCToE(truncatedMeandQdx);

    // get cut value
    std::string s;
    int rnd_trklen = std::round(track_length_mc);
    int cutvalue = 0;
    int j = 1;
    while (std::getline(infile,s)){
      if (j == rnd_trklen){
        cutvalue = std::atof(s.c_str());
        break;
      }
      j++;
    }

    infile.close();

    h_dEdx_length_total->Fill(dQdxInE*proton_dqdx_scaling, track_length_mc);

    if (std::abs(true_PDG) == 2212 && track_length_mc > track_length_low_cut && track_length_mc < track_length_high_cut){
      h_dEdx_length_proton->Fill(dQdxInE*proton_dqdx_scaling, track_length_mc);

      if (dQdxInE*proton_dqdx_scaling < cutvalue)
        h_muID_truep->Fill(dQdxInE*proton_dqdx_scaling);
      else h_pID_truep->Fill(dQdxInE*proton_dqdx_scaling);


    }
    else if (std::abs(true_PDG) == 13 && track_length_mc > track_length_low_cut && track_length_mc < track_length_high_cut){
      h_dEdx_length_muon->Fill(dQdxInE, track_length_mc);

      if (dQdxInE < cutvalue)
        h_muID_truemu->Fill(dQdxInE);
      else h_pID_truemu->Fill(dQdxInE);


    }
    else if (std::abs(true_PDG) == 211 && track_length_mc > track_length_low_cut && track_length_mc < track_length_high_cut){
      h_dEdx_length_pion->Fill(dQdxInE, track_length_mc);

      if (dQdxInE < cutvalue)
        h_muID_truepi->Fill(dQdxInE);
      else h_pID_truepi->Fill(dQdxInE);

    }
    else if (std::abs(true_PDG) == 321 && track_length_mc > track_length_low_cut && track_length_mc < track_length_high_cut){
      h_dEdx_length_kaon->Fill(dQdxInE, track_length_mc);

      if (dQdxInE < cutvalue)
        h_muID_truek->Fill(dQdxInE);
      else h_pID_truek->Fill(dQdxInE);

    }
    else if (track_length_mc > track_length_low_cut && track_length_mc < track_length_high_cut ){
      h_dEdx_length_other->Fill(dQdxInE, track_length_mc);

      if (dQdxInE < cutvalue)
        h_muID_trueother->Fill(dQdxInE);
      else h_pID_trueother->Fill(dQdxInE);

    }

  }

  /** loop onbeam */
  for (int i = 0; i < tree_data_onbeam->GetEntries() ; i++){

    std::fstream infile("dQdxSeparator2.txt");
    tree_data_onbeam->GetEntry(i);

    if (track_dEdx_data_onbeam_perhit_y->size() == 0)
      continue;

    /** convert dEdx values to dQdx values */
    std::vector<double> track_dQdx_data_onbeam_perhit_y = getdQdxVector(track_dEdx_data_onbeam_perhit_y);

    double truncatedMeandQdx = (double)CalcIterativeTruncMean(track_dQdx_data_onbeam_perhit_y, 1, 1, 0, 1, 0.1, 1.0, std::numeric_limits<double>::max());
    double dQdxInE = convertADCToE(truncatedMeandQdx);

    // get cut value
    std::string s;
    int rnd_trklen = std::round(track_length_data_onbeam);
    int cutvalue = 0;
    int j = 1;
    while (std::getline(infile,s)){
      if (j == rnd_trklen){
        cutvalue = std::atof(s.c_str());
        break;
      }
      j++;
    }

    infile.close();

    h_dEdx_length_data_onbeam->Fill(dQdxInE, track_length_data_onbeam);
    if ( track_length_data_onbeam > track_length_low_cut && track_length_data_onbeam < track_length_high_cut){
    if (dQdxInE < cutvalue)
      h_muID_data_onbeam->Fill(dQdxInE);
    else h_pID_data_onbeam->Fill(dQdxInE);
    }
  }

  /** loop offbeam */
  for (int i = 0; i < tree_data_offbeam->GetEntries() ; i++){

    std::fstream infile("dQdxSeparator2.txt");
    tree_data_offbeam->GetEntry(i);

    if (track_dEdx_data_offbeam_perhit_y->size() == 0)
      continue;

    /** convert dEdx values to dQdx values */
    std::vector<double> track_dQdx_data_offbeam_perhit_y = getdQdxVector(track_dEdx_data_offbeam_perhit_y);

    double truncatedMeandQdx = (double)CalcIterativeTruncMean(track_dQdx_data_offbeam_perhit_y, 1, 1, 0, 1, 0.1, 1.0, std::numeric_limits<double>::max());
    double dQdxInE = convertADCToE(truncatedMeandQdx);

    // get cut value
    std::string s;
    int rnd_trklen = std::round(track_length_data_offbeam);
    int cutvalue = 0;
    int j = 1;
    while (std::getline(infile,s)){
      if (j == rnd_trklen){
        cutvalue = std::atof(s.c_str());
        break;
      }
      j++;
    }

    infile.close();

    h_dEdx_length_data_offbeam->Fill(dQdxInE, track_length_data_offbeam);
    
    if ( track_length_data_onbeam > track_length_low_cut && track_length_data_onbeam < track_length_high_cut){

    if (dQdxInE < cutvalue)
      h_muID_data_offbeam->Fill(dQdxInE);
    else h_pID_data_offbeam->Fill(dQdxInE);
    }
  }

  h_dEdx_length_proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  h_dEdx_length_muon->SetMarkerColor(TColor::GetColor(8, 64, 129));
  h_dEdx_length_pion->SetMarkerColor(TColor::GetColor(166, 217, 106));
  h_dEdx_length_kaon->SetMarkerColor(TColor::GetColor(133, 1, 98));
  h_dEdx_length_other->SetMarkerColor(TColor::GetColor(197, 197, 197));

  h_dEdx_length_proton->Draw();
  h_dEdx_length_muon->Draw("same");
  h_dEdx_length_pion->Draw("same");
  h_dEdx_length_kaon->Draw("same");
  h_dEdx_length_other->Draw("same");

  TCanvas *c2 = new TCanvas();
  h_muID_truep->SetFillColor(TColor::GetColor(215, 48, 39));
  h_muID_truemu->SetFillColor(TColor::GetColor(8,64,129));
  h_muID_truepi->SetFillColor(TColor::GetColor(166,217,106));
  h_muID_truek->SetFillColor(TColor::GetColor(133, 1, 98));
  h_muID_trueother->SetFillColor(TColor::GetColor(197, 197, 197));
  h_muID_data_offbeam->SetFillColor(kBlack);
  h_muID_data_offbeam->SetFillStyle(3345);

  h_muID_truep->Scale(mc_scaling*proton_normalisation_scaling);
  h_muID_truemu->Scale(mc_scaling);
  h_muID_truepi->Scale(mc_scaling);
  h_muID_truek->Scale(mc_scaling);
  h_muID_trueother->Scale(mc_scaling);
  h_muID_data_offbeam->Scale(off_beam_scaling);

  THStack *hs = new THStack("hs", ";<dQ/dx>_{tr};");
  hs->Add(h_muID_data_offbeam);
  hs->Add(h_muID_truep);
  hs->Add(h_muID_truemu);
  hs->Add(h_muID_truepi);
  hs->Add(h_muID_truek);
  hs->Add(h_muID_trueother);

  hs->Draw();
  h_muID_data_onbeam->SetMarkerStyle(20);
  h_muID_data_onbeam->SetMarkerSize(0.6);
  h_muID_data_onbeam->Draw("psame");

  TCanvas *c3 = new TCanvas();
  h_pID_truep->SetFillColor(TColor::GetColor(215, 48, 39));
  h_pID_truemu->SetFillColor(TColor::GetColor(8,64,129));
  h_pID_truepi->SetFillColor(TColor::GetColor(166,217,106));
  h_pID_truek->SetFillColor(TColor::GetColor(133, 1, 98));
  h_pID_trueother->SetFillColor(TColor::GetColor(197, 197, 197));
  h_pID_data_offbeam->SetFillColor(kBlack);
  h_pID_data_offbeam->SetFillStyle(3345);

  h_pID_truep->Scale(mc_scaling*proton_normalisation_scaling);
  h_pID_truemu->Scale(mc_scaling);
  h_pID_truepi->Scale(mc_scaling);
  h_pID_truek->Scale(mc_scaling);
  h_pID_trueother->Scale(mc_scaling);
  h_pID_data_offbeam->Scale(off_beam_scaling);

  THStack *hs2 = new THStack("hs", ";<dQ/dx>_{tr};");
  hs2->Add(h_pID_data_offbeam);
  hs2->Add(h_pID_truep);
  hs2->Add(h_pID_truemu);
  hs2->Add(h_pID_truepi);
  hs2->Add(h_pID_truek);
  hs2->Add(h_pID_trueother);

  hs2->Draw("");
  h_pID_data_onbeam->SetMarkerStyle(20);
  h_pID_data_onbeam->SetMarkerSize(0.6);
  h_pID_data_onbeam->Draw("psame");

  std::cout << "number of protons: " << h_muID_truep->Integral() + h_pID_truep->Integral() << std::endl;
  std::cout << "proton Efficiency: " << h_pID_truep->Integral() / (h_muID_truep->Integral() + h_pID_truep->Integral()) << std::endl;
  std::cout << "proton purity:     " << h_pID_truep->Integral() / (h_pID_truemu->Integral() + h_pID_truepi->Integral() + h_pID_truek->Integral() + h_pID_trueother->Integral() + h_pID_truep->Integral()) << std::endl;

  std::cout << "number of muons: " << h_muID_truemu->Integral() + h_pID_truemu->Integral() << std::endl;
  std::cout << "muon Efficiency: " << h_muID_truemu->Integral() / (h_muID_truemu->Integral() + h_pID_truemu->Integral()) << std::endl;
  std::cout << "muon purity:     " << h_muID_truemu->Integral() / (h_muID_truemu->Integral() + h_muID_truepi->Integral() + h_muID_truek->Integral() + h_muID_trueother->Integral() + h_muID_truep->Integral()) << std::endl;

  TCanvas *c4 = new TCanvas();
  
  TH2D* hAdded = new TH2D("hAdded", ";<dQdx>_{tr} (e-/cm)", 50, -1, highval, 50, 0, 300);
  hAdded->Add(h_dEdx_length_muon);
  hAdded->Add(h_dEdx_length_pion);
  hAdded->Add(h_dEdx_length_kaon);
  hAdded->Add(h_dEdx_length_other);
  h_dEdx_length_proton->Scale(mc_scaling*proton_normalisation_scaling);
  hAdded->Add(h_dEdx_length_proton);

  hAdded->Scale(1./hAdded->Integral());
  h_dEdx_length_data_offbeam->Scale(off_beam_scaling);
  h_dEdx_length_data_onbeam->Add(h_dEdx_length_data_offbeam, -1.);
  h_dEdx_length_data_onbeam->Scale(1./h_dEdx_length_data_onbeam->Integral());

  h_dEdx_length_data_onbeam->Divide(hAdded);
  h_dEdx_length_data_onbeam->GetZaxis()->SetRangeUser(0,2);
  h_dEdx_length_data_onbeam->Draw("colz");

}
