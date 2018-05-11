struct treevars{
  // These are the variables that are filled directly from the tree
  int true_PDG=-9999;
  double track_neglogl_fwd_p=-9999;
  double track_neglogl_fwd_mu=-9999;
  double track_neglogl_fwd_pi=-9999;
  double track_neglogl_fwd_k=-9999;
  double track_neglogl_fwd_other=-9999;
  double track_neglogl_fwd_mip=-9999;
  double track_neglogl_bwd_p=-9999;
  double track_neglogl_bwd_mu=-9999;
  double track_neglogl_bwd_pi=-9999;
  double track_neglogl_bwd_k=-9999;
  double track_neglogl_bwd_other=-9999;
  double track_PIDA_mean=-9999;
  double track_PIDA_kde=-9999;
  double track_chi2mu=-9999;
  double track_chi2p=-9999;
  double track_chi2pi=-9999;
  double track_chi2k=-9999;

  // These are derived quantities - derived from the values above in CalcPIDvars
  double track_neglogl_p=-9999;
  double track_neglogl_mu=-9999;
  double track_neglogl_pi=-9999;
  double track_neglogl_k=-9999;
  double track_neglogl_mip=-9999;
  double track_neglogl_minmumip=-9999;
  double track_neglogl_muminusp=-9999;
  double track_neglogl_mipminusp=-9999;
  double track_neglogl_minmumipminusp=-9999;
  double track_chi2_muminusp=-9999;

  // These variables are for when "neglogl" is actually not a log, but just a likelihood (then we want the maximum likelihood, and likelihood ratio instead of subtraction)
  double track_max_neglogl_p=-9999;
  double track_max_neglogl_mu=-9999;
  double track_max_neglogl_pi=-9999;
  double track_max_neglogl_k=-9999;
  double track_max_neglogl_mip=-9999;
  double track_max_neglogl_maxmumip=-9999;
  double track_neglogl_muoverp=-9999;
  double track_neglogl_mipoverp=-9999;
  double track_neglogl_minmumipoverp=-9999;
  double track_Lmu_0to1=-9999;
  double track_Lmip_0to1=-9999;
  double track_Lpi_0to1=-9999;
  double track_Lp_0to1=-9999;
  double track_Lmumip_0to1=-9999;
  double track_Lmumippi_0to1=-9999;
};

void settreevars(TTree *intree, treevars *varstoset){
  intree->SetBranchAddress("true_PDG"              , &(varstoset->true_PDG));
  intree->SetBranchAddress("track_neglogl_fwd_p"   , &(varstoset->track_neglogl_fwd_p));
  intree->SetBranchAddress("track_neglogl_fwd_mu"  , &(varstoset->track_neglogl_fwd_mu));
  intree->SetBranchAddress("track_neglogl_fwd_pi"  , &(varstoset->track_neglogl_fwd_pi));
  intree->SetBranchAddress("track_neglogl_fwd_k"   , &(varstoset->track_neglogl_fwd_k));
  intree->SetBranchAddress("track_neglogl_fwd_mip" , &(varstoset->track_neglogl_fwd_mip));
  intree->SetBranchAddress("track_neglogl_bwd_p"   , &(varstoset->track_neglogl_bwd_p));
  intree->SetBranchAddress("track_neglogl_bwd_mu"  , &(varstoset->track_neglogl_bwd_mu));
  intree->SetBranchAddress("track_neglogl_bwd_pi"  , &(varstoset->track_neglogl_bwd_pi));
  intree->SetBranchAddress("track_neglogl_bwd_k"   , &(varstoset->track_neglogl_bwd_k));
  intree->SetBranchAddress("track_PIDA_mean"       , &(varstoset->track_PIDA_mean));
  intree->SetBranchAddress("track_PIDA_kde"        , &(varstoset->track_PIDA_kde));
  intree->SetBranchAddress("track_Chi2Muon", &(varstoset->track_chi2mu));
  intree->SetBranchAddress("track_Chi2Proton", &(varstoset->track_chi2p));
  intree->SetBranchAddress("track_Chi2Pion", &(varstoset->track_chi2pi));
  intree->SetBranchAddress("track_Chi2Kaon", &(varstoset->track_chi2k));
}

void CalcPIDvars(treevars *vars){
  vars->track_neglogl_p     = std::min(vars->track_neglogl_fwd_p     , vars->track_neglogl_bwd_p);
  vars->track_neglogl_mu    = std::min(vars->track_neglogl_fwd_mu    , vars->track_neglogl_bwd_mu);
  vars->track_neglogl_pi    = std::min(vars->track_neglogl_fwd_pi    , vars->track_neglogl_bwd_pi);
  vars->track_neglogl_k     = std::min(vars->track_neglogl_fwd_k     , vars->track_neglogl_bwd_k);
  vars->track_neglogl_mip   = vars->track_neglogl_fwd_mip;

  vars->track_neglogl_minmumip = std::min(vars->track_neglogl_mu, vars->track_neglogl_mip);

  vars->track_neglogl_muminusp = vars->track_neglogl_mu - vars->track_neglogl_p;
  vars->track_neglogl_mipminusp = vars->track_neglogl_mip - vars->track_neglogl_p;
  vars->track_neglogl_minmumipminusp = vars->track_neglogl_minmumip - vars->track_neglogl_p;
  vars->track_chi2_muminusp = vars->track_chi2mu - vars->track_chi2p;

  // These variables are for when "neglogl" is actually not a log, but just a likelihood (then we want the maximum likelihood, and likelihood ratio instead of subtraction)
  vars->track_max_neglogl_p = std::max(vars->track_neglogl_fwd_p     , vars->track_neglogl_bwd_p);
  vars->track_max_neglogl_mu = std::max(vars->track_neglogl_fwd_mu    , vars->track_neglogl_bwd_mu);
  vars->track_max_neglogl_pi = std::max(vars->track_neglogl_fwd_pi    , vars->track_neglogl_bwd_pi);
  vars->track_max_neglogl_k = std::max(vars->track_neglogl_fwd_k     , vars->track_neglogl_bwd_k);
  vars->track_max_neglogl_mip = vars->track_neglogl_fwd_mip;
  vars->track_max_neglogl_maxmumip = std::max(vars->track_neglogl_mu, vars->track_neglogl_mip);
  vars->track_neglogl_muoverp = vars->track_neglogl_mu / vars->track_neglogl_p;
  vars->track_neglogl_mipoverp = vars->track_neglogl_mip / vars->track_neglogl_p;
  vars->track_neglogl_minmumipoverp = vars->track_neglogl_minmumip / vars->track_neglogl_p;

  vars->track_Lmu_0to1 = vars->track_neglogl_mu/(vars->track_neglogl_p+vars->track_neglogl_mu+vars->track_neglogl_k+vars->track_neglogl_p+vars->track_neglogl_mip);
  vars->track_Lmip_0to1 = vars->track_neglogl_mip/(vars->track_neglogl_p+vars->track_neglogl_mu+vars->track_neglogl_k+vars->track_neglogl_p+vars->track_neglogl_mip);
  vars->track_Lpi_0to1 = vars->track_neglogl_pi/(vars->track_neglogl_p+vars->track_neglogl_mu+vars->track_neglogl_k+vars->track_neglogl_p+vars->track_neglogl_mip);
  vars->track_Lp_0to1 = vars->track_neglogl_p/(vars->track_neglogl_p+vars->track_neglogl_mu+vars->track_neglogl_k+vars->track_neglogl_p+vars->track_neglogl_mip);

  vars->track_Lmumip_0to1 = (vars->track_neglogl_mu+vars->track_neglogl_mip)/(vars->track_neglogl_p+vars->track_neglogl_mu+vars->track_neglogl_k+vars->track_neglogl_p+vars->track_neglogl_mip);
  vars->track_Lmumippi_0to1 = (vars->track_neglogl_mu+vars->track_neglogl_mip+vars->track_neglogl_pi)/(vars->track_neglogl_p+vars->track_neglogl_mu+vars->track_neglogl_k+vars->track_neglogl_p+vars->track_neglogl_mip);
}


// --------------------------------------------------- //
// This struct contains the histograms and all functions related to them

struct hist1D{
  TH1D *h_mu;
  TH1D *h_p;
  TH1D *h_pi;
  TH1D *h_k;
  TH1D *h_other;
  TH1D *h_all;

  TLegend *l;

  // Constructor for this struct of hists
  hist1D(std::string name, std::string title, double nbins, double binlow, double binhigh){
    h_mu = new TH1D(std::string(name+"_mu").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_p = new TH1D(std::string(name+"_p").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_pi = new TH1D(std::string(name+"_pi").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_k = new TH1D(std::string(name+"_k").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_other = new TH1D(std::string(name+"_other").c_str(),title.c_str(),nbins,binlow,binhigh);

    h_all = new TH1D(std::string(name+"_all").c_str(),title.c_str(),nbins,binlow,binhigh);

    h_mu->SetFillColor(TColor::GetColor(8,64,129));
    h_p->SetFillColor(TColor::GetColor(215, 48, 39));
    h_pi->SetFillColor(TColor::GetColor(166,217,106));
    h_k->SetFillColor(TColor::GetColor(133,1,98));
    h_other->SetFillColor(TColor::GetColor(197,197,197));
    // "All" styling for MC (data styling is set in DrawData)
    h_all->SetFillStyle(3345);
    h_all->SetFillColor(kGray+2);
    h_all->SetMarkerSize(0.); // bad hack because root keeps drawing markers and I can't make it stop

    l = new TLegend(0.59,0.64,0.81,0.87);
    l->AddEntry(h_p,"True proton","f");
    l->AddEntry(h_mu,"True muon","f");
    l->AddEntry(h_pi,"True pion","f");
    l->AddEntry(h_k,"True kaon","f");
    l->AddEntry(h_other,"True other","f");
  }
};

  void FillHist(hist1D *hists, double value, int pdg){
    // Fill "all" histogram for every entry
    hists->h_all->Fill(value);

    // Now fill histograms by particle type
    if (TMath::Abs(pdg)==13){ // muon
      hists->h_mu->Fill(value);
    }
    else if (TMath::Abs(pdg)==2212){ // proton
      hists->h_p->Fill(value);
    }
    else if (TMath::Abs(pdg)==211){ // pion
      hists->h_pi->Fill(value);
    }
    else if (TMath::Abs(pdg)==321){ // kaon
      hists->h_k->Fill(value);
    }
    else{ // other
      hists->h_other->Fill(value);
    }
  }

  void DrawMC(hist1D *hists, double POTScaling){
    if (POTScaling == 0.){ // area normalise
      POTScaling = 1./hists->h_all->Integral();
      hists->h_all->GetYaxis()->SetTitle("No. tracks (area normalised)");
    }
    else hists->h_all->GetYaxis()->SetTitle("No. tracks (POT normalised)");

    hists->h_mu->Sumw2();
    hists->h_p->Sumw2();
    hists->h_pi->Sumw2();
    hists->h_k->Sumw2();
    hists->h_other->Sumw2();
    hists->h_all->Sumw2();

    hists->h_mu->Scale(POTScaling);
    hists->h_p->Scale(POTScaling);
    hists->h_pi->Scale(POTScaling);
    hists->h_k->Scale(POTScaling);
    hists->h_other->Scale(POTScaling);
    hists->h_all->Scale(POTScaling);

    THStack *hs = new THStack("hs","hs");
    hs->Add(hists->h_p);
    hs->Add(hists->h_mu);
    hs->Add(hists->h_pi);
    hs->Add(hists->h_k);
    hs->Add(hists->h_other);

    hists->h_all->SetMaximum(hists->h_all->GetMaximum()*1.2);
    hists->h_all->SetMinimum(0);
    hists->h_all->Draw("hist"); // Draw this one first because it knows about the axis titles
    hs->Draw("same hist");
    hists->h_all->Draw("same E2"); // Draw it again so errors are on top
    hists->l->Draw();
  }

void OverlayData(TCanvas *c, hist1D *onbeam, hist1D *offbeam, double OffBeamScaling, double POTScaling){
  TH1D *h_onminusoff = (TH1D*)onbeam->h_all->Clone("h_onminusoff");

  h_onminusoff->Sumw2();
  offbeam->h_all->Sumw2();

  h_onminusoff->Add(offbeam->h_all,-1.0*OffBeamScaling);
  if (POTScaling==0){
    h_onminusoff->Scale(1.0/h_onminusoff->Integral());
  }

  h_onminusoff->SetMarkerStyle(20);
  h_onminusoff->SetMarkerSize(0.6);

  c->cd();
  h_onminusoff->Draw("same p E1");

  TLegend *l = (TLegend*)c->GetPrimitive("TPave");
  l->AddEntry(h_onminusoff,"Data (on-off beam)","lp");
}


void DrawMCEffPur(TCanvas *c, hist1D *hists, bool MIPlow){
  std::vector<TH1D*> histstoeval = {
                    hists->h_mu,
                    hists->h_pi,
                    hists->h_p
                    };

  std::vector<std::string> histtitles = {
                    "True muons",
                    "True pions",
                    "True protons"
                    };

  c->Divide(2,2,0.0005,0.0005);

  TLegend *l = new TLegend(0.2,0.2,0.8,0.8);
  l->SetTextFont(132);
  l->SetLineColor(kWhite);
  l->SetFillColor(kWhite);

  for (int i_h=0; i_h<histstoeval.size(); i_h++){
    TH1D *heff = (TH1D*)hists->h_all->Clone("heff");
    TH1D *hpur = (TH1D*)hists->h_all->Clone("hpur");
    TH1D *heffpur = (TH1D*)hists->h_all->Clone("heffpur");

    heff->Clear();
    hpur->Clear();
    heffpur->Clear();

    heff->SetTitle(histtitles.at(i_h).c_str());

    for (int i_bin=1; i_bin <= histstoeval.at(i_h)->GetXaxis()->GetNbins(); i_bin++){

      double eff, pur;//, efferr, purerr;
      if ((MIPlow && i_h < 2) || (!MIPlow && i_h ==2)){ // integrate from the bottom
        double selected_i = histstoeval.at(i_h)->Integral(1,i_bin);
        double total_i = histstoeval.at(i_h)->Integral();
        double selected_all = hists->h_all->Integral(1,i_bin);

        eff = selected_i/total_i;
        pur = selected_i/selected_all;

        if (selected_i==0 && selected_all==0) pur = 0;
      }
      else{ // integrate up to the top
        double selected_i = histstoeval.at(i_h)->Integral(i_bin,histstoeval.at(i_h)->GetXaxis()->GetNbins());
        double total_i = histstoeval.at(i_h)->Integral();
        double selected_all = hists->h_all->Integral(i_bin,histstoeval.at(i_h)->GetXaxis()->GetNbins());

        eff = selected_i/total_i;
        pur = selected_i/selected_all;

        if (selected_i==0 && selected_all==0) pur = 0;
      }

      double effpur = eff*pur;
      heff->SetBinContent(i_bin,eff);
      hpur->SetBinContent(i_bin,pur);
      heffpur->SetBinContent(i_bin,effpur);
    }

    heff->SetLineColor(kRed);
    heff->SetMarkerColor(kRed);
    heff->SetMarkerStyle(20);
    heff->SetMarkerSize(.3
    );
    hpur->SetLineColor(kBlue);
    hpur->SetMarkerColor(kBlue);
    hpur->SetMarkerStyle(20);
    hpur->SetMarkerSize(.3
    );
    heffpur->SetLineColor(kBlack);
    heffpur->SetMarkerColor(kBlack);
    heffpur->SetMarkerStyle(20);
    heffpur->SetMarkerSize(.3
    );

    heff->GetYaxis()->SetRangeUser(0,1.1);

    c->cd(i_h+1);
    heff->Draw("p");
    hpur->Draw("same p");
    heffpur->Draw("same p");

    if (i_h==0){
      l->AddEntry(heff,"Efficiency","p");
      l->AddEntry(hpur,"Purity","p");
      l->AddEntry(heffpur,"Efficiency #times Purity","p");
    }

    c->cd(histstoeval.size()+1);
    l->Draw();
  }

}
