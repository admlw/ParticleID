void plotNewDataMcFromTree(){


  int nbins = 40;
  int binlow = 0;
  int binhigh = 15;

  //double protonScaling = 1.563;
  //double muonScaling = 1.563; 

  double protonScaling = 1.0;
  double muonScaling = 1.0;

  //TFile *f_bnbcos  = new TFile("/uboone/data/users/alister1/particleID/180423-ParticleId/pid_bnbcos.root", "read");
  //TFile *f_bnbcos  = new TFile("pidtest.root", "read");
  TFile *f_bnbcos  = new TFile("/uboone/data/users/alister1/particleID/180423-ParticleId/pid_bnbcos_newdists.root", "read");
  TFile *f_offbeam  = new TFile("pid_bnbonbeam.root", "read");
  TFile *f_onbeam  = new TFile("pid_bnboffbeam.root", "read");

  //TFile *f_offbeam = new TFile("/uboone/data/users/alister1/particleID/180420-ParticleId/pid_offbeam.root", "read");
  //TFile *f_onbeam  = new TFile("/uboone/data/users/alister1/particleID/180420-ParticleId/pid_onbeam.root", "read");

  TTree *t_bnbcos  = (TTree*)f_bnbcos->Get("pidvalid/pidTree"); 
  TTree *t_onbeam = (TTree*)f_onbeam->Get("pidvalid/pidTree"); 
  TTree *t_offbeam = (TTree*)f_offbeam->Get("pidvalid/pidTree"); 


  //TTree *t_offbeam = (TTree*)f_offbeam->Get("pidvalid/pidTree"); 
  //TTree *t_onbeam  = (TTree*)f_onbeam->Get("pidvalid/pidTree"); 

  std::vector<TString> plotNames = {
    "track_neglogl_p",
    "track_neglogl_mu",
    "track_neglogl_pi",
    "track_neglogl_k",
    "track_neglogl_mip"};

    int true_PDG;
    double bnbcos_track_neglogl_fwd_p;
    double bnbcos_track_neglogl_fwd_mu;
    double bnbcos_track_neglogl_fwd_pi;
    double bnbcos_track_neglogl_fwd_k;
    double bnbcos_track_neglogl_fwd_other;
    double bnbcos_track_neglogl_fwd_mip;
    double bnbcos_track_neglogl_bwd_p;
    double bnbcos_track_neglogl_bwd_mu;
    double bnbcos_track_neglogl_bwd_pi;
    double bnbcos_track_neglogl_bwd_k;
    double bnbcos_track_neglogl_bwd_other;
    double bnbcos_track_PIDA_mean;
    double bnbcos_track_PIDA_kde;
    double offbeam_track_neglogl_fwd_p;
    double offbeam_track_neglogl_fwd_mu;
    double offbeam_track_neglogl_fwd_pi;
    double offbeam_track_neglogl_fwd_k;
    double offbeam_track_neglogl_fwd_other;
    double offbeam_track_neglogl_fwd_mip;
    double offbeam_track_neglogl_bwd_p;
    double offbeam_track_neglogl_bwd_mu;
    double offbeam_track_neglogl_bwd_pi;
    double offbeam_track_neglogl_bwd_k;
    double offbeam_track_neglogl_bwd_other;
    double offbeam_track_PIDA_mean;
    double offbeam_track_PIDA_kde;
    double onbeam_track_neglogl_fwd_p;
    double onbeam_track_neglogl_fwd_mu;
    double onbeam_track_neglogl_fwd_pi;
    double onbeam_track_neglogl_fwd_k;
    double onbeam_track_neglogl_fwd_other;
    double onbeam_track_neglogl_fwd_mip;
    double onbeam_track_neglogl_bwd_p;
    double onbeam_track_neglogl_bwd_mu;
    double onbeam_track_neglogl_bwd_pi;
    double onbeam_track_neglogl_bwd_k;
    double onbeam_track_neglogl_bwd_other;
    double onbeam_track_PIDA_mean;
    double onbeam_track_PIDA_kde;

    double track_PIDA_mean;
    double track_PIDA_kde;

    t_bnbcos->SetBranchAddress("true_PDG"              , &true_PDG);
    t_bnbcos->SetBranchAddress("track_neglogl_fwd_p"   , &bnbcos_track_neglogl_fwd_p);
    t_bnbcos->SetBranchAddress("track_neglogl_fwd_mu"  , &bnbcos_track_neglogl_fwd_mu);
    t_bnbcos->SetBranchAddress("track_neglogl_fwd_pi"  , &bnbcos_track_neglogl_fwd_pi);
    t_bnbcos->SetBranchAddress("track_neglogl_fwd_k"   , &bnbcos_track_neglogl_fwd_k);
    t_bnbcos->SetBranchAddress("track_neglogl_fwd_mip" , &bnbcos_track_neglogl_fwd_mip);
    t_bnbcos->SetBranchAddress("track_neglogl_bwd_p"   , &bnbcos_track_neglogl_bwd_p);
    t_bnbcos->SetBranchAddress("track_neglogl_bwd_mu"  , &bnbcos_track_neglogl_bwd_mu);
    t_bnbcos->SetBranchAddress("track_neglogl_bwd_pi"  , &bnbcos_track_neglogl_bwd_pi);
    t_bnbcos->SetBranchAddress("track_neglogl_bwd_k"   , &bnbcos_track_neglogl_bwd_k);
    t_bnbcos->SetBranchAddress("track_PIDA_mean"       , &bnbcos_track_PIDA_mean);
    t_bnbcos->SetBranchAddress("track_PIDA_kde"        , &bnbcos_track_PIDA_kde);

    t_onbeam->SetBranchAddress("track_neglogl_fwd_p"   , &onbeam_track_neglogl_fwd_p);
    t_onbeam->SetBranchAddress("track_neglogl_fwd_mu"  , &onbeam_track_neglogl_fwd_mu);
    t_onbeam->SetBranchAddress("track_neglogl_fwd_pi"  , &onbeam_track_neglogl_fwd_pi);
    t_onbeam->SetBranchAddress("track_neglogl_fwd_k"   , &onbeam_track_neglogl_fwd_k);
    t_onbeam->SetBranchAddress("track_neglogl_fwd_mip" , &onbeam_track_neglogl_fwd_mip);
    t_onbeam->SetBranchAddress("track_neglogl_bwd_p"   , &onbeam_track_neglogl_bwd_p);
    t_onbeam->SetBranchAddress("track_neglogl_bwd_mu"  , &onbeam_track_neglogl_bwd_mu);
    t_onbeam->SetBranchAddress("track_neglogl_bwd_pi"  , &onbeam_track_neglogl_bwd_pi);
    t_onbeam->SetBranchAddress("track_neglogl_bwd_k"   , &onbeam_track_neglogl_bwd_k);
    t_onbeam->SetBranchAddress("track_PIDA_mean"       , &onbeam_track_PIDA_mean);
    t_onbeam->SetBranchAddress("track_PIDA_kde"        , &onbeam_track_PIDA_kde);

    t_offbeam->SetBranchAddress("track_neglogl_fwd_p"   , &offbeam_track_neglogl_fwd_p);
    t_offbeam->SetBranchAddress("track_neglogl_fwd_mu"  , &offbeam_track_neglogl_fwd_mu);
    t_offbeam->SetBranchAddress("track_neglogl_fwd_pi"  , &offbeam_track_neglogl_fwd_pi);
    t_offbeam->SetBranchAddress("track_neglogl_fwd_k"   , &offbeam_track_neglogl_fwd_k);
    t_offbeam->SetBranchAddress("track_neglogl_fwd_mip" , &offbeam_track_neglogl_fwd_mip);
    t_offbeam->SetBranchAddress("track_neglogl_bwd_p"   , &offbeam_track_neglogl_bwd_p);
    t_offbeam->SetBranchAddress("track_neglogl_bwd_mu"  , &offbeam_track_neglogl_bwd_mu);
    t_offbeam->SetBranchAddress("track_neglogl_bwd_pi"  , &offbeam_track_neglogl_bwd_pi);
    t_offbeam->SetBranchAddress("track_neglogl_bwd_k"   , &offbeam_track_neglogl_bwd_k);
    t_offbeam->SetBranchAddress("track_PIDA_mean"       , &offbeam_track_PIDA_mean);
    t_offbeam->SetBranchAddress("track_PIDA_kde"        , &offbeam_track_PIDA_kde);

  for (int j = 0; j < plotNames.size(); j++){

    TH1D* h_bnbcos_p  = new TH1D("h_bnbcos_p", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_mu = new TH1D("h_bnbcos_mu", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_pi = new TH1D("h_bnbcos_pi", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_k  = new TH1D("h_bnbcos_k", ";-2NegLL_p;", nbins, binlow, binhigh);
    TH1D* h_bnbcos_other = new TH1D("h_bnbcos_other", ";-2NegLL_p;", nbins, binlow, binhigh);

    TH1D* h_offbeam = new TH1D("h_offbeam", ";;", nbins, binlow, binhigh);
    TH1D* h_onbeam = new TH1D("h_onbeam", ";;", nbins, binlow, binhigh);

    for (int i = 0; i < t_bnbcos->GetEntries(); i++){

      t_bnbcos->GetEntry(i);

      double track_neglogl_p     = std::min(bnbcos_track_neglogl_fwd_p     , bnbcos_track_neglogl_bwd_p);
      double track_neglogl_mu    = std::min(bnbcos_track_neglogl_fwd_mu    , bnbcos_track_neglogl_bwd_mu);
      double track_neglogl_pi    = std::min(bnbcos_track_neglogl_fwd_pi    , bnbcos_track_neglogl_bwd_pi);
      double track_neglogl_k     = std::min(bnbcos_track_neglogl_fwd_k     , bnbcos_track_neglogl_bwd_k);
      double track_neglogl_other = std::min(bnbcos_track_neglogl_fwd_other , bnbcos_track_neglogl_bwd_other);
      double track_neglogl_mip   = bnbcos_track_neglogl_fwd_mip;

        std::vector<double> minLogLikelihoods = {
         track_neglogl_p,
         track_neglogl_mu,
         track_neglogl_pi,
         track_neglogl_k,
         track_neglogl_mip
        };


      if (j < 5){ 
        if (std::abs(true_PDG) == 2212)
          h_bnbcos_p->Fill(minLogLikelihoods.at(j)*1./protonScaling);
        else if (std::abs(true_PDG) == 13)
          h_bnbcos_mu->Fill(minLogLikelihoods.at(j)*1./muonScaling);
        else if (std::abs(true_PDG) == 211)
          h_bnbcos_pi->Fill(minLogLikelihoods.at(j)*1./muonScaling);
        else if (std::abs(true_PDG) == 321)
          h_bnbcos_k->Fill(minLogLikelihoods.at(j)*1./muonScaling);
        else
          h_bnbcos_other->Fill(minLogLikelihoods.at(j)*1./muonScaling); 
      }
    }

    for (int i = 0; i < t_onbeam->GetEntries(); i++){
      t_onbeam->GetEntry(i);

      double track_neglogl_p     = std::min(onbeam_track_neglogl_fwd_p     , onbeam_track_neglogl_bwd_p);
      double track_neglogl_mu    = std::min(onbeam_track_neglogl_fwd_mu    , onbeam_track_neglogl_bwd_mu);
      double track_neglogl_pi    = std::min(onbeam_track_neglogl_fwd_pi    , onbeam_track_neglogl_bwd_pi);
      double track_neglogl_k     = std::min(onbeam_track_neglogl_fwd_k     , onbeam_track_neglogl_bwd_k);
      double track_neglogl_other = std::min(onbeam_track_neglogl_fwd_other , onbeam_track_neglogl_bwd_other);
      double track_neglogl_mip   = onbeam_track_neglogl_fwd_mip;

      std::vector<double> minLogLikelihoods = {
         track_neglogl_p,
         track_neglogl_mu,
         track_neglogl_pi,
         track_neglogl_k,
         track_neglogl_mip
        };

      h_onbeam->Fill(minLogLikelihoods.at(j));

    }
    for (int i = 0; i < t_offbeam->GetEntries(); i++){
       
      t_offbeam->GetEntry(i);
      
      double track_neglogl_p     = std::min(offbeam_track_neglogl_fwd_p     , offbeam_track_neglogl_bwd_p);
      double track_neglogl_mu    = std::min(offbeam_track_neglogl_fwd_mu    , offbeam_track_neglogl_bwd_mu);
      double track_neglogl_pi    = std::min(offbeam_track_neglogl_fwd_pi    , offbeam_track_neglogl_bwd_pi);
      double track_neglogl_k     = std::min(offbeam_track_neglogl_fwd_k     , offbeam_track_neglogl_bwd_k);
      double track_neglogl_other = std::min(offbeam_track_neglogl_fwd_other , offbeam_track_neglogl_bwd_other);
      double track_neglogl_mip   = offbeam_track_neglogl_fwd_mip;

      std::vector<double> minLogLikelihoods = {
         track_neglogl_p,
         track_neglogl_mu,
         track_neglogl_pi,
         track_neglogl_k,
         track_neglogl_mip
        };

      h_offbeam->Fill(minLogLikelihoods.at(j));
    
    }


    h_bnbcos_p->SetLineWidth(0);
    h_bnbcos_mu->SetLineWidth(0);
    h_bnbcos_pi->SetLineWidth(0);
    h_bnbcos_k->SetLineWidth(0);
    h_bnbcos_other->SetLineWidth(0);
    h_bnbcos_p->SetFillColor(TColor::GetColor(215, 48, 39));
    h_bnbcos_mu->SetFillColor(TColor::GetColor(8,64,129));
    h_bnbcos_pi->SetFillColor(TColor::GetColor(166,217,106));
    h_bnbcos_k->SetFillColor(TColor::GetColor(133,1,98));
    h_bnbcos_other->SetFillColor(TColor::GetColor(197,197,197));
    h_bnbcos_p->SetMarkerColor(TColor::GetColor(215, 48, 39));
    h_bnbcos_mu->SetMarkerColor(TColor::GetColor(8,64,129));
    h_bnbcos_pi->SetMarkerColor(TColor::GetColor(166,217,106));
    h_bnbcos_k->SetMarkerColor(TColor::GetColor(133,1,98));
    h_bnbcos_other->SetMarkerColor(TColor::GetColor(197,197,197));

    /*
       h_bnbcos_p->Sumw2();
       h_bnbcos_mu->Sumw2();
       h_bnbcos_pi->Sumw2();
       h_bnbcos_k->Sumw2();
       h_bnbcos_other->Sumw2();
       */

    TH1D* h_total = (TH1D*)h_bnbcos_p->Clone("h_total");
    h_total->Add(h_bnbcos_mu);
    h_total->Add(h_bnbcos_pi);
    h_total->Add(h_bnbcos_k);
    h_total->Add(h_bnbcos_other);  

    h_total->Sumw2();

    h_bnbcos_p->Scale(1./h_total->Integral());
    h_bnbcos_mu->Scale(1./h_total->Integral());
    h_bnbcos_pi->Scale(1./h_total->Integral());
    h_bnbcos_k->Scale(1./h_total->Integral());
    h_bnbcos_other->Scale(1./h_total->Integral());
    h_total->Scale(1./h_total->Integral());

    THStack *hs = new THStack("hs", "hs");
    hs->Add(h_bnbcos_p);
    hs->Add(h_bnbcos_mu);
    hs->Add(h_bnbcos_pi);
    hs->Add(h_bnbcos_k);
    hs->Add(h_bnbcos_other);

    hs->SetTitle(";"+plotNames.at(j)+";");
    
       
       h_offbeam->Scale(0.78);
       TH1D* h_onminusoff = (TH1D*)h_onbeam->Clone("h_onminusoff");

       h_onminusoff->Add(h_offbeam, -1);
       h_onminusoff->Sumw2();
       h_onminusoff->Scale(1./h_onminusoff->Integral());
   
    TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);

    c1->cd();
    hs->SetMaximum(std::max(h_onminusoff->GetMaximum(), h_total->GetMaximum())*1.15);
    hs->Draw();
    h_total->SetFillStyle(3345);
    h_total->SetFillColor(kGray+2);
    h_total->Draw("sameE2");
    h_onminusoff->SetMarkerStyle(20);
    h_onminusoff->SetMarkerSize(0.6);
    h_onminusoff->Draw("samepE1");

    TLegend *leg = new TLegend(0.55, 0.64, 0.89, 0.89);
    if (plotNames.at(j) == "track_neglogl_p"){
      std::cout << "gotya" << std::endl;
      leg->SetX1(0.15);
      leg->SetX2(0.45);
    }
    leg->AddEntry(h_bnbcos_p, "True proton");
    leg->AddEntry(h_bnbcos_mu, "True muon");
    leg->AddEntry(h_bnbcos_pi, "True pion");
    leg->AddEntry(h_bnbcos_k, "True kaon");
    leg->AddEntry(h_bnbcos_other, "True other");
    leg->AddEntry(h_onminusoff, "Data (on-off)");
    leg->SetLineWidth(0);
    leg->SetFillStyle(0);
    leg->Draw("same");

    c1->SaveAs(plotNames.at(j)+".png");

    h_bnbcos_p->Delete();
    h_bnbcos_mu->Delete();
    h_bnbcos_pi->Delete();
    h_bnbcos_k->Delete();
    h_bnbcos_other->Delete();
    h_offbeam->Delete();
    h_onbeam->Delete();

  }
}
