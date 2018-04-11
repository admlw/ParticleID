// Macro to make summary/overlaid plots from the output of the
// ParticleIDValidationPlots module. For now it only makes comparisons
// between true muons and true protons (because that's what we care about most)
// but true pion and kaon information is also available so those plots could be made.
//
// Run using:
// ~ root -b
// root[0] .L PIDValidationPlots_MakePretty.C
// root[1] PIDValidationPlots_MakePretty("particleIDMeta.root")
//
// The output is 1) a set of .pdf plots produced in the current directory, and
// 2) a root file PIDValidationPlots_out.root containing the same plots, also
// produced in the current directory. The percentage of particles that are
// correctly identified (by which we mean: for which the lowest neg2LL is the
// one corresponding to the correct particle type) is also printed to the screen.
//
// Kirsty Duffy, Fermilab, 11th April 2018

void PIDValidationPlots_MakePretty(std::string inputfile){

  std::vector<std::string> compare_1dplots = { "neglogl_mu",
                                            "neglogl_p",
                                            "neglogl_pi",
                                            "neglogl_K",
                                            "neglogl_muoverp",
                                            "neglogl_muminusp",
                                            "pull_mu",
                                            "pull_p",
                                            "pull_pi",
                                            "pull_K",
                                            "pull_MIP",
                                            "PIDA" };

  std::vector<std::string> compare_2dplots = { "neglogl_muvsp",
                                              "neglogl_muvspi",
                                              "dEdxtr_len" };

  std::vector<std::string> categories = {"TrueBragg", "All"};

  // --------------------------------------------------------------------- //
  // Now the actual plotting stuff starts!
  // --------------------------------------------------------------------- //

  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);

  TFile *fin = new TFile(inputfile.c_str(),"open");
  TFile *fout = new TFile("PIDValidationPlots_out.root","recreate");

  // --------------------------------------------------------------------- //
  // Comparison plots: 1D
  // --------------------------------------------------------------------- //

  // 1D comparison plots: look through all names in the compare_1dplots vector and overlay true muons with true protons
  // Remember to do this both for TrueBragg_truemu_blah and All_truemu_blah
  for (int i_name=0; i_name<compare_1dplots.size(); i_name++){
    for (int i_cat=0; i_cat < categories.size(); i_cat++){

      std::string savename = categories.at(i_cat)+"_"+compare_1dplots.at(i_name);
      std::cout << "Making 1D overlay plots for " << savename << std::endl;

      std::string name_mu = "pidvalid/"+categories.at(i_cat)+"_truemu_"+compare_1dplots.at(i_name);
      std::string name_p = "pidvalid/"+categories.at(i_cat)+"_truep_"+compare_1dplots.at(i_name);

      // Get plots
      fin->cd();
      TH1F *hmu = (TH1F*)fin->Get(name_mu.c_str());
      TH1F *hp = (TH1F*)fin->Get(name_p.c_str());

      // Style
      hmu->GetYaxis()->SetRangeUser(0,std::max(hmu->GetMaximum(),hp->GetMaximum())*1.15);
      hmu->SetLineColor(kBlue);
      hmu->SetFillColor(kBlue);
      hmu->SetFillStyle(3144);
      hp->SetLineColor(kRed);
      hp->SetFillColor(kRed);
      hp->SetFillStyle(3144);

      // Legend
      TLegend *l = new TLegend(0.6,0.6,0.87,0.87);
      l->SetLineColor(kWhite);
      l->SetFillStyle(0);
      l->SetTextFont(132);
      l->AddEntry(hmu,"True Muons","lf");
      l->AddEntry(hp,"True Protons","lf");

      // Draw
      TCanvas *c1 = new TCanvas();
      hmu->Draw();
      hp->Draw("same");
      l->Draw();

      // Save
      fout->cd();
      c1->Write(savename.c_str());
      c1->Print(std::string(savename+".pdf").c_str());

    } // loop over i_cat in categories
  } // loop over i_name in compare_1dplots



  // --------------------------------------------------------------------- //
  // Comparison plots: 2D
  // --------------------------------------------------------------------- //

  // 2D comparison plots: look through all names in the compare_2dplots vector and overlay true muons with true protons
  // Remember to do this both for TrueBragg_truemu_blah and All_truemu_blah
  for (int i_name=0; i_name<compare_2dplots.size(); i_name++){
    for (int i_cat=0; i_cat < categories.size(); i_cat++){

      std::string savename = categories.at(i_cat)+"_"+compare_2dplots.at(i_name);
      std::cout << "Making 2D overlay plots for " << savename << std::endl;

      std::string name_mu = "pidvalid/"+categories.at(i_cat)+"_truemu_"+compare_2dplots.at(i_name);
      std::string name_p = "pidvalid/"+categories.at(i_cat)+"_truep_"+compare_2dplots.at(i_name);

      // Get plots
      fin->cd();
      TH2F *hmu = (TH2F*)fin->Get(name_mu.c_str());
      TH2F *hp = (TH2F*)fin->Get(name_p.c_str());

      // Style
      //hmu->GetYaxis()->SetRangeUser(0,std::max(hmu->GetMaximum(),hp->GetMaximum())*1.15);
      hmu->SetMarkerColor(kBlue);
      hmu->SetMarkerStyle(7);
      hp->SetMarkerColor(kRed);
      hp->SetMarkerStyle(7);

      // Legend
      TLegend *l = new TLegend(0.6,0.6,0.87,0.87);
      l->SetLineColor(kWhite);
      l->SetFillStyle(0);
      l->SetTextFont(132);
      l->AddEntry(hmu,"True Muons","p");
      l->AddEntry(hp,"True Protons","p");

      // Draw
      TCanvas *c1 = new TCanvas();
      hmu->Draw();
      hp->Draw("same");
      l->Draw();

      // Save
      fout->cd();
      c1->Write(savename.c_str());
      c1->Print(std::string(savename+".pdf").c_str());

    } // loop over i_cat in categories
  } // loop over i_name in compare_1dplots



  // --------------------------------------------------------------------- //
  // Print out how often we got the PID right
  // --------------------------------------------------------------------- //

  // Remember to do this both for TrueBragg_truemu_blah and All_truemu_blah

  // TrueBragg
  std::cout << " --------- Tracks for which MCParticles have true p=0 at the end (true Bragg peaks) --------- " << std::endl;

  fin->cd();
  TH1F *h = (TH1F*)fin->Get("pidvalid/TrueBragg_truemu_smallest_neglogl");
  std::cout << "--- True muons: " << h->GetBinContent(1)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl << std::endl;
  h = (TH1F*)fin->Get("pidvalid/TrueBragg_truep_smallest_neglogl");
  std::cout << "--- True protons: " << h->GetBinContent(2)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl << std::endl;
  h = (TH1F*)fin->Get("pidvalid/TrueBragg_truepi_smallest_neglogl");
  std::cout << "--- True pions: " << h->GetBinContent(3)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl << std::endl;
  h = (TH1F*)fin->Get("pidvalid/TrueBragg_trueK_smallest_neglogl");
  std::cout << "--- True kaons: " << h->GetBinContent(4)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl << std::endl;


  // All particles
  std::cout << " --------- All tracks --------- " << std::endl;

  fin->cd();
  h = (TH1F*)fin->Get("pidvalid/All_truemu_smallest_neglogl");
  std::cout << "--- True muons: " << h->GetBinContent(1)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl << std::endl;
  h = (TH1F*)fin->Get("pidvalid/All_truep_smallest_neglogl");
  std::cout << "--- True protons: " << h->GetBinContent(2)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl << std::endl;
  h = (TH1F*)fin->Get("pidvalid/All_truepi_smallest_neglogl");
  std::cout << "--- True pions: " << h->GetBinContent(3)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl << std::endl;
  h = (TH1F*)fin->Get("pidvalid/All_trueK_smallest_neglogl");
  std::cout << "--- True kaons: " << h->GetBinContent(4)/h->Integral()*100.0 << "% identified correctly" << std::endl << std::endl
            << "         ID'd as muon: " << h->GetBinContent(1)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as proton: " << h->GetBinContent(2)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as pion: " << h->GetBinContent(3)/h->Integral()*100.0 << "%" << std::endl
            << "         ID'd as kaon: " << h->GetBinContent(4)/h->Integral()*100.0 << "%" << std::endl << std::endl;


}
