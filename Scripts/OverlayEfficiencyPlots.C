std::vector<int> marker_styles = {20,24,34};
std::vector<double> marker_sizes = {0.4,0.4,0.5};

void OverlayEfficiencyPlots(std::vector<std::string> inputfiles, std::vector<std::string> filedescriptions={}){

  gStyle->SetOptStat(0);

  // Check size of inputfiles vector
  std::cout << "Overlaying efficiency/purity curves from " << inputfiles.size() << " files." << std::endl;
  if (inputfiles.size() == 0){
    std::cout << "You need to give me some files! Exiting." << std::endl;\
    return;
  }

  TFile *file0 = new TFile(inputfiles.at(0).c_str(),"open");

  // Get list of plots
  // We want one entry in this list for each plane/variable
  // We *do not* want separate entries for efficiency, purity, and efficiency*purity, and we do not want separate entries for muons, pions, and protons
  // So: let's just ask for unique names that start with heff_ and end with _mu
  std::vector<std::string> plotnames;
  TIter next(file0->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())){
    // Get TH1s only
    if (!TString(key->GetClassName()).Contains("TH1")) continue;
    TString name = TString(key->GetName());
    if (!name.Contains("heff_")) continue;
    if (!name.EndsWith("_mu")) continue;
    name.Remove(0,5);
    name.Remove(TString::kTrailing,'u');
    name.Remove(TString::kTrailing,'m');
    name.Remove(TString::kTrailing,'_');
    plotnames.push_back(std::string(name.Data()));
  }

  // for (int i=0; i<plotnames.size(); i++){
  //   std::cout << plotnames.at(i) << std::endl;
  // }

  // Now we have a list of all the plots we want to overlay, let's get going!
  std::vector<std::string> particles = {"_mu","_pi","_p"};
  // Loop over plots
  for (size_t i_plot=0; i_plot<plotnames.size(); i_plot++){

    TCanvas *c1 = new TCanvas();
    c1->Divide(2,2,0.0005,0.0005);

    TLegend *l = new TLegend(0.2,0.2,0.8,0.8);
    l->SetTextFont(132);
    l->SetLineColor(kWhite);
    l->SetFillColor(kWhite);

    // Loop over input files
    for (size_t i_file=0; i_file < inputfiles.size(); i_file++){

      TFile *ftmp = new TFile(inputfiles.at(i_file).c_str(),"open");

      // Loop over particles (mu, pi, p)
      for (size_t i_p=0; i_p<particles.size(); i_p++){

        TH1D *heff = (TH1D*)ftmp->Get(TString::Format("heff_%s%s",plotnames.at(i_plot).c_str(),particles.at(i_p).c_str()).Data());
        TH1D *hpur = (TH1D*)ftmp->Get(TString::Format("hpur_%s%s",plotnames.at(i_plot).c_str(),particles.at(i_p).c_str()).Data());
        TH1D *heffpur = (TH1D*)ftmp->Get(TString::Format("heffpur_%s%s",plotnames.at(i_plot).c_str(),particles.at(i_p).c_str()).Data());

        if (!heff) std::cout << "Did not find heff in file " << inputfiles.at(i_file) << std::endl;
        if (!hpur) std::cout << "Did not find hpur in file " << inputfiles.at(i_file) << std::endl;
        if (!heffpur) std::cout << "Did not find heffpur in file " << inputfiles.at(i_file) << std::endl;

        int markerstyle = 20;
        if (marker_styles.size() < i_file) markerstyle = marker_styles.at(marker_styles.size()-1)+(i_file - marker_styles.size());
        else markerstyle = marker_styles.at(i_file);

        double markersize = 0.3;
        if (marker_sizes.size() >= i_file) markersize = marker_sizes.at(i_file);

        heff->SetFillStyle(0);
        heff->SetLineColor(kRed);
        heff->SetMarkerColor(kRed);
        heff->SetMarkerStyle(markerstyle);
        heff->SetMarkerSize(markersize);

        hpur->SetFillStyle(0);
        hpur->SetLineColor(kBlue);
        hpur->SetMarkerColor(kBlue);
        hpur->SetMarkerStyle(markerstyle);
        hpur->SetMarkerSize(markersize);

        heffpur->SetFillStyle(0);
        heffpur->SetLineColor(kBlack);
        heffpur->SetMarkerColor(kBlack);
        heffpur->SetMarkerStyle(markerstyle);
        heffpur->SetMarkerSize(markersize);

        heff->GetYaxis()->SetRangeUser(0,1.1);

        c1->cd(i_p+1);
        if (i_file==0){
          heff->Draw("lp");
          hpur->Draw("same lp");
          heffpur->Draw("same lp");
        }
        else{
          heff->Draw("same p");
          hpur->Draw("same p");
          heffpur->Draw("same p");
        }

        if (i_p==0 && i_file==0){
          l->AddEntry(heff,"Efficiency","l");
          l->AddEntry(hpur,"Purity","l");
          l->AddEntry(heffpur,"Efficiency #times Purity","l");
        }
        if (i_p==0){
          if (i_file==0){
            if (filedescriptions.size()==0){
              l->AddEntry(heffpur,TString::Format("File %d",(int)i_file).Data(),"pl");
            }
            else{
              l->AddEntry(heffpur,filedescriptions.at(i_file).c_str(),"pl");
            }
          }
          else{
            if (filedescriptions.size()==0){
              l->AddEntry(heffpur,TString::Format("File %d",(int)i_file).Data(),"p");
            }
            else{
              l->AddEntry(heffpur,filedescriptions.at(i_file).c_str(),"p");
            }
          }
        }

      } // end loop over particles (i_p)

      //ftmp->Close();

    } // end loop over input files (i_file)

    c1->cd(particles.size()+1);
    l->Draw();
    c1->Print(TString::Format("%s.png",plotnames.at(i_plot).c_str()).Data());
    delete c1;
  } // end loop over plots (i_plot)

}