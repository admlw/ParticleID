void makeTheoryCurvePlot(){

  TGraph* gr_proton = (TGraph*)_file0->Get("gr_proton");
  TGraph* gr_muon = (TGraph*)_file0->Get("gr_muon");
  TGraph* gr_pion = (TGraph*)_file0->Get("gr_pion");
  TGraph* gr_kaon = (TGraph*)_file0->Get("gr_kaon");

  gr_proton->SetLineColor(TColor::GetColor(215, 48, 39));
  gr_muon->SetLineColor(TColor::GetColor(8,64,129));
  gr_pion->SetLineColor(TColor::GetColor(166,217,106));
  gr_kaon->SetLineColor(TColor::GetColor(133,1,98));
  gr_proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  gr_muon->SetMarkerColor(TColor::GetColor(8,64,129));
  gr_pion->SetMarkerColor(TColor::GetColor(166,217,106));
  gr_kaon->SetMarkerColor(TColor::GetColor(133,1,98));
  gr_proton->SetLineWidth(3);
  gr_muon->SetLineWidth(3);
  gr_pion->SetLineWidth(3);
  gr_kaon->SetLineWidth(3);

  gr_proton->SetTitle(";Residual Range (cm);dE/dx (MeV/cm)");
  gr_proton->GetXaxis()->SetRangeUser(0,30);
  gr_proton->GetYaxis()->SetTitleOffset(0.4);

  TCanvas *c1 = new TCanvas("c1", "", 1500, 500);

  c1->SetRightMargin (0.07);
  c1->SetLeftMargin(0.07);

  gr_proton->Draw();
  gr_muon->Draw("same");
  gr_pion->Draw("same");
  gr_kaon->Draw("same");

  TLegend *leg = new TLegend(0.85, 0.69, 0.93, 0.89);
  leg->AddEntry(gr_proton, "Proton");
  leg->AddEntry(gr_kaon, "Kaon");
  leg->AddEntry(gr_pion, "Pion");
  leg->AddEntry(gr_muon, "Muon");

  leg->Draw("same");

  c1->SaveAs("theoryCurves.png");

}
