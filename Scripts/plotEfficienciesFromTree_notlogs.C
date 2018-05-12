#include "plotFromTreeHeader.h"

// What variables do we want these plots as a function of?
std::vector<double> GetPIDvarstoplot(treevars *vars){
  std::vector<double> varstoplot = {
                      vars->track_max_neglogl_p,
                      vars->track_max_neglogl_mu,
                      vars->track_max_neglogl_pi,
                      vars->track_max_neglogl_k,
                      vars->track_max_neglogl_mip,
                      vars->track_max_neglogl_maxmumip,
                      vars->track_chi2mu,
                      vars->track_chi2p,
                      vars->track_chi2pi,
                      vars->track_chi2k,
                      vars->track_PIDA_kde,
                      vars->track_PIDA_mean,
                      vars->track_neglogl_muoverp,
                      vars->track_neglogl_mipoverp,
                      vars->track_neglogl_minmumipoverp,
                      vars->track_chi2_muminusp,
                      vars->track_Lmu_0to1,
                      vars->track_Lmip_0to1,
                      vars->track_Lpi_0to1,
                      vars->track_Lp_0to1,
                      vars->track_Lmumip_0to1,
                      vars->track_Lmumippi_0to1,
                      vars->track_depE_minus_rangeE_mu,
                      vars->track_depE_minus_rangeE_p
                      };
  return varstoplot;
};

// Binning (nbins, binlow, binhigh) in the same order as the vector above
std::vector<std::vector<double>> bins = {
                    {20,0,0.6}, // track_max_neglogl_p
                    {40,0,0.6}, // track_max_neglogl_mu
                    {40,0,0.6}, // track_max_neglogl_pi
                    {40,0,0.4}, // track_max_neglogl_k
                    {40,0,0.6}, // track_max_neglogl_mip
                    {40,0,0.6}, // track_max_neglogl_minmumip
                    {25,0,125}, // track_chi2mu
                    {30,0,300}, // track_chi2p
                    {25,0,125}, // track_chi2pi
                    {30,0,300}, // track_chi2k
                    {40,0,30}, // track_PIDA_kde
                    {40,0,30}, // track_PIDA_mean
                    {60,0,60}, // track_neglogl_muoverp
                    {60,0,60}, // track_neglogl_mipoverp
                    {60,0,60}, // track_neglogl_minmumipoverp
                    {50,-400,100}, // track_chi2_muminusp
                    {50,0,1}, // track_Lmu_0to1
                    {50,0,1}, // track_Lmip_0to1
                    {50,0,3}, // track_Lpi_0to1
                    {50,0,1}, // track_Lp_0to1
                    {50,0,1}, // track_Lmumip_0to1
                    {50,0,3}, // track_Lmumippi_0to1
                    {100,-100,100}, // track_depE_minus_rangeE_mu
                    {100,-100,100} // track_depE_minus_rangeE_p
                    };

// Histogram titles in the same order as the vector above
std::vector<std::string> histtitles = {
                    ";L_{p};",
                    ";L_{#mu};",
                    ";L_{#pi};",
                    ";L_{K};",
                    ";L_{MIP};",
                    ";L_{#mu/MIP};",
                    ";#chi^{2}_{#mu};",
                    ";#chi^{2}_{p};",
                    ";#chi^{2}_{#pi};",
                    ";#chi^{2}_{K};",
                    ";PIDa (by KDE);",
                    ";PIDa (by mean);",
                    ";(L_{#mu})/(L_{p});",
                    ";(L_{MIP})/(L_{p});",
                    ";(L_{#mu/MIP})/(L_{p});",
                    ";#chi^{2}_{#mu}-#chi^{2}_{p};",
                    ";L_{#mu}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{MIP}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{#pi}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{p}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP}+L_{#pi})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";Dep. E - E by range (muon assumption) [MeV];",
                    ";Dep. E - E by range (proton assumption) [MeV];"
                  };

// What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
                  "effpur_Lp",
                  "effpur_Lmu",
                  "effpur_Lpi",
                  "effpur_Lk",
                  "effpur_Lmip",
                  "effpur_Lmumip",
                  "effpur_chi2mu",
                  "effpur_chi2p",
                  "effpur_chi2pi",
                  "effpur_chi2k",
                  "effpur_pida_kde",
                  "effpur_pida_mean",
                  "effpur_Lmuoverp",
                  "effpur_Lmipoverp",
                  "effpur_Lmumipoverp",
                  "effpur_chi2muminusp",
                  "effpur_Lmu0to1",
                  "effpur_Lmip0to1",
                  "effpur_Lpi0to1",
                  "effpur_Lp0to1",
                  "effpur_Lmumip0to1",
                  "effpur_Lmumippi0to1",
                  "effpur_depErangeEmu",
                  "effpur_depErangeEp"
                };

// For efficiency/purity we need to know whether MIPs are supposed to be low or high. In the same order as the vector above
std::vector<bool> MIPlow = {
                    true, // track_neglogl_p
                    false, // track_neglogl_mu
                    false, // track_neglogl_pi
                    true, // track_neglogl_k
                    false, // track_neglogl_mip
                    false, // track_neglogl_minmumip
                    true, // track_chi2mu
                    false, // track_chi2p
                    true, // track_chi2pi
                    false, // track_chi2k
                    true, // track_PIDA_kde
                    true, // track_PIDA_mean
                    false, // track_neglogl_muminusp
                    false, // track_neglogl_mipminusp
                    false, // track_neglogl_minmumipminusp
                    true, // track_chi2_muminusp
                    false, // track_Lmu_0to1
                    false, // track_Lmip_0to1
                    false, // track_Lpi_0to1
                    true, // track_Lp_0to1
                    false, // track_Lmumip_0to1
                    false, // track_Lmumippi_0to1
                    true, // track_depE_minus_rangeE_mu
                    true // track_depE_minus_rangeE_p
                    };

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void plotEfficienciesFromTree(std::string mcfile){

  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0.);

  TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
  TTree *t_bnbcos = (TTree*)f_bnbcos->Get("pidvalid/pidTree");
  treevars mc_vars;
  settreevars(t_bnbcos,&mc_vars);

  // Sanity check: the plot vectors should be the same size
  std::vector<double> PIDvarstoplot_dummy = GetPIDvarstoplot(&mc_vars);
  // if (PIDvarstoplot_dummy.size() != bins.size()) std::cout << "WARNING PIDvarstoplot_dummy.size() = " << PIDvarstoplot_dummy.size() << "and bins.size() = " << bins.size() << ". This is going to cause you problems!" << std::endl;
  std::cout << "PIDvarstoplot.size() = " << PIDvarstoplot_dummy.size() << std::endl;
  std::cout << "bins.size() = " << bins.size() << std::endl;
  std::cout << "histtitles.size() = " << histtitles.size() << std::endl;
  std::cout << "histnames.size() = " << histnames.size() << std::endl;
  std::cout << "MIPlow.size() = " << MIPlow.size() << std::endl;

  // ----------------- MC

  // Make histograms to fill
  const size_t nplots = PIDvarstoplot_dummy.size();
  hist1D *mc_hists[nplots];
  for (int i_h=0; i_h<nplots; i_h++){
    mc_hists[i_h] = new hist1D(std::string("h_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
  }

  // Loop through MC tree and fill plots
  for (int i = 0; i < t_bnbcos->GetEntries(); i++){
    t_bnbcos->GetEntry(i);
    CalcPIDvars(&mc_vars);
    std::vector<double> PIDvarstoplot = GetPIDvarstoplot(&mc_vars);

    for (size_t i_h = 0; i_h < nplots; i_h++){
      FillHist(mc_hists[i_h],PIDvarstoplot.at(i_h),mc_vars.true_PDG);
    }

  } // end loop over entries in tree


  // -------------------- Now make all the plots

  for (size_t i_h=0; i_h < nplots; i_h++){
    TCanvas *c1 = new TCanvas();
    DrawMCEffPur(c1, mc_hists[i_h],MIPlow.at(i_h));
    c1->Print(std::string(histnames[i_h]+".png").c_str());
  }
}
