#include "plotFromTreeHeader.h"

// What variables do we want these plots as a function of?
std::vector<double> GetPIDvarstoplot(treevars *vars){
  std::vector<double> varstoplot = {
                      vars->track_neglogl_p,
                      vars->track_neglogl_mu,
                      vars->track_neglogl_pi,
                      vars->track_neglogl_k,
                      vars->track_neglogl_mip,
                      vars->track_neglogl_minmumip,
                      vars->track_chi2mu,
                      vars->track_chi2p,
                      vars->track_chi2pi,
                      vars->track_chi2k,
                      vars->track_PIDA_kde,
                      vars->track_PIDA_mean,
                      vars->track_neglogl_muminusp,
                      vars->track_neglogl_mipminusp,
                      vars->track_neglogl_minmumipminusp,
                      vars->track_chi2_muminusp,
                      vars->track_depE_minus_rangeE_mu,
                      vars->track_depE_minus_rangeE_p
                      };
  return varstoplot;
};

// Binning (nbins, binlow, binhigh) in the same order as the vector above
std::vector<std::vector<int>> bins = {
                    {40,0,15}, // track_neglogl_p
                    {40,0,15}, // track_neglogl_mu
                    {40,0,15}, // track_neglogl_pi
                    {40,0,15}, // track_neglogl_k
                    {40,0,15}, // track_neglogl_mip
                    {40,0,15}, // track_neglogl_minmumip
                    {25,0,125}, // track_chi2mu
                    {30,0,300}, // track_chi2p
                    {25,0,125}, // track_chi2pi
                    {30,0,300}, // track_chi2k
                    {40,0,30}, // track_PIDA_kde
                    {40,0,30}, // track_PIDA_mean
                    {40,-20,7}, // track_neglogl_muminusp
                    {40,-20,15}, // track_neglogl_mipminusp
                    {50,-20,7}, // track_neglogl_minmumipminusp
                    {50,-400,100}, // track_chi2_muminusp
                    {100,-100,100}, // track_depE_minus_rangeE_mu
                    {100,-100,100} // track_depE_minus_rangeE_p
                    };

// Histogram titles in the same order as the vector above
std::vector<std::string> histtitles = {
                    ";-2LogL_{p};",
                    ";-2LogL_{#mu};",
                    ";-2LogL_{#pi};",
                    ";-2LogL_{K};",
                    ";-2LogL_{MIP};",
                    ";-2LogL_{#mu/MIP};",
                    ";#chi^{2}_{#mu};",
                    ";#chi^{2}_{p};",
                    ";#chi^{2}_{#pi};",
                    ";#chi^{2}_{K};",
                    ";PIDa (by KDE);",
                    ";PIDa (by mean);",
                    ";(-2LogL_{#mu})-(-2LogL_{p});",
                    ";(-2LogL_{MIP})-(-2LogL_{p});",
                    ";(-2LogL_{#mu/MIP})-(-2LogL_{p});",
                    ";#chi^{2}_{#mu}-#chi^{2}_{p};",
                    ";Dep. E - E by range (muon assumption) [MeV];",
                    ";Dep. E - E by range (proton assumption) [MeV];"
                  };

// What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
                  "effpur_neg2LLp",
                  "effpur_neg2LLmu",
                  "effpur_neg2LLpi",
                  "effpur_neg2LLk",
                  "effpur_neg2LLmip",
                  "effpur_neg2LLmumip",
                  "effpur_chi2mu",
                  "effpur_chi2p",
                  "effpur_chi2pi",
                  "effpur_chi2k",
                  "effpur_pida_kde",
                  "effpur_pida_mean",
                  "effpur_neglogl_muminusp",
                  "effpur_neglogl_mipminusp",
                  "effpur_neglogl_mumipminusp",
                  "effpur_chi2muminusp",
                  "effpur_depErangeEmu",
                  "effpur_depErangeEp"
                };

// For efficiency/purity we need to know whether MIPs are supposed to be low or high. In the same order as the vector above
std::vector<bool> MIPlow = {
                    false, // track_neglogl_p
                    true, // track_neglogl_mu
                    true, // track_neglogl_pi
                    false, // track_neglogl_k
                    true, // track_neglogl_mip
                    true, // track_neglogl_minmumip
                    true, // track_chi2mu
                    false, // track_chi2p
                    true, // track_chi2pi
                    false, // track_chi2k
                    true, // track_PIDA_kde
                    true, // track_PIDA_mean
                    true, // track_neglogl_muminusp
                    true, // track_neglogl_mipminusp
                    true, // track_neglogl_minmumipminusp
                    true, // track_chi2_muminusp
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
  if (PIDvarstoplot_dummy.size() != bins.size()) std::cout << "WARNING PIDvarstoplot_dummy.size() = " << PIDvarstoplot_dummy.size() << "and bins.size() = " << bins.size() << ". This is going to cause you problems!" << std::endl;

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
