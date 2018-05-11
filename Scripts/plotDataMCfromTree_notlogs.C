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
                      vars->track_Lmumippi_0to1
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
                    {50,0,3} // track_Lmumippi_0to1
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
                    ";(L_{#mu}+L_{MIP}+L_{#pi})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});"
                  };

// What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
                  "Lp",
                  "Lmu",
                  "Lpi",
                  "Lk",
                  "Lmip",
                  "Lmumip",
                  "chi2mu",
                  "chi2p",
                  "chi2pi",
                  "chi2k",
                  "pida_kde",
                  "pida_mean",
                  "Lmuoverp",
                  "Lmipoverp",
                  "Lmumipoverp",
                  "chi2muminusp",
                  "Lmu0to1",
                  "Lmip0to1",
                  "Lpi0to1",
                  "Lp0to1",
                  "Lmumip0to1",
                  "Lmumippi0to1"
                };

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void plotDataMCFromTree(std::string mcfile, double POTscaling=0., std::string onbeamdatafile="", std::string offbeamdatafile="", double offbeamscaling=0.){

  TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
  TTree *t_bnbcos = (TTree*)f_bnbcos->Get("pidvalid/pidTree");
  treevars mc_vars;
  settreevars(t_bnbcos,&mc_vars);

  TFile *f_onbeam=nullptr;
  TTree *t_onbeam=nullptr;
  treevars onbeam_vars;
  if (onbeamdatafile!=""){
    std::cout << "Making data-MC comparisons" << std::endl;
    f_onbeam = new TFile(onbeamdatafile.c_str(), "read");
    t_onbeam = (TTree*)f_onbeam->Get("pidvalid/pidTree");
    settreevars(t_onbeam,&onbeam_vars);
  }

  TFile *f_offbeam=nullptr;
  TTree *t_offbeam=nullptr;
  treevars offbeam_vars;
  if (offbeamdatafile!=""){
    f_offbeam = new TFile(offbeamdatafile.c_str(), "read");
    t_offbeam = (TTree*)f_offbeam->Get("pidvalid/pidTree");
    settreevars(t_offbeam,&offbeam_vars);
  }


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

  // ----------------- On-beam data
  hist1D *onb_hists[nplots];
  if (t_onbeam){
    // Make histograms to fill
    for (int i_h=0; i_h<nplots; i_h++){
      onb_hists[i_h] = new hist1D(std::string("h_ondat_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
    }

    // Loop through MC tree and fill plots
    for (int i = 0; i < t_onbeam->GetEntries(); i++){
      t_onbeam->GetEntry(i);
      CalcPIDvars(&onbeam_vars);
      std::vector<double> PIDvarstoplot = GetPIDvarstoplot(&onbeam_vars);

      for (size_t i_h = 0; i_h < nplots; i_h++){
        FillHist(onb_hists[i_h],PIDvarstoplot.at(i_h),0); // 0 because there is no "true PDG" for data
      }
    }
  }

    // ----------------- Off-beam data
  hist1D *offb_hists[nplots];
  if (t_offbeam){
    // Make histograms to fill
    for (int i_h=0; i_h<nplots; i_h++){
      offb_hists[i_h] = new hist1D(std::string("h_offdat_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
    }

    // Loop through tree and fill plots
    for (int i = 0; i < t_offbeam->GetEntries(); i++){
      t_offbeam->GetEntry(i);
      CalcPIDvars(&offbeam_vars);
      std::vector<double> PIDvarstoplot = GetPIDvarstoplot(&offbeam_vars);

      for (size_t i_h = 0; i_h < nplots; i_h++){
        FillHist(offb_hists[i_h],PIDvarstoplot.at(i_h),0.); // 0 because there is no "true PDG" for data
      }
    } // end loop over entries in tree
  }


  // -------------------- Now make all the plots

  for (size_t i_h=0; i_h < nplots; i_h++){
    TCanvas *c1 = new TCanvas();
    DrawMC(mc_hists[i_h],POTscaling);
    if (f_onbeam && f_offbeam){
      OverlayData(c1,onb_hists[i_h],offb_hists[i_h],offbeamscaling,POTscaling);
    }
    c1->Print(std::string(histnames[i_h]+".png").c_str());
  }
}
