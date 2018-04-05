// Header file for algorithm that calculates likelihood from comparison of data
// dEdx vs residual range to predicted Bragg peak
//
// Takes in dE/dx and residual range vectors, and a particle species, and returns
// -2logL, where the L is likelihood for the data to have been produced for that
// given particle species
//
// Likelihood is calculated by comparing the measured dE/dx at a given residual range
// to the expected dE/dx at that residual range. Assumes that dE/dx is distributed according to a Landau convoluted with a Gaussian
// with a peak value given by the theoreticl prediction in the Theory_dEdx_resrange class
// (calculated by Bruce Baller), and two widths that are user-configurable but have
// a default value of: landau_width = 0.09 for protons and 0.07 for muons/pions/kaons
// and gaus_width = 0.4 for protons and 0.1 for muons/pions/kaons.
//
// Can be run for the following particle species: muon, pion, Kaon, proton
//
// Tracks are fit in both directions (forwards and backwards) to account for reconstruction
// getting the track direction wrong, and the highest likelihood is returned

#ifndef BRAGG_NEGLOGL_H
#define BRAGG_NEGLOGL_H

#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"

#include <vector>
#include <iostream>

#include "Theory_dEdx_resrange.cxx"
#include "langaufun.h"

namespace particleid{

  class Bragg_negLogL_Estimator{

  public:
    double getNegLogL(std::vector<double> dEdx, std::vector<double> resRange, int particlehypothesis, bool forward, double *shift=nullptr);
    void setLanWidthMu(double lan_width){lan_width_mu = lan_width;}
    void setLanWidthPi(double lan_width){lan_width_pi = lan_width;}
    void setLanWidthK (double lan_width){lan_width_k  = lan_width;}
    void setLanWidthP (double lan_width){lan_width_p  = lan_width;}
    void setGausWidthMu(double gaus_width){gaus_width_mu = gaus_width;}
    void setGausWidthPi(double gaus_width){gaus_width_pi = gaus_width;}
    void setGausWidthK (double gaus_width){gaus_width_k  = gaus_width;}
    void setGausWidthP (double gaus_width){gaus_width_p  = gaus_width;}
    void setResRangeMaxShift (double rr_max_shift_input){rr_max_shift = rr_max_shift_input;}
    void setResRangeMinShift (double rr_min_shift_input){rr_min_shift = rr_min_shift_input;}

  private:
    double lan_width_mu  = 0.07;
    double lan_width_pi  = 0.07;
    double lan_width_k   = 0.07;
    double lan_width_p   = 0.09;
    double gaus_width_mu = 0.1;
    double gaus_width_pi = 0.1;
    double gaus_width_k  = 0.1;
    double gaus_width_p  = 0.4;
    double rr_max_shift  = 5.0;
    double rr_min_shift  = -5.0;
  };

}

#endif
