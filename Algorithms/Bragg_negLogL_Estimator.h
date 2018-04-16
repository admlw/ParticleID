// Implementation file for algorithm that calculates likelihood from comparison of data
// dEdx vs residual range to predicted Bragg peak
//
// Takes in dE/dx and residual range vectors, and a particle species, and returns
// -2logL, where the L is likelihood for the data to have been produced for that
// given particle species
//
// Likelihood is calculated by comparing the measured dE/dx at a given residual range
// to the expected dE/dx at that residual range. Assumes that dE/dx is Landau distributed
// with a mean given by the theoreticl prediction in the Theory_dEdx_resrange class
// (calculated by Bruce Baller), and a width that is user-configurable but has
// a default value of 0.2 for protons and 0.1 for muons/pions/kaons.
//
// Can be run for the following particle species: muon, pion, Kaon, proton
//
// Tracks are fit in both directions (forwards and backwards) to account for reconstruction
// getting the track direction wrong, and the highest likelihood is returned

#ifndef BRAGG_NEGLOGL_H
#define BRAGG_NEGLOGL_H

// cpp
#include <vector>
#include <iostream>

// ROOT
#include "TMath.h"
#include "TGraph.h"
#include "TF1.h"

// art
#include "fhiclcpp/ParameterSet.h"

// local
#include "Theory_dEdx_resrange.h"
#include "landauGaussian.h"

namespace particleid{

  class Bragg_negLogL_Estimator{

  public:
    void configure(fhicl::ParameterSet const &p);
    void printConfiguration();
    double getNegLogL(std::vector<double> dEdx, std::vector<double> resRange, int particlehypothesis, bool forward);

  private:
    double gausWidth_mu;
    double gausWidth_pi;
    double gausWidth_k;
    double gausWidth_p;
    double landauWidth_mu;
    double landauWidth_pi;
    double landauWidth_k;
    double landauWidth_p;
    int nHitsToDrop;
    double endPointFloatShort;
    double endPointFloatLong;
    double endPointFloatStepSize;

  };
  
}

#endif
