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

#ifndef BRAGG_NEGLOGL_CXX
#define BRAGG_NEGLOGL_CXX

#include "Bragg_negLogL_Estimator.h"

namespace particleid{

  double Bragg_negLogL_Estimator::getNegLogL(std::vector<double> dEdx, std::vector<double> resRange, int particlehypothesis, bool forward)
  {

    // Get theoretical prediction for given particle hypothesis
    // (This gives the MPV of a Landau distribution)
    // This is only available for muon, proton, kaon, and pion
    // Return an error if user tries to request a different particle type
    Theory_dEdx_resrange theory;
    TGraph *theorypred;
    double width;
    int absph = TMath::Abs(particlehypothesis);
    switch(absph){
      case 13: // muon
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to muon hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Muon;
        width = width_mu;
        break;
      case 2212: // proton
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to proton hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Proton;
        width = width_p;
        break;
      case 211: // pion
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to pion hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Pion;
        width = width_pi;
        break;
      case 321: // kaon
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to kaon hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Kaon;
        width = width_k;
        break;
      default:
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] ERROR: cannot calculate theoretical prediction for given particle hypothesis: " << particlehypothesis << ". Theoretical predictions are only available for charged muons (+/-13), pions (+/-211), kaons (+/-321), and protons (2212)" << std::endl;
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Exiting." << std::endl;
        throw;
    } // switch

    TF1 *landau = new TF1("landau", "landau", 0, 100);

    // Make neg2LogLikelihood
    double neg2LogL = 0.;

    // Now loop through hits (entries in dEdx and resRange vectors), compare to
    // theoretical prediction, and calculate likelihood
    int n_hits_used = 0;
    for (size_t i_hit=0; i_hit < resRange.size(); i_hit++){

      size_t rr_index;
      if (forward){ // Fit tracks "forward" (i.e. in the direction they already have)
        rr_index = i_hit;
      }
      else{ // Fit tracks "backward"
        rr_index = (resRange.size()-1)-i_hit;
      }

      double resrg_i = resRange.at(rr_index);
      double dEdx_i  = dEdx.at(i_hit);

      // Set theoretical Landau distribution for given residual range
      landau->SetParameters(1,theorypred->Eval(resrg_i,0,"S"),width);

      // Evaluate likelihood
      double neg2LogL_i = 0.;
      if (landau->Eval(dEdx_i) == 0){
        continue;
        neg2LogL_i = -2*std::log(1e-20);
      }
      else{
        neg2LogL_i = -2*std::log(landau->Eval(dEdx_i));
        n_hits_used++;
      }

      neg2LogL += neg2LogL_i;
/*
      std::cout << "Algorithm --- " << std::endl
        << "   theorypred->Eval(" << resrg_i << ",0,S) = " << theorypred->Eval(resrg_i,0,"S") << std::endl
        << "   theorypred->Eval(" << resrg_i << ") = " << theorypred->Eval(resrg_i) << std::endl
        << "   landau->Eval(" << dEdx_i << ") = " << landau->Eval(dEdx_i) << std::endl
        << "   neg2logL(i) = " << neg2LogL_i << std::endl
        << "   neg2logL total = " << neg2LogL << std::endl;
*/
    } // Loop over i_hit

    if (n_hits_used == 0)
      return 9999999;
    else return neg2LogL/n_hits_used;

  } // getneglogL


} // particleid namespace

#endif