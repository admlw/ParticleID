// Implementation file for algorithm that calculates likelihood from comparison of data
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

#ifndef BRAGG_NEGLOGL_CXX
#define BRAGG_NEGLOGL_CXX

#include "Bragg_negLogL_Estimator_new.h"

namespace particleid{

  double Bragg_negLogL_Estimator::getNegLogL(std::vector<double> dEdx, std::vector<double> resRange, int particlehypothesis, bool forward, double *shift)
  {

    // Get theoretical prediction for given particle hypothesis
    // (This gives the MPV of a Landau distribution)
    // This is only available for muon, proton, kaon, and pion
    // Return an error if user tries to request a different particle type
    Theory_dEdx_resrange theory;
    TGraph *theorypred;
    double lan_width;
    double gaus_width;
    int absph = TMath::Abs(particlehypothesis);
    switch(absph){
      case 13: // muon
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to muon hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Muon;
        lan_width = lan_width_mu;
        gaus_width = gaus_width_mu;
        break;
      case 2212: // proton
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to proton hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Proton;
        lan_width = lan_width_p;
        gaus_width = gaus_width_p;
        break;
      case 211: // pion
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to pion hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Pion;
        lan_width = lan_width_pi;
        gaus_width = gaus_width_pi;
        break;
      case 321: // kaon
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to kaon hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Kaon;
        lan_width = lan_width_k;
        gaus_width = gaus_width_k;
        break;
      default:
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] ERROR: cannot calculate theoretical prediction for given particle hypothesis: " << particlehypothesis << ". Theoretical predictions are only available for charged muons (+/-13), pions (+/-211), kaons (+/-321), and protons (2212)" << std::endl;
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Exiting." << std::endl;
        throw;
    } // switch

    // Landau convoluted with a Gaussian. Has four parameters:
    // par[0]: Landau width (determined from fitting width in MC and data dEdx)
    // par[1]: Landau peak (determined by theory prediction)
    // par[2]: Normalisation constant (doesn't matter; set to 1)
    // par[3]: Gaussian width (determined from fitting width in MC and data dEdx)
    TF1 *langaus = new TF1("langaus",langaufun, 0, 100, 4);

    // Now loop through hits (entries in dEdx and resRange vectors), compare to
    // theoretical prediction, and calculate likelihood
    // rr_shift allows use to shift the residual range so that we
    // can account for end point resolution
    double minLLNdf = 9999999;
    int n_hits_used_total = 0;
    double shift_final = 0.;

    // Determine range that we allow the residual range to shift
    // Negative can't be more than the track length
    // Positive can be larger because it might just be not at the Bragg peak
    rr_min_shift = -1.0*std::min(TMath::Abs(rr_min_shift), *std::max_element(resRange.begin(), resRange.end()));

    std::cout << "rr_min_shift = " << rr_min_shift << std::endl;
    std::cout << "rr_max_shift = " << rr_max_shift << std::endl;

    for (double rr_shift = rr_min_shift; rr_shift < rr_max_shift; rr_shift = rr_shift+0.05){

      // Make neg2LogLikelihood
      double neg2LogL = 0.;
      int n_hits_used = 0;

      for (size_t i_hit=0; i_hit < resRange.size(); i_hit++){

        // ignore first and last hit
        if (i_hit == 0) continue;
        if (i_hit == resRange.size()-1) continue;

        size_t rr_index;
        if (forward){ // Fit tracks "forward" (i.e. in the direction they already have)
          rr_index = i_hit;
        }
        else{ // Fit tracks "backward"
          rr_index = (resRange.size()-1)-i_hit;
        }

        double resrg_i = resRange.at(rr_index)+rr_shift;
        double dEdx_i  = dEdx.at(i_hit);

        // Set theoretical Landau/Gaussian distribution for given residual range
        // par[0]: Landau width (determined from fitting width in MC and data dEdx)
        // par[1]: Landau peak (determined by theory prediction)
        // par[2]: Normalisation constant (doesn't matter; set to 1)
        // par[3]: Gaussian width (determined from fitting width in MC and data dEdx)
        langaus->SetParameters(lan_width,theorypred->Eval(resrg_i,0,"S"),1,gaus_width);

        // Evaluate likelihood
        double neg2LogL_i = 0.;
        if (langaus->Eval(dEdx_i) == 0){
          continue;
          neg2LogL_i = -2*std::log(1e-20);
        }
        else{
          neg2LogL_i = -2*std::log(langaus->Eval(dEdx_i));
          n_hits_used++;
        }

        neg2LogL += neg2LogL_i;
        /*
           std::cout << "Algorithm --- " << std::endl
           << "   theorypred->Eval(" << resrg_i << ",0,S) = " << theorypred->Eval(resrg_i,0,"S") << std::endl
           << "   theorypred->Eval(" << resrg_i << ") = " << theorypred->Eval(resrg_i) << std::endl
           << "   langaus->Eval(" << dEdx_i << ") = " << langaus->Eval(dEdx_i) << std::endl
           << "   neg2logL(i) = " << neg2LogL_i << std::endl
           << "   neg2logL total = " << neg2LogL << std::endl;
           */
      } // Loop over i_hit

      if (neg2LogL/n_hits_used < minLLNdf){
        minLLNdf = neg2LogL/n_hits_used;
        n_hits_used_total = n_hits_used;
        shift_final = rr_shift;
      }

    } // residual range shift

    // save shift
    if (shift)
      *shift = shift_final;

    if (n_hits_used_total == 0)
      return 9999999;
    else return minLLNdf;

  } // getneglogL

} // end particleid namespace

#endif
