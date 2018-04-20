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

  void Bragg_negLogL_Estimator::configure(fhicl::ParameterSet const &p){

    gausWidth_p  = p.get<double>("dEdxGausWidthP"  , 0.4);
    gausWidth_mu = p.get<double>("dEdxGausWidthMu" , 0.1);
    gausWidth_pi = p.get<double>("dEdxGausWidthPi" , 0.1);
    gausWidth_k  = p.get<double>("dEdxGausWidthK"  , 0.1);
    landauWidth_p  = p.get<double>("dEdxLandauWidthP"  , 0.09);
    landauWidth_mu = p.get<double>("dEdxLandauWidthMu" , 0.07);
    landauWidth_pi = p.get<double>("dEdxLandauWidthPi" , 0.07);
    landauWidth_k  = p.get<double>("dEdxLandauWidthK"  , 0.07);
    nHitsToDrop    = p.get<int>("NHitsToDrop", 1);
    endPointFloatShort    = p.get<double>("EndPointFloatShort", -1.0);
    endPointFloatLong     = p.get<double>("EndPointFloatLong" , 1.0);
    endPointFloatStepSize = p.get<double>("EndPointFloatStepSize", 0.05);

  }

  void Bragg_negLogL_Estimator::printConfiguration(){

    std::cout << "[ParticleID::Bragg_negLogL_Estimator] PRINTING CONFIGURATION: " << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Proton dE/dx gaus width : " << gausWidth_p  << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Muon dE/dx gaus width   : " << gausWidth_mu << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Pion dE/dx gaus width   : " << gausWidth_pi << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Kaon dE/dx gaus width   : " << gausWidth_k  << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Proton dE/dx landau width : " << landauWidth_p  << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Muon dE/dx landau width   : " << landauWidth_mu << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Pion dE/dx landau width   : " << landauWidth_pi << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Kaon dE/dx landau width   : " << landauWidth_k  << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] Number of Hits to Drop: " << nHitsToDrop << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] End-point float long  : " << endPointFloatLong  << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] End-point float short : " << endPointFloatShort  << std::endl;
    std::cout << "[ParticleID::Bragg_negLogL_Estimator] End-point step size   : " << endPointFloatStepSize  << std::endl;

  }

  double Bragg_negLogL_Estimator::getNegLogL(std::vector<double> dEdx, std::vector<double> resRange, int particlehypothesis, bool forward)
  {

    /**
     * Get theoretical prediction for given particle hypothesis
     * (This gives the MPV of a Landau distribution)
     * This is only available for muon, proton, kaon, and pion
     * Return an error if user tries to request a different particle type
     */
    
    Theory_dEdx_resrange theory;
    TGraph *theorypred;
    double gausWidth;
    double landauWidth;
    int absph = TMath::Abs(particlehypothesis);
   
    switch(absph){
      case 13: // muon
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to muon hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Muon;
        gausWidth = gausWidth_mu;
        landauWidth = landauWidth_mu;
        break;
      case 2212: // proton
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to proton hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Proton;
        gausWidth = gausWidth_p;
        landauWidth = landauWidth_p;
        break;
      case 211: // pion
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to pion hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Pion;
        gausWidth = gausWidth_pi;
        landauWidth = landauWidth_pi;
        break;
      case 321: // kaon
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood with respect to kaon hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Kaon;
        gausWidth = gausWidth_k;
        landauWidth = landauWidth_k;
        break;
      case 0: // special case: fit to MIP region of muon prediction with no Bragg peak
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Calculating likelihood for non-Bragg MIP-like hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_MuonNoBragg;
        gausWidth = gausWidth_mu;
        landauWidth = landauWidth_mu;
        break;
      default:
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] ERROR: cannot calculate theoretical prediction for given particle hypothesis: " << particlehypothesis << ". Theoretical predictions are only available for charged muons (+/-13), pions (+/-211), kaons (+/-321), protons (2212), and non-Bragg MIP region (0)" << std::endl;
        std::cout << "[ParticleID::Bragg_negLogL_Estimator] Exiting." << std::endl;
        throw;
    } // switch

    TF1 *langaus = new TF1("langaus", landauGaussian, 0, 100, 4);

    /**
     * Now loop through hits (entries in dEdx and resRange vectors), compare to
     * theoretical prediction, and calculate likelihood
     * rr_shift allows us to shift the residual range so that we
     * can account for end point resolution
     */
    
    double minLLNdf = 9999999;
    int n_hits_used_total = 0;

    for (double rr_shift = endPointFloatShort; rr_shift < endPointFloatLong; rr_shift = rr_shift+endPointFloatStepSize){

      // Make neg2LogLikelihood
      double neg2LogL = 0.;
      int n_hits_used = 0;

      for (size_t i_hit=0; i_hit < resRange.size(); i_hit++){

        size_t rr_index;
        // Fit tracks "forward" (i.e. in the direction they already have)
        if (forward){           
          rr_index = i_hit;
          if ((int)rr_index >= (int)resRange.size() - nHitsToDrop){
            continue;
          }
        }
        // Fit tracks "backward"
        else{          
          rr_index = (resRange.size()-1)-i_hit;
          if ((int)i_hit < nHitsToDrop){
            continue;
          }
        }

        double resrg_i = resRange.at(rr_index)+rr_shift;
        double dEdx_i  = dEdx.at(i_hit);

        /**
         * Theory values are only defined up to 30 cm residual Range so we
         * can't compare beyond that
         */
        if (resrg_i > 30.0) continue;

        // Set theoretical Landau distribution for given residual range
        langaus->SetParameters(landauWidth,theorypred->Eval(resrg_i,0,"S"),1, gausWidth);

        // Evaluate likelihood
        double neg2LogL_i = 0.;
        if (langaus->Eval(dEdx_i) == 0){
          continue;
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
      }

    } // residual range shift

    if (n_hits_used_total == 0)
      return 9999999;
    else return minLLNdf;

  } // getneglogL


} // particleid namespace

#endif
