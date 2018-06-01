// Implementation file for algorithm that calculates likelihood from comparison of data
// dEdx vs residual range to predicted Bragg peak
//
// Takes in dE/dx and residual range vectors, and a particle species, and returns
// the likelihood, L, for the data to have been produced for that
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

#ifndef BRAGG_LIKELIHOOD_CXX
#define BRAGG_LIKELIHOOD_CXX

#include "Bragg_Likelihood_Estimator.h"

namespace particleid{

  void Bragg_Likelihood_Estimator::configure(fhicl::ParameterSet const &p){

    gausWidth_p   = p.get<std::vector<double>>("dEdxGausWidthP"  );
    gausWidth_mu  = p.get<std::vector<double>>("dEdxGausWidthMu" );
    gausWidth_pi  = p.get<std::vector<double>>("dEdxGausWidthPi" );
    gausWidth_k   = p.get<std::vector<double>>("dEdxGausWidthK"  );
    gausWidth_mip = p.get<std::vector<double>>("dEdxGausWidthMIP");
    landauWidth_p   = p.get<std::vector<double>>("dEdxLandauWidthP"  );
    landauWidth_mu  = p.get<std::vector<double>>("dEdxLandauWidthMu" );
    landauWidth_pi  = p.get<std::vector<double>>("dEdxLandauWidthPi" );
    landauWidth_k   = p.get<std::vector<double>>("dEdxLandauWidthK"  );
    landauWidth_mip = p.get<std::vector<double>>("dEdxLandauWidthMIP");

    offset_p       = p.get<double>("PeakOffsetP"  , 0);
    offset_mu      = p.get<double>("PeakOffsetMu" , 0);
    offset_pi      = p.get<double>("PeakOffsetPi" , 0);
    offset_k       = p.get<double>("PeakOffsetK"  , 0);
    offset_mip     = p.get<double>("PeakOffsetMIP", 0);

    nHitsToDrop    = p.get<int>("NHitsToDrop", 1);
    endPointFloatShort    = p.get<double>("EndPointFloatShort", -1.0);
    endPointFloatLong     = p.get<double>("EndPointFloatLong" , 1.0);
    endPointFloatStepSize = p.get<double>("EndPointFloatStepSize", 0.05);

  }

  void Bragg_Likelihood_Estimator::printConfiguration(){

    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] PRINTING CONFIGURATION: " << std::endl;
    for (int i = 0; i < 3; i++){
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Proton dE/dx gaus width : " << gausWidth_p.at(i)  << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Muon dE/dx gaus width   : " << gausWidth_mu.at(i) << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Pion dE/dx gaus width   : " << gausWidth_pi.at(i) << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Kaon dE/dx gaus width   : " << gausWidth_k.at(i)  << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Kaon dE/dx gaus width   : " << gausWidth_mip.at(i)  << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Proton dE/dx landau width : " << landauWidth_p.at(i)  << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Muon dE/dx landau width   : " << landauWidth_mu.at(i) << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Pion dE/dx landau width   : " << landauWidth_pi.at(i) << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " Kaon dE/dx landau width   : " << landauWidth_k.at(i)  << std::endl;
      std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Plane " << i << " MIP dE/dx landau width    : " << landauWidth_mip.at(i) << std::endl;
    }
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Proton MPV Offset : " << offset_p   << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Muon MPV Offset   : " << offset_mu  << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Pion MPV Offset   : " << offset_pi  << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Kaon MPV Offset   : " << offset_k   << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> MIP MPV Offset    : " << offset_mip << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> Number of Hits to Drop : " << nHitsToDrop << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> End-point float long   : " << endPointFloatLong  << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> End-point float short  : " << endPointFloatShort  << std::endl;
    std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> End-point step size    : " << endPointFloatStepSize  << std::endl;

  }

  double Bragg_Likelihood_Estimator::getLikelihood(std::vector<double> dEdx, std::vector<double> resRange, int particlehypothesis, bool forward, int planenum)
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
    double offset;
    int absph = TMath::Abs(particlehypothesis);

    switch(absph){
      case 13: // muon
        //std::cout << "[ParticleID::Bragg_Likelihood_Estimator] Calculating likelihood with respect to muon hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Muon;
        gausWidth = gausWidth_mu.at(planenum);
        landauWidth = landauWidth_mu.at(planenum);
        offset = offset_mu;
        break;
      case 2212: // proton
        //std::cout << "[ParticleID::Bragg_Likelihood_Estimator] Calculating likelihood with respect to proton hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Proton;
        gausWidth = gausWidth_p.at(planenum);
        landauWidth = landauWidth_p.at(planenum);
        offset = offset_p;
        break;
      case 211: // pion
        //std::cout << "[ParticleID::Bragg_Likelihood_Estimator] Calculating likelihood with respect to pion hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Pion;
        gausWidth = gausWidth_pi.at(planenum);
        landauWidth = landauWidth_pi.at(planenum);
        offset = offset_pi;
        break;
      case 321: // kaon
        //std::cout << "[ParticleID::Bragg_Likelihood_Estimator] Calculating likelihood with respect to kaon hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_Kaon;
        gausWidth = gausWidth_k.at(planenum);
        landauWidth = landauWidth_k.at(planenum);
        offset = offset_k;
        break;
      case 0: // special case: fit to MIP region of muon prediction with no Bragg peak
        //std::cout << "[ParticleID::Bragg_Likelihood_Estimator] Calculating likelihood for non-Bragg MIP-like hypothesis" << std::endl;
        theorypred = theory.g_ThdEdxRR_MuonNoBragg;
        gausWidth = gausWidth_mip.at(planenum);
        landauWidth = landauWidth_mip.at(planenum);
        offset = offset_mip;
        break;
      default:
        std::cout << "[ParticleID::Bragg_Likelihood_Estimator] ERROR: cannot calculate theoretical prediction for given particle hypothesis: " << particlehypothesis << ". Theoretical predictions are only available for charged muons (+/-13), pions (+/-211), kaons (+/-321), protons (2212), and non-Bragg MIP region (0)" << std::endl;
        std::cout << "[ParticleID::Bragg_Likelihood_Estimator] Exiting." << std::endl;
        throw;
    } // switch

    TF1 *langaus = new TF1("langaus", landauGaussian, 0, 100, 4);

    // create langaus of correct width, and calculate offset between MPV and mean
    langaus->SetParameters(landauWidth, 10, 1, gausWidth);
    double langaus_mean = langaus->Mean(0, 100);
    double langaus_mpv  = langaus->GetMaximumX();
    double langaus_mean_mpv_offset = langaus_mean - langaus_mpv;

    /*
       std::cout << "[ParticleID::Bragg_Likelihood_Estimator] SHIFT SUMMARY" << std::endl;
       std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> langaus mean : " << langaus_mean << std::endl;
       std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> langaus mpv  : " << langaus_mpv << std::endl;
       std::cout << "[ParticleID::Bragg_Likelihood_Estimator] >> offset value : " << langaus_mean_mpv_offset << std::endl;
       */
    /**
     * Now loop through hits (entries in dEdx and resRange vectors), compare to
     * theoretical prediction, and calculate likelihood
     * rr_shift allows us to shift the residual range so that we
     * can account for end point resolution
     */

    double likelihoodNdf = 0;
    int n_hits_used_total = 0;

    for (double rr_shift = endPointFloatShort; rr_shift < endPointFloatLong; rr_shift = rr_shift+endPointFloatStepSize){

      // Make likelihoodikelihood
      double likelihood = 0.;
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
        langaus->SetParameters(landauWidth,theorypred->Eval(resrg_i,0,"S")-langaus_mean_mpv_offset+offset,1, gausWidth);

        // Evaluate likelihood
        double likelihood_i = 0.;
        if (langaus->Eval(dEdx_i) == 0){
          continue;
        }
        else{
          likelihood_i = langaus->Eval(dEdx_i);
          n_hits_used++;
        }

        likelihood += likelihood_i;
        /*
           std::cout << "Algorithm --- " << std::endl
           << "   theorypred->Eval(" << resrg_i << ",0,S) = " << theorypred->Eval(resrg_i,0,"S") << std::endl
           << "   theorypred->Eval(" << resrg_i << ") = " << theorypred->Eval(resrg_i) << std::endl
           << "   langaus->Eval(" << dEdx_i << ") = " << langaus->Eval(dEdx_i) << std::endl
           << "   Likelihood(i) = " << likelihood_i << std::endl
           << "   Likelihood total = " << likelihood << std::endl;
           */
      } // Loop over i_hit

      if (likelihood/n_hits_used > likelihoodNdf){
        likelihoodNdf = likelihood/n_hits_used;
        n_hits_used_total = n_hits_used;
      }

    } // residual range shift

    if (n_hits_used_total == 0)
      return -9999999;
    else return likelihoodNdf;

  } // getLikelihood


} // particleid namespace

#endif
