#include "PIDA.h"

namespace particleid{

  double PIDA::getPida(std::vector<double> dEdx, std::vector<double> resRange, std::string method){

    kde::KernelDensityEstimator kde;
    double pida = -1.0;

    std::vector<double> pidaValues;
    for (size_t i = 0; i < resRange.size(); i++){

      pidaValues.push_back(dEdx.at(i)*std::pow(resRange.at(i),0.42));

    }

    if (method == "mean") 
      pida = TMath::Mean(pidaValues.begin(), pidaValues.end());

    else if (method == "median")
      pida = TMath::Median(pidaValues.size(), &pidaValues[0]);

    else if (method == "kde")
      pida = kde.getKernelDensityMpv(pidaValues);
  
    return pida;

  }

}
