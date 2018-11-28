#ifndef PIDA_H
#define PIDA_H

#include "ubana/ParticleID/Algorithms/KernelDensityEstimator.h"
#include "TMath.h"
#include <vector>

namespace particleid{

  class PIDA{

    public:
      
<<<<<<< HEAD
      double getPida(std::vector<double> dEdx, std::vector<double> resRange, std::string method);
=======
      double getPida(std::vector<float> dEdx, std::vector<float> resRange, std::string method);
>>>>>>> temp

  };

}

#endif
