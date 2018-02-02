#ifndef KERNELDENSITYESTIMATOR_H
#define KERNELDENSITYESTIMATOR_H

#include <vector>

namespace kde{

  class KernelDensityEstimator{

    public:

      double getKernelDensityMpv(std::vector<double> pidaVals);

  };

}

#endif
