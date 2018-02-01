#ifndef FIDUCIALVOLUME_H
#define FIDUCIALVOLUME_H

#include "fhiclcpp/ParameterSet.h"

namespace fidvol{

  class fiducialVolume{

    public:
      
      void setFiducialVolume(double xl, double xh, double yl, double yh, double zl, double zh, fhicl::ParameterSet const & p);

      void printFiducialVolume(double xl, double xh, double yl, double yh, double zl, double zh);

  };

}

#endif

