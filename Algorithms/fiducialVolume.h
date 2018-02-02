#ifndef FIDUCIALVOLUME_H
#define FIDUCIALVOLUME_H

#include "fhiclcpp/ParameterSet.h"
#include "TVector3.h"

namespace fidvol{

  class fiducialVolume{

    public:
      
      void setFiducialVolume(double xl, double xh, double yl, double yh, double zl, double zh, std::vector<double> fv, fhicl::ParameterSet const & p);

      void printFiducialVolume(double xl, double xh, double yl, double yh, double zl, double zh);

      bool isInFiducialVolume(TVector3 xyz, std::vector<double> fv);

      std::vector<double> getTpcDimensions();

  };

}

#endif

