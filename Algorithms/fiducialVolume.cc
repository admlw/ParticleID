#include "fiducialVolume.h"

namespace fidvol{

  void fiducialVolume::setFiducialVolume(double xl, double xh, double yl, double yh, double zl, double zh, fhicl::ParameterSet const & p){

    xl = p.get< double > ("X_LOW");
    xh = p.get< double > ("X_HIGH");
    yl = p.get< double > ("Y_LOW");
    yh = p.get< double > ("Y_HIGH");
    zl = p.get< double > ("Z_LOW");
    zh = p.get< double > ("Z_HIGH");

  }

  void fiducialVolume::printFiducialVolume(double xl, double xh, double yl, double yh, double zl, double zh){

    std::cout << "----- PRINTING FIDUCIAL VOLUME INFORMATION -----" << std::endl;
    std::cout << "X_LOW: " << xl << std::endl;
    std::cout << "X_HIGH: " << xh << std::endl;
    std::cout << "Y_LOW: " << yl << std::endl;
    std::cout << "Y_HIGH: " << yh << std::endl;
    std::cout << "Z_LOW: " << zl << std::endl;
    std::cout << "Z_HIGH: " << zh << std::endl;
    std::cout << "------------------------------------------------" << std::endl;

  }

}
