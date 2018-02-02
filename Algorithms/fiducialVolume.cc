#include "fiducialVolume.h"

namespace fidvol{

  void fiducialVolume::setFiducialVolume(double xl, double xh, double yl, double yh, double zl, double zh, std::vector<double> fv, fhicl::ParameterSet const & p){

    xl = p.get< double > ("X_LOW");
    xh = p.get< double > ("X_HIGH");
    yl = p.get< double > ("Y_LOW");
    yh = p.get< double > ("Y_HIGH");
    zl = p.get< double > ("Z_LOW");
    zh = p.get< double > ("Z_HIGH");

    fv.push_back(xl);
    fv.push_back(xh);
    fv.push_back(yl);
    fv.push_back(yh);
    fv.push_back(zl);
    fv.push_back(zh);

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

  bool fiducialVolume::isInFiducialVolume(TVector3 vec, std::vector<double> fv){

    double x = vec.X();
    double y = vec.Y();
    double z = vec.Z();

    std::vector<double> tpcd = fiducialVolume::getTpcDimensions();

    if (x > (tpcd.at(0)+fv.at(0)) && x < (tpcd.at(1)-fv.at(1)) && 
        y > (tpcd.at(2)+fv.at(2)) && y < (tpcd.at(3)-fv.at(3)) && 
        z > (tpcd.at(4)+fv.at(4)) && z < (tpcd.at(5)-fv.at(5))) 
      return true;

    else return false;

  }

  std::vector<double> fiducialVolume::getTpcDimensions(){

    std::vector<double> tpcd;
    double tpc_xl = 0.0;
    double tpc_xh = 256.0;
    double tpc_yl = -116.5;
    double tpc_yh = 116.5;
    double tpc_zl = 0.0;
    double tpc_zh = 1036.0;

    tpcd.push_back(tpc_xl);
    tpcd.push_back(tpc_xh);
    tpcd.push_back(tpc_yl);
    tpcd.push_back(tpc_yh);
    tpcd.push_back(tpc_zl);
    tpcd.push_back(tpc_zh);

    return tpcd;
  }

}
