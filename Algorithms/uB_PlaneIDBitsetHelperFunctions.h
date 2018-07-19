// MicroBooNE-specific helper functions to translate bitset planeIDs used in anab::ParticleID

#include <bitset>

namespace UBPID{

  inline int uB_getSinglePlane(std::bitset<5> planeid){
    if ((planeid[0]+planeid[1]+planeid[2]+planeid[3]+planeid[4]>1) || (planeid[0]+planeid[1]+planeid[2]+planeid[3]+planeid[4]==0) || planeid[3]==1 || planeid[4]==1){
      std::cout << "[ParticleID::PlaneIDBitsetHelperFunctions] Cannot return a single MicroBooNE plane for bitset " << planeid << ". Returning -1 (invalid planeID). This may cause you problems!" << std::endl;
      return -1;
    }
    else if (planeid[0]==1) return 2;
    else if (planeid[1]==1) return 1;
    else if (planeid[2]==1) return 0;

    // Default: invalid return
    return -1;
  }

  inline std::bitset<5> uB_SinglePlaneGetBitset(int planeid){
    if (planeid<0 || planeid>2){
      std::cout << "[ParticleID::PlaneIDBitsetHelperFunctions] Cannot return a bitset for MicroBooNE planeid " << planeid << ". Returning 00000 (invalid planeID). This may cause you problems!" << std::endl;
      return std::bitset<5>("00000");
    }
    else if (planeid==2) return std::bitset<5>("10000");
    else if (planeid==1) return std::bitset<5>("01000");
    else if (planeid==0) return std::bitset<5>("00100");

    // Default: invalid return
    return std::bitset<5>("00000");
  }

}
