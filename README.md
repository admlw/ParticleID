# Particle ID module

This package is intended to be used for analyses on MicroBooNE. It's currently under active development. The intention is to have many algorithms produce log-likelihood values for each particle type, and by combining these we will be able to infer the particle type.

## Dependencies

The producer module has a single dependency: lardataobj feature branch `feature/feature/kduffy_pidrefactor_v1_11_00_04`. This contains an extension to the anab::ParticleID class, in the form of a struct: 

```
struct sParticleIDAlgScores {

  std::string fAlgName;
  kVariableType fVariableType;
  int fAssumedPdg;
  float fValue;
  geo::PlaneID fPlaneID;

}

```

which holds the output of any generic PID algorithm.

## Algorithm details

The producer module only makes use of a single algorithm right now. Extension of this is currently being investigated.

__Bragg_negLogL_Estimator__
This algorithm uses reconstructed hits in the Bragg peak of each track, along with the dE/dx width of each particle species to estimate the likelihood of each data point belonging to each particle species, and then combines those values to get a single likelihood value for the track.

