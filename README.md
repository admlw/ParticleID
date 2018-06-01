# Particle ID module

This package is intended to be used for analyses on MicroBooNE. It's currently under active development. The intention of this framework is to centralise the PID efforts, and make it easy to create and implement additional algorithms in an experiment agnostic way. In addition to this, we have implemented a number of new particle ID algorithms in addition to the older, more tested algorithms.

## Algorithm details

__Bragg_Likelihood Estimator__
This algorithm uses B. Baller's theory, along with reconstructed hits in the Bragg peak of each track to estimate the likelihood of each data point belonging to each particle species, and then combines those values to get a single likelihood value for the track.

This can be thought of as an extension to the chi-square method. The underlying dE/dx distribution is assumed to be a Landau-Gaussian with widths measured in simulation and data rather than the Gaussian assumption of the chi-square method. There are a number of other additions which are expected to help give better separation.

__PIDA__
This makes use of B. Baller's home-brewed PIDA method which was used to great success on ArgoNEUT. There are thee implementations available here:
1. Mean (B. Baller)
2. Median (T. Yang + V. Meddage)
3. Kernel density estimator (A. Lister)

__Truncated Mean dE/dx Versus Track Length__
Makes use of D. Caratelli's "Truncated Mean" algorithm to plot the mean dE/dx of a track versus its length.

__Deposited Energy vs Energy By Range__
Calculates deposited energy, and compares against the energy by range under different assumptions for separation.

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

which holds the output of any generic PID algorithm. Here, `kVariableType` is an enum which can take the following values:

```
enum kVariableType{
  kGOF,
  kLikelihood,
  kLikelihood_fwd,
  kLikelihood_bwd,
  kLogL,
  kLogL_fwd,
  kLogL_bwd,
  kScore,
  kPIDA,
  kdEdxtruncmean,
  kdQdxtruncmean,
  kTrackLength,
  kEdeposited,
  kEbyRange,
  kNotSet
}
```


