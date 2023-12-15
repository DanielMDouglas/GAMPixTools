detectormodel ReadMe
updated 5/10/23.  TS.

This package is a collection of tools to simulate the detector response of liquid noble TPCs, being developed in support of measuring multiply scattering MeV gamma rays for GammaTPC, and now for the GAMPix readout of DUNE.  Both of these are single phase detectors.

There are currently two mostly separate area of functionality.
1. Readout of charge from high fidelity electron tracks simulated by PENELOPE in a geometry-less infinite fluid, with several different charge readout schemes implemented.   
2. Vectorized treatment of data sets of gamma ray events in a full G4 detector geometry simulated by the Cosima tool in the MEGALib package for gamma ray astrophysics.  These tools include a framework for treating in detail the production and readout of charge and light signals and recombination fluctuations, and how these affect energy resolution in the MeV energy range.  

Notionally, the results of (1) will be applied in a parameterized way to the tracks from gamma rays in (2), though a version where (1) is fully implemented in (2) is possible. 

Examples are provided as a means of documenting the code.