# GAMPixTools

This package is a collection of tools to simulate the detector response of liquid noble TPCs, being developed in support of measuring multiply scattering MeV gamma rays for GammaTPC, and now for the GAMPix readout of DUNE.  Both of these are single phase detectors.

There are currently two mostly separate area of functionality.
1. Readout of charge from high fidelity electron tracks simulated by PENELOPE in a geometry-less infinite fluid, with several different charge readout schemes implemented.   
2. Vectorized treatment of data sets of gamma ray events in a full G4 detector geometry simulated by the Cosima tool in the MEGALib package for gamma ray astrophysics.  These tools include a framework for treating in detail the production and readout of charge and light signals and recombination fluctuations, and how these affect energy resolution in the MeV energy range.  

Notionally, the results of (1) will be applied in a parameterized way to the tracks from gamma rays in (2), though a version where (1) is fully implemented in (2) is possible. 

Examples are provided as a means of documenting the code.

## Getting started on S3DF

S3DF's shared data space (`/sdf/group/neutrino`) contains a userspace with some previous studies and inputs with various primary particles and different momentum distributions.  Take a look in `/sdf/group/neutrino/gampix/FD_size_samples` for more details.  Each sub-directory contains a raw edep-sim output file (in Geant4's RooTracker format), and an hdf5 dump (labelled "dumpTree"), which is a bit simpler to access in python.

An example is included (`chargeEff.py`) which demonstrates how to load one of these files and simulate the GAMPix charge readout for each track in the file.  The example then varies the depth of the track with respect to the anode position and compares the measured charge hits as a function of this depth.

