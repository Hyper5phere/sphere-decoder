## TODO-list for the Report and Documentation

### Input files
##### Settings file examples (7 in total)
- SISO Alamouti 8-PAM ✔
- MIMO Alamouti 4-PAM
- MIMO MIDO ✔
- SISO Golden code 8-PAM or 16-PAM
- MIMO Golden code 8-PAM ✔

##### Input matrix examples
- Alamouti ✔
- Golden code ✔
- MIDO ✔
- Wiretap examples

### Third chapter
- Finish mathematical background (how sphdec works)

##### Sphere decoder implementation
- explain use of finite signal sets (restricted to PAM)
- accepts complex matrices as input
- complex matrices converted to equivalent real representation internally

### Fourth chapter
- Remove technical information, add graphs and compare simulation results

##### Program implementation
- simulation implements siso and mimo channels
- fixes noise variance to 1
- simulates different SNRs by varying fading coefficients of the channel matrix

##### Simulation examples
- Compare orthonormal golden code (SISO) to Z^8
- run two simulations, plot and compare results
- Graphs for MIMO and Wiretap (no comparison)
- Wiretap: Bob's (base) lattice Z^4, index 16 (determinant of sublattice), integer entries.
- describe lattices and channel model

### Conclusion
- repeat yourself
- this program implements efficient sphdec
- written in C++, is fast, uses linear algebra library
- general importantness for lattice code research
- what could it be used for
- easy to modify and implement new channel models
- test new lattice codes and compare with old ones
- facilitates threoretical simulations to evaluate new code constructions

### Program User Guide
- leave LLL reduction documentation out (doesn't work properly)
- include settings examples with basis matrix file embedded as TEXT
- leave comments out of settings file and explain them in detail (refer to line numbers?)
- input examples: ⋅Alamouti, SISO, Wiretap, Spherical shaping


