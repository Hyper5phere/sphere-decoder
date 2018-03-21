Planewalker - C++ Sphere decoder simulator for space-time lattice codes
=======================================================================

Description
-----------
Planewalker is a sphere decoder implementation with diverse lattice code simulation capabilitites built around it. Its implementation was carried out as a special assignment by Pasi Pyrrö, advised by Oliver Gnilke, Marcus Greferath, and Camilla Hollanti, Department of Mathematics and Systems Analysis, Aalto University, Finland. 

How to cite
-----------
If you use or adapt this software in your research, please cite it as follows:
---
P. Pyrrö, O. Gnilke, C. Hollanti, and M. Greferath, “*Planewalker* sphere decoder implementation,” 2018. [Online]. Available: https://version.aalto.fi/gitlab/pasi.pyrro/sphere-decoder/
---
Or with Bibtex:
```tex
@misc{PGHG18,
    Author = {P. Pyrr\"o and O. Gnilke and C. Hollanti and M. Greferath},
    Title = {\emph{Planewalker} sphere decoder implementation},
    Url = {https://version.aalto.fi/gitlab/pasi.pyrro/sphere-decoder/},
    Year = {2018}
}
```
![example simulation results](https://version.aalto.fi/gitlab/pasi.pyrro/sphere-decoder/raw/master/docs/src/alamouti_cosets12345_8-PAM.png "example simulation results")

Key features
------------
- Comprehensive generator matrix and simulation parameters input (using text based configuration files)
- Complex matrix support
- Spherical constellation (codebook shaping) support with radius estimator for desired codebook size
- [Schnorr-Euchner smart implementation of the sphere decoder algorithm](https://doi.org/10.1109/TIT.2003.817444) with additional optimizations
- Native handling of q-PAM signaling sets
- Supports wiretap simulations (coset encoding)
- Comprehensive console (stdout), log.txt and csv simulation output
- Direct output plotting utilities via [gnuplot](http://www.gnuplot.info/)
- Possibility to program own custom simulations with the help of predefined, reusable functions (API)

Quick setup guide on Linux
--------------------------
- Download the whole repository (as a zip or using git clone) somewhere and open terminal there
- Install all the required packages (e.g. with apt-get) and the [Armadillo C++ linear algebra library](http://arma.sourceforge.net/download.html)
- If you wish to use plotting directly from the program you need to install C++ Boost library (e.g. *libboost1.58-all-dev* from the package manager)
- Compile the program with minimal features in terminal by typing: **make**
- Make also accepts argument string **with=plotting+gpu** where the value is a list of features to be installed (for more details see the documentation)
- Open the *settings.ini* in the */settings/* folder and change the variables there to setup the program, it should look something like this:

```ini
// configuration settings and simulation parameters for the sphere decoder program //

basis_file=alamouti.txt         // Text file containing the basis matrices (located in the /bases/ folder)
output_file=                    // Optionally specify the output csv filename (located in the /output/ folder)
coset_file=                     // Optionally specify the coset encoding sublattice basis matrix text file (located in the /bases/ folder)
error_file=                     // Optionally specify a csv file containing error requirements for the SNR simulations. (located in the /settings/ folder)
channel_model=mimo              // Define the channel model for the simulation (either 'mimo' or 'siso')
x-PAM=4                         // The size of the PAM signaling set (even positive integer)
energy_estimation_samples=-1    // Number of samples to make the code energy estimation (-1 = sample all)
no_of_matrices=4                // Number of basis matrices (dimension of the data vectors)
matrix_coefficient=1.0          // Multiply all basis matrices by this constant
time_slots=2                    // Number of time slots used in the code
no_of_transmit_antennas=2       // Number of transmit antennas
no_of_receiver_antennas=2       // Number of receiver antennas
snr_min=-6                      // Minimum value for signal-to-noise ratio
snr_max=20                      // Maximum value for signal-to-noise ratio
snr_step=2                      // Increase SNR by this value per each iteration
simulation_rounds=10000         // Number of simulation rounds to run
required_errors=-1              // Demand at minimum this many errors before the simulation ends
plot_results=-1                 // Draw plots? (1 = yes, -1 = no)
stat_display_interval=-1        // Defines after each how many rounds to display the current simulation stats (-1 = disabled)
spherical_shaping_max_power=-1  // Defines the maximum squared distance from origin for codebook elements (-1 = unbounded)
codebook_size_exponent=-1       // The codebook will have 2^s codewords where s is this parameter (overrides above parameter)
radius_search_density=100       // Defines how accurate the codebook squared radius estimation will be (shortest vector of generator matrix is divided by this)
```

- There should be a file called *alamouti.txt* in */bases/* folder that contains your basis matrices (for alamouti example code)
- You can edit this file or create a new one to fit your simulation needs 
- Basis matrices can be inputted in many formats, but Mathematica format (the one used in example basis files) is preferred for complex matrices (IMPORTANT: DO NOT USE WHITE SPACES AS SEPARATORS FOR MATRIX ELEMENTS)
- Run the program with: **./pwalk** or **make run** (the latter runs with the default settings file)
- If you configured the program correctly it should now run the simulation
- You can have multiple settings files (in the */settings/* folder) and use them in the simulation by giving their name as an command line argument for the program like: **./pwalk alamouti_settings.ini**
- Program output should be found at */output/* folder

(Naturally in order for this to work you need a C++ compiler (g++) that supports C++11 standard installed on your system, should be no problem on Aalto computers)

Known issues
------------
- Random generation of spherical codewords seems to have some selection bias, that shifts the average energy of the codebook slightly
- Spherical shaping radius estimation does not work 100% of the time (used with *codebook_size_exponent* parameter)
- There is a weird issue with openmp parallel computing library that causes the initial radius for the sphere decoder sometimes to be too small when using spherical shaping (quick workarounds: either disable parallelism or use infinite initial radius)

License
-------
Copyright (c) 2018 Pasi Pyrrö

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
