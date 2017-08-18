C++11 Sphere decoder simulation for space-time lattice codes
============================================================

Key features
------------
- Basis matrix input
- Complex matrix support
- Spherical constellation (codebook shaping) support with radius estimator for desired codebook size
- Schnorr-Euchnerr smart implementation of the sphere decoder algorithm with additional optimizations
- Supports wiretap simulations (coset encoding)
- Comprehensive console (stdout), log.txt and csv simulation output
- Direct output plotting utilities via [gnuplot](http://www.gnuplot.info/)

Quick setup guide on Linux
--------------------------
- Download the whole repository (as a zip or using git clone) somewhere and open terminal there
- Install all the required packages (e.g. with apt-get) and the [Armadillo C++ linear algebra library](http://arma.sourceforge.net/download.html)
- If you wish to use plotting directly from the program you need to install C++ Boost library (e.g. *libboost1.58-all-dev* from the package manager)
- Otherwise comment out/remove the line that defines PLOTTING in the *src/misc.hpp* file
- Compile the program in terminal by typing: **make**
- Open the *settings.ini* in the */settings/* folder and change the variables there to setup the program, it should look something like this:

```ini
// configuration settings and simulation parameters for the sphere decoder program

basis_file=bases.txt            // Text file containing the basis matrices
output_file=                    // Optionally spesify the output filename
coset_file=                     // Optionally specify the coset encoding sublattice basis matrix file
error_file=                     // Optionally spesify the file containing error requirements for the SNR simulations.
channel_model=mimo              // Define the channel model for the simulation (either 'mimo' or 'siso')
x-PAM=4                         // The size of the PAM signaling set (even positive integer)
energy_estimation_samples=-1    // Number of samples to make the code energy estimation (-1 = sample all)
no_of_matrices=4                // Number of basis matrices (dimension of the data vectors)
matrix_coefficient=1.0          // Multiply all basis matrices by this constant
time_slots=2                    // Number of time slots used in the code
no_of_transmit_antennas=2       // Number of transmit antennas
no_of_receiver_antennas=2       // Number of receiver antennas
snr_min=-40                     // Minimum value for signal-to-noise ratio
snr_max=20                      // Maximum value for signal-to-noise ratio
snr_step=2                      // Increase SNR by this value per each iteration
simulation_rounds=1000000       // Number of simulation rounds to run
required_errors=-1              // Demand at minimum this many errors before the simulation ends
plot_results=1                  // Draw plots? (1 = yes, -1 = no)
stat_display_interval=100000    // Defines after each how many rounds to display the current simulation stats (-1 = disabled)
spherical_shaping_max_power=-1  // Defines the maximum distance from origin for codebook elements (-1 = unbounded)
codebook_size_exponent=-1       // The codebook will have 2^s codewords where s is this parameter (overrides above parameter)
radius_search_density=-1        // Defines how accurate the codebook squared radius estimation will be
```

- Create file called *bases.txt* in */bases/* folder and put your code basis matrices there (e.g. in Mathematica format)
- Run the program with: **./sphdec** or **make run** (the latter runs with the default settings file)
- If you configured the program correctly it should now run the simulation
- You can have multiple settings files (in the */settings/* folder) and use them in the simulation by giving their name as an command line argument for the program like: **./sphdec alamouti_settings.ini**
- Program output should be found at */output/* folder

(Naturally in order for this to work you need a C++ compiler (g++) that supports C++11 standard installed on your system, should be no problem on Aalto computers)