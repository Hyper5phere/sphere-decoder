C++14 Sphere decoder simulation for space time block codes
==========================================================

Key features to implement
-------------------------
- Matrix input
- Complex matrices
- Spherical constellation
- Euchner-Schnorr enumeration smart implementation of the algorithm

How to setup and run the program on Linux
-----------------------------------------
- Download the whole repository (as a zip or using git clone) somewhere and open terminal there
- Install all the required packages (e.g. with apt-get) and the [Armadillo C++ linear algebra library](http://arma.sourceforge.net/download.html)
- If you wish to use plotting directly from the program you need to install C++ Boost library (e.g. *libboost1.58-all-dev* from the package manager)
- Otherwise comment out/remove the line that defines PLOTTING in the *src/misc.hpp* file
- Compile the program in terminal by typing: **make**
- Open the *settings.ini* in the */settings/* folder to setup the program, it should look something like this:

```ini
// configuration settings and simulation parameters for the sphere decoder program //

basis_file=alamouti.txt         // Text file containing the basis matrices
output_file=auto                // Optionally spesify the output filename (auto = automatic)
x-PAM=4                         // The size of the PAM signaling set
energy_estimation_samples=-1    // Number of samples to make the code energy estimation (-1 = sample all)
no_of_matrices=4                // Number of basis matrices
no_of_transmit_antennas=2       // Number of broadcast antennas
no_of_receiver_antennas=2       // Number of receiver antennas
snr_min=6                       // Minimum value for signal-to-noise ratio
snr_max=20                      // Maximum value for signal-to-noise ratio
snr_step=2		                // Increase SNR by this value per each iteration
simulation_rounds=10000         // Number of simulation rounds to run
required_errors=-1              // Demand at minimum this many errors before the simulation ends
time_slots=2		            // Number of time slots used in the code
plot_results=1                  // Draw plots? (1 = yes, -1 = no)
spherical_shaping_max_power=-1  // Defines the maximum distance from origin for codebook elements (-1 = unbounded)
```

- Create file called *bases.txt* in */bases/* folder and put your code basis matrices there (e.g. in Mathematica format)
- Run the program with: **./sphdec** or **make run** (the latter runs with the default settings file)
- If you configured the program correctly it should now run the simulation
- You can have multiple settings files (in the */settings/* folder) and use them in the simulation by giving their name as an command line argument for the program like: **./sphdec alamouti_settings.ini**
- Program output should be found at */output/* folder

(Naturally in order for this to work you need a C++ compiler (g++) that supports C++14 standard installed on your system, should be no problem on Aalto computers)
