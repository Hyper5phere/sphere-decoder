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
- Otherwise comment out/remove the line that defines PLOTTING in the *main.cpp*
- Compile the program in terminal by typing: **make**
- Open the *settings.ini* in the */settings/* folder to setup the program, it should look something like this:

```ini
// configuration settings and simulation parameters for the sphere decoder program //

basis_file=bases.txt          // Text file containing the basis matrices
x-PAM=2                       // The size of the PAM signaling set
energy_estimation_samples=10  // Number of samples to make the code energy estimation (-1 = sample all)
no_of_matrices=2              // Number of basis matrices (dimension of the data vectors)
time_slots=2                  // Number of time slots used in the code
no_of_transmit_antennas=2     // Number of transmit antennas
no_of_receiver_antennas=2     // Number of receiver antennas
snr_min=6                     // Minimum value for signal-to-noise ratio
snr_max=12                    // Maximum value for signal-to-noise ratio
snr_step=2                    // Increase SNR by this value per each iteration
simulation_rounds=100000      // Number of simulation rounds to run
required_errors=500           // Demand at minimum this many errors before the simulation ends
```

- Create file called *bases.txt* in */bases/* folder and put your code basis matrices there (e.g. in Mathematica format)
- Run the program with: **./sphdec** or **make run** (the latter runs with the default settings file)
- If you configured the program correctly it should now run the simulation
- You can have multiple settings files (in the */settings/* folder) and use them in the simulation by giving their name as an command line argument for the program like: **./sphdec alamouti_settings.ini**
- Program output should be found at */output/* folder

(Naturally in order for this to work you need a C++ compiler (g++) that supports C++14 standard installed on your system, should be no problem on Aalto computers)