C++14 Sphere decoder simulation for space time block codes
==========================================================

Key features to implement
-------------------------
- Matrix input
- Complex matrices
- Spherical constellation
- Euchner-Schnorr enumeration smart implementation of the algorithm

How to run the program on Linux
-------------------------------
1. Download the *sphdec.cpp* file in some folder and open terminal there
2. Install all the required packages (e.g. with apt-get) and the [Armadillo C++ linear algebra library](http://arma.sourceforge.net/download.html)
2. Compile the program in terminal by typing: **g++ sphdec.cpp -o sphdec -O2 -larmadillo -llapack -lblas -std=c++14**
3. Run the program with: **./sphdec**
4. The program should exit and have genererated *settings.ini* file in the same folder
5. Open the settings.ini to setup the program, it should look something like this:
```
// configuration settings and simulation parameters for the sphere decoder program //

basis_file=bases.txt          // Text file containing the basis matrices
output_file=output.txt        // Text file used for simulation output
x-PAM=2                       // The size of the PAM signaling set
energy_estimation_samples=10  // Number of samples to make the code energy estimation (-1 = sample all)
no_of_matrices=2              // Number of basis matrices (dimension of the data vectors)
time_slots=2                  // Number of time slots used in the code
no_of_transmit_antennas=2     // Number of transmit antennas
no_of_receiver_antennas=2     // Number of receiver antennas
snr_min=6                     // Minimum value for signal-to-noise ratio
snr_max=12                    // Maximum value for signal-to-noise ratio
required_errors=500           // Demand at minimum this many errors before the simulation ends
```
6. Create file called *bases.txt* and put your code basis matrices there (e.g. in Mathematica format)
7. If you configured the program correctly it should now run the simulation with: **./sphdec**
8. You can have multiple settings files and use them in the simulation by giving their name as an command line argument for the program like: **./sphdec alamouti_settings.ini**
9. program output should be found at *output.txt* file

(Naturally you need a C++ compiler (g++) that supports C++14 standard installed on your system, should be no problem on Aalto computers)