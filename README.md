This is the version (10/2026) of FAR3d that includes the helical couplings that are needed for
the nonlinear simulation of stellarators and other 3D configurations (i.e. tokamaks with ripple,
RMPs, etc.). This version uses an input file (farin) based upon namelist formatting. To compile
follow these steps:
create clean build directory
cd into build
type: cmake ../
Then type: make

Input files and an example Slurm batch script are included in the Example directory.
farin is the overall input control file, Wisc_reactor_eq is the magnetic 3d equilibrium input,
and Wisc_reactor_profiles.txt provides the 1D profiles needed by FAR3d).
Example directory. This is a nonlinear run using 5 MPI processes and set up 
for 50 time steps (the maxstp parameter in farin determines the number of time steps).
It took about 10 minutes to run. Longer or shorter runs can be created by changing maxstp.
Note that a full nonlinear run to saturation for this case may take 30,000 to 50,000 time steps.
