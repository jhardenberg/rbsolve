# RBsolve
## Pseudospectral Rayleigh-Bénard Solver
---------------------
A fast parallel Fortran code which integrates the Rayleigh-Bénard convection equations in the Boussinesq approximation for an horizontally periodic geometry. Free-slip and no-slip, Dirichlet and von Neumann boundary conditions are supported. Supports also a second active scalar allowingt o solve doubly diffusive problems. The code can also be used for DNS simulations of the Navier-Stokes equations.

(c) 2003-2007 Rayleigh-Bénard, Double diffusion by J. von Hardenberg (ISAC-CNR) 
(c) 2003 MPI version by J. von Hardenberg (ISAC-CNR) and G. Passoni (POLIMI)
(c) 2002 Navier-Stokes version by G. Passoni (POLIMI)

---------------------

# Very short code description 
The code is pseudospectral in the horizontal and finite differences in the vertical using a grid which is denser close to the vertical boundaries. Time advancement uses a third-order fractional step method (Kim and Moin, J. Comp. Phys. 59 (2), 1985)
MPI parallelization uses a vertical partitioning in slices.

# Configuration 

- config.h
The file config.h contains all options (as CPP define keys) which select the physical structure of the model. So the type of boundary conditions, if the problem is Rayleigh-Benard or Salt Fingering, of if you have the FFTW library is selected here. If you want a 2D problem select ONLY2D.

config.h also contains a line
#define NOMPI 
uncomment this line if you would like to compile a scalar ( not MPI ) version of the code, but leave it for compiling the tools which are all serial.

At the bottom notice the definitions
```
#define VDIFF (Pr)
#define BUOYT (Ra*Pr)
#define TDIFF (1.d0)
```
which imply that the problem is solved using the normalization where in the equations: `Pr` is in front of the Laplacian in the momentum equation, `Ra*Pr` is in front of the buoyancy term and the normal thermal diffusivity is 1. In theory one could change these definitions and solve directly the equations using a different normalization convention. Suggestion: do not touch.

- `param.h`
The file param.h is used to select the resolution of the model. 
Just comment or uncomment one of the many example lines.
You will need to set Nx and Nz (the horizontal resolution) and Ny (the vertical resolution) - The code follows an engineering convention by which the vertical is y.
There is a second set of resolutions to define, the spectral resolutions. These should take into account de-aliasing. For this code a 4/5 dealiasing (softer than the classical 2/3) has been found to work well. So, set kx=4/5 Nx and kz=4/5 Nx (rounded at the closest EVEN integer).

If you are compiling for MPI you will also need to set the two parameter variables (at the bottom):
```
      parameter (Nylmem=17)
      parameter  (NPROC=8)
```
where NPROC is the number of cores you plan to use and Nylmem is how may vertical layers to allocate for each 'slice'   associated with each core. This should be slightly in excess of `Ny/NPROC`.

- `param0`
Set up also the file param0 with all the physical parameters for the experiment (Rayleigh number, domain size etc).

# Compiling 

- Makefile
Edit the makefile to set up you compiler and its options some examples are included. For the tools you need to specify in the Makefile that you are compiling without MPI leaving the line "NOMPI=1" uncommented.
For compiling the main code with MPI you need to comment this line first.

To make sure that everything is compiled correctly better do a "make clean" first.
Compile the main code with "make"


# Running 
We suggest to create a separate 'run' directory and to copy there the following files:

* `param0` (The main physical parameter file)
* `param1` (parameters for timestepping - leave this one alone)
* `cpuweights` (a silly file with a lot of 1 values - This one is used only if you have a cluster of very different machines and you want to distribute the load differently among them - just use it as it is)
* rb  (our executable)

You can create initial conditions with "inicond" which will create a  randomly perturbed (on T) linear conductive solution. This will create initial files for t,u,v,w and a nrec.d file (containing the number 0).
The file nrec.d contains always the timestep of the latest save. When the code starts it checks for this file and restarts from the timestep it indicates.

# Tools 

The subdirectory 'tools' contains some useful tools, such as tools for extracting vertical profiles (prof) , horizontal and vertical slices (sectionh and sectionv) or spectra. To compile, if not already there, copy the relevant source .F file in the same directory as the Makefile. The command alone gives a short explanation on its usage.

# Out-of-the-box examples 

The example code is set up for a Ra=1e7 run, starting from an initial conductive solution with RB standard Dirichlet BCs. The code is set up for a MPI run with 8 cores.

After setting up the Makefile correctly for your compiler (see above) you can try the following:
 
1. Create some tools:
- Edit Makefile and leave active the line "NOMPI = 1"
```
make inicond   
make sectionh 
make sectionv
make prof
```

2. Compile and run the code

- Edit Makefile and comment out the line "NOMPI = 1"
```
make clean
make
./inicond
mpiexec -n 8 ./rb
```
3. Make a vertical profile at time 2000:
```
./prof t 2000
```

4. Make a vertical section (at y=96, step=2000)
```
./sectionv t 96 2000
```

# References to the code 

- J. von Hardenberg, A.Parodi, G. Passoni, A. Provenzale, E.A Spiegel, 2008: Large-scale patterns in Rayleigh-Benard convection. Physics Letters A, 372, 2223-2229.
- A. Parodi, J. von Hardenberg, G. Passoni, A. Provenzale and E.A Spiegel. Clustering of plumes in turbulent convection. Phys. Rev. Lett. 92 (19), 194503 (2004).
- G. Passoni, G. Alfonsi, and M. Galbiati, Int. J. Numer. Methods Fluids 38, 1069 (2002).


