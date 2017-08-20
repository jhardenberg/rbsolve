# RBsolve
## Pseudospectral Rayleigh-Bénard Solver


(c) 2003-2007 Rayleigh-Bénard, Double diffusion by J. von Hardenberg (ISAC-CNR) 

(c) 2003 MPI version by J. von Hardenberg (ISAC-CNR) and G. Passoni (POLIMI)

(c) 2002 Navier-Stokes version by G. Passoni (POLIMI)


Very short code description 
---------------------------

RBsolve is a CFD code for Rayleigh-Benard (natural convection) problem simulation with Boussinesq approximation.
The code is written in FORTRAN and solves a system of coupled equations for momentum and temperature in the incompressibility hypothesis.
Introduction of a second active scalar (for double diffusive problem) or the usage for DNS simulations are additional available options.


NUMERICAL FEATURES:

TIME ADVANCEMENT is performed using a third-order Runge-Kutta method (Kim and Moin, J. Comp. Phys. 59 (2), 1985). 
SPATIAL INTEGRATION is based on a PSEUDOSPECTRAL method for the horizontal directions and a FINITE DIFFERENCES method for the vertical dimension, using a grid which is denser close to the vertical boundaries.
Horizontal periodic geometry is assumed.
MPI parallelization is available and performed using a vertical partitioning in slices.




** Configuration **
---------------------
- config.h

The file config.h contains all options which select the physical structure of the model.

_If you want a 2D problem select ONLY2D.  This option reduces the 3D domain to a 2D using x and y directions. 	

_Boundary conditions for TOP and BOTTOM:

SPEED = FREE SLIP / NO SLIP (NO SLIP default mode).

TEMPERATURE = Dirichelet conditions (FIXED VALUES) or Neumann conditions (FLUX / NO FLUX) are available.

_Pressure = a mean pressure gradient is defined. 

At the bottom of config.h notice the following definitions (for Rayleigh-Benard problem):

    #define VDIFF (Pr)
    #define BUOYT (Ra*Pr)
    #define TDIFF (1.d0)

which imply that the problem is solved using the normalization where in the equations: Pr is in front of the Laplacian in the momentum equation, Ra*Pr is in front of the buoyancy term and the normal THERMAL DIFFUSIVITY IS 1.
(Similar definitions are mentioned for Salt Fingering problem).
In theory these definitions could be change and the equations could be directly solve using a different normalization convention. Suggestion: do not touch.



- param.h

The file param.h is used to select the RESOLUTION of the model. 
You can choose the resolution just comment or uncomment one of the many example lines.

_Spatial resolution

To set the spatial resolution, uncomment a triple of Nx, Nz (the horizontal resolution) and Ny (the vertical resolution).

The code follows an engineering convention by which the vertical is y.

_Spectral resolution

Also the spectral resolution has to be defined. These should take into account de-aliasing. 
For this code a 4/5 dealiasing (softer than the classical 2/3) has been found to work well. 
So, set kx=4/5 Nx and kz=4/5 Nx (rounded at the closest EVEN integer).

If you are compiling for MPI you will also need to set the two parameter variables (at the bottom of the file):

    parameter (Nylmem = ...)
    parameter  (NPROC = ...)

where NPROC is the number of cores you plan to use and Nylmem is how may vertical layers to allocate for each 'slice' associated with each core. This should be slightly in excess of Ny/NPROC.
EXAMPLE:
 
   parameter (Nylmem=17)
   parameter  (NPROC=8)
   
for the Ny=129 case.



- param0

Set up also the file param0 with all the physical parameters for the experiment:

    Lx & Lz = horizontal domain sizes
    DD = vertical domain size
    dt = time step width
    ndt = tot # of time steps
    saveevery = writing time (example: dt=2e-7  saveevery=500  write at= )
    Ra = Rayleigh number
    Pr = Prandtl number
    Le = Lewis number
    Re = Reynolds number
    Omega = 2*angular speed
    Rs = Saline Rayleigh number
    Rt = Thermal Rayleigh number
    Kscalar =
    ptan =
    qtan =




** Compiling **
----------------
- Makefile

Edit the makefile to set up you compiler and its options some examples are included. 

To compile the tools you need to specify in the Makefile that you are compiling without MPI leaving the line "NOMPI=1" uncommented.
To compile the main code with MPI you need to comment this line first.

To make sure that everything is compiled correctly better do a "make clean" first.
To compile the main code use the command "make".

config.h also contains a line
#define NOMPI 
Uncomment this line if you would like to compile a scalar (not MPI) version of the code, but leave it for compiling the tools which are all serial.




** Running **
--------------
It's highly recommended to create a separate 'run' directory and to copy there the following files:

* param0 (The main physical parameter file)
* param1 (parameters for timestepping - leave this one alone)
* cpuweights (a silly file with a lot of 1 values - This one is used only if you have a cluster of very different machines and you want to distribute the load differently among them - just use it as it is)
* rb  (our executable)

You can create initial conditions with "inicond" which will create a randomly perturbed (on T) linear conductive solution. 
This will create initial files for t,u,v,w and a nrec.d file (containing the number 0).
The file nrec.d contains always the timestep of the latest save. When the code starts it checks for this file and restarts from the timestep it indicates.



** Tools **
-----------
The subdirectory 'tools' contains some useful tools.
The command alone gives a short explanation on its usage. Below a short list with some examples:

_prof : extracts vertical temperature profiles and save it in a file (*.prf).

Usage: ./prof t 2000

_profvar : extracts variance value and save it in a file (*.var).

Usage: ./profvar u 2000

_sectionh/sectionv : extracts horizontal and vertical slices.

Usage: ./sectionv t 96 2000 (component=t, vertical section at y=96, step=2000).

Usage: ./sectionh w 15 400 (component=w, horizontal section at y=15, step=400)

_spectrumh : compute horizontal spectrum (for a specified component) and save it in a file (*.spe).

Usage: ./spectrumh u 50 4000 (component=u, y=50, step=4000)

_nusselt : shows the nusselt number  Nu = 1 + <wT>(is this the assumed formula?) to estimates heat flux.

Usage: ./nusselt 2000 (nusselt step=2000)

_cfl : shows the cfl in the selected range. 

Usage: ./cfl 2000 3000 40 

To compile, if not already there, copy the relevant source .F file in the same directory as the Makefile. 
Tools have to be compiled necessarily BEFORE the compilation of the code (see out-of-the-box examples).




** Out-of-the-box examples **
----------------------------------
The example code is set up for a Ra=1e7 run, starting from an initial conductive solution with RB standard Dirichlet BCs. 
The code is set up for a MPI run with 8 cores.

After setting up the Makefile correctly for your compiler (see above) you can try the following:
 
1. Create some tools:

  Edit Makefile and leave active the line "NOMPI = 1"

  make inicond   

  make sectionh 

  make sectionv

  make prof

2. Compile and run the code:

  Edit Makefile and comment out the line "NOMPI = 1"

  make clean

  make

  ./inicond

  mpiexec -n 8 ./rb

3. Make a vertical profile (at time 2000):

  ./prof t 2000

4. Make a vertical section (at y=96, step=2000):

  ./sectionv t 96 2000



** References to the code **
--------------------------
- J. von Hardenberg, A.Parodi, G. Passoni, A. Provenzale, E.A Spiegel, 2008: Large-scale patterns in Rayleigh-Benard convection. Physics Letters A, 372, 2223-2229.
- A. Parodi, J. von Hardenberg, G. Passoni, A. Provenzale and E.A Spiegel. Clustering of plumes in turbulent convection. Phys. Rev. Lett. 92 (19), 194503 (2004).
- G. Passoni, G. Alfonsi, and M. Galbiati, Int. J. Numer. Methods Fluids 38, 1069 (2002).

