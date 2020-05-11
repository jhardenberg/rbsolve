# RBsolve
## Pseudospectral Rayleigh-Bénard Solver

(c) 2020      Non-homogeneous boundaries by J.von Hardenberg (PoliTO)
(c) 2003-2013 Rayleigh-Bénard, Double diffusion by J. von Hardenberg (ISAC-CNR) 
(c) 2003      MPI version by J. von Hardenberg (ISAC-CNR) and G. Passoni (POLIMI)
(c) 2002      Navier-Stokes version by G. Passoni (POLIMI)

## Very short code description ##

RBsolve is a CFD code for Rayleigh-Benard (natural convection) problem simulation with Boussinesq approximation.
The code is written in FORTRAN and solves a system of coupled equations for momentum and temperature in the incompressibility hypothesis.
Introduction of a second active scalar (for double diffusive problem) or the usage for DNS simulations are additional available options.


## Numerical features ##

Time advancement is performed using a third-order Runge-Kutta method (Kim and Moin, J. Comp. Phys. 59 (2), 1985). 
spatial integration is based on a pseudospectral method for the horizontal directions and a finite differences method for the vertical dimension, using a grid which is denser close to the vertical boundaries.
Horizontal periodic geometry is assumed.
MPI parallelization is available and performed using a vertical partitioning in slices.

The convention in is this code is to call _x_ and _z_ the horizontal coordinates and _y_ the vertical. Correspondingly the velocity components are _(u,w,v)_

## Configuration ##

### config.h ### 
The file `config.h` contains options which define the physical structure of the model.

* `ONLY2D`:  If you want a 2D problem.  This option reduces the 3D domain to a 2D using x and y directions. Remember to set also `kzmax=1` in `param.h`	

* `TEMPERATURE`: if defined, a convective problem is solved (RB). Else the Navier-Stokes equations are solved.

* `FINGER, SALINITY`: to run a salt fingering experiment.

* `FREE_SLIP_TOP`etc. :  Choose _free slip_ or _no slip_ separately for top and bottom boundaries (no slip is default mode).
  
* `FLUXT_BOTTOM` etc. : used to define fixed flux boundary conditions. 

* `TEMPERATURE_BOTTOM`etc: Dirichlet BCs for top and bottom.

* `PRESSURE_GRADIENT`: defines a mean pressure gradient (useful for NS experiments)

At the bottom of config.h notice the following definitions (for Rayleigh-Benard problem):

    #define VDIFF (Pr)
    #define BUOYT (Ra*Pr)
    #define TDIFF (1.d0)

which imply that the problem is solved using the normalization where in the equations: Pr is in front of the Laplacian in the momentum equation, Ra*Pr is in front of the buoyancy term and thermal diffusivity is 1. That is, the equations are non-dimensionalized respect to te diffusive time scale. 

### param.h ###

The file `param.h` is used to select the resolution of the model. 

* _Spatial resolution_

To set the spatial resolution, uncomment a triple of Nx, Nz (the horizontal resolution) and Ny (the vertical resolution).
The code follows an engineering convention by which the vertical is y.

* _Spectral resolution_

Also the spectral resolution has to be defined. These should take into account de-aliasing. 
For this code a 4/5 dealiasing (softer than the classical 2/3) has sometimes been found to work well, even if a 2/3 dealiasing is the safe option to avoid aliasing artifacts. 
So, set kx=4/5 Nx and kz=4/5 Nz (rounded to the closest even integer) or kx=2/3 Nx and kz=2/3 Nz.

If you are compiling for MPI you will also need to set the two parameter variables (at the bottom of the file):

    parameter (Nylmem = ...)
    parameter  (NPROC = ...)

where NPROC is the number of cores you plan to use and Nylmem is how may vertical layers to allocate for each 'slice' associated with each core. This should be slightly in excess of `Ny/NPROC`.

_Example:_
 
   parameter (Nylmem=17)
   parameter  (NPROC=8)
   
for the Ny=129 case.


### param_namelist ###

We recently changed the configuration file to a standard Fortran namelist. It is still possible to use the old format (param0) by commenting `#define NAMELIST_PARAMS`in `config.h`.

Set up  the file param_namelist with all the physical parameters for the experiment (this is only a subset of the possible parameters):

    &NUMERIC
      dt=1e-6     # Timestep
      ttot=50000  # Total number of steps
      nsave=250   # How often to save
    &END

    &PHYS
      Lx=6.28318530717958   # Domain length
      Lz=6.28318530717958   # Domain width
      Ra=1e7                # Rayleigh number
      Pr=0.71d0             # Prandtl number
      Omegay=0.d0           # Rotation component y direction
      Omegaz=0.d0           # rotation component z direction
      lapse=0.d0            # Lapse rate
      Re=300.d0             # Reynolds number
    &END

Not all paramters are used for all configurations. For example `lapse`and `Re` are ignored for the standard Rayleigh-Bénard problem.

## Compiling ##

Edit the `Makefile` to set up you compiler and its options some examples are included. 

To compile the tools you need to specify in the Makefile that you are compiling without MPI leaving the line "NOMPI=1" uncommented.
To compile the main code with MPI you need to comment this line first.

To make sure that everything is compiled correctly better do a "make clean" first.
To compile the main code use the command "make".

config.h also contains a line
      #define NOMPI 
Uncomment this line if you would like to compile a scalar (not MPI) version of the code, but leave it for compiling the tools which are all serial. (This is actually redundant with the NOMPI in the makefile and will be removed soon)

## Running ##

It's highly recommended to create a separate 'run' directory and to copy there the following files:

* `param_namelist` (The main  parameter file)
* `rb`  (our executable)

You can create initial conditions with the tool `inicond` which will create a randomly perturbed (on T) linear conductive solution. 
This will create initial files for t,u,v,w and a nrec.d file (containing the number 0).
The file `nrec.d` contains always the timestep of the latest save. When the code starts it checks for this file and restarts from the timestep it indicates.

## Tools ##

The subdirectory 'tools' contains some useful tools.
In order to compile them you will need to copy them to the same source directory as the rest of the code.
The command alone gives a short explanation on its usage. Below a short list with some examples:

_prof_ : extracts vertical temperature profiles and save it in a file (*.prf).

Usage example: `./prof t 2000` 

_sectionh/sectionv_ : extract horizontal and vertical slices.

Usage examples: `./sectionv t 96 2000` (component=t, vertical section at y=96, step=2000).

Usage examples: `./sectionh w 15 400` (component=w, horizontal section at y=15, step=400)

If you specify a negative position for the vertical slice, sectioning is in the y-z plane.

_spectrumh_ : compute horizontal spectrum (for a specified component) and save it in a file (*.spe).

Usage example: `./spectrumh u 50 4000` (component=u, y=50, step=4000)

_nusselt_ : computes the nusselt number  Nu = 1 + <wT> to estimates heat flux.

Usage example: `./nusselt 2000` (nusselt step=2000)

_cfl_ : computes the CFL number in the selected range. 

Usage: `./cfl 2000 3000 40 `

To compile, if not already there, copy the relevant source .F file in the same directory as the Makefile. 
Tools have to be compiled necessarily BEFORE the compilation of the code (see out-of-the-box examples).

## Usage examples ##

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

     `./prof t 2000`

4. Make a vertical section (at y=96, step=2000):

    `./sectionv t 96 2000`


## References ##

[1] J. von Hardenberg; D. Goluskin; A. Provenzale; E. A. Spiegel. Generation of Large-Scale Winds in Horizontally Anisotropic Convection. Physical Review Letters, Volume 115, Issue 13, id.134501 (2015).

[2] J. von Hardenberg, A.Parodi, G. Passoni, A. Provenzale, E.A Spiegel. Large-scale patterns in Rayleigh-Benard convection. Physics Letters A, 372, 2223-2229 (2008).

[3] A. Parodi, J. von Hardenberg, G. Passoni, A. Provenzale and E.A Spiegel. Clustering of plumes in turbulent convection. Phys. Rev. Lett. 92 (19), 194503 (2004).

[4] G. Passoni, G. Alfonsi, and M. Galbiati, Int. J. Numer. Methods Fluids 38, 1069 (2002).

