Runtime Parameters
==================

The physics and numerics of a given `rbsolve` simulation are controlled entirely at runtime by the text file `param_namelist`. This allows the exact same binary executable to run entirely different regimes without recompiling.

The `param_namelist` File
-------------------------

The file is divided into two logical sections: `&NUMERIC` for integration parameters and `&PHYS` for the physical properties of the fluid and the domain geometry.

&NUMERIC Options
^^^^^^^^^^^^^^^^

.. code-block:: fortran

    &NUMERIC
          ptan=1.7    ! Geometric grid stretching factor (boundary layer)
          qtan=1.5    ! Internal grid mapping
          dt=1e-7     ! Timestep 
          ttot=4000   ! Absolute final timestep step to reach
          nsave=250   ! How often to save snapshot fields (in steps)
          nsnc=1000   ! Switch to a new NetCDF output file every nsnc steps
          alfa=0.0d0, ...  ! Timestepping coefficients 
    &END

* ``ttot`` sets the final absolute step number where the run will stop.
* ``nsave`` is the cadence. Snapshot records will be generated every ``nsave`` steps.
* ``nsnc`` controls file chunking. For example, if ``nsnc=1000``, the code groups saves such that records generated in the interval `[0, 1000)` go to ``rb_0000000.nc``, steps `[1000, 2000)` go to ``rb_0001000.nc``, etc.

&PHYS Options
^^^^^^^^^^^^^

.. code-block:: fortran

    &PHYS
          DD=1.d0                 ! Vertical domain extension
          Lx=6.28318530717958     ! Domain length in X
          Lz=6.28318530717958     ! Domain length in Z
          Re=300.d0               ! Reynolds number (Prandtl, NS)
          Ra=1e7                  ! Rayleigh number (RB)
          Pr=0.71d0               ! Prandtl number
    &END

1. **Domain Size**: The aspect ratio is strictly given by :math:`Lx / DD` and :math:`Lz / DD`. In the example above, an aspect ratio of :math:`2\pi \times 2\pi` is used.
2. **Rayleigh Number (Ra)**: Controls the strength of buoyancy forces versus viscous forces. Higher Ra creates more intense turbulence and a thinner thermal boundary layer.
3. **Prandtl Number (Pr)**: Ratio of momentum diffusivity to thermal diffusivity. A value of `0.71` corresponds approximately to air.

Grid Resolution Modification
----------------------------

The spatial resolution (the number of grid points or Fourier modes) is strictly defined in `src/param.h`. Changing these requires a full re-compilation.

* ``Nx``, ``Nz``: Number of points in the X and Z physical directions.
* ``Ny``: Number of vertical grid levels.
* ``nwax``, ``nwaz``: The spectral cutoff modes in X and Z.
