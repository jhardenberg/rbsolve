Compilation & Configuration
===========================

``rbsolve`` leverages conditional compilation to enable or disable major model components and optimizations without branching the code. This is managed via the ``src/config.h`` file.

Basic Build
-----------

To compile the model on macOS or Linux:

.. code-block:: bash

    cd src/
    make clean
    make inicond
    make

This produces two executables:
1. ``inicond``: Generates the initial condition files.
2. ``rb``: The main simulation binary.

config.h Options
----------------

The ``config.h`` file contains C-preprocessor macros that fundamentally alter the equations or I/O methods. Once you change these, you **must** recompile the code.

Key definitions:

* ``#define TEMPERATURE``: Includes the active temperature scalar field, necessary for Rayleigh-Bénard convection.
* ``#define SALINITY``: Enables a second active scaler (e.g., for double-diffusive convection).
* ``#define NETCDF_OUTPUT``: Overrides the legacy unformatted binary grid outputs, saving standard ``.nc`` NetCDF4 files in physical space.
* ``#define SAVE_PRESSURE``: Instructs the code to explicitly write the pressure field `p` to disk alongside velocity and temperature. Note that pressure can also be diagnostically recovered from the velocity field.
* ``#define PRINF``: Enables the calculation of the pressure influence matrix (usually for specialized boundary conditions).
* ``#define LINSTABDBG``: Outputs extensive diagnostic information on linear stability to unit `20`.

Example Configuration
---------------------

To run a classical Rayleigh-Bénard setup with standard NetCDF output:

.. code-block:: C

    // Inside config.h
    #define TEMPERATURE
    #undef SALINITY
    #define NETCDF_OUTPUT
    #undef SAVE_PRESSURE
    #undef PRINF

After editing, run ``make clean`` and ``make`` to apply.
