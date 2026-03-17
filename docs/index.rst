.. rbsolve documentation master file, created by
   sphinx-quickstart.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to rbsolve's documentation!
===================================

``rbsolve`` is a high-performance Fortran code for solving Rayleigh-Bénard turbulence in a fully 3D convective setup. It simulates the fluid dynamics by combining a pseudo-spectral formulation in the periodic horizontal directions with a finite difference scheme in the vertical direction.

Recently, the code has been updated to support optional **NetCDF4** output in physical space, making post-processing and visual analysis significantly more straightforward while maintaining exact restart capabilities.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   physics
   compilation
   parameters
   usage

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
