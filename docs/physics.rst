Physics of the Model
====================

``rbsolve`` implements a complete fully-resolved 3D Rayleigh-Bénard convection model. The system is governed by the non-dimensional incompressible Boussinesq equations for a fluid heated from below and cooled from above.

The Equations
-------------

The non-dimensional governing equations are:

.. math::

    \nabla \cdot \mathbf{u} &= 0 \\
    \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} &= -\nabla P + \sqrt{\frac{Pr}{Ra}} \nabla^2 \mathbf{u} + T \mathbf{\hat{z}} \\
    \frac{\partial T}{\partial t} + (\mathbf{u} \cdot \nabla)T &= \frac{1}{\sqrt{Ra Pr}} \nabla^2 T

Where:

* :math:`\mathbf{u} = (u,v,w)` is the velocity vector.
* :math:`T` is the temperature field.
* :math:`P` is the kinematic pressure.
* :math:`Ra` is the Rayleigh number, characterizing the strength of the thermal driving force against viscous dissipation and heat conduction.
* :math:`Pr` is the Prandtl number, representing the ratio of momentum diffusivity (kinematic viscosity) to thermal diffusivity.

Plume Clustering (PRL 92, 194503)
---------------------------------

As detailed by Parodi et al. in *Physical Review Letters 92, 194503 (2004)*, turbulent Rayleigh-Bénard convection produces intense coherent structures in the form of updrafts (warm plumes) and downdrafts (cold plumes). These plumes dominate the vertical heat transport.

The computational setup resolved by ``rbsolve`` was instrumental in proving that these convective thermal plumes are not randomly distributed. Instead, they spontaneously organize into macroscopic **clusters**. 

Key findings captured by the model dynamics:
1. Plumes have horizontal scales of the order of the thermal boundary layer thickness.
2. The clusters they form grow in size dynamically and can eventually span the entire domain.
3. This massive organization drives a large-scale mean circulation or wind which dramatically alters the macroscopic heat transport scaling.

The model solves these fields using a spectral decomposition in the horizontal :math:`(x,z)` plane and finite differences in the vertical :math:`y` direction, allowing extremely precise resolution of the boundary layers driving this plume clustering behavior.
