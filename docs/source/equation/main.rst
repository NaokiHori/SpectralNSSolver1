###################
Governing equations
###################

I consider the incompressibility constraint:

.. math::

   \der{u_i}{x_i}
   =
   0

and the momentum balance:

.. math::

   \der{u_i}{t}
   =
   -
   u_j \der{u_i}{x_j}
   -
   \der{p}{x_i}
   +
   \frac{1}{Re}
   \der{}{x_j} \der{u_i}{x_j}
   +
   a_i

to describe the motion of the fluid, where :math:`Re` is the Reynolds number.
Also a passive scalar field :math:`T` is transported, which is governed by the advection-diffusion equation:

.. math::

   \der{T}{t}
   =
   -
   u_j \der{T}{x_j}
   +
   \frac{1}{Re Sc}
   \der{}{x_j} \der{T}{x_j},

where :math:`Sc` is the Schmidt number (the ratio of the fluid diffusivity to the scalar diffusivity).

For later convenience, I consider the advective terms in the divergence form:

.. math::

   \der{u_j q}{x_j},

where the incompressibility constraint is used:

.. math::

   \der{u_j q}{x_j}
   \equiv
   q \der{u_j}{x_j}
   +
   u_j \der{q}{x_j}.

.. toctree::
   :maxdepth: 1

   galerkin
   spectral

