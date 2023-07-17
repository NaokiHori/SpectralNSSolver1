###############
Advective terms
###############

********
Overview
********

To evaluate the convolution sums efficiently, instead of computing

.. math::

   \sum_{\ix_0^{\prime} + \ix_1^{\prime} = \ix}
   \sum_{\iy_0^{\prime} + \iy_1^{\prime} = \iy}
   \wav{p}{\ix_0^{\prime} \iy_0^{\prime}}
   \wav{q}{\ix_1^{\prime} \iy_1^{\prime}},

where :math:`p` and :math:`q` are the arguments and whose computational cost is :math:`\mathcal{O} \left( N_x^2 \times N_y^2 \right)`, I consider their representations in the physical domain, i.e. computing the product:

.. math::

   p q

directly, whose cost is :math:`\mathcal{O} \left( N_x \log N_x \times N_y \log N_y \right)`, which is known as the transform method.

**************
Implementation
**************

#. ``src/fluid/physical.c``

   This file contains a function to compute the description of each flow field in the physical place by performing the multi-dimensional inverse Fourier transform.

#. ``src/fluid/slope.c``

   This file contains several functions which evaluate the right-hand-side terms of the Runge-Kutta scheme.
   In particular, ``convolute`` computes the convolution sum (product of the given two arrays), while ``compute_adv`` computes the advective terms by multiplying pre-factors corresponding to the spatial derivative in each direction and the pre-factors.
   Finally ``project_velocity`` is necessary for the velocity field to satisfy the divergence-free condition (see the :ref:`equations <spectral>`, where two terms are involved in each direction to describe the advection).

