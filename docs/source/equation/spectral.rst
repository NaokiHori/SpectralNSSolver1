
.. _spectral:

################################
Equations in the spectral domain
################################

I apply the Fourier series expansion to the governing equations to describe them in the spectral domain.

*****************
Incompressibility
*****************

.. math::

   \lx \ly
   \left(
      I \kx \frac{2 \pi}{\lx} \wav{\ux}{\kx \ky}
      +
      I \ky \frac{2 \pi}{\ly} \wav{\uy}{\kx \ky}
   \right)
   =
   0.

******************
Temporal evolution
******************

.. math::

   \lx \ly
   \der{\wav{\ux}{\kx \ky}}{t},

   \lx \ly
   \der{\wav{\uy}{\kx \ky}}{t},

   \lx \ly
   \der{\wav{T}{\kx \ky}}{t}.

***************
Advective terms
***************

The external forcing terms :math:`\wav{a_x}{\kx \ky}` and :math:`\wav{a_y}{\kx \ky}` are included.

.. math::

   \lx \ly
   \left[
      -
      I \kx \frac{2 \pi}{\lx}
      \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
      \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
      \wav{\ux}{\kx_0^{\prime} \ky_0^{\prime}}
      \wav{\ux}{\kx_1^{\prime} \ky_1^{\prime}}
      -
      I \ky \frac{2 \pi}{\ly}
      \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
      \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
      \wav{\uy}{\kx_0^{\prime} \ky_0^{\prime}}
      \wav{\ux}{\kx_1^{\prime} \ky_1^{\prime}}
      +
      \wav{a_x}{\kx \ky}
   \right],

   \lx \ly
   \left[
      -
      I \kx \frac{2 \pi}{\lx}
      \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
      \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
      \wav{\ux}{\kx_0^{\prime} \ky_0^{\prime}}
      \wav{\uy}{\kx_1^{\prime} \ky_1^{\prime}}
      -
      I \ky \frac{2 \pi}{\ly}
      \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
      \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
      \wav{\uy}{\kx_0^{\prime} \ky_0^{\prime}}
      \wav{\uy}{\kx_1^{\prime} \ky_1^{\prime}}
      +
      \wav{a_y}{\kx \ky}
   \right],

   \lx \ly
   \left[
      -
      I \kx \frac{2 \pi}{\lx}
      \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
      \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
      \wav{\ux}{\kx_0^{\prime} \ky_0^{\prime}}
      \wav{  T}{\kx_1^{\prime} \ky_1^{\prime}}
      -
      I \ky \frac{2 \pi}{\ly}
      \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
      \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
      \wav{\uy}{\kx_0^{\prime} \ky_0^{\prime}}
      \wav{  T}{\kx_1^{\prime} \ky_1^{\prime}}
   \right].

***********************
Pressure-gradient terms
***********************

.. math::

   - \lx \ly I \kx \frac{2 \pi}{\lx} \wav{p}{\kx \ky},

   - \lx \ly I \ky \frac{2 \pi}{\ly} \wav{p}{\kx \ky}.

***************
Diffusive terms
***************

.. math::

   - \frac{1}{Re}    \lx \ly \left[ \left( \kx \frac{2 \pi}{\lx} \right)^2 + \left( \ky \frac{2 \pi}{\ly} \right)^2 \right] \wav{\ux}{\kx \ky},

   - \frac{1}{Re}    \lx \ly \left[ \left( \kx \frac{2 \pi}{\lx} \right)^2 + \left( \ky \frac{2 \pi}{\ly} \right)^2 \right] \wav{\uy}{\kx \ky},

   - \frac{1}{Re Sc} \lx \ly \left[ \left( \kx \frac{2 \pi}{\lx} \right)^2 + \left( \ky \frac{2 \pi}{\ly} \right)^2 \right] \wav{  T}{\kx \ky}.

*************
Forcing terms
*************

They are embedded in the advective terms.

*******
Summary
*******

.. math::

   I \kx \frac{2 \pi}{\lx} \wav{\ux}{\kx \ky}
   +
   I \ky \frac{2 \pi}{\ly} \wav{\uy}{\kx \ky}
   =
   0.

.. math::

   \der{\wav{\ux}{\kx \ky}}{t}
   =
   \wav{h_x}{\kx \ky}
   -
   I \kx \frac{2 \pi}{\lx} \wav{p}{\kx \ky}
   -
   \frac{1}{Re} \left[ \left( \kx \frac{2 \pi}{\lx} \right)^2 + \left( \ky \frac{2 \pi}{\ly} \right)^2 \right] \wav{\ux}{\kx \ky},

   \der{\wav{\uy}{\kx \ky}}{t}
   =
   \wav{h_y}{\kx \ky}
   -
   I \ky \frac{2 \pi}{\ly} \wav{p}{\kx \ky}
   -
   \frac{1}{Re} \left[ \left( \kx \frac{2 \pi}{\lx} \right)^2 + \left( \ky \frac{2 \pi}{\ly} \right)^2 \right] \wav{\uy}{\kx \ky},

   \der{\wav{T}{\kx \ky}}{t}
   =
   \wav{g}{\kx \ky}
   -
   \frac{1}{Re Sc} \left[ \left( \kx \frac{2 \pi}{\lx} \right)^2 + \left( \ky \frac{2 \pi}{\ly} \right)^2 \right] \wav{T}{\kx \ky},

where

.. math::

   \wav{h_x}{\kx \ky}
   \equiv
   -
   I \kx \frac{2 \pi}{\lx}
   \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
   \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
   \wav{\ux}{\kx_0^{\prime} \ky_0^{\prime}}
   \wav{\ux}{\kx_1^{\prime} \ky_1^{\prime}}
   -
   I \ky \frac{2 \pi}{\ly}
   \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
   \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
   \wav{\uy}{\kx_0^{\prime} \ky_0^{\prime}}
   \wav{\ux}{\kx_1^{\prime} \ky_1^{\prime}}
   +
   \wav{a_x}{\kx \ky},

   \wav{h_y}{\kx \ky}
   \equiv
   -
   I \kx \frac{2 \pi}{\lx}
   \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
   \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
   \wav{\ux}{\kx_0^{\prime} \ky_0^{\prime}}
   \wav{\uy}{\kx_1^{\prime} \ky_1^{\prime}}
   -
   I \ky \frac{2 \pi}{\ly}
   \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
   \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
   \wav{\uy}{\kx_0^{\prime} \ky_0^{\prime}}
   \wav{\uy}{\kx_1^{\prime} \ky_1^{\prime}}
   +
   \wav{a_y}{\kx \ky},

   \wav{g}{\kx \ky}
   \equiv
   -
   I \kx \frac{2 \pi}{\lx}
   \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
   \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
   \wav{\ux}{\kx_0^{\prime} \ky_0^{\prime}}
   \wav{  T}{\kx_1^{\prime} \ky_1^{\prime}}
   -
   I \ky \frac{2 \pi}{\ly}
   \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
   \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
   \wav{\uy}{\kx_0^{\prime} \ky_0^{\prime}}
   \wav{  T}{\kx_1^{\prime} \ky_1^{\prime}}.

To eliminate the pressure from the momentum equation, I consider the inner product of the wave vector and the momentum balance, namely the sum of

.. math::

   I \kx \frac{2 \pi}{\lx}
   \der{\wav{\ux}{\kx \ky}}{t}
   =
   I \kx \frac{2 \pi}{\lx}
   \wav{h_x}{\kx \ky}
   +
   \left( \kx \frac{2 \pi}{\lx} \right)^2 \wav{p}{\kx \ky}
   -
   \frac{1}{Re} \left[
      \left( \kx \frac{2 \pi}{\lx} \right)^2
      +
      \left( \ky \frac{2 \pi}{\ly} \right)^2
   \right]
   I \kx \frac{2 \pi}{\lx} \wav{\ux}{\kx \ky}

and

.. math::

   I \ky \frac{2 \pi}{\ly}
   \der{\wav{\uy}{\kx \ky}}{t}
   =
   I \ky \frac{2 \pi}{\ly}
   \wav{h_y}{\kx \ky}
   +
   \left( \ky \frac{2 \pi}{\ly} \right)^2 \wav{p}{\kx \ky}
   -
   \frac{1}{Re} \left[
      \left( \kx \frac{2 \pi}{\lx} \right)^2
      +
      \left( \ky \frac{2 \pi}{\ly} \right)^2
   \right]
   I \ky \frac{2 \pi}{\ly} \wav{\uy}{\kx \ky}.

Since the temporal derivative terms and the diffusive terms are zero because of the incompressibility, I obtain

.. math::

   -
   I \kx \frac{2 \pi}{\lx}
   \wav{h_x}{\kx \ky}
   -
   I \ky \frac{2 \pi}{\ly}
   \wav{h_y}{\kx \ky}
   =
   \left[
      \left( \kx \frac{2 \pi}{\lx} \right)^2
      +
      \left( \ky \frac{2 \pi}{\ly} \right)^2
   \right]
   \wav{p}{\kx \ky},

which is used to eliminate the pressure, giving

.. math::

   \der{\wav{\ux}{\mkx \mky}}{t}
   =
   \wav{h_x}{\mkx \mky}
   -
   \frac{\mkx}{\mkx^2 + \mky^2}
   \left(
      \mkx \wav{h_x}{\mkx \mky}
      +
      \mky \wav{h_y}{\mkx \mky}
   \right)
   -
   \frac{1}{Re} \left( \mkx^2 + \mky^2 \right) \wav{\ux}{\mkx \mky},

   \der{\wav{\uy}{\mkx \mky}}{t}
   =
   \wav{h_y}{\mkx \mky}
   -
   \frac{\mky}{\mkx^2 + \mky^2}
   \left(
      \mkx \wav{h_x}{\mkx \mky}
      +
      \mky \wav{h_y}{\mkx \mky}
   \right)
   -
   \frac{1}{Re} \left( \mkx^2 + \mky^2 \right) \wav{\uy}{\mkx \mky}.

Recall that

.. math::

   \mkx
   \equiv
   \kx \frac{2 \pi}{\lx},

   \mky
   \equiv
   \ky \frac{2 \pi}{\ly}.

