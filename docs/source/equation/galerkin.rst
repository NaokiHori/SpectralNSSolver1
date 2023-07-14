#######################
Fourier-Galerkin method
#######################

.. note::

   For simplicity I consider a two-dimensional space in this part.

**************
Fourier series
**************

I consider a fully-periodic domain :math:`\left[ 0, \lx \right) \times \left[ 0, \ly \right)`, where a scalar field :math:`q \in \mathbb{R}^2` (velocity, pressure, passive scalar in the above equations) can be expanded using the Fourier series (trigonometric functions as the trial functions):

.. math::

   q \left( t, x, y \right)
   =
   \sum_{\kx^{\prime}}
   \sum_{\ky^{\prime}}
   \wav{q}{\kx^{\prime} \ky^{\prime}} \left( t \right)
   \exp \left( I \kx^{\prime} \frac{2 \pi}{\lx} x \right)
   \exp \left( I \ky^{\prime} \frac{2 \pi}{\ly} y \right),

where :math:`\kx^{\prime}, \ky^{\prime} \in \mathbb{Z}` and take :math:`\left[ - N / 2, - N / 2 + 1, \cdots, N / 2 - 1 \right]`, where :math:`N` is the degree of freedom in the direction.
Note that :math:`\wav{q}{\kx \ky}` are complex numbers :math:`\in \mathbb{C}`.
Hereafter the temporal dependency is dropped for notational simplicity.

Integrating this equation in the whole domain with the trigonometric functions as the weighting function

.. math::

   \int_{0}^{\ly}
   \int_{0}^{\lx}
   q \left( x, y \right)
   \exp \left( - I \kx \frac{2 \pi}{\lx} x \right)
   \exp \left( - I \ky \frac{2 \pi}{\ly} y \right)
   dx
   dy

yields

.. math::

   &
   \int_{0}^{\ly}
   \int_{0}^{\lx}
   \sum_{\kx^{\prime}}
   \sum_{\ky^{\prime}}
   \wav{q}{\kx^{\prime} \ky^{\prime}}
   \exp \left( I \kx^{\prime} \frac{2 \pi}{\lx} x \right)
   \exp \left( I \ky^{\prime} \frac{2 \pi}{\ly} y \right)
   \exp \left( - I \kx \frac{2 \pi}{\lx} x \right)
   \exp \left( - I \ky \frac{2 \pi}{\ly} y \right)
   dx
   dy \\
   &
   =
   \lx \ly \wav{q}{\kx \ky},

where the orthogonality is used to reach the final relation.

******************
Spatial derivative
******************

For the spatial derivative of this quantity, I have

.. math::

   \der{q}{x}
   & =
   \der{}{x}
   \sum_{\kx^{\prime}}
   \sum_{\ky^{\prime}}
   \wav{q}{\kx^{\prime} \ky^{\prime}}
   \exp \left( I \kx^{\prime} \frac{2 \pi}{\lx} x \right)
   \exp \left( I \ky^{\prime} \frac{2 \pi}{\ly} y \right) \\
   & =
   \sum_{\kx^{\prime}}
   \sum_{\ky^{\prime}}
   I \kx^{\prime} \frac{2 \pi}{\lx}
   \wav{q}{\kx^{\prime} \ky^{\prime}}
   \exp \left( I \kx^{\prime} \frac{2 \pi}{\lx} x \right)
   \exp \left( I \ky^{\prime} \frac{2 \pi}{\ly} y \right),

and thus the integral

.. math::

   \int_{0}^{\ly}
   \int_{0}^{\lx}
   \der{q}{x}
   \exp \left( - I \kx \frac{2 \pi}{\lx} x \right)
   \exp \left( - I \ky \frac{2 \pi}{\ly} y \right)
   dx
   dy

yields

.. math::

   \lx \ly I \kx \frac{2 \pi}{\lx} \wav{q}{\kx \ky}.

***************
Convolution sum
***************

Integrand:

.. math::

   p
   q
   =
   \sum_{\kx_0^{\prime}} \sum_{\ky_0^{\prime}}
   \wav{p}{\kx_0^{\prime} \ky_0^{\prime}}
   \exp \left( I \kx_0^{\prime} \frac{2 \pi}{\lx} x \right)
   \exp \left( I \ky_0^{\prime} \frac{2 \pi}{\ly} y \right)
   \sum_{\kx_1^{\prime}} \sum_{\ky_1^{\prime}}
   \wav{q}{\kx_1^{\prime} \ky_1^{\prime}}
   \exp \left( I \kx_1^{\prime} \frac{2 \pi}{\lx} x \right)
   \exp \left( I \ky_1^{\prime} \frac{2 \pi}{\ly} y \right),

Integral:

.. math::

   &
   \int_{0}^{\ly}
   \int_{0}^{\lx}
   pq
   \exp \left( - I \kx \frac{2 \pi}{\lx} x \right)
   \exp \left( - I \ky \frac{2 \pi}{\ly} y \right)
   dx
   dy \\
   =
   &
   \int_{0}^{\ly}
   \int_{0}^{\lx}
   \sum_{\kx_0^{\prime}}
   \sum_{\ky_0^{\prime}}
   \sum_{\kx_1^{\prime}}
   \sum_{\ky_1^{\prime}}
   \wav{p}{\kx_0^{\prime} \ky_0^{\prime}}
   \wav{q}{\kx_1^{\prime} \ky_1^{\prime}}
   \exp \left\{ I \left( \kx_0^{\prime} + \kx_1^{\prime} - \kx \right) \frac{2 \pi}{\lx} x \right\}
   \exp \left\{ I \left( \ky_0^{\prime} + \ky_1^{\prime} - \ky \right) \frac{2 \pi}{\ly} y \right\}
   dx
   dy \\
   =
   &
   \lx \ly
   \sum_{\kx_0^{\prime} + \kx_1^{\prime} = \kx}
   \sum_{\ky_0^{\prime} + \ky_1^{\prime} = \ky}
   \wav{p}{\kx_0^{\prime} \ky_0^{\prime}}
   \wav{q}{\kx_1^{\prime} \ky_1^{\prime}}.

