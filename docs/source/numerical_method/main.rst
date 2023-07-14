################
Numerical method
################

I integrate the following three equations in time:

.. math::

   \der{\wav{\ux}{\ix \iy}}{t}
   =
   \wav{h_x}{\ix \iy}
   +
   \frac{\mkx}{\mkx^2 + \mky^2}
   \left(
      \mkx \wav{h_x}{\ix \iy}
      +
      \mky \wav{h_y}{\ix \iy}
   \right)
   -
   \frac{1}{Re} \left( \mkx^2 + \mky^2 \right) \wav{\ux}{\ix \iy},

.. math::

   \der{\wav{\uy}{\ix \iy}}{t}
   =
   \wav{h_y}{\ix \iy}
   +
   \frac{\mky}{\mkx^2 + \mky^2}
   \left(
      \mkx \wav{h_x}{\ix \iy}
      +
      \mky \wav{h_y}{\ix \iy}
   \right)
   -
   \frac{1}{Re} \left( \mkx^2 + \mky^2 \right) \wav{\uy}{\ix \iy},

.. math::

   \der{\wav{T}{\ix \iy}}{t}
   =
   \wav{g}{\ix \iy}
   -
   \frac{1}{Re Sc} \left( \mkx^2 + \mky^2 \right) \wav{T}{\ix \iy},

where the non-linear terms are given by

.. math::

   \wav{h_x}{\ix \iy}
   \equiv
   -
   I \mkx
   \sum_{\ix_0^{\prime} + \ix_1^{\prime} = \ix}
   \sum_{\iy_0^{\prime} + \iy_1^{\prime} = \iy}
   \wav{\ux}{\ix_0^{\prime} \iy_0^{\prime}}
   \wav{\ux}{\ix_1^{\prime} \iy_1^{\prime}}
   -
   I \mky
   \sum_{\ix_0^{\prime} + \ix_1^{\prime} = \ix}
   \sum_{\iy_0^{\prime} + \iy_1^{\prime} = \iy}
   \wav{\uy}{\ix_0^{\prime} \iy_0^{\prime}}
   \wav{\ux}{\ix_1^{\prime} \iy_1^{\prime}}
   +
   \wav{a_x}{\ix \iy},

.. math::

   \wav{h_y}{\ix \iy}
   \equiv
   -
   I \mkx
   \sum_{\ix_0^{\prime} + \ix_1^{\prime} = \ix}
   \sum_{\iy_0^{\prime} + \iy_1^{\prime} = \iy}
   \wav{\ux}{\ix_0^{\prime} \iy_0^{\prime}}
   \wav{\uy}{\ix_1^{\prime} \iy_1^{\prime}}
   -
   I \mky
   \sum_{\ix_0^{\prime} + \ix_1^{\prime} = \ix}
   \sum_{\iy_0^{\prime} + \iy_1^{\prime} = \iy}
   \wav{\uy}{\ix_0^{\prime} \iy_0^{\prime}}
   \wav{\uy}{\ix_1^{\prime} \iy_1^{\prime}}
   +
   \wav{a_y}{\ix \iy},

.. math::

   \wav{g}{\ix \iy}
   \equiv
   -
   I \mkx
   \sum_{\ix_0^{\prime} + \ix_1^{\prime} = \ix}
   \sum_{\iy_0^{\prime} + \iy_1^{\prime} = \iy}
   \wav{\ux}{\ix_0^{\prime} \iy_0^{\prime}}
   \wav{  T}{\ix_1^{\prime} \iy_1^{\prime}}
   -
   I \mky
   \sum_{\ix_0^{\prime} + \ix_1^{\prime} = \ix}
   \sum_{\iy_0^{\prime} + \iy_1^{\prime} = \iy}
   \wav{\uy}{\ix_0^{\prime} \iy_0^{\prime}}
   \wav{  T}{\ix_1^{\prime} \iy_1^{\prime}}.

.. note::

   For now the external forcing terms :math:`a_i` are omitted.

.. toctree::
   :maxdepth: 1

   domain
   wave
   adv
   time

