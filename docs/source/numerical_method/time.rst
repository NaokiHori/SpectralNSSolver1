############
Time marcher
############

****************************
Integrating-factor technique
****************************

I would like to integrate the following equations in time:

.. math::

   \der{\wav{\ux}{\ix \iy}}{t}
   +
   \frac{1}{Re} \left( \mkx^2 + \mky^2 \right) \wav{\ux}{\ix \iy}
   =
   \wav{f_x}{\ix \iy},

.. math::

   \der{\wav{\uy}{\ix \iy}}{t}
   +
   \frac{1}{Re} \left( \mkx^2 + \mky^2 \right) \wav{\uy}{\ix \iy}
   =
   \wav{f_y}{\ix \iy},

.. math::

   \der{\wav{T}{\ix \iy}}{t}
   +
   \frac{1}{Re Sc} \left( \mkx^2 + \mky^2 \right) \wav{T}{\ix \iy}
   =
   \wav{g}{\ix \iy},

where :math:`f_x`, :math:`f_y`, and :math:`g` are the other terms such as the non-linear terms.

In general, these equations have the following form:

.. math::

   \der{p}{t}
   +
   C p
   =
   q,

where :math:`C` is a constant.

Here I introduce an integrating factor

.. math::

   \exp{
      \left( C t \right)
   }

and multiply it to the equation:

.. math::

   \exp{
      \left( C t \right)
   }
   \left(
      \der{p}{t}
      +
      C p
   \right)
   =
   \exp{
      \left( C t \right)
   }
   q.

The left-hand-side terms are re-arranged as

.. math::

   \der{}{t}
   \left\{
      p
      \exp{
         \left( C t \right)
      }
   \right\},

and thus I obtain the resulting differential equation:

.. math::

   \der{P}{t}
   =
   Q,

where

.. math::

   P
   \equiv
   p
   \exp{
      \left( C t \right)
   }

and

.. math::

   Q
   \equiv
   q
   \exp{
      \left( C t \right)
   }.

Up to this point, no approximations have been employed, ensuring that all relations are analytically accurate.
From this point onward, I will introduce an approximation for the right-hand-side integral in order to handle it numerically.
Specifically, I will utilise the classical fourth-order Runge-Kutta scheme in conjunction with the integrating-factor technique.

************************************
General explicit Runge-Kutta schemes
************************************

Here I aim to numerically solve

.. math::

   \der{P}{t}
   =
   Q,

where

.. math::

   p
   \equiv
   P
   \exp{
      \left( C t \right)
   }

and

.. math::

   Q
   \equiv
   q
   \exp{
      \left( C t \right)
   },

where :math:`\exp{\left( C t \right)}` is the integrating factor.

A general explicit Runge-Kutta scheme reads

.. math::

   & P^1 = P^n + Q^0 a_{10} \Delta t \\
   & P^2 = P^n + Q^0 a_{20} \Delta t + Q^1 a_{21} \Delta t \\
   & P^3 = P^n + Q^0 a_{30} \Delta t + Q^1 a_{31} \Delta t + Q^2 a_{32} \Delta t \\
   & \vdots

The :math:`k`-th row is written as

.. math::

   P^k = P^n + \sum_{l = 0}^{k - 1} Q^l a_{kl} \Delta t.

Note that

.. math::

   t^0 = t^n,
   P^0 = P^n.

******************
Integrating factor
******************

For notational simplicity, I define

.. math::

   E \left( x \right)
   \equiv
   \exp{
      \left(
         C x
      \right)
   }.

Assigning :math:`P` and :math:`Q` to the above relation yields

.. math::

   p^k
   E \left(
      t^n
      +
      c_k \Delta t
   \right)
   =
   p^n
   E \left(
      t^n
   \right)
   +
   \sum_{l = 0}^{k - 1}
   q^l
   E \left(
      t^n
      +
      c_l \Delta t
   \right)
   a_{kl} \Delta t.

Dividing the equation by :math:`E \left( t^n \right)` leads to

.. math::

   p^k
   E \left(
      c_k \Delta t
   \right)
   =
   p^n
   +
   \sum_{l = 0}^{k - 1}
   q^l
   E \left(
      c_l \Delta t
   \right)
   a_{kl} \Delta t.

*******************************
Fourth-order Runge-Kutta scheme
*******************************

The Butcher tableau of the classical fourth-order scheme

.. math::

   \begin{array}{c|cccc}
      c_0 &        &        &        &        \\
      c_1 & a_{10} &        &        &        \\
      c_2 &        & a_{21} &        &        \\
      c_3 &        &        & a_{32} &        \\
      \hline
      c_4 & a_{40} & a_{41} & a_{42} & a_{43} \\
   \end{array}

is

.. math::

   \begin{array}{c|cccc}
        0 &     &     &     &     \\
      1/2 & 1/2 &     &     &     \\
      1/2 &     & 1/2 &     &     \\
        1 &     &     & 1   &     \\
      \hline
        1 & 1/6 & 1/3 & 1/3 & 1/6 \\
   \end{array}

The whole process to update a field from
:math:`p^n` to :math:`p^{n+1}` is as follows:

.. math::

   p^0 = p^n.

.. math::
   \newcommand{\va}{1}
   \newcommand{\vb}{0}
   p^{\va}
   E \left(
      c_{\va} \Delta t
   \right)
   =
   p^n
   +
   q^{\vb}
   E \left(
      c_{\vb} \Delta t
   \right)
   a_{\va \vb} \Delta t.

.. math::
   \newcommand{\va}{2}
   \newcommand{\vb}{1}
   p^{\va}
   E \left(
      c_{\va} \Delta t
   \right)
   =
   p^n
   +
   q^{\vb}
   E \left(
      c_{\vb} \Delta t
   \right)
   a_{\va \vb} \Delta t.

.. math::
   \newcommand{\va}{3}
   \newcommand{\vb}{2}
   p^{\va}
   E \left(
      c_{\va} \Delta t
   \right)
   =
   p^n
   +
   q^{\vb}
   E \left(
      c_{\vb} \Delta t
   \right)
   a_{\va \vb} \Delta t.

.. math::
   p^{4}
   E \left(
      c_{4} \Delta t
   \right)
   & =
   p^n \\
   & +
   q^{0}
   E \left(
      c_{0} \Delta t
   \right)
   a_{40} \Delta t \\
   & +
   q^{1}
   E \left(
      c_{1} \Delta t
   \right)
   a_{41} \Delta t \\
   & +
   q^{2}
   E \left(
      c_{2} \Delta t
   \right)
   a_{42} \Delta t \\
   & +
   q^{3}
   E \left(
      c_{3} \Delta t
   \right)
   a_{43} \Delta t.

.. math::

   p^{n+1}
   =
   p^4.

Here, :math:`p^k` is stored to a buffer which is allocated to store the intermediate field.
Since :math:`a_{kl}` is all zero except :math:`k - l = 1`, only one additional buffer is needed for this purpose.

On the other hand, the derivatives :math:`q^k` are all to be stored since the classical scheme is not a low-storage scheme.

