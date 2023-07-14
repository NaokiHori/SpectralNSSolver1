####################
Domain decomposition
####################

`Pencil-like domain decomposition <https://github.com/NaokiHori/SimpleDecomp>`_ is adopted to parallelise the domain.

By default the domain is decomposed in the :math:`y` direction (:math:`x1` pencil), in which all variables in the spectral domain

.. math::

   \wav{\ux}{\ix \iy \iz},
   \wav{\uy}{\ix \iy \iz},
   \wav{\uz}{\ix \iy \iz},
   \wav{  T}{\ix \iy \iz}

are defined.
Note that they are all complex numbers whose data type is ``fftw_complex``:

.. myliteralinclude:: /../../include/fluid.h
   :language: c
   :tag: array in spectral domain, x1 pencil

Here the prefix ``s_`` and ``x1_`` denote the variables are in the spectral domain and stored as the ``x1`` pencils.

To evaluate the non-linear terms in the physical domain based on the transform method, multi-dimensional discrete Fourier transforms are to be performed, whose results are stored as the :math:`y1` pencil (or :math:`z1` pencil for three-dimensional domains).
Since they are real numbers in theory, I define them as ``double`` to save the storage and to roughly halve the number of operations needed:

.. myliteralinclude:: /../../include/fluid.h
   :language: c
   :tag: array in physical domain, y1 pencil

.. myliteralinclude:: /../../include/fluid.h
   :language: c
   :tag: array in physical domain, z1 pencil

