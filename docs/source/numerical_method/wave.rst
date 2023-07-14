#################################
Wave number and angular frequency
#################################

Since the `FFTW3 <https://www.fftw.org>`_ stores the spectra from the zero-th wave number, I need to map the indices to the wave numbers :math:`k` as follows:

.. myliteralinclude:: /../../src/domain.c
   :language: c
   :tag: compute wave numbers

Also for convenience, the angular frequency

.. math::

   \frac{2 \pi}{L} k

are computed and stored:

.. myliteralinclude:: /../../src/domain.c
   :language: c
   :tag: compute angular frequency

