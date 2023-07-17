####################
Spectral NS Solver 1
####################

|License|_ |CI|_ |DOCS|_ |LastCommit|_

.. |License| image:: https://img.shields.io/github/license/NaokiHori/SpectralNSSolver1
.. _License: https://opensource.org/license/MIT

.. |CI| image:: https://github.com/NaokiHori/SpectralNSSolver1/actions/workflows/ci.yml/badge.svg
.. _CI: https://github.com/NaokiHori/SpectralNSSolver1/actions/workflows/ci.yml

.. |DOCS| image:: https://github.com/NaokiHori/SpectralNSSolver1/actions/workflows/documentation.yml/badge.svg
.. _DOCS: https://naokihori.github.io/SpectralNSSolver1/

.. |LastCommit| image:: https://img.shields.io/github/last-commit/NaokiHori/SpectralNSSolver1/main
.. _LastCommit: https://github.com/NaokiHori/SpectralNSSolver1/commits/main

.. image:: https://github.com/NaokiHori/SpectralNSSolver1/blob/main/docs/source/thumbnail.jpg
   :target: https://youtu.be/123d48J34eo
   :width: 800

.. image:: https://github.com/NaokiHori/SpectralNSSolver1/blob/main/docs/source/thumbnail2.jpg
   :target: https://youtu.be/eVjTrzJ4mMY
   :width: 800

.. image:: https://github.com/NaokiHori/SpectralNSSolver1/blob/main/docs/source/thumbnail3.jpg
   :width: 800

.. contents::
   :depth: 1

********
Overview
********

This library numerically solves the incompressible Navier-Stokes equations with a scalar field in two- and three-dimensional Cartesian domains using the spectral method.
This is developed for a self-study purpose to experience by myself the difference between `the finite-difference methods <https://github.com/NaokiHori/SimpleNSSolver>`_ and the spectral methods.

*******
Feature
*******

* Fourier-Galerkin method.
* `Orszag–Patterson algorithm <https://doi.org/10.1063/1.1692445>`_.
* 2/3 de-aliasing.
* `Pencil-based MPI parallelisation <https://github.com/NaokiHori/SimpleDecomp>`_ for more than 10^4 process.
* Fourth-order Runge-Kutta method (non-linear terms) combined with the integrating-factor technique (linear terms) for temporal integration.

Please refer to the `documentation <https://naokihori.github.io/SpectralNSSolver1/>`_ for details (under construction).

**********
Dependency
**********

* `C compiler <https://gcc.gnu.org>`_
* `GNU Make <https://www.gnu.org/software/make/>`_
* `MPI <https://www.open-mpi.org>`_
* `FFTW3 <https://www.fftw.org>`_

To initialise the flow field easily, I recommend

* `Python <https://www.python.org>`_ with `NumPy <https://numpy.org>`_ (see ``initial_condition/main.py``)

***********
Quick start
***********

#. Prepare workplace

   .. code-block:: console

      mkdir -p /path/to/your/directory
      cd       /path/to/your/directory

#. Get source

   .. code-block:: console

      git clone --recurse-submodules https://github.com/NaokiHori/SpectralNSSolver1
      cd SpectralNSSolver1

#. Set initial condition

   Although the scalar field can be arbitrary, the velocity field should be solenoidal.
   ``main.py`` offers several examples.

   .. code-block:: console

      cd initial_condition
      python3 main.py 0
      cd ..

#. Build

   .. code-block:: console

      make clean
      make output
      make all

#. Execute

   Parameters are defined in ``exec.sh``.

   .. code-block:: console

      bash exec.sh

   This may take a few minutes, depending on your machine spec.

The flow fields are saved under ``output/save/`` in `NPY <https://numpy.org/devdocs/reference/generated/numpy.lib.format.html>`_ format.
Note that these velocities are in the spectral domain; you need to perform the inverse Fourier transform (and the normalisation) to recover the velocities in the physical domain.

If proper Python libraries are installed, you can visualise the flow fields by

.. code-block:: console

   python3 visualise/2d.py

*************
3D simulation
*************

Please checkout ``3d`` branch and re-compile the whole source files.
You also need to re-generate the initial flow field.

*********
Reference
*********

* Canuto et al., Spectral Methods - Fundamentals in Single Domains, Springer

* Canuto et al., Spectral Methods - Evolution to Complex Geometries and Applications to Fluid Dynamics, Springer

***************
Acknowledgement
***************

I would like to thank `Dr. Chris Howland <https://chowland.github.io>`_ for fruitful discussions.

