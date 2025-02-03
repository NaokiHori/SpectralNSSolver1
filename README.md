# Spectral NS Solver 1

[![License](https://img.shields.io/github/license/NaokiHori/SpectralNSSolver1)](https://opensource.org/license/MIT)
[![CI](https://github.com/NaokiHori/SpectralNSSolver1/actions/workflows/ci.yml/badge.svg)](https://github.com/NaokiHori/SpectralNSSolver1/actions/workflows/ci.yml)
[![Docs](https://github.com/NaokiHori/SpectralNSSolver1/actions/workflows/documentation.yml/badge.svg)](https://naokihori.github.io/SpectralNSSolver1/)
[![Last Commit](https://img.shields.io/github/last-commit/NaokiHori/SpectralNSSolver1/main)](https://github.com/NaokiHori/SpectralNSSolver1/commits/main)

[![Thumbnail 1](https://github.com/NaokiHori/SpectralNSSolver1/blob/main/docs/source/thumbnail.jpg)](https://youtu.be/123d48J34eo)
[![Thumbnail 2](https://github.com/NaokiHori/SpectralNSSolver1/blob/main/docs/source/thumbnail2.jpg)](https://youtu.be/eVjTrzJ4mMY)
![Thumbnail 3](https://github.com/NaokiHori/SpectralNSSolver1/blob/main/docs/source/thumbnail3.jpg)

## Overview

This library numerically solves the incompressible Navier-Stokes equations with a scalar field in two- and three-dimensional Cartesian domains using the spectral method. It was developed as a self-study project to explore the differences between [finite-difference methods](https://github.com/NaokiHori/SimpleNSSolver) and spectral methods.

## Features

- Fourier-Galerkin method
- [Orszag–Patterson algorithm](https://doi.org/10.1063/1.1692445)
- 2/3 de-aliasing
- [Pencil-based MPI parallelization](https://github.com/NaokiHori/SimpleDecomp) for scaling up to 10⁴ processes
- Fourth-order Runge-Kutta method (for nonlinear terms) combined with the integrating-factor technique (for linear terms) for temporal integration

Refer to the [documentation](https://naokihori.github.io/SpectralNSSolver1/) for details (currently under construction).

## Dependencies

- [C compiler](https://gcc.gnu.org)
- [GNU Make](https://www.gnu.org/software/make/)
- [MPI](https://www.open-mpi.org)
- [FFTW3](https://www.fftw.org)

To easily initialize the flow field, it is recommended to use:

- [Python](https://www.python.org) with [NumPy](https://numpy.org) (see `initial_condition/main.py`)

## Quick Start

1. **Set up your workspace**

   ```console
   mkdir -p /path/to/your/directory
   cd /path/to/your/directory
   ```

2. **Clone the repository**

   ```console
   git clone --recurse-submodules https://github.com/NaokiHori/SpectralNSSolver1
   cd SpectralNSSolver1
   ```

3. **Set the initial condition**

   The velocity field must be solenoidal, while the scalar field can be arbitrary. `main.py` provides several example initial conditions.

   ```console
   cd initial_condition
   python3 main.py 0
   cd ..
   ```

4. **Build the solver**

   ```console
   make clean
   make output
   make all
   ```

5. **Run the simulation**

   Execution parameters are defined in `exec.sh`.

   ```console
   bash exec.sh
   ```

   The runtime depends on your system specifications.

6. **Output and Visualization**

   The flow fields are stored in `output/save/` as [NPY files](https://numpy.org/devdocs/reference/generated/numpy.lib.format.html). These velocities are in the spectral domain, so an inverse Fourier transform (with normalization) is needed to obtain physical velocities.

   If the necessary Python libraries are installed, you can visualize the results with:

   ```console
   python3 visualise/2d.py
   ```

## 3D Simulation

To run a 3D simulation, switch to the `3d` branch and recompile all source files. You will also need to regenerate the initial flow field.

## Reference

- Canuto et al., *Spectral Methods - Fundamentals in Single Domains*, Springer
- Canuto et al., *Spectral Methods - Evolution to Complex Geometries and Applications to Fluid Dynamics*, Springer

## Acknowledgement

I would like to thank [Dr. Chris Howland](https://chowland.github.io) for valuable discussions.

