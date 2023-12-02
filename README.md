# lbmtau1

lbmtau1 is the high-performance implementation of a recently proposed formulation of the Lattice Boltzmann Method (LBM) with a fixed relaxation time τ = 1. An important advantage of the presented algorithm is that its performance (measured in fluid lattice updates per second, FLUPS) is high (up to over four times higher than standard LBM for D3Q19 scheme) and almost constant regardless of lattice arrangements whereas, for standard LBM, the performance drops with the number of lattice links. For the double-precision D3Q19 scheme, we achieved 497 MFLUPS on the system based on Ryzen 9 7945HX processor with 45 GB/s real memory bandwidth, resulting in 11 MFLUPS per GB/s of available memory bandwidth, whereas the best implementations for standard LBM achieve up to 3 MFLUPS per GB/s.

## Usage

lbmtau1 is designed as a single C++ header file (``LBMTau1.hpp``).
It requires `clang` with OpenMP support (tested on AMD clang version 14.0.6).
Additionally, the example ``main.cpp`` file is provided showing simple simulation setup (channel with randomly placed spherical obstacles) and the appropriate ``Makefile``.
File ``run_measurements.sh`` is a simple bash script for performance measurements (complete run takes many hours).

### Compilation

Simply run

  `make`

You may need to tune clang options inside ``Makefile`` for your system.

#### Different kernel versions

lbmtau1 contains four versions of kernels:
- manually optimised vector version resulting in the highest performance,
- automatically vectorized ``for`` loop (slightly slower than manual version, but significantly simpler),
- simple scalar reference code (very slow, for reference only),
- memory copy kernel without computations used for performance comparison.

The version of kernel can be set at the beginning of ``Makefile`` by setting `OPTIONS` to one of the available values.
The manually optimized kernel is used by default.

### Launching

After successfull compilation, run the ``main`` file

  `./main`

The simulation result will be saved to ``.vtk`` files that can be processed with [ParaView](https://www.paraview.org/).

#### Parameters

Up to four parameters can be passed to main:
1. The number of time steps computed per single iteration. After computing these time steps, the performance is measured and reported.
2. The number of measures. Each measure computes the number of time steps given as the first parameter.
3. The number of LBM nodes per width and depth of the simulated channel. The length of the channel is twice as long. ***WARNING! The number of nodes must be divisible by 4!***
4. The number of spherical obstacles placed in the channel.

For example, 1000 time steps (with performance reports after each 100 steps) for channel containing 160x80x80 nodes and 30 obstacles can be run with
  
  ``./main 100 10 80 30``
  
