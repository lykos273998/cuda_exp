# Study Experiments on CUDA

** Requires CUDA installation in order to compile, actually version > 5.0 with `nvcc`

All experiments are tested on my system, they are not optimized, they are study projects done for fun.

System: Intel Core i7 6700HQ + 8GB DDR4 + NVIDA 950M 4GB

Not the best but not the worst

## cudaMandelbrot 

Implementation of Mandelbrot Like set, compile it with `make`, or if you want code profiling use `make profile`. Then run the executable.

Outputs a `.ppm` greyscale image. 

Usage: `./cuM [OUTPUT FILE NAME] [SCALE] [CENTER-x] [CENTER-y] [HEIGHT] [WIDTH]`

## heat

Heat equation parabolic PDE numerical solution, very raw one. The goal was to gain confidence with `eigen` library for linear algebra. 

Euler method was used as integrator, in the current implementation as expected the method is unstable. The animaiton is made using `SFML` library. More info on it [here](https://www.sfml-dev.org/)
