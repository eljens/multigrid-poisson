# Multigrid Methods for the Poisson Equation
This project contains a framework for finite difference multigrid methods. The header files in `include` contain class templates for a number of functionalities, the `src` directory contains the source files, and the `drivers` directory contains examples of how to use the solvers.

## Compilation
To compile one of the examples, you may type
```{bash}
make GPU_ARCH=<sm_70,sm_80> APP=<minimal_jacobi,minimal_multigrid>
```
and an example can for instance be executed with
```{bash}
./drivers/minimal_multigrid -x 513 -y 257 -z 129 -l 6 -maxiter 10
```
To compile only the archive, type
```{bash}
make archive
```
to make the archibe `libpoisson.a`. To use this archive in your project, compile with `-Iinlcude -Llib -lpoisson`.

## The Array Class and the Device Array Class
The most important class templates in this project are probably `array.h` and the derived template class `devicearray.h`. They contain functions for indexing 3D arrays that abstracts away the halos. Futher, they contain functionality for mapping arrays to GPUs and saving arrays with or without the halos. See for instance the following figure:

Saved With Halo            | Saved Without Halo
:-------------------------:|:-------------------------:
![Halo](figures/halo.png)  | ![Ignoring Halo](figures/no_halo.png)

## Verifying the Order of Accuracy
From theory, the scheme should be second order convergent. That means that halving the grid spacing should result in a four times smaller maximal absolute error. As it can be seen, the multigrid algorithm has the expected order of convergence. The test was performed for the mixed boundary problem.
<p align="center">
  <img src="./figures/mg_convergence.png" width="400" title="Convergence as function of grid spacing.">
</p>
