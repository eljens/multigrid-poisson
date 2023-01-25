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
to make the archibe `libpoisson.a`. To use this archive in your project, compile with `-Iinlcude -Llib -lpoisson`.<br><br>
One can select between three test problems if one wants to use one of the drivers. The test problems can be selected by compiling with `PROBLEM=1`, `PROBLEM=2`, or, `PROBLEM=3`. 

## The Array Class and the Device Array Class
The most important class templates in this project are probably `array.h` and the derived template class `devicearray.h`. They contain functions for indexing 3D arrays that abstracts away the halos. Futher, they contain functionality for mapping arrays to GPUs and saving arrays with or without the halos. See for instance the following figure:

Saved With Halo            | Saved Without Halo
:-------------------------:|:-------------------------:
![Halo](figures/halo.png)  | ![Ignoring Halo](figures/no_halo.png)

## Test Problems
To generate a number of test problems, we may consider some function $u$ and find the Laplacian $f=\Delta u$ and the first order derivatives at the boundaries $\frac{\partial u}{\partial x}$ and $\frac{\partial u}{\partial y}$
### Test Problem I
We may consider some trigonometric function
```math
u(x,y,x) = \sin(k_x x)\sin(k_y y)\sin(k_z z)
```
for which it is easy to derive the analytic Laplacian
```math
f = \Delta u = -(k_x^2+k_y^2+k_z^2)u(x,y,x)
```
and the Neumann boundary conditions are given by
```math
\frac{\partial u}{\partial x} = k_x\cos(k_xx)\sin(k_yy)\sin(k_zz)
```
```math
\frac{\partial u}{\partial y} = k_y\sin(k_xx)\cos(k_yy)\sin(k_zz).
```
True Solution              | Numerical Solution
:-------------------------:|:-------------------------:
![Problem I reference solution](figures/problem_1_true.png)  | ![Problem I numerical solution](figures/problem_1.png)

### Test Problem II
The polynomial
```math
u(x,y,x) =  x^3y^2z
```
has the Laplacian
```math
f = \Delta u = 2x^3z + 6xy^2z.
```
The Neumann boundary conditions are given by
```math
\frac{\partial u}{\partial x} = 3x^2y^2z
```
```math
\frac{\partial u}{\partial y} = 2x^3yz
```
True Solution              | Numerical Solution
:-------------------------:|:-------------------------:
![Problem II reference solution](figures/problem_2_true.png)  | ![Problem II numerical solution](figures/problem_2.png)

### Test Problem III
The trigonometric funtion
```math
u(x,y,x) = \cos(xz^2)\sin(y^3)
```
has the Laplacian
```math
f = \Delta u = ((-4x^2z^2 - 9y^4 - z^4)\sin(y^3) + 6y\cos(y^3))\cos(xz^2) - 2x\sin(xz^2)sin(y^3).
```
The Neumann boundary conditions are given by
```math
\frac{\partial u}{\partial x} = -z^2\sin(xz^2)\sin(y^3)
```
```math
\frac{\partial u}{\partial y} = 3y^2\cos(xz^2)\cos(y^3)
```
True Solution              | Numerical Solution
:-------------------------:|:-------------------------:
![Problem III reference solution](figures/problem_3_true.png)  | ![Problem III numerical solution](figures/problem_3.png)

## Verifying the Order of Accuracy
From theory, the scheme should be second order convergent. That means that halving the grid spacing should result in a four times smaller maximal absolute error. As it can be seen, the multigrid algorithm has the expected order of convergence for the three test problems. The files used to generate the convergence tests are `convergence.sh` and `plot_convergence.m`.
<p align="center">
  <img src="./figures/mg_convergence.png" width="400" title="Convergence as function of grid spacing.">
</p>
