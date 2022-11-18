# Multigrid Solvers for the Poisson Equation in 3D
This project contains Matlab and C++ implementations of multigrid solvers for Dirichlet boundary condititions and mixed boundary conditions.

## Overview
The `matlab` subdirectory contains three implementations of V-cycles based on Jacobi iterations. `matlab/dirichlet_2d` and `matlab/dirichlet_3d` contains simple implementations for solving the Poisson equation in two and three spatial dimensions respectively with Dirichlet condition on all boundary points.<br>
The subdirectory `matlab/mixed_3d` contains a more advanced example with Dirichlet conditions on horizontal boundaries and Neumann conditions on all vertical boundaries. 

## The Poisson Equation
The Poisson equation is an elliptic partial differential equation, and the left-hand side is given by the Laplacian operator,
```math
\Delta u = f
```
It is assumed that it is a linear PDE, meaning that the right-hand side is only a function of the spatial variables.<br>
In the discretization it is assumed that the number of cells in each dimension is even, or equivalently that the number of lattices is odd.

## Test Problem
We may consider some trigonometric function
```math
u(x,y,x) = \sin(k_x x)\sin(k_y y)\sin(k_z z)
```
for which it is easy to derive the analytic Laplacian
```math
f = \Delta u = -(kx^2+ky^2+kz^2)u(x,y,x).
```