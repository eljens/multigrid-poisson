# Matlab Implementations of Geometric Multigrid Poisson Solvers
This subdirectory contains three different implementations of multigrid poisson solvers.

## Pure Dirichlet Conditions
There are two different implementations for two and three spatial dimensions in `dirichlet_2d` and `dirichlet_3d` respectively. <br>
Note that when all boundaries are of the same type as in this case, it results in a symmetric discretization of the Laplacian operator. Thus the coarsest level is solved with an LDL<sup>T</sup> factorization as
```math
    (LDL^T)u=f.
```

## Mixed Boundary Conditions
The model in `mixed_3d` has Neumann bondary conditions on all vertical boundaries and Dirichlet type boundary conditions on all horizontal boundaries. This results in a more complicated algorithm since the Neumann conditions for the defect/residual are not the same on each level like in the pure Dirichlet case.
Further, the discrete Laplacian is asymmetric and the LDL<sup>T</sup> cannot be used. Instead, we use the LU factorization,
```math
    (LU)u=f.
```

## General Comments
The multigrid algorithms are based on the following properties/assumptions:
- The multigrid solvers are meant for linear elliptic PDEs only.
- The domain is discretized into an equidistant cartesian grid.
- Mixed boundary conditions are considered in the implementation in `mixed_3d`.
    - Inhomogeneous Neumann conditions are imposed on vertical boundaries.
    - Inhomogeneous Dirichlet conditions are imposed on horizontal boundaries.
- The ghost points have been eliminated in the discretization.
- The multigrid scheme is based on V-cycles.
- Injection is used as the restriction operator.
- Bi- or trilinear interpolation is used as interpolation operator in 2d and 3d respectively.

## Help
Matlab function headers has been written following the [standard](https://se.mathworks.com/help/simulink/mdl_gd/maab/na_0025matlabfunctionheader.html) such that `help <function>.m` will print documentation in a formatted way that is easier to read than the comments in the Matlab editor.