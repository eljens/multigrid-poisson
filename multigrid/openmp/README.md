# Multigrid Methods for the Poisson Equation
This project contains a framework for finite difference multigrid methods. The header files in `include` contain class templates for a number of functionalities. 

## The Array Class and the Device Array Class
The most important class templates in this project are probably `array.h` and the derived template class `devicearray.h`. They contain functions for indexing 3D arrays that abstracts away the halos. Futher, they contain functionality for mapping arrays to GPUs and saving arrays with or without the halos. See for instance the following figure:

Saved With Halo            |  Saved Without Halo
:-------------------------:|:-------------------------:
![Halo](figures/halo.png)  |  ![Ignoring Halo](figures/no_halo.png)