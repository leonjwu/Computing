# Computing Projects

## Credibility: Some projects below were part of Year 3 Imperial College Mathematics High-Performance Computing Module, in which I ranked 1st in the year


---
# [Project 1: Simulating Fluid Flow PDEs with Finite Differences and Simulating Weakly-Coupled Osicalltors on a network (Parallel Fortran (OMP and MPI)](https://github.com/leonwu4951/Computing/blob/master/Fluid-Oscillators/)

## Overview
- Simulated blood flow through a deformed artery
- Implemented finite difference methods in Python
- Designed fast sparse matrix solvers in parallel Fortran code (OMP) to speedup the bottleneck in the finite-difference schemes by 100x
- Simulated Weakly-Coupled Osicalltors on a network using 1D and 2D decompositions for MPI parallel Computing in Fortran.
- Setup communication effectively between multiple processes on a grid using appropriate MPI directives
 

Fortran vs Python Speedup  |  Number of Threads Speedup
:-------------------------:|:-------------------------:
![](https://github.com/leonwu4951/Computing/blob/master/Fluid-Oscillators/speedup.png)  |  ![](https://github.com/leonwu4951/Computing/blob/master/Fluid-Oscillators/threads.png)
---



---
# [Project 3: Simulating Bacteria Interaction Dynamics in Python and Parallel Fortran (OMP)](https://github.com/leonwu4951/Computing/blob/master/Bacteria/)

## Overview
- Developed vectorised Python code using NumPy and SciPy for bacteria simulations
- Analysed dynamics of the simulations with various plots
- Developed parallel Fortran code using OMP for the simulations, resulting in a 100x speedup over the vectorised Python code through better memory management and parallelisation

Parallel Fortran Speedup  |  Speedup Python vs Fortran 
:-------------------------:|:-------------------------:
![](https://github.com/leonwu4951/Computing/blob/master/Bacteria/speedup.png)  |  ![](https://github.com/leonwu4951/Computing/blob/master/Bacteria/speedup2.png)
---




