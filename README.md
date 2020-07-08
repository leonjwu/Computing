# Computing Projects

## Credibility: Some projects below were submitted as part of the Year 3 Imperial College Mathematics High-Performance Computing Module, in which I ranked 1st in the year



# [Project 1: Computational PDEs: 1D and 2D Wave Equations](https://leonjwu.github.io/comp-pdes/)

## Overview
- Implemented finite difference schemes for 1D and 2D Wave Equations
- Designed boundary conditions to allow for reflecting and non-reflecting boundaries (Neumann and Dirichlet)
- Created visualisations of 1D and 2D waves using Matplotlib Animations


2D Wave - Single  |  2D Wave - Multi
:-------------------------:|:-------------------------:
![](https://github.com/leonjwu/comp-pdes/blob/master/2D_Wave_Symetric_Single.gif)  |  ![](https://github.com/leonjwu/comp-pdes/blob/master/2D_Wave_Symmetric_Multi.gif)
---



---
# [Project 2: Simulating Fluid Flow PDEs with Finite Differences and Simulating Weakly-Coupled Oscillattors on a network (Parallel Fortran (OMP and MPI))](https://github.com/leonjwu/Computing/blob/master/Fluid-Oscillators/)

## Overview
- Simulated blood flow through a deformed artery
- Designed fast sparse matrix solvers in parallel Fortran code (OMP) to speedup the bottleneck in the finite-difference schemes by 100x
- Simulated Weakly-Coupled Osicalltors on a network using 1D and 2D decompositions for MPI parallel Computing in Fortran.
- Setup communication effectively between multiple processes on a grid using appropriate MPI directives
- Implemented finite difference methods in Python

 

Fortran vs Python Speedup  |  Number of Threads Speedup
:-------------------------:|:-------------------------:
![](https://github.com/leonjwu/Computing/blob/master/Fluid-Oscillators/speedup.png)  |  ![](https://github.com/leonjwu/Computing/blob/master/Fluid-Oscillators/threads.png)
---









# [Project 3: Scientific Computing: Multidimensional SVD, Image Repair Algorithm, Chaotic Systems and Path-Finding Algorithms (Python)](https://github.com/leonjwu/Computing/blob/master/Scientific-Computing/)

## Overview
- Designed multidimensional SVD techniques to reduced dataset storage by 95%
- Optimised code for repairing images through vectorisation
- Reformulated matrix computations for most efficient use of Python vectorisation
- Simulated Chaotic Bacteria Interactions through PDEs and Finite Differences
- Implemented banded matrix solvers using SciPy to massively speed up simulations 
- Designed fast path-finding algorithms utilising dequeues, hash tables and binary heaps
- Modified BFS, DFS and Dijkstra algorithms for specific path-finding problems



Chaotic System Contours  |  Repaired Image
:-------------------------:|:-------------------------:
![](https://github.com/leonjwu/Computing/blob/master/Scientific-Computing/Chaos.png)  |  ![](https://github.com/leonjwu/Computing/blob/master/Scientific-Computing/Repaired_Image.png)
---




---


