# Data-Driven Diffusion Mechanics 

## Context
This repository contains Matlab implementation of data-driven mechanics for steady-state and transient mass diffusion problems formulated in terms of chemical potential field. A specific class of data-driven approach is adopted here called *distance-minimizing data-driven solvers*. For more details on the subject the user of the code is directed to the papers in the literature.

- Poineer work in the data-driven mechanics https://doi.org/10.1016/j.cma.2016.02.001
- Extension to the dynamic problems https://doi.org/10.1002/nme.5716
- Variational framework of data-driven problems https://doi.org/10.1016/j.cma.2020.112898

## Features
This code provides

- A very clean implementation of distance-minimizing data-driven solver. 
- A simple 2D finite element implementation for scalar fields (which can easily be extended to 3D finite elements as well). 
- Pre- and post-processing parsing functions for Gmsh and Paraview.

Author is fully aware of the fact that this code is a very naive approach towards a very complex subject (finite-elements) in hand and there exist much better and robust implementations and open source packages e.g., [deal.II](https://www.dealii.org/) and [FEniCS Project](https://fenicsproject.org/). However, reporting of any comment, found bugs, improvements and dicussions on the theoratical and implementation aspects are highly appreciated by the users.
You can contact the author at  engineerabdullah@ymail.com . 

## Disclaimer
This code is purely for educational purposes. All rights are preseved, however, author shall not be liable in any event caused by the use of the code.
