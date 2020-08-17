This is a program written in MATLAB(TM) able to perform atomistic-to-continuum multiscale simulations of mechanical problems.
It uses the two-dimensional, finite-temperature version of the Coupled Atomistic and Discrete Dislocation (CADD) method. Thus, the atomistic region uses molecular dynamics and the continuum region uses the (linear) finite element method.

One example (a two-dimensional tensile test) at various temperatures is included, which can be run by running the script "run_Tensile_Test.m".

The program
* can model the atomistic region under constant energy conditions (NVE) or constant temperature (NVT) using a Nose-Hoover or a Langevin thermostat.
* allows the user to choose from two interatomic potentials (Lennard-Jones and EAM).
* Allows the user to choose from a quasi-static, dynamic or hybrid (a complementary superposition of a quasi-static and a dynamic subproblem) linear finite element model
* supports input from GMSH (www.gmsh.info)
* supports .VTK output, to be displayed e.g. in ParaView(TM).
* supports the detection of dislocation, the passing of dislocations to the continuum as discrete dislocations and the evolution of the discrete dislocations
(Although these capabilities are limited to a single glide system. A generalization is possible but was not needed in the project to which this program was applied)

I wrote everything in the code that is concerned with atoms and the coupling between the atomic and continuum regions from scratch. 
To simulate the continuum, I heavily modified a finite element code (soofeaM) developed at our institute by former colleagues and myself (see below).

soofeaM - Software for Object Oriented Finite Element Analysis with Matlab
(c) 2008-2013: Michael Hammer
(c) 2017-2020: Benedikt Weger, Patrick Wurm

soofeaM is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Here is a short summary of the differences between the current version of soofeaM and the standard soofeaM.
This version:
* is much more optimized and is therefore not well suited as a learning tool
* uses velocity verlet explicit time integration with a lumped mass
  matrix for extra efficiency
* is a linear finite element version with a vectorized version of the material routine/internal node force 
  evaluation 
* does only support displacements as boundary conditions (it would be straightforward
  to include nodal forces, volume forces, ... just like in standard soofeaM,
  but they were simply not needed for the project and therefore not included)
  
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
(c) Patrick Wurm

