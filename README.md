This is a program written in MATLAB(TM) able to perform atomistic-to-continuum multiscale simulations of mechanical problems.
It uses the two-dimensional, finite-temperature version of the Coupled Atomistic and Discrete Dislocation (CADD) method. Thus, the atomistic region uses molecular dynamics and the continuum region uses the (linear) finite element method.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

Details of the program:

Everything in the code that is concerned with atoms and the coupling between the atomic and continuum regions was written from scratch by me. 
The program can model the atomistic region under constant energy conditions (NVE) or constant temperature (NVT).
To simulate constant temperature conditions, the Nose-Hoover and Langevin thermostat are included.

The program supports input from GMSH (www.gmsh.info) and .VTK output, to be displayed e.g. in ParaView(TM).

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

Limitations:
* The capabilities for dislocation detection, passing and evolution of discrete dislocation are limited:
  The CADD method uses a very sophisticated and complex approach to transfer dislocations between the atomistic and continuum region.
  This version includes the detection of dislocations in the atomistic region, a mechanism to pass the dislocation to the continuum region as a discrete dislocation, and the evolution of the discrete dislocations in the continuum.
  These capabilities were however only of secondary interest in the projects in which this program was used and only a single glide direction was studied. Thus, the detection, passing and evolution is currently limited to a single glide direction
  in a hexagonal material with c-axis normal to the plane. A generalization is possible.    

(c) Patrick Wurm

