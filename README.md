HyPar - Hyperbolic-Parabolic Partial Differential Equations Solver
------------------------------------------------------------------

`HyPar` is a finite-difference algorithm to solve hyperbolic-parabolic partial differential
equations (with source terms) on Cartesian grids. It is a unified framework that can handle 
systems of PDEs with arbitrary number of spatial dimensions and solution components. It 
provides the spatial discretization and time integration functions, functions to read and 
write solutions from/to files, as well as functions required to solve the system on parallel 
(MPI) platforms. The physical models define the physics-specific functions such as the exact 
forms of the hyperbolic flux, parabolic flux, source terms, upwinding functions, etc.

This fork is part of my master's dissertation project at the University of Warwick titled **Converting HyPar to OP2 eDSL** under the supervision of (https://warwick.ac.uk/fac/sci/dcs/people/Gihan_Mudalige/)[Gihan Mudalige].

Documentation
-------------

See documentation on how to compile and run along with examples:
http://hypar.github.io/


