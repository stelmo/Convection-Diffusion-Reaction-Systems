Convection-Diffusion-Reaction-Systems
=====================================

### The analysis and implementation of coupled and uncoupled convection diffusion reaction systems with the finite element method.
All the programs were written in the Julia Language.

# Uncoupled System
The reactor is assumed to be isothermal.

## One Dimensional Analysis
The elements are comprised of piecewise linear basis functions (roof functions). The spatial problem was solved using the finite element method. The unsteady (temporal) problem was solved using a combination of the FEM and the Crank-Nicholson method (to solve the system of resultant ODEs).

### Steady State System
The Aris-Taylor approximation is used to model a reacting system with radial and axial diffusion, and axial convection. First order kinetics are used with a constant rate.

### Unsteady System
The unsteady case was an extension of the steady case. The initial condition approximates an inert reactor.

## Two Dimensional Analysis
Piecewise bilinear basis functions (rectangular elements) are used as elements in the spatial dimensions. Again the Crank-Nicholson method was used to solve the system of resultant ODEs for the unsteady problem.

### Steady State System
Axial and radial diffusion with axial convection was considered. First order kinetics are used. Both plug flow and laminar flow were implemented.

### Unsteady System
The unsteady case was again just an extension of the steady case. An inert reactor was approximated for the initial condition. Both plug flow and laminar flow were implemented.

# Coupled Systems
The isothermal assumption is dropped. The same systems as above, but this time incorporating a non-constant rate term which depends on temperature were modelled. Only the one dimensional case was modelled.
