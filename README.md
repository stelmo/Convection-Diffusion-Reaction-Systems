Convection-Diffusion-Reaction-Systems
=====================================

## The analysis and implementation of coupled convection diffusion reaction systems with the finite element method.

The implementation of the project is in Julia 0.3. The following packages are required:
1. Gadfly -> For static visualisation
2. Interact -> For interactive visualisation
3. IJulia -> To use the IPython Notebook interface (makes everything look pretty)
4. PyPlot -> For 3D visualisations (not as nice as the Gadfly visualisation implemented in the notebook)

To use:

Launch the IJulia Notebook interface and run the notebook "Convection Diffusion Reaction Uncoupled". 
It is possible to change the forcing function and all the initial conditions but, if you don't modify the
auxiliary modules, you can only modify the parameters shown in the notebook directly. 

Done: Implement an uncoupled linear convection diffusion reaction parabolic equation with constant coefficients 
and Von Neumann boundary conditions.

To do: Implement the FEM solution of a system of coupled nonlinear convection diffusion reaction equations. 


