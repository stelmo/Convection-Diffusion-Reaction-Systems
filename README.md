Convection-Diffusion-Reaction-Systems
=====================================

### The analysis and implementation of coupled and uncoupled convection diffusion reaction systems with the finite element method.
All the programs were written in the Julia Language.

The PyPlot package is required for the graphical output. The packages: Images and ImageView are required to view the PDF mass matrices as seen in the report.

## Guide to the programs
All programs are stand alone. There is some room for improvement in programming style...

1. aris_fem_us_1d.jl = One dimensional isothermal system unsteady.

2. 2d_fem_rf_us.jl = Two dimensional isothermal system unsteady.

3. adi_1d_us.jl = Adiabatic one dimensional system unsteady.

4. coupled_1d_us.jl = Coupled isothermal one dimensional system unsteady.

Note: If a program's name includes "us" it implies it is the unsteady solution. The same name without the "us" is the steady state solution.

Note: If a program's name includes "rf" it implies it is the laminar flow profile solution. The same name without the "rf" is the plug flow solution. This applies only to the 2 dimensional problems.

Note: The other programs were used to generate the convergence results etc. but consist only of the programs mentioned above with some formatting.
