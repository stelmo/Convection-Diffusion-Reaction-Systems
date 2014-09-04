# Convection Diffusion Reaction Partial Differential Equation Model

using PyPlot
plt.ioff()

include("helper.jl")
include("fem1D.jl")
include("fdm1D.jl")

# One Dimensional Case

a = 5.0 # convection coefficient
b = 1.0 # diffusion coefficient
c = 0.0 # Dissipation

x0 = 0.0 # space boundary
xL = 10.0 # space boundary

t0 = 0.0 # time boundary
tL = 1.0 # time boundary

Nx = 100 # Number of space nodes
Nt = 250 # Number of time nodes

femImp = Helper.Model(x0, xL, t0, tL, Nx, Nt, a, b, c, "FEM Implicit")
Helper.enforceBCs!(femImp)
Fem1Du.solveImplicit!(femImp)

femExp = Helper.Model(x0, xL, t0, tL, Nx, Nt, a, b, c, "FEM Explicit")
Helper.enforceBCs!(femExp)
Fem1Du.solveExplicit!(femExp)

fdmFTCS = Helper.Model(x0, xL, t0, tL, Nx, Nt, a, b, c, "FDM FTCS")
Helper.enforceBCs!(fdmFTCS)
Fdm1Du.solveFTCS!(fdmFTCS)

fdmFTCSup = Helper.Model(x0, xL, t0, tL, Nx, Nt, a, b, c, "FDM Upwind FTCS")
Helper.enforceBCs!(fdmFTCSup)
Fdm1Du.solveFTCSup!(fdmFTCSup)


## 3D plotting
p1 = plot_wireframe(femImp.ts, femImp.xs, femImp.us, color="b")
p2 = plot_wireframe(femExp.ts, femExp.xs, femExp.us, color="r")
p3 = plot_wireframe(fdmFTCS.ts, fdmFTCS.xs, fdmFTCS.us, color="g")
p4 = plot_wireframe(fdmFTCSup.ts, fdmFTCSup.xs, fdmFTCSup.us, color="y")

legend([p1, p2, p3, p4], ["Implicit FEM","Explicit FEM", "FTCS", "Upwind FTCS"])
xlabel("Time")
ylabel("Space")
plt.show()



