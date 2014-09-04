# Convection Diffusion Reaction Partial Differential Equation Model

using PyPlot
plt.ioff()

include("helper.jl")
include("fem1D.jl")
include("fdm1D.jl")

# One Dimensional Case

a = 5.0 # convection coefficient
b = 1.0 # diffusion coefficient
c = 0.0 # dispersion coefficient

x0 = 0.0 # space boundary
xL = 10.0 # space boundary

t0 = 0.0 # time boundary
tL = 1.0 # time boundary

Nx = 100 # Number of space nodes
Nt = 250 # Number of time nodes

femImp = Helper_u.Model(x0, xL, t0, tL, Nx, Nt, a, b, c, "FEM Implicit")
Helper_u.enforceBCs!(femImp)
Fem1D_u.solveImplicit!(femImp)

femExp = Helper_u.Model(x0, xL, t0, tL, Nx, Nt, a, b, c, "FEM Explicit")
Helper_u.enforceBCs!(femExp)
Fem1D_u.solveExplicit!(femExp)

fdmFTCS = Helper_u.Model(x0, xL, t0, tL, Nx, Nt, a, b, c, "FDM FTCS")
Helper_u.enforceBCs!(fdmFTCS)
Fdm1D_u.solveFTCS!(fdmFTCS)

fdmFTCSup = Helper_u.Model(x0, xL, t0, tL, Nx, Nt, a, b, c, "FDM Upwind FTCS")
Helper_u.enforceBCs!(fdmFTCSup)
Fdm1D_u.solveFTCSup!(fdmFTCSup)


## 3D plotting
p1 = plot_wireframe(femImp.ts, femImp.xs, femImp.us, color="b")
p2 = plot_wireframe(femExp.ts, femExp.xs, femExp.us, color="r")
p3 = plot_wireframe(fdmFTCS.ts, fdmFTCS.xs, fdmFTCS.us, color="g")
p4 = plot_wireframe(fdmFTCSup.ts, fdmFTCSup.xs, fdmFTCSup.us, color="y")

legend([p1, p2, p3, p4], ["Implicit FEM","Explicit FEM", "FTCS", "Upwind FTCS"])
xlabel("Time")
ylabel("Space")
plt.show()



