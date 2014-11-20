# Aris Taylor FEM approximation
# f_xx - a*f_x + a*b*f = 0
# f(x=0) = 1 and f_x(x=1) = 0
# The equations are dimensionless

using PyPlot

# Inputs
a = 1. # 7.5 # Peclet number
b = 1. #5.129 #  L/U

N = 5 # number of nodes
xs = linspace(0., 1., N)
h = xs[2] - xs[1]

# Functions
kT(x) = 233253*exp(-4881*(1./(300+100*(1-x))))
kTd(x) = kT(x)*(-4811*((400-100*x)^(-2))*100)

# Guess a solution
u0 = zeros(N)
u0[1] = 1.
u0[2] = 0.6


# Symmetric diffusion matrices
dm = fill(2./h, N)
dm[1] = 1./h
dm[end] = 1./h

dl = fill(-1./h, N-1)
du = fill(-1./h, N-1)

K1 = diagm(dl,-1) + diagm(dm, 0) + diagm(du, 1)
# Now modify K1 to conform
K1 = K1[2:end, :]
K3 = copy(K1)
K1 = -1.*K1

# Anti-symmetric convection matrices
dm2 = fill(0., N)
dm2[1] = -0.5
dm2[end] = 0.5

dl2 = fill(-0.5, N-1)
du2 = fill(0.5, N-1)

K2 = diagm(dl2,-1) + diagm(dm2, 0) + diagm(du2, 1)
# Now modify K1 to conform
K2 = K2[2:end, :]
K4 = copy(K2)
K4 = a*K4
K2 = -a*K2

# Non-linear reaction matrices
dm3 = fill(2.*h/3., N)
dm3[1] = h/3.
dm3[end] = h/3.

dl3 = fill(h/6., N-1)
du3 = fill(h/6., N-1)

K5 = diagm(dl3,-1) + diagm(dm3, 0) + diagm(du3, 1)
K5 = K5[2:end, :]




















# plt.show()
