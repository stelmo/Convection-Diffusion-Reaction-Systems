# Aris Taylor FEM approximation
# f_xx - a*f_x + a*b*f = 0
# f(x=0) = 1 and f_x(x=1) = 0
# The equations are dimensionless

using PyPlot

# Inputs
a = 7.5 # Peclet number
b = 5.129 #  L/U

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

for newton=1:4

  dm = fill(2./h, N)
  dm[1] = 1./h  - a/2.
  dm[end] = 1./h + a/2.

  dl = fill(-1./h - a/2., N-1)
  du = fill(-1./h + a/2., N-1)

  K1 = diagm(dl,-1) + diagm(dm, 0) + diagm(du, 1) #what a pity we can't use Tridiagonal like this :(
  # Now modify K1 to conform
  K1 = -1.*K1
  K1 = K1[2:end,:]
  B1 = K1[:,1]
  K1 = K1[:,2:end]

  dm2 = fill(2./3.*h, N)
  dm2[1] = 1./3.*h
  dm2[end] = 1./3.*h

  dl2 = fill(1./6.*h, N-1)
  du2 = fill(1./6.*h, N-1)

  K2 = diagm(dl2,-1) + diagm(dm2, 0) + diagm(du2, 1)
  K3 = diagm(dl2,-1) + diagm(dm2, 0) + diagm(du2, 1)

  for j=1:N
    K2[:,j] = -a*b*kT(u0[j])
    K3[:,j] = -a*b*kTd(u0[j])*u0[j]
  end

  K2 = K2[2:end,:]
  B2 = K2[:,1]
  K2 = K2[:,2:end]

  K3 = K3[2:end,:]
  B3 = K3[:,1]
  K3 = K3[:,2:end]

  K = K1 + K2 + K3
  Bleft = B1 + B2 + B3


# plt.show()
