# Aris Taylor Analysis
# Unsteady Problem
# Uncoupled One Dimensional

using PyPlot

# Inputs
L = 6.36
R = 0.05
DAB = 7.6E-05
k = 0.25
U = 1.24
Da = DAB + U^2*R^2/(48*DAB)

a = U*L/Da
b = k*L/U


N = 100
xs = linspace(0., 1., N+1)
h = xs[2] - xs[1]


d = fill(2./h, N)
d[end] = 1./h + a/2.

dl = fill(-1./h - a/2., N-1)
du = fill(-1./h + a/2., N-1)

K1 = Tridiagonal(dl, d, du)

d2 = fill(2./3.*h, N)
d2[end] = 1./3.*h

dl2 = fill(1./6.*h, N-1)
du2 = fill(1./6.*h, N-1)

K2 = Tridiagonal(dl2, d2, du2)

K = (-1.)*K1 - b*K2

B = zeros(N)
B[1] = (-(-1./h - a/2.) - b*h/6.)

m = fill((2./3.)*h, N)
m[end] = (1./3.)*h

md = fill((1./6.)*h, N-1)

M = Tridiagonal(md, m, md)
M = a*M
# Now start the time series
Nt = 50
ts = linspace(0., 1., Nt) # this is dimensionless time... It also takes long...
dt = ts[2] - ts[1]

us = zeros(N+1, Nt)
finit(x) = exp(-x*10)#4.0.*x.^2 -4.0.*x +1 # # Initial condition
us[:,1] = finit(xs) # Initial condition
us[1, :] = 1. # Boundary condition

KK1 = M - (dt/2.)*K
KK2 = M + (dt/2.)*K


for (tind, t) in enumerate(ts[2:end])

  us[2:end, tind+1] = \(KK1, KK2*us[2:end, tind] + dt*B)

end


surf(xs, ts, us', linewidth=0.0)
xlabel("Space")
ylabel("Time")
zlabel("Concentration")

# plot(xs, us[:, end])

plt.show()
