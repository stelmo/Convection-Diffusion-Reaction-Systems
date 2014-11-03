# Aris Taylor FEM approximation
# f_xx - a*f_x + a*b*f = 0
# f(x=0) = 1 and f_x(x=1) = 0
# The equations are dimensionless

using PyPlot

# Inputs
a = 7.5 # Peclet number
b = 1.29 # Damkohler number

N = 4
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

K = (-1.)*K1 - a*b*K2

B = zeros(N)
B[1] = -(-(-1./h - a/2.) - a*b*h/6.)

fs = \(K, B)

plot(xs, [1.,fs])

plt.show()
