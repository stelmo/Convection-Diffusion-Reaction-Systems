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

c12 = 1./h + (1./6.)*a*b*h
c11 = -2./h + (2./3.)*a*b*h
c13 = -1./h + (1./3.)*a*b*h

# Kx = 0
d = fill(c11, N)
d[end] = c13-a/2.
dl = fill(c12+a/2., N-1)
du = fill(c12-a/2., N-1)

K = Tridiagonal(dl, d, du)
B = zeros(N)
B[1] = c12+a/2.

fs = \(K, B)

plot(xs, [1.,fs])

plt.show()
