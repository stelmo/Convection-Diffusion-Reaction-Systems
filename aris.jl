# Aris Taylor Approximation
# f_xx - a*f_x + a*b*f = 0
# f(x=0) = 1 and f_x(x=1) = 0
# The equations are dimensionless

using PyPlot

# Inputs
a = 7.5 # Peclet number
b = 1.29 # Damkohler number

r1 = 0.5*(a + sqrt(a^2 + 4.0*a*b))
r2 = 0.5*(a - sqrt(a^2 + 4.0*a*b))

if typeof(r1) != Float64 || typeof(r2) != Float64
  println("Error: Aris Taylor system has complex roots!")
end

A = [1 1;r1*exp(r1) r2*exp(r2)]
B = [1; 0]
cs = \(A, B) # the constants

f(x) = cs[1].*exp(r1.*x) + cs[2].*exp(r2.*x)

xs = linspace(0.,1.,100)
fs = f(xs)
plot(xs, fs)

plt.show()
