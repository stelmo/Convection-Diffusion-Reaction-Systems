# Convergence results of problem 1
# Aris Taylor Analysis
# Unsteady Problem
# Uncoupled One Dimensional

using PyPlot

function p1_1d(N, Nt)

  # Inputs
  a = 7.5 # Peclet number
  b = 1.29 # Damkohler number
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
  B[1] = (-(-1./h - a/2.) - a*b*h/6.)

  m = fill((2./3.)*h, N)
  m[end] = (1./3.)*h

  md = fill((1./6.)*h, N-1)

  M = Tridiagonal(md, m, md)
  M = a*M
  # Now start the time series
  ts = linspace(0., 3., Nt) # this is dimensionless time... It also takes long...
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
  return us[end, :]
end

Nt = 250
ts = linspace(0., 3., Nt)
elements = [4:1:15]
elen = length(elements)
change = zeros(elen-1, Nt)
for k=1:elen-1
  change[k,:] = abs(p1_1d(elements[k+1], Nt) - p1_1d(elements[k], Nt))
end

contourf(ts,elements[2:end], change, 20)
xlabel("Time")
ylabel("Elements")
colorbar(label="Absolute Difference")
