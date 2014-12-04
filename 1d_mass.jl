using Images, ImageView

a = 7.5 # Peclet number
b = 1.29 # Damkohler number

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

K = (-1.)*K1 - a*b*K2
K = full(K)

for k=1:length(K)
  if K[k] != 0.0
    K[k] = 0.0
  else
    K[k] = 1.0
  end
end
imwrite(K, "p1_mass.pdf")
