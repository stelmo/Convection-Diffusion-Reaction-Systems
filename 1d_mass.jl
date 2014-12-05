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

n = int(sqrt(length(K)))
b0 = zeros(3,3)
b1 = diagm([1,1,1.])

newK = zeros(3*n, 3*n)
counter = 1

for j=1:3:3*n
  for i=1:3:3*n
    if K[counter] != 0.
      newK[i:i+2, j:j+2] = b1
    else
      newK[i:i+2, j:j+2] = b0
    end
    counter += 1
  end
end
for k=1:length(newK)
  if newK[k] == 0.0
    newK[k] = 1.
  else
    newK[k] = 0.
  end
end
imwrite(newK, "p1_mass.pdf")
