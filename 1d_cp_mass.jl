using Images, ImageView
# Inputs
a = 7.5 # Peclet number
Ca0 = 1.0
Cb0 = 0.5
b1 = 6.36^2 * Cb0 * 0.25 / 1.05
b2 = 6.36^2 * Ca0 * 0.25 / 1.05


N = 100 # number of nodes - need many of these
xs = linspace(0., 1., N)
h = xs[2] - xs[1]

# Guess a solution
u0 = zeros(N)
u0[1] = 1.0
for i=2:N
  u0[i] =u0[i-1]*0.8
end

w0 = zeros(N)
w0[1] = 1.0
for i=2:N
  w0[i] =w0[i-1]*0.8
end

# Symmetric diffusion matrices
dm = fill(2./h, N)
dm[1] = 1./h
dm[end] = 1./h

dl = fill(-1./h, N-1)
du = fill(-1./h, N-1)

K1 = diagm(dl,-1) + diagm(dm, 0) + diagm(du, 1)
K11 = copy(-K1[2:end, :])
K21 = copy(-K1[2:end, :])
B11 = copy(K1[2:end, :])
B21 = copy(K1[2:end, :])

# Anti-symmetric convection matrices
dm2 = fill(0., N)
dm2[1] = -0.5
dm2[end] = 0.5

dl2 = fill(-0.5, N-1)
du2 = fill(0.5, N-1)

K2 = diagm(dl2,-1) + diagm(dm2, 0) + diagm(du2, 1)
K12 = copy(-a*K2[2:end, :])
K22 = copy(-a*K2[2:end, :])
B12 = copy(a*K2[2:end, :])
B22 = copy(a*K2[2:end, :])

K0 = zeros(N-1, N)

newtonFlag = true
counter = 0
alpha = 0.8 # slow down NR - dampening
maxIter = 1000


# Semi-symmetric reaction matrices - depends on the interation number!
dl3 = zeros(N-1)
for i=1:N-1
  dl3[i] = (h/12.)*(u0[i]+u0[i+1])
end

du3 = zeros(N-1)
for i=1:N-1
  du3[i] = (h/12.)*(u0[i]+u0[i+1])
end

dm3 = zeros(N)
dm3[1] = -9. # it gets cut away so w/e
for i=2:N-1
  dm3[i] = (u0[i-1]*h/12. + u0[i]*h/2. + u0[i+1]*h/12.)
end
dm3[end] = (u0[end-1]*h/12. + u0[end]*h/4.)

K3 = diagm(dl3,-1) + diagm(dm3, 0) + diagm(du3, 1)
K13 = copy(-b1*K3[2:end, :])
K23 = copy(-b2*K3[2:end, :])

dl4 = zeros(N-1)
for i=1:N-1
  dl4[i] = (h/12.)*(w0[i]+w0[i+1])
end

du4 = zeros(N-1)
for i=1:N-1
  du4[i] = (h/12.)*(w0[i]+w0[i+1])
end

dm4 = zeros(N)
dm4[1] = -9. # it gets cut away so w/e
for i=2:N-1
  dm4[i] = (w0[i-1]*h/12. + w0[i]*h/2. + w0[i+1]*h/12.)
end
dm4[end] = (w0[end-1]*h/12. + w0[end]*h/4.)

K4 = diagm(dl4,-1) + diagm(dm4, 0) + diagm(du4, 1)
K14 = copy(-b1*K4[2:end, :])
K24 = copy(-b2*K4[2:end, :])
B14 = copy(b1*K4[2:end, :])
B24 = copy(b2*K4[2:end, :])

K = [K11+K12+K13 K14;K23 K21+K22+K24]
Ks = [K[:,2:N] K[:, N+2:end]]

n = int(sqrt(length(Ks)))
b0 = zeros(3,3)
b1 = diagm([1,1,1.])

newK = zeros(3*n, 3*n)
counter = 1

for j=1:3:3*n
  for i=1:3:3*n
    if Ks[counter] != 0.
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
imwrite(newK, "p4_mass.pdf")
