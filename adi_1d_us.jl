using PyPlot

# Inputs
a = 7.5 # Peclet number
b = 5.129 #  L/U

N = 100 # number of nodes - need many of these
xs = linspace(0., 1., N)
h = xs[2] - xs[1]

Nt = 100
ts = linspace(0., 1., Nt) # this is dimensionless time... It also takes long...
dt = ts[2] - ts[1]

us = zeros(N, Nt)

# Rate functions
# kT(x) = 0.25
# kT(x) = -1.75*x + 2.
# kT(x) = -3.*x+3.2
# kT(x) = (3.)*x.^2 -6.*x + 3.2
kT(x) = 1633253.0*exp(-4881.*(1./(300.+100.*(1.-x)))) # should be 2 910 122.0
kTd(x) = (kT(x+0.00001)-kT(x))/0.00001 # derivative approximation

# Initial Condition
finit(x) = exp(-x*100) # 4.0.*x.^2 -4.0.*x +1 # ## Initial condition
us[:,1] = finit(xs) # Initial condition

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


# Conform the constant matrices
K1 = K1[:, 2:end]
K2 = K2[:, 2:end]

# Start time loop
for tloop=2:Nt

  newtonFlag = true
  counter = 0
  alpha = 0.8 # slow down NR - dampening
  maxIter = 1000

  # Guess a solution
  uguess = us[:, tloop-1] # guess the previous solution

    while newtonFlag

    # Assemble the non-linear non-constant matrices
    dm3 = fill(2.*h/3., N)
    dm3[1] = h/3.
    dm3[end] = h/3.

    dl3 = fill(h/6., N-1)
    du3 = fill(h/6., N-1)

    K5 = diagm(dl3,-1) + diagm(dm3, 0) + diagm(du3, 1)
    K5 = K5[2:end, :]
    K6 = copy(K5)
    M = copy(K5)
    M = a*M
    K6 = a*b*K6
    K5 = -a*b*K5

    for row=1:N-1
      K6[row, :] = kT(uguess[row])*K6[row, :]
      K5[row, :] = kT(uguess[row])*K5[row, :]
    end

    dl4 = zeros(N-1)
    for i=1:N-1
      dl4[i] = kTd(uguess[i+1])*(h/12.)*(uguess[i]+uguess[i+1])
    end

    du4 = zeros(N-1)
    for i=1:N-1
      du4[i] = kTd(uguess[i])*(h/12.)*(uguess[i]+uguess[i+1])
    end

    dm4 = zeros(N)
    dm4[1] = -9. # it gets cut away so w/e
    for i=2:N-1
      dm4[i] = kTd(uguess[i])*(uguess[i-1]*h/12. + uguess[i]*h/2. + uguess[i+1]*h/12.)
    end
    dm4[end] = kTd(uguess[end])*(uguess[end-1]*h/12. + uguess[end]*h/4.)

    K7 = diagm(dl4,-1) + diagm(dm4, 0) + diagm(du4, 1)
    K7 = K7[2:end, :]
    K7 = -a*b*K7

    B1 = (K3 + K4 + K6)*uguess + (1./dt)*M*(uguess - us[:,tloop-1])

    K5 = K5[:, 2:end]
    K7 = K7[:, 2:end]
    # du0 = 0 because we know u0 = 1 = u

    KK = 0.5*(K1 + K2 + K5 + K7) - M[:,2:end]*(1./dt)
    du = \(KK, B1)
    uguess = uguess + alpha*[0., du]

    (abs(maximum(du)) < 1.E-05 || counter > maxIter) ? newtonFlag = false : newtonFlag = true
    counter += 1
  end

  us[:,tloop] = uguess

end

mesh(xs, ts, us')
xlabel("Space")
ylabel("Time")

# plot(xs, us[:, end])

plt.show()
