# One dimensional finite element method for the convection diffusion reaction scheme.

# Model form:
# u_t = b*u_xx - a*u_x - c*u + f(u, t)
# u_x(0,t) = g(t)
# u_x(L, t) = h(t)
# u(x, 0) = i(t)
# with x0 < x < xL
# and t0 < t < tL
# The interpolant is used to approximate f(u, t)

module Fem1Du

using Helper: Model, forcingFunction, vnBC0, vnBCL

function solveImplicit!(model::Model)
# FEM approximation will be of the form:
# Mu' = (-b*K1-c*M-a*K2)u + K3 + b*u_x(L,t)[0,0,...,1] - b*u_x(0,t)*[1,0,...0]
# Crank-Nicolson FDM is used to solve the FEM approximation

   K1 = createK1(model)
   K2 = createK2(model)
   M = createM(model)

   gamma = (-model.b*K1 - model.c*M - model.a*K2)

   # time evaluation loop
   for i=2:model.Nt

      Kspast = createK3(model, model.ts[i-1]) + createBoundaryVectors(model, model.ts[i-1])
      Ksnow = createK3(model, model.ts[i]) + createBoundaryVectors(model, model.ts[i])
      A = M - (model.k/2.0)*gamma
      B = (M + (model.k/2.0)*gamma)*model.us[:,i-1] + (model.k/2.0)*(Kspast + Ksnow)
      model.us[:,i] = \(A, B)

   end

end

function solveExplicit!(model::Model)
# FEM approximation will be of the form:
# Mu' = (-b*K1-c*M-a*K2)u + K3 + b*u_x(L,t)[0,0,...,1] - b*u_x(0,t)*[1,0,...0]
# Mass lumping and forward difference FDM is used to solve the FEM approximation

   K1 = createK1(model)
   K2 = createK2(model)
   M = createM(model)
   Dinv = inv(createD(model)) # perform this only once

   gamma = Dinv*(-model.b*K1 - model.c*M - model.a*K2) # perform only once

   #time evaluation loop
   for i=2:model.Nt

      Kspast = createK3(model, model.ts[i-1]) + createBoundaryVectors(model, model.ts[i-1])
      model.us[:,i] = model.us[:,i-1] + model.k*(gamma*model.us[:,i-1] + Dinv*Kspast)

   end

end

function createD(model::Model)

   d1 = model.h*ones(Float64, model.Nx)
   d1[1] = model.h/2.0
   d1[model.Nx] = model.h/2.0
   D = Diagonal(d1)
   return D

end

function createM(model::Model)

   m1 = (2.0/3.0)*model.h*ones(Float64, model.Nx) # diagonal
   m1[1] = model.h/3.0 # fix boundary case
   m1[model.Nx] = model.h/3.0 # fix boundary case
   m2 = (model.h/6.0)*ones(Float64, model.Nx-1) # off diagonal
   M = Tridiagonal(m2,m1,m2)

   return M
end

function createK1(model::Model)

   k1 = (2.0/model.h)*ones(Float64, model.Nx) # diagonal - extend trivially to incorporate pure terms
   k1[1,1] = 1.0/model.h # fix boundary case
   k1[model.Nx] = 1.0/model.h # fix boundary case
   k2 = (-1.0/model.h)*ones(Float64, model.Nx-1)
   K = Tridiagonal(k2, k1, k2)

   return K
end

function createK2(model::Model)

   n2 = (1.0/2.0)*ones(Float64, model.Nx-1)
   n1 = zeros(Float64, model.Nx)
   n1[1] = -1.0/2.0 # fix boundary case
   n1[model.Nx] = 1.0/2.0 # fix boundary case
   N = Tridiagonal(-n2, n1, n2)

   return N
end

function createK3(model::Model, timeState::Float64)

   fx = zeros(Float64, model.Nx)
   for k=1:model.Nx
      fx[k] = forcingFunction(model.xs[k], timeState)
   end
   M = createM(model)
   K3 = M*fx

   return K3

end

function createBoundaryVectors(model::Model, timeState::Float64)

   bvecL = zeros(Float64, model.Nx)
   bvecL[model.Nx] = 1.0
   bvec0 = zeros(Float64, model.Nx)
   bvec0[1] = 1.0

   bvecs = model.b * vnBCL(timeState) * bvecL
   bvecs = bvecs - model.b * vnBC0(timeState) * bvec0

   return bvecs
end

end # module





