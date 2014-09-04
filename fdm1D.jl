# One dimensional finite difference method for the convection diffusion reaction scheme.

# Model form:
# u_t = b*u_xx - a*u_x - c*u + f(u, t)
# u_x(0,t) = g(t)
# u_x(L, t) = h(t)
# u(x, 0) = i(t)
# with x0 < x < xL
# and t0 < t < tL

module Fdm1Du

using Helper: Model, forcingFunction, vnBC0, vnBCL

function solveFTCS!(model::Model)
# Forward time central space

gamma = 1.0/(1.0+model.c*0.5*model.k)
c1 = model.b*model.k/(model.h^2)
c2 = model.a*model.k*0.5/model.h

for tIndex = 1:model.Nt-1 # solve for next step

   # x=0 boundary case
   vsub1 = model.us[2, tIndex] - 2.0*model.h*vnBC0(model.ts[tIndex])
   k1 = gamma*(c1+c2)*vsub1
   k2 = gamma*(1.0 - 2.0*c1 - model.c*0.5*model.k)*model.us[1, tIndex]
   k3 = gamma*(c1 - c2)*model.us[2, tIndex]
   k4 = model.k*0.5*gamma*(forcingFunction(model.xs[1], model.ts[tIndex+1]) + forcingFunction(model.xs[1], model.ts[tIndex]))
   model.us[1, tIndex+1] = k1 + k2 + k3 + k4

   # 0<x<L loop
   for xIndex=2:model.Nx-1
      k1 = gamma*(c1+c2)*model.us[xIndex-1, tIndex]
      k2 = gamma*(1.0 - 2.0*c1 - model.c*0.5*model.k)*model.us[xIndex, tIndex]
      k3 = gamma*(c1 - c2)*model.us[xIndex+1, tIndex]
      k4 = model.k*0.5*gamma*(forcingFunction(model.xs[xIndex], model.ts[tIndex+1]) + forcingFunction(model.xs[xIndex], model.ts[tIndex]))
      model.us[xIndex, tIndex+1] = k1 + k2 + k3 + k4
   end

   # x=L boundary case
   vplus1 = model.us[model.Nx-1, tIndex] + 2.0*model.h*vnBCL(model.ts[tIndex])
   k1 = gamma*(c1+c2)*model.us[model.Nx-1, tIndex]
   k2 = gamma*(1.0 - 2.0*c1 - model.c*0.5*model.k)*model.us[model.Nx, tIndex]
   k3 = gamma*(c1 - c2)*vplus1
   k4 = model.k*0.5*gamma*(forcingFunction(model.xs[model.Nx], model.ts[tIndex+1]) + forcingFunction(model.xs[model.Nx], model.ts[tIndex]))
   model.us[model.Nx, tIndex+1] = k1 + k2 + k3 + k4

end


end

function solveFTCSup!(model::Model)
# Forward time central space with upwind differencing for the convection term

gamma = 1.0/(1.0+model.c*0.5*model.k)
c1 = model.b*model.k/(model.h^2)
c2 = model.a*model.k/model.h

for tIndex = 1:model.Nt-1 # solve for next step

   # x=0 boundary case
   vsub1 = model.us[2, tIndex] - 2.0*model.h*vnBC0(model.ts[tIndex])
   k1 = gamma*(c1+c2)*vsub1
   k2 = gamma*(1.0 - 2.0*c1 - model.c*0.5*model.k - c2)*model.us[1, tIndex]
   k3 = gamma*c1*model.us[2, tIndex]
   k4 = model.k*0.5*gamma*(forcingFunction(model.xs[1], model.ts[tIndex+1]) + forcingFunction(model.xs[1], model.ts[tIndex]))
   model.us[1, tIndex+1] = k1 + k2 + k3 + k4

   # 0<x<L loop
   for xIndex=2:model.Nx-1
      k1 = gamma*(c1+c2)*model.us[xIndex-1, tIndex]
      k2 = gamma*(1.0 - 2.0*c1 - model.c*0.5*model.k - c2)*model.us[xIndex, tIndex]
      k3 = gamma*c1*model.us[xIndex+1, tIndex]
      k4 = model.k*0.5*gamma*(forcingFunction(model.xs[xIndex], model.ts[tIndex+1]) + forcingFunction(model.xs[xIndex], model.ts[tIndex]))
      model.us[xIndex, tIndex+1] = k1 + k2 + k3 + k4
   end

   # x=L boundary case
   vplus1 = model.us[model.Nx-1, tIndex] + 2.0*model.h*vnBCL(model.ts[tIndex])
   k1 = gamma*(c1+c2)*model.us[model.Nx-1, tIndex]
   k2 = gamma*(1.0 - 2.0*c1 - model.c*0.5*model.k - c2)*model.us[model.Nx, tIndex]
   k3 = gamma*c1*vplus1
   k4 = model.k*0.5*gamma*(forcingFunction(model.xs[model.Nx], model.ts[tIndex+1]) + forcingFunction(model.xs[model.Nx], model.ts[tIndex]))
   model.us[model.Nx, tIndex+1] = k1 + k2 + k3 + k4

end




end

end # module





