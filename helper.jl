# Helper functions - basically the globally common functions/structs
module Helper

type Model

   xs :: Array{Float64, 1}
   ts :: Array{Float64, 1}
   h :: Float64
   k :: Float64
   us :: Array{Float64, 2}
   a :: Float64
   b :: Float64
   c :: Float64
   Nx :: Int64
   Nt :: Int64
   name:: String

   function Model(x0, xL, t0, tL, Nx, Nt, a, b, c, name)
      xs = linspace(x0, xL, Nx)
      ts = linspace(t0, tL, Nt)
      h = xs[2] - xs[1] # space discretisation step
      k = ts[2] - ts[1] # time discretisation step
      us = zeros(Float64, (Nx, Nt)) # solution matrix
      new(xs, ts, h, k, us, a, b, c, Nx, Nt, name)
   end
end

function enforceBCs!(model::Model)

   # u(x, 0) condition
   for i=1:model.Nx
      if 2.0 <= model.xs[i] <= 2.0 + pi
         model.us[i, 1] = sin(model.xs[i]-2.0) # arbitrary function
      else
         model.us[i, 1] = 0.0
      end
   end

end

function vnBC0(time::Float64)
   # Von Neumann boundary condition at x=0
   return 0.0
end

function vnBCL(time::Float64)
   # Von Neumann boundary condition at x=L
   return 0.0
end

function forcingFunction(space::Float64, time::Float64)
   return 0.0
end

function combine(models::Model...)
   space = models[1].xs
   values = models[1].us
   solvers = fill(models[1].name, models[1].Nx)
   for model in models[2:end]
      space = [space, model.xs]
      values = [values, model.us]
      solvers = [solvers, fill(model.name, model.Nx)]
   end
   return space, values, solvers
end

end # module
