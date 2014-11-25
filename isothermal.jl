# Comparison of the isothermal systems

using PyPlot

# Inputs
k = 0.25
L = 6.36
R = 0.05
De = 7.6E-05
U = 1.24


N = 50
Nt = 50
xs = linspace(0., 1., N+1)
ts = linspace(0., 5., Nt)

function fem_1d(k, L, R, De, U, N, Nt)

  Da = De + (U^2 * R^2 )/(48*De)
  a = U*L/Da
  b = L*k/U

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
  ts = linspace(0., 5., Nt) # this is dimensionless time... It also takes long...
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


  return us
end

function fem_2d(kreact, L, R, De, U, N, Nt)

  xN = N+1 # because for some reason I have + 1 in the 1D case
  xs = linspace(0., 1., xN)
  dx = xs[2] - xs[1]

  yN = N # Number of elements - 1
  ys = linspace(0., 1., yN)
  dy = ys[2] - ys[1]

  tN = Nt
  ts = linspace(0., 5., tN) # it takes long to reach SS. Try 5.0 for the end time.
  dt = ts[2] - ts[1]

  reaction_constant = kreact # should be 0.25
  Dab = De
  U0 = U

  function getRecs(node, yN, xN)
    # Returns all the admissible rectangles adjacent to the node

    total = yN*xN # total number of nodes
    recs = Array{Int64,1}[]

    tR = false
    bR = false
    lC = false
    rC = false

    # Check if in the top row
    (node in [1:yN:total]) ? tR = true : tR = false
    # Check if in the bottom row
    (node in [yN:yN:total]) ? bR = true : bR = false
    # Check if in the left most column
    (node in [1:1:yN]) ? lC = true : lC = false
    # Check if in the right most column
    (node in [total-yN+1:1:total]) ? rC = true : rC = false

    # Check all possible types
    if tR # top row
      if lC # Top left
        push!(recs, [node, node+1, node+yN, node+yN+1])
      elseif rC # Top right
        push!(recs, [node-yN, node-yN+1, node, node+1])
      else
          push!(recs, [node-yN, node-yN+1, node, node+1])
          push!(recs, [node, node+1, node+yN, node+yN+1])
      end
    elseif bR # bottom row
      if lC # Bottom left
        push!(recs, [node-1, node, node+yN-1, node+yN])
      elseif rC # Bottom right
        push!(recs, [node-yN-1, node-yN, node-1, node])
      else
        push!(recs, [node-yN-1, node-yN, node-1, node])
        push!(recs, [node-1, node, node+yN-1, node+yN])
      end
    elseif lC && !tR && !bR
      push!(recs, [node-1, node, node+yN-1, node+yN])
      push!(recs, [node, node+1, node+yN, node+yN+1])
    elseif rC && !tR && !bR
      push!(recs, [node-yN-1, node-yN, node-1, node])
      push!(recs, [node-yN, node-yN+1, node, node+1])
    else # center
      push!(recs, [node, node+1, node+yN, node+yN+1])
      push!(recs, [node-yN-1, node-yN, node-1, node])
      push!(recs, [node-1, node, node+yN-1, node+yN])
      push!(recs, [node-yN, node-yN+1, node, node+1])
    end

    # Just in case -  don't uncomment!
    # push!(recs, [node-yN-1, node-yN, node-1, node]) # NW
    # push!(recs, [node-1, node, node+yN-1, node+yN]) # NE
    # push!(recs, [node-yN, node-yN+1, node, node+1]) # SW
    # push!(recs, [node, node+1, node+yN, node+yN+1]) # SE

    return recs
  end

  function getIndex(ind)

    if ind == 4
      return 2
    elseif ind == 3
      return 3
    elseif ind == 2
      return 1
    else
      return 4
    end
  end

  function getTheta(rec, yN, dy)

    topvertex = rec[1]

    flag = false
    topvertex <= yN ? flag = false : flag = true

    while flag
      topvertex = topvertex - yN
      (topvertex <= yN) ? flag = false : flag = true
    end

    m = (1.-dy)/(1.-yN)
    c = dy/2. - yN*m
    theta = m*topvertex + c

    return theta
  end

  function getK4(row_node, col_node, yN, xN, dy, dx)
    # Get each entry for the K4 matrix

    # M matrix - its symmetric hence only these entries
    m11 = (4./36.)*dx*dy
    m12 = (2./36.)*dx*dy
    m13 = (1./36.)*dx*dy
    m14 = (2./36.)*dx*dy
    m23 = (2./36.)*dx*dy
    m24 = (1./36.)*dx*dy
    m34 = (2./36.)*dx*dy

    M = zeros(4,4)
    M[1,1] = m11
    M[2,2] = m11
    M[3,3] = m11
    M[4,4] = m11
    M[1,2] = m12
    M[2,1] = m12
    M[1,3] = m13
    M[3,1] = m13
    M[1,4] = m14
    M[4,1] = m14
    M[2,3] = m23
    M[3,2] = m23
    M[2,4] = m24
    M[4,2] = m24
    M[3,4] = m34
    M[4,3] = m34

    rn_recs = getRecs(row_node, yN, xN)
    cn_recs = getRecs(col_node, yN, xN)

    entry = 0.0

    for rn_rec in rn_recs
      for cn_rec in cn_recs
        if rn_rec == cn_rec
          ind_r = getIndex(findfirst(rn_rec, row_node))
          ind_c = getIndex(findfirst(cn_rec, col_node))
          # println(ind_r, ind_c) # for fault checking
          entry = entry + M[ind_r, ind_c]
        end
      end
    end
    # println("****") # for fault checking
    return entry
  end

  function getK1(row_node, col_node, yN, xN, dy, dx)
    # Get each entry for the K1 matrix

    # K1 matrix - its symmetric hence only these entries
    d = 2.*dx/dy
    s1 = dx/dy
    s2 = -2.*dx/dy

    m11 = 2.*d
    m12 = 2.*s1
    m13 = -d
    m14 = 2.*s2
    m23 = 2.*s2
    m24 = -d
    m34 = 2.*s1

    M = zeros(4,4) # it's symmetric
    M[1,1] = m11
    M[2,2] = m11
    M[3,3] = m11
    M[4,4] = m11
    M[1,2] = m12
    M[2,1] = m12
    M[1,3] = m13
    M[3,1] = m13
    M[1,4] = m14
    M[4,1] = m14
    M[2,3] = m23
    M[3,2] = m23
    M[2,4] = m24
    M[4,2] = m24
    M[3,4] = m34
    M[4,3] = m34

    M = 1./12.*M

    rn_recs = getRecs(row_node, yN, xN)
    cn_recs = getRecs(col_node, yN, xN)

    entry = 0.0

    for rn_rec in rn_recs
      for cn_rec in cn_recs
        if rn_rec == cn_rec
          ind_r = getIndex(findfirst(rn_rec, row_node))
          ind_c = getIndex(findfirst(cn_rec, col_node))
          # println(ind_r, ind_c) # for fault checking
          entry = entry + M[ind_r, ind_c]
        end
      end
    end
    # println("****") # for fault checking
    return entry
  end

  function getK2(row_node, col_node, yN, xN, dy, dx)
    # Get each entry for the K2 matrix

    # K1 matrix - its symmetric hence only these entries
    d = 2.*dy/dx
    s1 = -2.*dy/dx
    s2 = dy/dx

    m11 = 2.*d
    m12 = 2.*s1
    m13 = -d
    m14 = 2.*s2
    m23 = 2.*s2
    m24 = -d
    m34 = 2.*s1

    M = zeros(4,4) # it's symmetric
    M[1,1] = m11
    M[2,2] = m11
    M[3,3] = m11
    M[4,4] = m11
    M[1,2] = m12
    M[2,1] = m12
    M[1,3] = m13
    M[3,1] = m13
    M[1,4] = m14
    M[4,1] = m14
    M[2,3] = m23
    M[3,2] = m23
    M[2,4] = m24
    M[4,2] = m24
    M[3,4] = m34
    M[4,3] = m34

    M = 1./12.*M

    rn_recs = getRecs(row_node, yN, xN)
    cn_recs = getRecs(col_node, yN, xN)

    entry = 0.0

    for rn_rec in rn_recs
      for cn_rec in cn_recs
        if rn_rec == cn_rec
          ind_r = getIndex(findfirst(rn_rec, row_node))
          ind_c = getIndex(findfirst(cn_rec, col_node))
          # println(ind_r, ind_c) # for fault checking
          entry = entry + M[ind_r, ind_c]
        end
      end
    end
    # println("****") # for fault checking
    return entry
  end

  function getK3(row_node, col_node, yN, xN, dy, dx)
    # Get each entry for the K3 matrix

    M = zeros(4,4) # it's not symmetric
    M[1,1] = -2.
    M[1,2] = 2.
    M[1,3] = 1.
    M[1,4] = -1.

    M[2,1] = -2.
    M[2,2] = 2.
    M[2,3] = 1.
    M[2,4] = -1.

    M[3,1] = -1.
    M[3,2] = 1.
    M[3,3] = 2.
    M[3,4] = -2.

    M[4,1] = -1.
    M[4,2] = 1.
    M[4,3] = 2.
    M[4,4] = -2.

    M = (dy/12)*M

    rn_recs = getRecs(row_node, yN, xN)
    cn_recs = getRecs(col_node, yN, xN)

    entry = 0.0

    for rn_rec in rn_recs
      for cn_rec in cn_recs
        if rn_rec == cn_rec
          ind_r = getIndex(findfirst(rn_rec, row_node))
          ind_c = getIndex(findfirst(cn_rec, col_node))
          # println(ind_r, ind_c) # for fault checking
          entry = entry + (1.-getTheta(rn_rec,yN, dy)^2)*M[ind_r, ind_c]
        end
      end
    end
    # println("****") # for fault checking
    return entry
  end

  # Create the general K1 matrix
  K1 = zeros(yN*xN, xN*yN) # total # nodes * total # nodes
  for j = [1:xN*yN]
    for i = [1:xN*yN]
      K1[i,j] = getK1(i, j, yN, xN, dy, dx)
    end
  end

  # Create the general K2 matrix
  K2 = zeros(yN*xN, xN*yN) # total # nodes * total # nodes
  for j = [1:xN*yN]
    for i = [1:xN*yN]
      K2[i,j] = getK2(i, j, yN, xN, dy, dx)
    end
  end

  # Create the general K4 matrix
  K3 = zeros(yN*xN, xN*yN) # total # nodes * total # nodes
  for j = [1:xN*yN]
    for i = [1:xN*yN]
      K3[i,j] = getK3(i, j, yN, xN, dy, dx)
    end
  end

  # Create the general K4 matrix
  K4 = zeros(yN*xN, xN*yN) # total # nodes * total # nodes
  for j = [1:xN*yN]
    for i = [1:xN*yN]
      K4[i,j] = getK4(i, j, yN, xN, dy, dx)
    end
  end

  K5 = zeros(yN*xN, xN*yN) # total # nodes * total # nodes
  for j = [1:xN*yN]
    for i = [1:xN*yN]
      K5[i,j] = getK4(i, j, yN, xN, dy, dx)
    end
  end

  # Strip out unnecessary rows for K1
  K1 = (-Dab/(R^2))*K1[yN+1:end, :]
  K1_const = K1[:, 1:yN]
  K1 = K1[:, yN+1:end]

  # Strip out unnecessary rows for K2
  K2 = (-Dab/(L^2))*K2[yN+1:end, :]
  K2_const = K2[:, 1:yN]
  K2 = K2[:, yN+1:end]

  # Strip out unnecessary rows for K4
  K3 = -(2.*U0/L)*K3[yN+1:end, :]
  K3_const = K3[:, 1:yN]
  K3 = K3[:, yN+1:end]

  # Strip out unnecessary rows for K4
  K4 = -reaction_constant*K4[yN+1:end, :]
  K4_const = K4[:, 1:yN]
  K4 = K4[:, yN+1:end]

  # Strip out unnecessary rows for K4
  K5 = (U0/((L^2)*R))*K5[yN+1:end, :]
  K5 = K5[:, yN+1:end]

  # Calculate the magic!
  KK = K1 + K2 + K3 + K4
  BB = zeros(xN*yN-yN, 1)
  u1 = zeros(yN)

  # Left boundary condition
  for i=1:yN
    BB[:] = BB[:] +  (K1_const[:,i] + K2_const[:,i] + K3_const[:,i] + K4_const[:,i])
    u1[i] = 1.0
  end

  # Initial condition
  us = zeros(yN, xN)
  us[:,1] = u1[:]
  for i=2:xN
    us[:,i] = us[:, i-1]*exp(-10.*xs[i])
  end
  ## Weird conditions for fun :)
  ## Weird one
  # ymid = int(yN/2)
  # xmid = int(xN/2)
  # us[ymid-2:ymid+2, xmid-2:xmid+2] = 1.
  ## Weird two
  # us[:, 2:end] = rand(yN,xN-1)

  # Normal code follows

  u0 = reshape(us[:,2:end], ((xN-1)*yN, 1))
  u_profile = zeros((xN-1)*yN, tN)
  u_profile[:,1] = u0

  for t=2:tN
    A = K5 - (dt/2.)*KK
    B = (K5 + (dt/2)*KK)*u_profile[:,t-1] + dt*BB
    u_profile[:, t] = \(A, B)
  end

  uave = zeros(xN, tN)
  for time=1:tN
    u = reshape(u_profile[:, time], (yN, xN-1))
    u = [u1 u]
    # Get average concentration
    aveconv = zeros(xN)
    for j=1:xN
      for i=1:yN
        aveconv[j] += 2.*dy*ys[i]*u[i,j]
      end
    end
    uave[:, time] = aveconv[:]
  end


  return uave
end



us_1d = fem_1d(k, L, R, De, U, N, Nt)
us_2d = fem_2d(k, L, R, De, U, 30, 30)




# # Start
# figure(1)
# u = reshape(u_profile[:, 1], (yN, xN-1))
# u = [u1 u]
# contourf(xs, reverse(ys), u, 20)
# colorbar()
# xlabel("Axial Distance")
# ylabel("Radial Distance")
#
# # Middle
# figure(2)
# mid = int((tN/2))
# u = reshape(u_profile[:, mid], (yN, xN-1))
# u = [u1 u]
# contourf(xs, reverse(ys), u, 20)
# colorbar()
# xlabel("Axial Distance")
# ylabel("Radial Distance")
#
# # End
# figure(3)
# u = reshape(u_profile[:, end], (yN, xN-1))
# u = [u1 u]
# contourf(xs, reverse(ys), u, 20)
# colorbar()
# xlabel("Axial Distance")
# ylabel("Radial Distance")
#
# plt.show()
#
# # Get average concentration
# aveconv = 0.0
# for i=1:yN
#   aveconv += 2.*dy*ys[i]*u[i,end]
# end
# println(1.-aveconv)



plt.show()
