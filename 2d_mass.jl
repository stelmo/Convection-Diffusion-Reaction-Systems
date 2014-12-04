using Images, ImageView

xN = 10 # Number of elements - 1
xs = linspace(0., 1., xN)
dx = xs[2] - xs[1]

yN = 10 # Number of elements - 1
ys = linspace(0., 1., yN)
dy = ys[2] - ys[1]

reaction_constant = 0.25 #0.25
Dab = 7.6E-05 # should be 7.6E-05
L = 6.36
R = 0.05 # should be 0.05
U0 = 1.24

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

# Strip out unnecessary rows for K1
K1 = (-Dab/(R^2))*K1[yN+1:end, :]
K1_const = K1[:, 1:yN]
K1 = K1[:, yN+1:end]

# Strip out unnecessary rows for K2
K2 = (-Dab/(L^2))*K2[yN+1:end, :]
K2_const = K2[:, 1:yN]
K2 = K2[:, yN+1:end]

# Strip out unnecessary rows for K4
K3 = -2.*(U0/L)*K3[yN+1:end, :]
K3_const = K3[:, 1:yN]
K3 = K3[:, yN+1:end]

# Strip out unnecessary rows for K4
K4 = -reaction_constant*K4[yN+1:end, :]
K4_const = K4[:, 1:yN]
K4 = K4[:, yN+1:end]


# Calculate the magic!
KK = K1 + K2 + K3 + K4

for k=1:length(KK)
  if KK[k] != 0.0
    KK[k] = 0.0
  else
    KK[k] = 1.0
  end
end

imwrite(KK, "p2_mass.pdf")
