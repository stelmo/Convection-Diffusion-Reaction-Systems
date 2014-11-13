# Two dimensional problem

xN = 3 # Number of elements - 1
xs = linspace(0., 1., xN)
dx = xs[2] - xs[1]

yN = 3 # Number of elements - 1
ys = linspace(0., 1., yN)
dy = ys[2] - ys[1]

us = zeros(yN, xN)

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

# Create the general matrix
K4 = zeros(yN*xN, xN*yN) # total # nodes * total # nodes
for j = [1:xN*yN]
  for i = [1:xN*yN]
    K4[i,j] = getK4(i, j, yN, xN, dy, dx)
  end
end

# Strip out unnecessary rows
K4 = K4[yN+1:end, :]
K4_const = K4[:, 1:yN]
K4 = K4[:, yN+1:end]
