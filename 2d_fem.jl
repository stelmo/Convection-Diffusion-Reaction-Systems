# Two dimensional problem

xN = 5 # Number of elements - 1
xs = linspace(0., 1., xN)
dx = xs[2] - xs[1]

yN = 4 # Number of elements - 1
ys = linspace(0., 1., yN)
dy = ys[2] - ys[1]

us = zeros(yN, xN)

function getRecs(node, yN, xN)
  # Returns all the admissable rectangles adjacent to the node

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



function get_K4(K4, dy, dx)

  # M matrix - its symmetric hence only these entries
  m11 = (4./36.)*dx*dy
  m12 = (2./36.)*dx*dy
  m13 = (1./36.)*dx*dy
  m14 = (2./36.)*dx*dy
  m23 = (2./36.)*dx*dy
  m24 = (1./36.)*dx*dy
  m34 = (2./36.)*dx*dy

  # Construct K4
  (yN, xN) = size(K4)
  for i=1:xN
    for j=1:yN
      # K4[j,i] = (j, i, yN, xN, dy, dx)
    end
  end

end
