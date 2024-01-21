function SandT(x1, y1, x2, y2, x3, y3)
    global LOCALEDGENODES

    area = .5abs(det([1. x1 y1
                      1. x2 y2
                      1. x3 y3]))
    temp = inv([x1 x2 x3
                y1 y2 y3
                1. 1. 1.])
    b = temp[:,1]
    c = temp[:,2]
    a = temp[:,3]
    ∇λ = [[b[1], c[1]],
          [b[2], c[2]],
          [b[3], c[3]]]
    φ = zeros(3, 3) # φij = ∇λi ⋅ ∇λj = (bi⋅ci + bj⋅cj)
    v = zeros(3, 3) # vij = ∇λi × ∇λj = (bi⋅cj - bj⋅ci)ẑ
    for ii = 1:3
       for jj = 1:3
           φ[ii, jj] = b[ii] * b[jj] + c[ii] * c[jj]
           v[ii, jj] = b[ii] * c[jj] - b[jj] * c[ii]
       end 
    end

    M = [2. 1. 1.
         1. 2. 1.
         1. 1. 2.] / 12.0

    # Compute S and T
    S = zeros(3, 3)
    T = zeros(3, 3)
    for ii = 1:3
       for jj = 1:3
         i1 = LOCALEDGENODES[ii,1]
         i2 = LOCALEDGENODES[ii,2]
         j1 = LOCALEDGENODES[jj,1]
         j2 = LOCALEDGENODES[jj,2]
         S[ii,jj] = 4area*v[i1, i2] * v[j1, j2]
         T[ii,jj] = area*(φ[i2, j2] * M[i1, j1] +
                        - φ[i2, j1] * M[i1, j2] +
                        - φ[i1, j2] * M[i2, j1] +
                        + φ[i1, j1] * M[i2, j2])
       end 
    end
    
    return S, T
end

function simplex2D(elem_num, xp, yp)
    global ELEMENTS
    global NODE_COORD 

    trinodes = ELEMENTS[elem_num, :]; 
    x1 = NODE_COORD[trinodes[1],1];
    y1 = NODE_COORD[trinodes[1],2];
    x2 = NODE_COORD[trinodes[2],1];
    y2 = NODE_COORD[trinodes[2],2];
    x3 = NODE_COORD[trinodes[3],1];
    y3 = NODE_COORD[trinodes[3],2];

    σ0 = det([1. x1 y1; 1. x2 y2; 1. x3 y3]);
    σ1 = det([1. xp yp; 1. x2 y2; 1. x3 y3]);
    σ2 = det([1. x1 y1; 1. xp yp; 1. x3 y3]);
    σ3 = det([1. x1 y1; 1. x2 y2; 1. xp yp]);
    λ  = [σ1/σ0, σ2/σ0, σ3/σ0]
    return λ
end

DOF_NONE = 0
DOF_PEC  = 1

function dof_type(a, b)
    global NUM_EDGES
    global EDGES
    global NODE_COORD
    dof_flag = zeros(Int64, NUM_EDGES)
    
    for i_edge = 1:NUM_EDGES
       if count(==(i_edge), ELEMENT_EDGES) == 1
         dof_flag[i_edge] = DOF_PEC
       end
    end
    return dof_flag
end

function dof_renumber!(dof, Γ)
    global NUM_EDGES
    global NUM_DOFS
    last = 0
    for i_edge = 1:NUM_EDGES
       if Γ[i_edge] == 0
         dof[i_edge] = last + 1
         last += 1
       else
         dof[i_edge] = 0
       end
    end
    NUM_DOFS = last
    return dof
end

function assemble()
  # ASSUMPTION: Waveguide is homogenous
  # Assemble stiffness and mass matrices
  S = zeros(NUM_DOFS, NUM_DOFS) # 1/μr ∫( ∇Ni × ∇Nj )dΩ
  T = zeros(NUM_DOFS, NUM_DOFS) # μ0 ε ∫( ∇Ni ⋅ ∇Nj )dΩ
  
  # NODE_COORD - współrzędne wszystkich punktów z których zbudowane są trójkąty (a być może po prosty figura?)
  # ELEMENTS - numery punktów z których zbudowane są trójkąty
  # NUM_ELEMS - liczba trójkątów
  
  for ielem = 1:NUM_ELEMS # Assemble by elements
    trinodes = ELEMENTS[ielem, :]
    Se, Te = SandT(NODE_COORD[trinodes[1],1], NODE_COORD[trinodes[1],2],
                   NODE_COORD[trinodes[2],1], NODE_COORD[trinodes[2],2],
                   NODE_COORD[trinodes[3],1], NODE_COORD[trinodes[3],2])
    
    for jedge = 1:3
      jj = dof[ELEMENT_EDGES[ielem, jedge]]
      if jj == 0 continue end
      for kedge = 1:3
        kk = dof[ELEMENT_EDGES[ielem, kedge]]
        if kk == 0 continue end
        S[jj, kk] = S[jj, kk] + (1/μr) * Se[jedge, kedge]
        T[jj, kk] = T[jj, kk] + (μ0*ε) * Te[jedge, kedge]
      end
    end
  end
  return S, T
end
