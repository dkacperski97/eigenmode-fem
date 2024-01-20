using Printf
using Statistics
using LinearAlgebra
using VoronoiDelaunay
import VoronoiDelaunay: getx, gety

mutable struct IndexedPoint2D <: AbstractPoint2D
    _x::Float64
    _y::Float64
    _idx::Int64
    IndexedPoint2D(x, y, idx) = new(x, y, idx)
    IndexedPoint2D(x, y) = new(x, y, 0)
end

getx(p::IndexedPoint2D) = p._x
gety(p::IndexedPoint2D) = p._y
getidx(p::IndexedPoint2D) = p._idx

function l_initpoints( h0, bbox )
    # Create an initial grid point/vertex distribution.

    n_sdim = size(bbox,2)
    pinit  = Array{Any}(undef, n_sdim)
    for i=1:n_sdim
        if( i==2 )
            pinit[i] = bbox[1,i]:h0*sqrt(3)/2:bbox[2,i]
        else
            pinit[i] = bbox[1,i]:h0:bbox[2,i]
        end
    end
    pp = Array{Any}(undef, n_sdim)
    n  = length(pinit[1])
    m  = length(pinit[2])
    xh = reshape(pinit[1], n, 1)
    yh = reshape(pinit[2], 1, m)

    pp[1] = repeat(xh, 1, m)
    pp[2] = repeat(yh, n, 1)

    pp[1][:,2:2:end] = pp[1][:,2:2:end] .+ h0/2

    p = Array{Float64}(undef, prod(size(pp[1])),n_sdim)

    for i=1:n_sdim
        p[:,i] = pp[i][:]
    end

    return p
end

#------------------------------------------------------------------------------#
"return true for Points inside the Polygon"
function inpolygonv( x, y, xv, yv )
    # Octave inpolygon implementation.

    npol = length(xv)
    in = fill(false,length(x))
    on = fill(false,length(x))

    j = npol
    for i = 1:npol
        delta_xv = xv[j] - xv[i]
        delta_yv = yv[j] - yv[i]

        ## distance = [distance from (x,y) to edge] * length(edge)
        distance = delta_xv * (y .- yv[i]) - (x .- xv[i]) * delta_yv

        ## is y between the y-values of edge i,j AND (x,y) on the left of the edge?
        idx1 = ( (( (yv[i] .<= y) .& (y .< yv[j]) ) .| ( (yv[j] .<= y) .& (y .< yv[i]) ))
                 .& (0 .< distance*delta_yv) )
        in[idx1] = .! in[idx1]

        ## Check if (x,y) are actually on the boundary of the polygon.
        idx2 = ((( (yv[i] .<= y) .& (y .<= yv[j]) ) .| ( (yv[j] .<= y) .& (y .<= yv[i]) ))
                .& (( (xv[i] .<= x) .& (x .<= xv[j]) ) .| ( (xv[j] .<= x) .& (x .<= xv[i]) ))
                .& ( (0 .== distance) .| (delta_xv .== 0)))
        on[idx2] .= true

        j = i
    end

    return in .| on
    # return (in .| on, on)
end

function dpolygon( p, v )

    n_p = size(p,1)
    n_s = size(v,1)-1
    s = Vector{Float64}(undef, n_p)
    dist = Array{Float64}(undef, n_p,n_s)
    for i_s=1:n_s
        v_i = v[[i_s,i_s+1],:] # sprawdzana krawędz
        w   = v_i[2,:]-v_i[1,:] # wektor sprawdzanej krawędzi
        vp  = repeat(v_i[1,:]',n_p,1)-p # wektory pomiędzy punktami a sprawdzaną krawędzią
        w1  = repeat(w',n_p,1) # powielony o tyle ile jest punktów wektor sprawdzanej krawędzi
        s   = sum(w1.*vp,dims=2) # dot product dla kazdego punktu
        u   = -s./(w[1]^2+w[2]^2)
        u[u.<0.0] .= 0.0
        u[u.>1.0] .= 1.0
        h = w1.*[u u]+vp
        dist[:,i_s] .= sqrt.(sum(h.^2,dims=2))
    end
    # dist = (-1.0).^(inpolygonv(p[:,1],p[:,2],v[:,1],v[:,2])).*minimum(dist,2)
    dist = minimum(dist,dims=2)
    ind  = inpolygonv( p[:,1], p[:,2], v[:,1], v[:,2] )
    dist[ind] = -dist[ind]

    return dist
end


#------------------------------------------------------------------------------#
function l_deduplicate( a, s=[] )
    # Deduplicate array rows.

    if( isempty(a) )
        return (a, 0, 0)
    end
    if( isempty(s) )
        s = 1e-6*maximum(maximum(a)-minimum(a))
    end

    c = s*round.(a/s)

    ix = 1:size(c,1)
    c = [c ix[:]]
    c = sortslices( c, dims=1 )
    k = trunc.(Int,c[:,end])
    c = c[:,1:end-1]

    c_bool = c[1:size(c,1)-1,:] .!= c[2:size(c,1),:]
    ix = any(c_bool , dims=2)
    # dump(ix)
    j  = Vector{Int}(undef,length(k))
    j[k] = cumsum( [true;ix], dims=1 )
    tmp_A = (LinearIndices(ix))[findall(ix)]
    # dump(tmp_A)
    # tmp_B = getindex.(tmp_A, [1 2])
    i  = k[ [1;tmp_A[:,1] .+ 1] ]

    jj = sortperm(i)
    i  = i[jj]
    jjinv = Vector{Int}(undef,length(jj))
    jjinv[jj] = 1:length(jj)
    j  = vec(jjinv[j])

    b = a[i,:]

    return (b, i, j)
end
#------------------------------------------------------------------------------#
function l_simpvol( p, t )
    # SIMPVOL Triangle and tetrahedron volumes.

    n_sdim = size(p,2)
    if( n_sdim==2 )
        d12 = p[t[:,2],:] - p[t[:,1],:]
        d13 = p[t[:,3],:] - p[t[:,1],:]
        v = 0.5*( d12[:,1].*d13[:,2] - d12[:,2].*d13[:,1] )
    elseif( n_sdim==3 )
        d12 = p[t[:,2],:] - p[t[:,1],:]
        d13 = p[t[:,3],:] - p[t[:,1],:]
        d14 = p[t[:,4],:] - p[t[:,1],:]
        v = dot( cross(d12,d13,2), d14, 2 )/6
    else
        v = []
    end

    return v
end
#------------------------------------------------------------------------------#
function l_triangulate( p, fd::Function, e_fix=[], dtrm_tol=1e-4 )
    # Generate triangulation for vertices p.

    AV_TOL = eps()*1e1   # Minimum accepted absolute area/volume.
    # tmp_A = .!isnan.(p)
    # ind_keep = (LinearIndices(tmp_A))[findall( tmp_A )]
    # # rows, = ind2sub( size(p), ind_keep )
    # # i2s = CartesianIndices(p)
    # # rows = i2s[ind_keep]
    # rows=ind_keep
    # p  = p[unique(rows, dims=1),:]
    # p = filter((x) -> !isnan(x), p)
    rows = .!any(isnan.(p), dims=2)
    p = p[vec(rows), :]
    p, = l_deduplicate( p )

    # Generate triangulation for grid points p.
    # tic()
    t = l_delaunay_triangulation( p, e_fix )
    # td = toq()

    # Calculate simplex centers.
    n_sdim = size(p,2)
    pc = zeros(size(t,1),n_sdim)
    for i=1:n_sdim
        pc[:,i] = mean(reshape(p[t,i],size(t)),dims=2)
    end

    # Remove simplices with center outside region.
    dist = fd( pc )
    test=dist.<dtrm_tol
    t = t[vec(test),:]

    # # Reorient simplices.
    av = l_simpvol( p, t )
    ix_flip = vec(av.<0)
    t[ix_flip,[1 2]] = t[ix_flip,[2 1]]

    # Remove simplices with volume < AV_TOL.
    t = t[vec(any(abs.(av).>=AV_TOL, dims=2)),:]

    if( isempty(t) )
        # tic()
        t = l_delaunay_triangulation( p, e_fix )
        # td = td + toq()
    end

    return (p, t, 0)
end


#------------------------------------------------------------------------------#
function l_delaunay_triangulation( p, c )
    # Generate triangulation for vertices p.

    # Coordinate scaling.
    n_sdim = size(p,2)
    pp = Array{Any}(undef, n_sdim)
    for i=1:n_sdim
        pmin = minimum(p[:,i])
        pmax = maximum(p[:,i])
        pp[i] = ((p[:,i] .- pmin)/(pmax - pmin)*(max_coord - min_coord)) .+ min_coord
    end

    n_p = size(p,1)
    # dump(pp)
    ip = IndexedPoint2D[ IndexedPoint2D(pp[1][i],pp[2][i],i) for i in 1:n_p]
    tess = DelaunayTessellation2D{IndexedPoint2D}(n_p)

    push!( tess, ip )

    # Extract triangle indices.
    n_t = 0
    ind_t = Vector{Int64}(undef,length(tess._trigs))
    for i in 2:tess._last_trig_index
        isexternal(tess._trigs[i]) && continue
        n_t += 1
        ind_t[n_t] = i
    end
    t = Array{Int64}(undef,n_t,n_sdim+1)

    for i=1:n_t
        tri = tess._trigs[ind_t[i]]
        t[i,1] = tri._a._idx
        t[i,2] = tri._b._idx
        t[i,3] = tri._c._idx
    end

    return t
end


#------------------------------------------------------------------------------#
function l_uniqueperm( a )
    # Unique index permutation vector.

    ind = sortperm( a[:] )
    tmp = a[ind]
    msk = [ true; diff(tmp).!=0 ]

    ind1 = ind[ vec(msk) ]
    ind2 = Vector{Int}(undef,length(ind))
    ind2[ind] = cumsum(msk, dims=1)

    return (ind1, ind2)
end

function distmesh( fd=p->sqrt.(sum(p.^2,2))-1, fh=p->ones(size(p,1)), h0=0.1,
                   bbox=[-1 -1;1 1], p_fix=[], e_fix=[], it_max=1000, fid=1 )

    t0 = time_ns()
    #------------------------------------------------------------------------------#
    # Initialization and meshing parameters.
    #------------------------------------------------------------------------------#
    IT_MIN  = 20                # Minimum number of iterations.
    IT_MINC = 50                # Minimum number of iter. after which to call constraint function.
    IT_PRT  = 25                # Output every IT_PRT iterations.

    N_RECV  = 2                 # Number of recovery iteration steps to move points outside back to boundary.
    N_DCF   = 3000                # Frequency of density control checks.
    n_sdim  = size(bbox,2)
    # @printf "%i" n_sdim
    
    dp_tol   = -0.001*h0    # Abs point rejection tol (p(dist(p)>=dp0_tol) are rejected).
    dtrm_tol = -0.001*h0    # Abs dist tol for tri rejection (t(dist(p_tcent)>=dtrm_tol) are rejected).
    rt_tol   =  0.3         # Rel fraction of h0 to trigger retriangulation.
    F_scale  =  1.2         # Rel force scaling factor.
    F_DCF    =  2.0         # Fraction of L to L_target to allow.
    dp_scale =  0.2         # Rel fraction of computed new distance to move points in update step.

    dpc_tol = 0.001*h0          # Abs tol for grid point movements during convergence check.
    gradeps = sqrt(eps())*h0    # Gradient computation offset.
    #------------------------------------------------------------------------------#

    # Initial grid point distribution, p, confined to the bounding box.
    p = l_initpoints( h0, bbox )
    
    # Remove points outside the region and apply the rejection method.
    tmp_A = fd(p).<-dp_tol
    p = p[ (LinearIndices(tmp_A))[findall(tmp_A)], : ]
    if( isempty(p) )
        return (p, [], [])
    end

    r0 = fh( p )   # Probability to keep point.
    p  = p[ findall( rand(size(p,1)) .< minimum(r0)^n_sdim./r0.^n_sdim ), : ]
    p_fix,  = l_deduplicate( p_fix )
    n_p_fix = size(p_fix,1)
    if( !isempty(p_fix) )
        p, = l_deduplicate( [ p_fix; p ] )
    end
    n_p = size( p, 1 )

    # l_message( fid, "Grid generation (DistMesh):" )
    time_ns()
    if( it_max<=0 )
        t = l_delaunay_triangulation( p, e_fix )
        return (p, t, [])
    end
    t_tri = time_ns()
    it    = 0
    p0    = Inf
    n_tri = 0
    n_dcs = 0
    dist  = []
    delta_p = []
    is_converged = false

    ind = Vector{Int64}()
    edge_pairs = Matrix{Int64}[]

    while( it<it_max )
        it = it + 1

        # Retriangulate, if grid points have moved significantly.
        delta_p_max = maximum( sqrt.(sum((p.-p0).^2,dims=2)) )
        rt_h0 = rt_tol*h0
        if( rt_h0<delta_p_max )
            n_tri = n_tri + 1
            (p, t, td) = l_triangulate( p, fd, e_fix, dtrm_tol )
            # if( !isempty(e_fix) && it>IT_MINC )
            #     [p,t] = l_call_function( e_fix, p, t, n_sdim, 1:n_p_fix )
            # end
            p0  = p
            n_p = size(p,1)
            t_tri = t_tri + td

            # Describe each edge by a unique edge_pairs of nodes.
            e = [ t[:,[1,2]]; t[:,[2,3]]; t[:,[3,1]] ]

            e = sort(e,dims=2)
            e_max = maximum(e[:])
            if( e_max*(e_max+1)<floatmax() )
                ecomp = (e_max+1)*e[:,1] + e[:,2]
                ind, = l_uniqueperm( ecomp )
                edge_pairs = e[ind,:]
            else
                error( "row wise unique not supported" )
                # edge_pairs = unique( e, "rows" )
            end
        end


        # Move mesh points based on edge lengths L and forces F.
        p1 = p[edge_pairs[:,1],:]
        p2 = p[edge_pairs[:,2],:]
        bars = p1 - p2                  # Bar vectors.
        L = sqrt.(sum(bars.^2,dims=2))       # Bar lengths.
        hbars = fh( 0.5*( p1 + p2 ) )   # Rel bar mid point sizes.
        L_target = hbars*F_scale*(sum(L.^n_sdim)/sum(hbars.^n_sdim))^(1/n_sdim)   # Bar target lengths.

        # Density control, remove points that are too close to each other.
        if( mod(it,N_DCF)==0 && any(L_target.>F_DCF*L) )
            n_dcs = n_dcs + 1
            ind_del  = find( L_target .> F_DCF*L )
            ind_del  = setdiff(reshape(edge_pairs[ind,:],2*length(ind),1),1:n_p_fix)
            ind_keep = setdiff(1:n_p,ind_del)
            p   = p[ind_keep,:]
            n_p = size(p,1)
            p0  = Inf
            continue
        end

        # Compute grid point movements.
        F = L_target - L   # Scalar bar forces.
        F[vec(F.<0.0)] .= 0.0
        F_bar = F./L*ones(1,n_sdim).*bars
        delta_p = zeros(n_p,n_sdim)
        for i=1:n_sdim
            for j=1:size(edge_pairs,1)
                delta_p[edge_pairs[j,1],i] += F_bar[j,i]
                delta_p[edge_pairs[j,2],i] -= F_bar[j,i]
            end
        end
        delta_p[1:n_p_fix,:] .= 0.0
        delta_p = dp_scale * delta_p
        p = p + delta_p


        # Move grid points outside geometry back to the boundary.
        for jt=1:N_RECV

            dist = fd( p )
            ix = dist .> 0
            ix[1:n_p_fix] .= 0
            if( any(ix) )
                ixv = vec(ix)
                grad_dist = zeros(sum(ixv),n_sdim)
                for i=1:n_sdim
                    doff = zeros(1,n_sdim)
                    doff[i] = gradeps
                    p_ix =  p[ixv,:]
                    ones_ix = ones(sum(ix),1)
                    fd_i = (p_ix+ones_ix*doff)
                    dist_offset_i  = fd( fd_i )
                    grad_dist[:,i] = ( dist_offset_i - dist[ixv] )/gradeps
                end
                gradnm = sum( grad_dist.^2, dims=2 )
                p[ixv,:] -= dist[ixv]./gradnm*ones(1,n_sdim) .* grad_dist
            end
        end

        # break
        # Statistics/output.
        delta_p_max = abs( maximum( [sqrt.(sum(delta_p[vec(dist.<dp_tol),:].^2,dims=2)); -Inf] ) )
        # # if( mod(it,IT_PRT)==0 )
        # #     s = @sprintf( "Iteration %4i: %i vertices, %i cells, max(delta_p) = %g\n", it, size(p,1), size(t,1), delta_p_max )
        # #     l_message( fid, s )
        # # end


        # Check for convergence.
        if( (delta_p_max<dpc_tol) || size(t,1)<=2 || it>it_max )
            is_converged = delta_p_max<dpc_tol
            break
        end

    end

    
    return (p, t, [])
end

function plotgrid( p, t )
    # Plot triangulation.
    # dump(p)
    # dump(t)
    n_t = size(t,1)
    x = Array{Real}(undef,0);
    y = Array{Real}(undef,0);
    for i=1:n_t
        x = [x; NaN; p[t[i,:],1]]
        y = [y; NaN; p[t[i,:],2]]
    end

    Plots.plot( Plots.Shape(x,y), opacity=0.5 )
end

import Plots
pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
      1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5]
fd = p -> dpolygon( p, pv )
fh = p -> ones(size(p,1))
h0=0.1
bbox=[-1 -1;2 1]
p_fix=pv
e_fix=[]
it_max=2
(p, t) = distmesh( fd, fh, h0, bbox, p_fix, e_fix, it_max )
x = pv[:,1]
y = pv[:,2]
plotgrid( p, t )                                                        
# # y = fd(x)
Plots.plot!(x, y)
# scatter!(p[:,1], p[:,2],markersize=1)
# # scatter!(p[:,1], p[:,2])
# # pv[:,1]