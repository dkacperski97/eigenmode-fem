# helper functions
using LinearAlgebra
include("whitney.jl")
include("meshing.jl")
include("fem.jl")
include("distmesh-julia/test.jl")
include("getmejulia/GetmeJulia.jl")

# Parameters
ε0 = 1.
μ0 = 1.
εr = 1.
μr = 1.
σ  = 0.
ε  = ε0 * εr
μ  = μ0 * μr

a, b = 1.3, 0.9
n, m = 13, 9

# using PyPlot

a, b = 1.3, 0.9

pv = [0 0; a 0; a b; 0 b; 0 0]
fd = p -> dpolygon( p, pv )
fh = p -> ones(size(p,1))
h0=0.15
bbox=[0 0;a b]
p_fix=pv
e_fix=[]
it_max=10
(p, t) = distmesh( fd, fh, h0, bbox, p_fix, e_fix, it_max )
sort!(t, dims=2) # element nodes have to be sorted
x = pv[:,1]
y = pv[:,2]
plotgrid( p, t )

epsilon = 0.001
# bottom_row_points
p[p[:, 1] .< epsilon, 1] .= 0
# top_row_points
p[p[:, 2] .> b - epsilon, 2] .= b
# left_column_points
p[p[:, 2] .< epsilon, 2] .= 0
# right_column_points
p[p[:, 1] .> a - epsilon, 1] .= a


plotgrid( p, t )                                                        
# # y = fd(x)
Plots.plot!(x, y)

import .Smoothing

# function PolygonalMesh(
#   nodes::Vector{Mathematics.Vector2D}, 
#   polygons::Vector{Mathematics.Polygon}, 
#   fixedNodeIndices::Set{Int}=Set{Int}()
# 
nodes = Vector{Mathematics.Vector2D}(undef,size(p,1))
for index in 1:size(p,1)
  nodes[index] = Mathematics.Vector2D(p[index,1], p[index,2])
end
polygons = Vector{Mathematics.Polygon}(undef,size(t,1))
for index in 1:size(t,1)
  polygons[index] = Mathematics.Polygon([t[index,1], t[index,2], t[index,3]])
end
# TODO: fixedNodeIndices
mesh = Mesh.PolygonalMesh(nodes, polygons)

config = Smoothing.GetmeSequentialConfig(
    Mesh.getMaximalNumberOfPolygonNodes(mesh), 
    Smoothing.Generic
)

results = Smoothing.getmeSequential(mesh, config)
dump(results)