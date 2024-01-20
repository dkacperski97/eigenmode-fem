include( "DistMesh.jl" )
pv = [-0.4 -0.5;0.4 -0.2;0.4 -0.7;1.5 -0.4;0.9 0.1;
1.6 0.8;0.5 0.5;0.2 1;0.1 0.4;-0.7 0.7;-0.4 -0.5]
fd = p -> dpolygon( p, pv )
fh = p -> ones(size(p,1))
(p, t) = distmesh( fd, fh, 0.1, [-1 -1;2 1], pv )
plotgrid( p, t )