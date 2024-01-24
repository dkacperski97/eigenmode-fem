const PolygonalMeshKeyword = "planar_polygonal_mesh"
const NodesKeyword = "nodes"
const PolygonsKeyword = "polygons"
const FixedNodeIndicesKeyword = "fixed_node_indices"
const MeanRatioKeyword = "polygon_mean_ratio_quality_numbers"

# TODO: convert write methods to Julia

using Random

const PolygonalMeshKeyword = "planar_polygonal_mesh"
const NodesKeyword = "nodes"
const PolygonsKeyword = "polygons"
const FixedNodeIndicesKeyword = "fixed_node_indices"
const MeanRatioKeyword = "polygon_mean_ratio_quality_numbers"

function readPolygonalMeshHeader(infile::IOStream)
    line = readline(infile)
    containsMeshKeyword = occursin(PolygonalMeshKeyword, line)
    @assert containsMeshKeyword "Mesh type information not found."
end

function readMeshNodes(infile::IOStream)::Vector{Mathematics.Vector2D}
    keyword, numberOfNodes = split(readline(infile))
    @assert keyword == NodesKeyword "Nodes keyword expected but not found."
    nodes = Vector{Mathematics.Vector2D}(undef, parse(Int, numberOfNodes))
    for nodeIndex in 1:parse(Int, numberOfNodes)
        x, y = split(readline(infile))
        nodes[nodeIndex] = Mathematics.Vector2D(parse(Float64, x), parse(Float64, y))
    end
    return nodes
end

function readMeshPolygons(infile::IOStream)::Vector{Mathematics.Polygon}
    keyword, numberOfPolygons = split(readline(infile))
    @assert keyword == PolygonsKeyword "Polygons keyword expected but not found."
    polygons = Vector{Mathematics.Polygon}(undef, parse(Int, numberOfPolygons))
    for polygonIndex in 1:parse(Int, numberOfPolygons)
        nodeIndices = split(readline(infile))
        polygons[polygonIndex] = Polygon(parse.(Int, nodeIndices[2:end]))
    end
    return polygons
end

function readFixedNodeIndices(infile::IOStream)::Set{Int}
    keyword, numberOfEntries = split(readline(infile))
    @assert keyword == FixedNodeIndicesKeyword "Fixed node indices keyword expected but not found."
    fixedNodeIndices = Set{Int}()
    for index in 1:parse(Int, numberOfEntries)
        nodeIndex = parse(Int, readline(infile))
        push!(fixedNodeIndices, nodeIndex)
    end
    return fixedNodeIndices
end

function readMeshFile(infilePath::String)::PolygonalMesh
    @assert isfile(infilePath) "Did not find input file $infilePath."
    try
        open(infilePath, "r") do infile
            readPolygonalMeshHeader(infile)
            nodes = readMeshNodes(infile)
            polygons = readMeshPolygons(infile)
            fixedNodeIndices = readFixedNodeIndices(infile)
            return PolygonalMesh(nodes, polygons, fixedNodeIndices)
        end
    catch e
        throw(e)
    end
end

function computeMeanRatioQualityNumberOfPolygons(polygons::Vector{Mathematics.Polygon}, nodes::Vector{Mathematics.Vector2D}, meanRatioQualityNumbers::Vector{Float64})
    @assert length(polygons) == length(meanRatioQualityNumbers) "Mean ratio quality numbers vector size has to match number of polygons."
    for (i, polygon) in enumerate(polygons)
        meanRatioQualityNumbers[i] = Mathematics.getMeanRatio(polygon, nodes)
    end
end

function computeMeanRatioQualityNumberOfPolygons(polygons::Vector{Mathematics.Polygon}, nodes::Vector{Mathematics.Vector2D})::Vector{Float64}
    meanRatioNumbers = fill(-1.0, length(polygons))
    computeMeanRatioQualityNumberOfPolygons(polygons, nodes, meanRatioNumbers)
    return meanRatioNumbers
end

function areEqual(first::PolygonalMesh, second::PolygonalMesh, nodesEqualTolerance::Float64)::Bool
    return Mathematics.areEqual(first.nodes, second.nodes, nodesEqualTolerance) &&
           (first.polygons == second.polygons) &&
           (first.fixedNodeIndices == second.fixedNodeIndices)
end

function distortNodesLocally(mesh::PolygonalMesh, maxDistortionRadius::Float64)::PolygonalMesh
    newNodes = mesh.nodes
    for nodeIndex in mesh.nonFixedNodeIndices
        newNodes[nodeIndex] += Mathematics.getRandomVector(maxDistortionRadius)
    end
    return PolygonalMesh(newNodes, mesh.polygons, mesh.fixedNodeIndices)
end