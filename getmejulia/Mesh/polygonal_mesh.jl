# using Set

mutable struct PolygonalMesh
    nodes::Vector{Mathematics.Vector2D}
    polygons::Vector{Mathematics.Polygon}
    fixedNodeIndices::Set{Int}
    nonFixedNodeIndices::Vector{Int}
    areAllPolygonNodesFixed::Vector{Bool}
    indicesOfEdgeConnectedNodes::Vector{Set{Int}}
    attachedPolygonIndices::Vector{Set{Int}}
    indicesOfNeighborPolygons::Vector{Set{Int}}
    maximalNumberOfPolygonNodes::Int

    function PolygonalMesh(nodes::Vector{Mathematics.Vector2D}, polygons::Vector{Mathematics.Polygon}, fixedNodeIndices::Set{Int}=Set{Int}())
        newMesh = new(
            nodes, 
            polygons, 
            fixedNodeIndices, 
            Vector{Int}(), 
            fill(false, length(polygons)), 
            fill(Set{Int}(), length(nodes)), 
            fill(Set{Int}(), length(nodes)), 
            fill(Set{Int}(), length(polygons))
        )
        setNonFixedNodes(newMesh)
        setFixedPolygonAndNodeTopologyData(newMesh)
        setIndicesOfNeighborPolygons(newMesh)
        return newMesh
    end
end

getNodes(mesh::PolygonalMesh) = mesh.nodes
getMutableNodes(mesh::PolygonalMesh) = mesh.nodes
function setNodes(mesh::PolygonalMesh, newNodes::Vector{Mathematics.Vector2D})
    @assert length(mesh.nodes) == length(newNodes) "Non matching number of nodes."
    mesh.nodes = newNodes
end
getNumberOfNodes(mesh::PolygonalMesh) = length(mesh.nodes)
getPolygons(mesh::PolygonalMesh) = mesh.polygons
getFixedNodeIndices(mesh::PolygonalMesh) = mesh.fixedNodeIndices
getNonFixedNodeIndices(mesh::PolygonalMesh) = mesh.nonFixedNodeIndices
isFixedPolygon(mesh::PolygonalMesh, polygonIndex::Int) = mesh.areAllPolygonNodesFixed[polygonIndex]
getNumberOfPolygons(mesh::PolygonalMesh) = length(mesh.polygons)
getIndicesOfEdgeConnectedNodes(mesh::PolygonalMesh, nodeIndex::Int) = mesh.indicesOfEdgeConnectedNodes[nodeIndex]
getAttachedPolygonIndices(mesh::PolygonalMesh, nodeIndex::Int) = mesh.attachedPolygonIndices[nodeIndex]
getIndicesOfNeighborPolygons(mesh::PolygonalMesh, polygonIndex::Int) = mesh.indicesOfNeighborPolygons[polygonIndex]
getMaximalNumberOfPolygonNodes(mesh::PolygonalMesh) = mesh.maximalNumberOfPolygonNodes

function setNonFixedNodes(mesh::PolygonalMesh)
    numberOfMeshNodes = length(mesh.nodes)
    @assert length(mesh.fixedNodeIndices) <= numberOfMeshNodes "Number of fixed nodes cannot be larger than number of nodes."
    mesh.nonFixedNodeIndices = [nodeIndex for nodeIndex in 1:numberOfMeshNodes if !(nodeIndex in mesh.fixedNodeIndices)]
end

function setFixedPolygonAndNodeTopologyData(mesh::PolygonalMesh)
    isNodeFixedPredicate = nodeIndex -> nodeIndex in mesh.fixedNodeIndices

    for polygonIndex in 1:length(mesh.polygons)
        polygon = mesh.polygons[polygonIndex]
        mesh.areAllPolygonNodesFixed[polygonIndex] = all(isNodeFixedPredicate, Mathematics.getNodeIndices(polygon))
        numberOfPolygonNodes = Mathematics.getNumberOfNodes(polygon)
        if numberOfPolygonNodes > mesh.maximalNumberOfPolygonNodes
            mesh.maximalNumberOfPolygonNodes = numberOfPolygonNodes
        end
        for nodeNumber in 1:numberOfPolygonNodes
            predecessorNodeIndex = Mathematics.getPredecessorNodeIndex(polygon, nodeNumber)
            currentNodeIndex = Mathematics.getNodeIndex(polygon, nodeNumber)
            successorNodeIndex = Mathematics.getSuccessorNodeIndex(polygon, nodeNumber)

            @assert currentNodeIndex <= length(mesh.nodes) "Node index exceeds number of mesh nodes."

            push!(mesh.indicesOfEdgeConnectedNodes[currentNodeIndex], predecessorNodeIndex)
            push!(mesh.indicesOfEdgeConnectedNodes[currentNodeIndex], successorNodeIndex)

            push!(mesh.attachedPolygonIndices[currentNodeIndex], polygonIndex)
        end
    end
end

function setIndicesOfNeighborPolygons(mesh::PolygonalMesh)
    for polygonIndex in 1:length(mesh.polygons)
        for nodeIndex in Mathematics.getNodeIndices(mesh.polygons[polygonIndex])
            union!(mesh.indicesOfNeighborPolygons[polygonIndex], mesh.attachedPolygonIndices[nodeIndex])
        end
        delete!(mesh.indicesOfNeighborPolygons[polygonIndex], polygonIndex)
    end
end
