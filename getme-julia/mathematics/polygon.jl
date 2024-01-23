struct Polygon
    nodeIndices::Vector{Int}

    function Polygon(nodeIndices::Vector{Int})
        if length(nodeIndices) < 3
            throw(DomainError(nodeIndices, "Polygon must consist of at least 3 nodes."))
        end
        if length(nodeIndices) != length(unique(nodeIndices))
            throw(DomainError(nodeIndices, "Duplicate node indices are not allowed in polygonal elements."))
        end
        new(nodeIndices)
    end
end

function getNodeIndices(polygon::Polygon)
    return polygon.nodeIndices
end

function getNumberOfNodes(polygon::Polygon)
    return length(polygon.nodeIndices)
end

function getNodeIndex(polygon::Polygon, nodeNumber::Int)
    return polygon.nodeIndices[nodeNumber]
end

function getPredecessorNodeIndex(polygon::Polygon, nodeNumber::Int)
    return polygon.nodeIndices[nodeNumber == 1 ? end : nodeNumber - 1]
end

function getSuccessorNodeIndex(polygon::Polygon, nodeNumber::Int)
    return polygon.nodeIndices[nodeNumber == end ? 1 : nodeNumber + 1]
end

Base.:(==)(p1::Polygon, p2::Polygon) = p1.nodeIndices == p2.nodeIndices