struct Polygon
    nodeIndices::Vector{Int}

    function Polygon(nodeIndices::Vector{Int})
        if !(length(nodeIndices) >= 3)
            throw(ArgumentError("Polygon must consist of at least 3 nodes."))
        end

        uniqueNodeIndices = Set(nodeIndices)
        if !(length(nodeIndices) == length(uniqueNodeIndices))
            throw(ArgumentError("Duplicate node indices are not allowed in polygonal elements."))
        end

        new(nodeIndices)
    end
end