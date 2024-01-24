function getNodesOfRegularPolygon(numberOfNodes::Int)
    if numberOfNodes < 3
        throw(ArgumentError("At least three nodes required"))
    end
    nodes = Vector{Vector2D}(undef, numberOfNodes)
    for index in 1:numberOfNodes
        angle = index * 2.0 * π / numberOfNodes
        nodes[index] = Vector2D(cos(angle), sin(angle))
    end
    return nodes
end

function getMeanRatioSummand(polygon::Polygon, polygonNodeNumber::Int, nodes::Vector{Vector2D}, numberOfPolygonNodes::Int)
    predecessorNodeIndex = getPredecessorNodeIndex(polygon, polygonNodeNumber)
    centerNodeIndex = getNodeIndex(polygon, polygonNodeNumber)
    successorNodeIndex = getSuccessorNodeIndex(polygon, polygonNodeNumber)

    regularPolygonAngle = 2.0 * π / numberOfPolygonNodes
    a = cos(regularPolygonAngle) - 1.0
    b = sin(regularPolygonAngle)
    diffSuccessorCenter = nodes[successorNodeIndex] - nodes[centerNodeIndex]
    diffPredecessorCenter = nodes[predecessorNodeIndex] - nodes[centerNodeIndex]

    d11 = diffSuccessorCenter.x
    d12 = diffPredecessorCenter.x
    d21 = diffSuccessorCenter.y
    d22 = diffPredecessorCenter.y

    detS = (d12 * d21 - d11 * d22) / (2.0 * a * b)
    if detS < 0.0
        polygon.nodeIndices[2], polygon.nodeIndices[3] = polygon.nodeIndices[3], polygon.nodeIndices[2] # HACK: fix
        predecessorNodeIndex = getPredecessorNodeIndex(polygon, polygonNodeNumber)
        centerNodeIndex = getNodeIndex(polygon, polygonNodeNumber)
        successorNodeIndex = getSuccessorNodeIndex(polygon, polygonNodeNumber)

        regularPolygonAngle = 2.0 * π / numberOfPolygonNodes
        a = cos(regularPolygonAngle) - 1.0
        b = sin(regularPolygonAngle)
        diffSuccessorCenter = nodes[successorNodeIndex] - nodes[centerNodeIndex]
        diffPredecessorCenter = nodes[predecessorNodeIndex] - nodes[centerNodeIndex]

        d11 = diffSuccessorCenter.x
        d12 = diffPredecessorCenter.x
        d21 = diffSuccessorCenter.y
        d22 = diffPredecessorCenter.y

        detS = (d12 * d21 - d11 * d22) / (2.0 * a * b)
    end
    if detS < 0.0
        return -1.0
    end
    trace = ((d11 - d12) * (d11 - d12) + (d21 - d22) * (d21 - d22)) / (4.0 * b * b) + ((d11 + d12) * (d11 + d12) + (d21 + d22) * (d21 + d22)) / (4.0 * a * a)
    summand = detS / trace
    return summand
end

function getMeanRatio(polygon::Polygon, nodes::Vector{Vector2D})
    numberOfNodes = getNumberOfNodes(polygon)
    if numberOfNodes == 3
        summand = getMeanRatioSummand(polygon, 1, nodes, numberOfNodes)
        return summand < 0.0 ? -1.0 : min(1.0, 2.0 * summand)
    else
        sum = 0.0
        for nodeNumber in 1:numberOfNodes
            summand = getMeanRatioSummand(polygon, nodeNumber, nodes, numberOfNodes)
            if summand < 0.0
                return -1.0
            end
            sum += summand
        end
        return min(1.0, 2.0 * sum / numberOfNodes)
    end
end
