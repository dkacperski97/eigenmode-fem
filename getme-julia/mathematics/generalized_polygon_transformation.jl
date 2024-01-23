struct GeneralizedPolygonTransformation
    lambda::Float64
    theta::Float64
    c1::Float64
    c2::Float64
    c3::Float64

    function GeneralizedPolygonTransformation(lambda::Float64, theta::Float64)
        if !(lambda > 0.0 && lambda < 1.0)
            throw(DomainError(lambda, "Lambda has to be in (0,1)."))
        end
        if !(theta > 0.0 && theta < π / 2.0)
            throw(DomainError(theta, "Theta has to be in (0,π/2)."))
        end
        c1 = (1.0 - lambda) * tan(theta)
        c2 = lambda * (1.0 - lambda) - c1 * c1
        c3 = 1.0 - 2.0 * c2
        new(lambda, theta, c1, c2, c3)
    end

    function GeneralizedPolygonTransformation(numberOfPolygonNodes::Int)
        new(0.5, numberOfPolygonNodes < 3 ? π / 4.0 : π / numberOfPolygonNodes)
    end
end

getLambda(gpt::GeneralizedPolygonTransformation) = gpt.lambda

getTheta(gpt::GeneralizedPolygonTransformation) = gpt.theta


function getNodesOfTransformedPolygon(gpt::GeneralizedPolygonTransformation, polygon::Polygon, nodes::Vector{Vector2D})
    transformedNodes = Vector{Vector2D}(undef, getNumberOfNodes(polygon))
    for nodeNumber in 1:getNumberOfNodes(polygon)
        predecessorNode = nodes[getPredecessorNodeIndex(polygon, nodeNumber)]
        node = nodes[getNodeIndex(polygon, nodeNumber)]
        successorNode = nodes[getSuccessorNodeIndex(polygon, nodeNumber)]
        transformedNodes[nodeNumber] = 
            gpt.c1 * Vector2D(successorNode.y - predecessorNode.y, predecessorNode.x - successorNode.x) + 
            gpt.c2 * (predecessorNode + successorNode) + 
            gpt.c3 * node
    end
    return transformedNodes
end

function getEigenvalues(gpt::GeneralizedPolygonTransformation, numberOfPolygonNodes::Int)
    eigenvalues = zeros(numberOfPolygonNodes)
    w = complex(gpt.lambda, (1 - gpt.lambda) * tan(gpt.theta))
    r = exp(2im * π / numberOfPolygonNodes)
    for k in 1:numberOfPolygonNodes
        magnitude = abs(1.0 - conj(w) + r^k * conj(w))
        eigenvalues[k] = magnitude * magnitude
    end
    return eigenvalues
end

function isCounterclockwiseRegularizingTransformation(gpt::GeneralizedPolygonTransformation, numberOfPolygonNodes::Int)
    eigenvalues = getEigenvalues(gpt, numberOfPolygonNodes)
    expectedDominantEigenvalue = eigenvalues[2]
    isLessOrEqualToDominating(eigenvalue) = expectedDominantEigenvalue >= eigenvalue
    return all(isLessOrEqualToDominating, eigenvalues)
end
