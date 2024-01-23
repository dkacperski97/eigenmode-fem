function applyEdgeLengthScaling(
    polygon::Mathematics.Polygon,
    originalMeshNodes::Vector{Mathematics.Vector2D},
    transformedElementNodes::Vector{Mathematics.Vector2D})
    commonPolygonCentroid = Mathematics.Vector2D(0.0, 0.0)
    originalPolygonLength = 0.0
    transformedPolygonLength = 0.0
    polygonNodeIndices = Mathematics.getNodeIndices(polygon)
    previousMeshNodeIndex = polygonNodeIndices[end]
    previousNodeIndex = length(polygonNodeIndices)
    for nodeNumber in 1:length(polygonNodeIndices)
        meshNodeIndex = polygonNodeIndices[nodeNumber]
        commonPolygonCentroid += originalMeshNodes[meshNodeIndex]
        originalPolygonLength += Mathematics.getLength(originalMeshNodes[meshNodeIndex] - originalMeshNodes[previousMeshNodeIndex])
        previousMeshNodeIndex = meshNodeIndex
        transformedPolygonLength += Mathematics.getLength(transformedElementNodes[nodeNumber] - transformedElementNodes[previousNodeIndex])
        previousNodeIndex = nodeNumber
    end
    commonPolygonCentroid /= length(polygonNodeIndices)
    scalingFactor = originalPolygonLength / transformedPolygonLength
    oneMinusScalingFactor = 1.0 - scalingFactor
    for nodeNumber in 1:length(polygonNodeIndices)
        transformedElementNodes[nodeNumber] = oneMinusScalingFactor * commonPolygonCentroid + scalingFactor * transformedElementNodes[nodeNumber]
    end
end

function updateMaxSquaredNodeRelocationDistance(
    oldNodePosition::Mathematics.Vector2D,
    newNodePosition::Mathematics.Vector2D,
    maxSquaredNodeRelocationDistance::Float64)

    squaredNodeRelocationDistance = Mathematics.getLengthSquared(newNodePosition - oldNodePosition)
    if squaredNodeRelocationDistance > maxSquaredNodeRelocationDistance
        maxSquaredNodeRelocationDistance = squaredNodeRelocationDistance
    end
end

function transformAndScaleElement(
    gpt::Mathematics.GeneralizedPolygonTransformation,
    polygon::Mathematics.Polygon,
    meshNodes::Vector{Mathematics.Vector2D})

    transformedElementNodes = Mathematics.getNodesOfTransformedPolygon(gpt, polygon, meshNodes)
    applyEdgeLengthScaling(polygon, meshNodes, transformedElementNodes)
    return transformedElementNodes
end

function transformScaleAndRelaxElement(
    gpt::Mathematics.GeneralizedPolygonTransformation,
    relaxationFactorRho::Float64,
    polygon::Mathematics.Polygon,
    meshNodes::Vector{Mathematics.Vector2D})

    newElementNodes = transformAndScaleElement(gpt, polygon, meshNodes)
    if relaxationFactorRho != 1.0
        oneMinusRho = 1.0 - relaxationFactorRho
        for nodeNumber in 1:Mathematics.getNumberOfNodes(polygon)
            newElementNodes[nodeNumber] = oneMinusRho * meshNodes[Mathematics.getNodeIndex(polygon, nodeNumber)] + relaxationFactorRho * newElementNodes[nodeNumber]
        end
    end
    return newElementNodes
end

function getIndicesOfNodesToReset(
    polygonMeanRatioValues::Vector{Float64},
    mesh::Mesh.PolygonalMesh)
    indicesOfNodesToReset = Set{Int}()
    for polygonIndex in 1:Mesh.getNumberOfPolygons(mesh)
        if polygonMeanRatioValues[polygonIndex] <= 0.0
            polygonNodeIndices = Mathematics.getNodeIndices(Mesh.getPolygons(mesh)[polygonIndex])
            union!(indicesOfNodesToReset, polygonNodeIndices)
        end
    end
    return indicesOfNodesToReset
end

function resetNodesAndGetAffectedPolygonIndices(
    indicesOfNodesToReset::Set{Int},
    mesh::Mesh.PolygonalMesh,
    newNodePositions::Vector{Mathematics.Vector2D})
    indicesOfAffectedPolygons = Set{Int}()
    for nodeIndex in indicesOfNodesToReset
        newNodePositions[nodeIndex] = Mesh.getNodes(mesh)[nodeIndex]
        for attachedPolygonIndex in Mesh.getAttachedPolygonIndices(mesh, nodeIndex)
            push!(indicesOfAffectedPolygons, attachedPolygonIndex)
        end
    end
    return indicesOfAffectedPolygons
end

function iterativelyResetNodesResultingInInvalidElementsSetNewMeshNodesAndUpdateElementQualityNumbers(
    newNodePositions::Vector{Mathematics.Vector2D},
    polygonMeanRatioValues::Vector{Float64},
    mesh::Mesh.PolygonalMesh)
    Mesh.computeMeanRatioQualityNumberOfPolygons(Mesh.getPolygons(mesh), newNodePositions, polygonMeanRatioValues)
    while true
        indicesOfNodesToReset = getIndicesOfNodesToReset(polygonMeanRatioValues, mesh)
        if isempty(indicesOfNodesToReset)
            break
        end
        indicesOfAffectedPolygons = resetNodesAndGetAffectedPolygonIndices(indicesOfNodesToReset, mesh, newNodePositions)
        for polygonIndex in indicesOfAffectedPolygons
            polygonMeanRatioValues[polygonIndex] = Mathematics.getMeanRatio(Mesh.getPolygons(mesh)[polygonIndex], newNodePositions)
        end
    end
    Mesh.setNodes(mesh, newNodePositions)
    return Mesh.MeshQuality(polygonMeanRatioValues, false)
end

function checkTransformations(
    maxNumberOfPolygonNodes::Int,
    transformations::Vector{Mathematics.GeneralizedPolygonTransformation})
    @assert maxNumberOfPolygonNodes <= length(transformations) "Regularizing transformations must be given for all polygon types."
    for numberOfPolygonNodes in 3:maxNumberOfPolygonNodes
        @assert Mathematics.isCounterclockwiseRegularizingTransformation(transformations[numberOfPolygonNodes], numberOfPolygonNodes) "Transformation for $(numberOfPolygonNodes)-gons is not regularizing."
    end
end

function checkTransformations(
    mesh::Mesh.PolygonalMesh,
    transformations::Vector{Mathematics.GeneralizedPolygonTransformation})
    checkTransformations(Mesh.getMaximalNumberOfPolygonNodes(mesh), transformations)
end
