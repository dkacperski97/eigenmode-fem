
mutable struct LocalQualityResult
    areAllElementsValid::Bool
    transformedElementMeanRatioNumber::Float64
    neighborElementIndexAndMeanRatioNumber::Vector{Tuple{Int, Float64}}

    LocalQualityResult() = new(false, -1.0, Vector{Tuple{Int, Float64}}())
end

mutable struct GetmeSequential
    mesh::Mesh.PolygonalMesh
    config::GetmeSequentialConfig
    minHeap::PolygonQualityMinHeap
    isNodeFixed::Vector{Bool}
    temporaryNodes::Vector{Mesh.Node}
    smoothingTimeInSeconds::Float64
    iterationsApplied::Int
end

function GetmeSequential(mesh::Mesh.PolygonalMesh, config::GetmeSequentialConfig)
    minHeap = PolygonQualityMinHeap(mesh)
    getme = GetmeSequential(mesh, config, minHeap, Vector{Bool}(), Vector{Mesh.Node}(), 0.0, 0)
    checkInputData(getme)
    initHelperData(getme)
    applySmoothing(getme)
    @assert getme.minHeap.isConsistent(), "Inconsistent min heap."
    return getme
end

function getResult(getme::GetmeSequential)
    return SmoothingResult("GETMe sequential", getme.mesh, getme.smoothingTimeInSeconds, getme.iterationsApplied)
end

function checkInputData(getme::GetmeSequential)
    @assert !minHeap.containsAnInvalidPolygon(), "GETMe sequential can only be applied to valid initial meshes."
    @assert getme.config.qualityEvaluationCycleLength < getme.config.maxIterations, "Quality evaluation cycle length must be <= maximal number of iterations."
    checkTransformations(getme.mesh, getme.config.polygonTransformations)
end

function initHelperData(getme::GetmeSequential)
    getme.isNodeFixed = fill(false, Mesh.getNumberOfNodes(getme.mesh))
    for fixedNodeIndex in Mesh.getFixedNodeIndices(getme.mesh)
        getme.isNodeFixed[fixedNodeIndex] = true
    end
    getme.temporaryNodes = copy(Mesh.getNodes(getme.mesh)) # TODO: check if copy is necessary
end

function applySmoothing(getme::GetmeSequential)
    iteration = 0
    polygons = getme.mesh.polygons

    lastTransformedPolygonIndex = typemax(Int)

    lastQMinStar = getQMinStar(getme.minHeap)
    bestQMinStarValue = lastQMinStar
    bestQMinStarNodes = copy(getme.mesh.nodes) # TODO: check if copy is necessary
    numberOfConsecutiveNoImproveCycles = 0

    stopWatch = time()
    while true
        iteration += 1
        transformedPolygonIndex = getLowestQualityPolygonIndex(getme.minHeap)

        if lastTransformedPolygonIndex == transformedPolygonIndex
            addToPenaltySum(getme.minHeap, transformedPolygonIndex, getme.config.penaltyRepeated)
        end

        transformedPolygon = polygons[transformedPolygonIndex]
        transformPolygonAndSetTemporaryNodes(getme, transformedPolygon)
        localQualityInfo = assessLocalQuality(getme, transformedPolygonIndex)
        if !localQualityInfo.areAllElementsValid
            copyNodes(getme, transformedPolygonIndex, getme.mesh.nodes, getme.temporaryNodes)
            addToPenaltySum(getme.minHeap, transformedPolygonIndex, getme.config.penaltyInvalid)
        else
            copyNodes(getme, transformedPolygonIndex, getme.temporaryNodes, getme.mesh.nodes)
            updateMeanRatioNumberAndAddToPenaltySum(
                getme.minHeap,
                transformedPolygonIndex,
                localQualityInfo.transformedElementMeanRatioNumber,
                -getme.config.penaltySuccess)
            for (polygonIndex, newMeanRatioNumber) in localQualityInfo.neighborElementIndexAndMeanRatioNumber
                updateMeanRatioNumberIfNotFixedPolygon(getme.minHeap, polygonIndex, newMeanRatioNumber)
            end
        end
        lastTransformedPolygonIndex = transformedPolygonIndex

        if iteration % getme.config.qualityEvaluationCycleLength == 0
            qMinStar = getQMinStar(getme.minHeap)
            if qMinStar > bestQMinStarValue
                bestQMinStarValue = qMinStar
                bestQMinStarNodes = copy(getme.mesh.nodes)
                numberOfConsecutiveNoImproveCycles = 0
            else
                numberOfConsecutiveNoImproveCycles += 1
            end
        end
        if iteration == getme.config.maxIterations || numberOfConsecutiveNoImproveCycles == getme.config.maxNoImprovementCycles
            break
        end
    end
    stopWatch = time() - stopWatch

    getme.mesh.nodes = bestQMinStarNodes
    getme.iterationsApplied = iteration
    getme.smoothingTimeInSeconds = stopWatch
end

function transformPolygonAndSetTemporaryNodes(getme::GetmeSequential, polygon::Mathematics.Polygon)
    transformedNodes = transformScaleAndRelaxElement(
        getme.config.polygonTransformations[Mathematics.getNumberOfNodes(polygon)],
        getme.config.relaxationParameterRho, polygon, getme.mesh.nodes)
    for (nodeNumber, nodeIndex) in enumerate(polygon.nodeIndices)
        if !getme.isNodeFixed[nodeIndex]
            getme.temporaryNodes[nodeIndex] = transformedNodes[nodeNumber]
        end
    end
end

function assessLocalQuality(getme::GetmeSequential, transformedPolygonIndex::Int)
    polygons = getme.mesh.polygons
    result = LocalQualityResult()
    result.transformedElementMeanRatioNumber = Mathematics.getMeanRatio(polygons[transformedPolygonIndex], getme.temporaryNodes)
    if result.transformedElementMeanRatioNumber <= 0.0
        return result
    end
    for neighborPolygonIndex in getme.mesh.getIndicesOfNeighborPolygons(transformedPolygonIndex)
        neighborMeanRatioNumber = Mathematics.getMeanRatio(polygons[neighborPolygonIndex], getme.temporaryNodes)
        # Preliminary termination if one neighbor was invalidated.
        if neighborMeanRatioNumber <= 0.0
            return result
        end
        push!(result.neighborElementIndexAndMeanRatioNumber, (neighborPolygonIndex, neighborMeanRatioNumber))
    end
    result.areAllElementsValid = true
    return result
end

function copyNodes(getme::GetmeSequential, polygonIndex::Int, sourceNodes::Vector{Mathematics.Vector2D}, targetNodes::Vector{Mathematics.Vector2D})
    for nodeIndex in Mesh.getNodeIndices(Mesh.getPolygons(getme.mesh)[polygonIndex])
        targetNodes[nodeIndex] = sourceNodes[nodeIndex]
    end
end
