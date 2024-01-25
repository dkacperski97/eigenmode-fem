function basicGetmeSimultaneous(mesh::Mesh.PolygonalMesh, config::BasicGetmeSimultaneousConfig)
    checkTransformations(mesh, config.polygonTransformations)
    iteration = 0
    polygons = getPolygons(mesh)
    newNodePositions = [Vector2D(0.0, 0.0) for _ in 1:getNumberOfNodes(mesh)]

    stopWatch = StopWatch()
    while true
        for polygon in polygons
            numberOfPolygonNodes = getNumberOfNodes(polygon)
            transformedNodes = transformAndScaleElement(
                config.polygonTransformations[numberOfPolygonNodes], polygon,
                getNodes(mesh))
            for nodeNumber in 1:numberOfPolygonNodes
                newNodePositions[getNodeIndex(polygon, nodeNumber)] += transformedNodes[nodeNumber]
            end
        end
        maxSquaredNodeRelocationDistance = 0.0
        for nodeIndex in getNonFixedNodeIndices(mesh)
            newNodePosition = newNodePositions[nodeIndex] / length(getAttachedPolygonIndices(mesh, nodeIndex))
            maxSquaredNodeRelocationDistance = updateMaxSquaredNodeRelocationDistance(getNodes(mesh)[nodeIndex],
                                                    newNodePosition,
                                                    maxSquaredNodeRelocationDistance)
            getMutableNodes(mesh)[nodeIndex] = newNodePosition
        end
        if iteration >= config.maxIterations || maxSquaredNodeRelocationDistance <= config.maxSquaredNodeRelocationDistanceThreshold
            break
        end
        newNodePositions = [Vector2D(0.0, 0.0) for _ in 1:getNumberOfNodes(mesh)]
        iteration += 1
    end
    stopWatch.stop()
    return SmoothingResult("Basic GETMe simultaneous", mesh, getElapsedTimeInSeconds(stopWatch), iteration)
end

function getmeSimultaneous(mesh::Mesh.PolygonalMesh, config::GetmeSimultaneousConfig)
    checkTransformations(mesh, config.polygonTransformations)
    iteration = 0
    polygons = getPolygons(mesh)
    polygonMeanRatioValues = computeMeanRatioQualityNumberOfPolygons(polygons, getNodes(mesh))
    oldMeshQuality = MeshQuality(polygonMeanRatioValues, false)
    transformedNodeSums = [Vector2D(0.0, 0.0) for _ in 1:getNumberOfNodes(mesh)]
    newNodePositions = getNodes(mesh)
    nodeWeightSums = zeros(getNumberOfNodes(mesh))
    bestQMeanValue = getQMean(oldMeshQuality)
    bestQMeanNodes = getNodes(mesh)

    stopWatch = StopWatch()
    while true
        for polygonIndex in 1:length(polygons)
            polygon = polygons[polygonIndex]
            numberOfPolygonNodes = getNumberOfNodes(polygon)
            transformedNodes = transformScaleAndRelaxElement(
                config.polygonTransformations[numberOfPolygonNodes],
                config.relaxationParameterRho, polygon, getNodes(mesh))
            weight = config.weightExponentEta == 0.0 ? 1.0 :
                (1.0 - polygonMeanRatioValues[polygonIndex]) ^ config.weightExponentEta
            for nodeNumber in 1:numberOfPolygonNodes
                nodeIndex = getNodeIndex(polygon, nodeNumber)
                transformedNodeSums[nodeIndex] += weight * transformedNodes[nodeNumber]
                nodeWeightSums[nodeIndex] += weight
            end
        end
        for nodeIndex in getNonFixedNodeIndices(mesh)
            if nodeWeightSums[nodeIndex] > 0.0
                newNodePositions[nodeIndex] = transformedNodeSums[nodeIndex] / nodeWeightSums[nodeIndex]
            end
        end
        newMeshQuality = iterativelyResetNodesResultingInInvalidElementsSetNewMeshNodesAndUpdateElementQualityNumbers(
            newNodePositions, polygonMeanRatioValues, mesh)
        if getQMean(newMeshQuality) > bestQMeanValue
            bestQMeanValue = getQMean(newMeshQuality)
            bestQMeanNodes = getNodes(mesh)
        end
        qMeanImprovement = getQMean(newMeshQuality) - getQMean(oldMeshQuality)
        if iteration >= config.maxIterations || qMeanImprovement <= config.qMeanImprovementThreshold
            break
        end
        iteration += 1
        oldMeshQuality = newMeshQuality
        transformedNodeSums = [Vector2D(0.0, 0.0) for _ in 1:getNumberOfNodes(mesh)]
        nodeWeightSums = zeros(getNumberOfNodes(mesh))
        newNodePositions = getNodes(mesh)
    end
    stopWatch.stop()

    setNodes(mesh, bestQMeanNodes)
    return SmoothingResult("GETMe simultaneous", mesh, getElapsedTimeInSeconds(stopWatch), iteration)
end

function getmeSequential(mesh::Mesh.PolygonalMesh, config::GetmeSequentialConfig)
    algorithm = GetmeSequential(mesh, config)
    return getResult(algorithm)
end

function getme(mesh::Mesh.PolygonalMesh, config::GetmeConfig)
    getmeSimultaneousResult = getmeSimultaneous(mesh, config.getmeSimultaneousConfig)
    getmeSequentialResult = getmeSequential(getmeSimultaneousResult.mesh, config.getmeSequentialConfig)
    return GetmeResult(getmeSimultaneousResult, getmeSequentialResult)
end