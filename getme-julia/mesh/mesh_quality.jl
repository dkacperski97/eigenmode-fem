using Statistics

function computeQMinAndQMeanTerminateIfInvalid(polygonMeanRatioQualityNumbers::Vector{Float64})
    qMin = Inf
    sumOfMeanRatioNumbers = 0.0
    for meanRatioNumber in polygonMeanRatioQualityNumbers
        if meanRatioNumber <= 0.0
            return (-1.0, -1.0)
        end
        if meanRatioNumber < qMin
            qMin = meanRatioNumber
        end
        sumOfMeanRatioNumbers += meanRatioNumber
    end
    return (qMin, sumOfMeanRatioNumbers / length(polygonMeanRatioQualityNumbers))
end

function computeQMinQMinStarQMeanAndNumberOfInvalid(polygonMeanRatioQualityNumbers::Vector{Float64}, mesh::PolygonalMesh=nothing)
    qMin = Inf
    qMinStar = Inf
    qMinStarOptional = nothing
    sumOfMeanRatioNumbers = 0.0
    numberOfInvalid = 0
    for polygonIndex in 1:length(polygonMeanRatioQualityNumbers)
        meanRatioNumber = polygonMeanRatioQualityNumbers[polygonIndex]
        if meanRatioNumber <= 0.0
            numberOfInvalid += 1
        end
        sumOfMeanRatioNumbers += meanRatioNumber
        if meanRatioNumber < qMin
            qMin = meanRatioNumber
        end
        if mesh !== nothing && !isFixedPolygon(mesh, polygonIndex) && meanRatioNumber < qMinStar
            qMinStar = meanRatioNumber
        end
    end
    if numberOfInvalid > 0
        return (-1.0, qMinStarOptional, -1.0, numberOfInvalid)
    end
    if mesh !== nothing && qMinStar <= 1.0
        qMinStarOptional = qMinStar
    end
    return (qMin, qMinStarOptional, sumOfMeanRatioNumbers / length(polygonMeanRatioQualityNumbers), 0)
end

mutable struct MeshQuality
    qMin::Float64
    qMinStar::Union{Float64, Nothing}
    qMean::Float64
    numberOfInvalidElements::Union{Int, Nothing}

    function MeshQuality(mesh::PolygonalMesh)
        elementMeanRatioNumbers = computeMeanRatioQualityNumberOfPolygons(getPolygons(mesh), getNodes(mesh))
        newMeshQuality = new()
        newMeshQuality.qMin, newMeshQuality.qMinStar, newMeshQuality.qMean, newMeshQuality.numberOfInvalidElements = computeQMinQMinStarQMeanAndNumberOfInvalid(elementMeanRatioNumbers, mesh)
        return newMeshQuality
    end

    function MeshQuality(elementMeanRatioNumbers::Vector{Float64}, determineNumberOfInvalidElements::Bool)
        newMeshQuality = new()
        if determineNumberOfInvalidElements
            newMeshQuality.qMin, newMeshQuality.qMinStar, newMeshQuality.qMean, newMeshQuality.numberOfInvalidElements = computeQMinQMinStarQMeanAndNumberOfInvalid(elementMeanRatioNumbers)
        else
            newMeshQuality.qMin, newMeshQuality.qMean = computeQMinAndQMeanTerminateIfInvalid(elementMeanRatioNumbers)
        end
        return newMeshQuality
    end
end

getQMin(quality::MeshQuality) = quality.qMin
getQMinStar(quality::MeshQuality) = quality.qMinStar
getQMean(quality::MeshQuality) = quality.qMean
getNumberOfInvalidElements(quality::MeshQuality) = quality.numberOfInvalidElements
isValidMesh(quality::MeshQuality) = quality.qMin > 0.0
