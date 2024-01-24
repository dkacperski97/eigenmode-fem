
struct GetmeSimultaneousConfig
    weightExponentEta::Float64
    relaxationParameterRho::Float64
    qMeanImprovementThreshold::Float64
    maxIterations::Int
    polygonTransformations::Vector{Mathematics.GeneralizedPolygonTransformation}

    function GetmeSimultaneousConfig(
        maxNumberOfPolygonNodes::Int,
        transformationSet::PolygonTransformationSet = PolygonTransformationSet.GETMeBookExamples)
        new(
            0.0,
            1.0,
            qMeanImprovementThreshold,
            maxIterations,
            getRegularizingPolygonTransformations(maxNumberOfPolygonNodes, transformationSet)
        )
    end
end
