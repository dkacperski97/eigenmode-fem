struct BasicGetmeSimultaneousConfig
    maxSquaredNodeRelocationDistanceThreshold::Float64
    maxIterations::Int
    polygonTransformations::Vector{Mathematics.GeneralizedPolygonTransformation}

    function BasicGetmeSimultaneousConfig(
        maxNodeRelocationDistanceThreshold::Float64,
        maxNumberOfPolygonNodes::Int,
        transformationSet::PolygonTransformationSet = PolygonTransformationSet.GETMeBookExamples)
        new(
            maxNodeRelocationDistanceThreshold^2,
            maxIterations,
            getRegularizingPolygonTransformations(maxNumberOfPolygonNodes, transformationSet)
        )
    end
end
