
struct GetmeSequentialConfig
    relaxationParameterRho::Float64
    maxIterations::Int
    qualityEvaluationCycleLength::Int
    maxNoImprovementCycles::Int
    penaltyInvalid::Float64
    penaltyRepeated::Float64
    penaltySuccess::Float64
    polygonTransformations::Vector{Mathematics.GeneralizedPolygonTransformation}

    function GetmeSequentialConfig(
        maxNumberOfPolygonNodes::Int,
        transformationSet::PolygonTransformationSet = PolygonTransformationSet.GETMeBookExamples)
        new(
            0.01,
            1_000_000,
            100,
            20,
            1.0e-4,
            1.0e-5,
            1.0e-3,
            getRegularizingPolygonTransformations(maxNumberOfPolygonNodes, transformationSet)
        )
    end
end
