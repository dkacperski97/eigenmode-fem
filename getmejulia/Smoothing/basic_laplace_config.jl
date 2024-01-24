struct BasicLaplaceConfig
    maxSquaredNodeRelocationDistanceThreshold::Float64
    maxIterations::Int

    function BasicLaplaceConfig(maxNodeRelocationDistanceThreshold::Float64)
        new(
            maxNodeRelocationDistanceThreshold^2,
            maxIterations
        )
    end
end
