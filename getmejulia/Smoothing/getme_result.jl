struct GetmeResult
    smoothingResultBase::SmoothingResultBase
    getmeSimultaneousIterations::Int
    getmeSequentialIterations::Int

    function GetmeResult(
        getmeSimultaneousSmoothingResult::SmoothingResult,
        getmeSequentialSmoothingResult::SmoothingResult)
        new(
            SmoothingResultBase(
                "GETMe",
                getmeSequentialSmoothingResult.mesh,
                getmeSimultaneousSmoothingResult.smoothingWallClockTimeInSeconds
                    + getmeSequentialSmoothingResult.smoothingWallClockTimeInSeconds
            ),
            getmeSimultaneousSmoothingResult.iterations,
            getmeSequentialSmoothingResult.iterations
        )
    end
end