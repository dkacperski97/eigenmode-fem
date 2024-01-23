
struct GetmeConfig
    getmeSimultaneousConfig::GetmeSimultaneousConfig
    getmeSequentialConfig::GetmeSequentialConfig

    function GetmeConfig(
        maxNumberOfPolygonNodes::Int,
        transformationSet::PolygonTransformationSet = PolygonTransformationSet.GETMeBookExamples)
        new(
            GetmeSimultaneousConfig(maxNumberOfPolygonNodes, transformationSet),
            GetmeSequentialConfig(maxNumberOfPolygonNodes, transformationSet)
        )
    end

    function GetmeConfig(
        getmeSimultaneousConfig::GetmeSimultaneousConfig,
        getmeSequentialConfig::GetmeSequentialConfig)
        new(
            getmeSimultaneousConfig,
            getmeSequentialConfig
        )
    end
end