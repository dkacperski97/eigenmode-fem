
struct SmoothingResult
    smoothingResultBase::SmoothingResultBase
    iterations::Int

    function SmoothingResult(
        algorithmName::String,
        mesh::Mesh.PolygonalMesh,
        smoothingWallClockTimeInSeconds::Float64,
        iterations::Int)
        new(
            SmoothingResultBase(algorithmName, mesh, smoothingWallClockTimeInSeconds),
            iterations
        )
    end
end