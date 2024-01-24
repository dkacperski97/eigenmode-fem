import .Mesh: PolygonalMesh, MeshQuality

struct SmoothingResultBase
    algorithmName::String
    mesh::Mesh.PolygonalMesh
    meshQuality::Mesh.MeshQuality
    smoothingWallClockTimeInSeconds::Float64

    function SmoothingResultBase(
        algorithmName::String,
        mesh::Mesh.PolygonalMesh,
        smoothingWallClockTimeInSeconds::Float64)
        new(
            algorithmName,
            mesh,
            Mesh.MeshQuality(mesh),
            smoothingWallClockTimeInSeconds
        )
    end
end
