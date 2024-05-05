const expectedSampleMeshQMin = 6.298366572977736e-01
const expectedSampleMeshQMinStar = expectedSampleMeshQMin
const expectedSampleMeshQMean = 8.567191148264383e-01
const qualityTolerance = 1.0e-15
const invalidMeshMeanRatioNumbers = [0.1, 0.2, -1.0, -1.0, 0.3, 0.4, 0.5]

function getSampleMeshElementMeanRatioNumbers()
    mesh = Testdata.getMixedSampleMesh()
    return Mesh.computeMeanRatioQualityNumberOfPolygons(Mesh.getPolygons(mesh), Mesh.getNodes(mesh))
end

@testset "MeshQuality meshConstructorValidMesh" begin
    mesh = Testdata.getMixedSampleMesh()
    meshQuality = Mesh.MeshQuality(mesh)
    @test abs(expectedSampleMeshQMin - Mesh.getQMin(meshQuality)) <= qualityTolerance
    @test abs(expectedSampleMeshQMinStar - Mesh.getQMinStar(meshQuality)) <= qualityTolerance
    @test abs(expectedSampleMeshQMean - Mesh.getQMean(meshQuality)) <= qualityTolerance
    @test Mesh.getNumberOfInvalidElements(meshQuality) == 0
    @test Mesh.isValidMesh(meshQuality)
end

@testset "MeshQuality qualityConstructorValidMeshWithoutValidCount" begin
    determineNumberOfInvalidElements = false
    meshQuality = Mesh.MeshQuality(getSampleMeshElementMeanRatioNumbers(), determineNumberOfInvalidElements)
    @test abs(expectedSampleMeshQMin - Mesh.getQMin(meshQuality)) <= qualityTolerance
    @test Mesh.getQMinStar(meshQuality) == nothing
    @test abs(expectedSampleMeshQMean - Mesh.getQMean(meshQuality)) <= qualityTolerance
    @test Mesh.getNumberOfInvalidElements(meshQuality) == nothing
    @test Mesh.isValidMesh(meshQuality)
end

@testset "MeshQuality qualityConstructorValidMeshWithValidCount" begin
    determineNumberOfInvalidElements = true
    meshQuality = Mesh.MeshQuality(getSampleMeshElementMeanRatioNumbers(), determineNumberOfInvalidElements)
    @test abs(expectedSampleMeshQMin - meshQuality.getQMin()) <= qualityTolerance
    @test meshQuality.getQMinStar() == nothing
    @test abs(expectedSampleMeshQMean - meshQuality.getQMean()) <= qualityTolerance
    @test meshQuality.getNumberOfInvalidElements() == 0
    @test meshQuality.isValidMesh()
end

@testset "MeshQuality meshConstructorInvalidMesh" begin
    mesh = Testdata.getInvalidMixedSampleMesh()
    meshQuality = Mesh.MeshQuality(mesh)
    @test Mesh.getQMin(meshQuality) == -1.0
    @test Mesh.getQMinStar(meshQuality) == nothing
    @test Mesh.getQMean(meshQuality) == -1.0
    @test Mesh.getNumberOfInvalidElements(meshQuality) != nothing
    @test Mesh.getNumberOfInvalidElements(meshQuality) == 2
    @test !Mesh.isValidMesh(meshQuality)
end

@testset "MeshQuality qualityConstructorInvalidMeshWithoutValidCount" begin
    determineNumberOfInvalidElements = false
    meshQuality = Mesh.MeshQuality(invalidMeshMeanRatioNumbers, determineNumberOfInvalidElements)
    @test Mesh.getQMin(meshQuality) == -1.0
    @test Mesh.getQMinStar(meshQuality) == nothing
    @test Mesh.getQMean(meshQuality) == -1.0
    @test Mesh.getNumberOfInvalidElements(meshQuality) == nothing
    @test !Mesh.isValidMesh(meshQuality)
end

@testset "MeshQuality qualityConstructorInvalidMeshWithValidCount" begin
    determineNumberOfInvalidElements = true
    meshQuality = Mesh.MeshQuality(invalidMeshMeanRatioNumbers, determineNumberOfInvalidElements)
    @test Mesh.getQMin(meshQuality) == -1.0
    @test Mesh.getQMinStar(meshQuality) == nothing
    @test Mesh.getQMean(meshQuality) == -1.0
    @test Mesh.getNumberOfInvalidElements(meshQuality) != nothing
    @test Mesh.getNumberOfInvalidElements(meshQuality) == 2
    @test !Mesh.isValidMesh(meshQuality)
end

@testset "MeshQuality qMinAndQMinStarDiffer" begin
    nodes = Testdata.getMixedSampleMeshNodes()
    polygons = Testdata.getMixedSampleMeshPolygons()
    fixedNodeIndices = Testdata.getMixedSampleMeshFixedNodeIndices()
    push!(fixedNodeIndices, 10)  # Indices are 1-based in Julia
    push!(fixedNodeIndices, 11)  # Indices are 1-based in Julia
    deleteat!(fixedNodeIndices, findfirst(==(5), fixedNodeIndices))  # Indices are 1-based in Julia
    mesh = Mesh.PolygonalMesh(nodes, polygons, fixedNodeIndices)
    expectedMeanRatioQualityNumbers = getSampleMeshElementMeanRatioNumbers()
    expectedQminStar = expectedMeanRatioQualityNumbers[3]  # Indices are 1-based in Julia
    meshQuality = Mesh.MeshQuality(mesh)
    @test abs(expectedSampleMeshQMin - Mesh.getQMin(meshQuality)) <= qualityTolerance
    @test abs(expectedQminStar - Mesh.getQMinStar(meshQuality)) <= qualityTolerance
    @test abs(expectedSampleMeshQMean - Mesh.getQMean(meshQuality)) <= qualityTolerance
    @test Mesh.getNumberOfInvalidElements(meshQuality) == 0
    @test Mesh.isValidMesh(meshQuality)
end