
@testset "Smoothing applyEdgeLengthScaling" begin
    polygon = Mathematics.Polygon([1, 2, 3, 4, 5, 6])  # Indices are 1-based in Julia
    meshNodes = [
        Mathematics.Vector2D(3.0, 5.0),  
        Mathematics.Vector2D(-9.0, -2.0), 
        Mathematics.Vector2D(7.0, 3.0),
        Mathematics.Vector2D(9.0, -7.0), 
        Mathematics.Vector2D(4.0, 4.0),   
        Mathematics.Vector2D(5.0, -9.0)
    ]
    transformedNodes = [
        Mathematics.Vector2D(2.120029401627464, 3.050931073468095),
        Mathematics.Vector2D(-3.993368337253149, -0.308356785127020),
        Mathematics.Vector2D(3.629849854773023, -2.437418888363173),
        Mathematics.Vector2D(7.820578277134416, -2.594900517549698),
        Mathematics.Vector2D(4.845821607436489, 0.038824027685005),
        Mathematics.Vector2D(4.577089196281754, -3.749078910113208)
    ]
    expectedNodes = [
        Mathematics.Vector2D(0.706538394108473, 8.521742056117311),
        Mathematics.Vector2D(-13.663046103965184, 0.625712254157441),
        Mathematics.Vector2D(4.255382016765926, -4.378663234046956),
        Mathematics.Vector2D(14.105718727334228, -4.748824914040057),
        Mathematics.Vector2D(7.113532196663723, 1.441763203056728),
        Mathematics.Vector2D(6.481874769092831, -7.461729365244467)
    ]
    Smoothing.applyEdgeLengthScaling(polygon, meshNodes, transformedNodes)
    tolerance = 1.0e-14
    @test Mathematics.areEqual(expectedNodes, transformedNodes, tolerance)
end

@testset "Smoothing updateMaxSquaredNodeRelocationDistance" begin
    oldNode = Mathematics.Vector2D(3.0, -2.0)
    newNode = Mathematics.Vector2D(1.0, 2.0)
    maxSquaredNodeRelocationDistance = 0.0

    maxSquaredNodeRelocationDistance = Smoothing.updateMaxSquaredNodeRelocationDistance(oldNode, newNode, maxSquaredNodeRelocationDistance)
    @test maxSquaredNodeRelocationDistance == 20.0

    largerMaxValue = 100.0
    maxSquaredNodeRelocationDistance = largerMaxValue
    maxSquaredNodeRelocationDistance = Smoothing.updateMaxSquaredNodeRelocationDistance(oldNode, newNode, maxSquaredNodeRelocationDistance)
    @test maxSquaredNodeRelocationDistance == largerMaxValue
end

@testset "Smoothing transformAndScaleElement" begin
    lambda = 0.3
    theta = pi / 5.0
    transformation = Mathematics.GeneralizedPolygonTransformation(lambda, theta)
    polygon = Mathematics.Polygon([1, 2, 3, 4, 5])  # Indices are 1-based in Julia
    meshNodes = [
        Mathematics.Vector2D(6.0, -8.0), 
        Mathematics.Vector2D(8.0, -4.0), 
        Mathematics.Vector2D(-7.0, 1.0), 
        Mathematics.Vector2D(8.0, 9.0), 
        Mathematics.Vector2D(3.0, 9.0)
    ]
    expectedNodes = [
        Mathematics.Vector2D(+0.955412051458474, -6.836136432133172),
        Mathematics.Vector2D(+9.828471223820758, +2.138754261004771),
        Mathematics.Vector2D(+0.138965075543568, +1.053163587137887),
        Mathematics.Vector2D(+9.598106992441132, +3.244484494788414),
        Mathematics.Vector2D(-2.520955343263935, +7.399734089202100)
    ]
    transformedAndScaledNodes = Smoothing.transformAndScaleElement(transformation, polygon, meshNodes)
    tolerance = 1.0e-14
    @test Mathematics.areEqual(expectedNodes, transformedAndScaledNodes, tolerance)
end

@testset "Smoothing transformScaleAndRelaxElement" begin
    lambda = 0.3
    theta = pi / 5.0
    transformation = Mathematics.GeneralizedPolygonTransformation(lambda, theta)
    polygon = Mathematics.Polygon([1, 2, 3, 4, 5])  # Indices are 1-based in Julia
    meshNodes = [
        Mathematics.Vector2D(6.0, -8.0), 
        Mathematics.Vector2D(8.0, -4.0), 
        Mathematics.Vector2D(-7.0, 1.0), 
        Mathematics.Vector2D(8.0, 9.0), 
        Mathematics.Vector2D(3.0, 9.0)
    ]
    expectedNodes = [
        Mathematics.Vector2D(+2.468788436020932, -7.185295502493220),
        Mathematics.Vector2D(+9.279929856674531, +0.297127982703340),
        Mathematics.Vector2D(-2.002724447119503, +1.037214510996521),
        Mathematics.Vector2D(+9.118674894708793, +4.971139146351890),
        Mathematics.Vector2D(-0.864668740284754, +7.879813862441470)
    ]
    relaxationFactorRho = 0.7
    transformedAndScaledNodes = Smoothing.transformScaleAndRelaxElement(transformation, relaxationFactorRho, polygon, meshNodes)
    tolerance = 1.0e-14
    @test Mathematics.areEqual(expectedNodes, transformedAndScaledNodes, tolerance)
end

@testset "Smoothing iterativelyResetNodesResultingInInvalidElementsSetNewMeshNodesAndUpdateElementQualityNumbers" begin
    mesh = Testdata.getMixedSampleMesh()
    meanRatioNumbers = Mesh.computeMeanRatioQualityNumberOfPolygons(Mesh.getPolygons(mesh), Mesh.getNodes(mesh))
    newNodePositions = Mesh.getNodes(mesh)
    newNodePositions[10] = Mathematics.Vector2D(12.0, 3.0)  # Indices are 1-based in Julia
    newNodePositions[11] = Mathematics.Vector2D(4.0, 3.0)   # Indices are 1-based in Julia
    expectedMesh = copy(mesh)
    Mesh.getMutableNodes(expectedMesh)[11] = newNodePositions[11]  # Indices are 1-based in Julia
    expectedMeanRatioNumbers = Mesh.computeMeanRatioQualityNumberOfPolygons(Mesh.getPolygons(expectedMesh), Mesh.getNodes(expectedMesh))
    expectedMeshQuality = Mesh.MeshQuality(expectedMesh)
    meshQuality = Smoothing.iterativelyResetNodesResultingInInvalidElementsSetNewMeshNodesAndUpdateElementQualityNumbers(newNodePositions, meanRatioNumbers, mesh)
    @test expectedMeshQuality.qMin == meshQuality.qMin
    @test expectedMeshQuality.qMean == meshQuality.qMean
    @test !meshQuality.getNumberOfInvalidElements().has_value()
    @test !meshQuality.getQMinStar().has_value()
    @test Mesh.areEqual(expectedMesh, mesh, 0.0)
    @test expectedMeanRatioNumbers == meanRatioNumbers
    @test Mesh.getNodes(expectedMesh) == newNodePositions
end

@testset "Smoothing checkTransformations" begin
    maxNumberOfPolygonNodes = 5
    transformations = [Mathematics.GeneralizedPolygonTransformation(i) for i in 1:maxNumberOfPolygonNodes]
    @test try
        Smoothing.checkTransformations(maxNumberOfPolygonNodes, transformations)
        true
    catch
        false
    end
end