
@testset "MinHeapEntry" begin
    polygonIndex = 17
    polygonMeanRatioNumber = 0.74
    isFixedPolygon = true
    entry = Smoothing.MinHeapEntry(polygonIndex, polygonMeanRatioNumber, isFixedPolygon)

    @test Smoothing.getPolygonIndex(entry) == polygonIndex
    @test Smoothing.getPenaltyCorrectedMeanRatioNumber(entry) == polygonMeanRatioNumber
    @test Smoothing.getMeanRatioNumber(entry) == polygonMeanRatioNumber
    @test Smoothing.isAllFixedNodesPolygon(entry) == isFixedPolygon
end

@testset "MinHeapEntry" begin
    polygonIndex = 3
    polygonMeanRatioNumber = 0.4
    isFixedPolygon = false
    initialPenaltySum = 0.06
    entry = Smoothing.MinHeapEntry(polygonIndex, polygonMeanRatioNumber, isFixedPolygon, initialPenaltySum)

    @test Smoothing.getPolygonIndex(entry) == polygonIndex
    @test Smoothing.getPenaltyCorrectedMeanRatioNumber(entry) == polygonMeanRatioNumber + initialPenaltySum
    @test Smoothing.getMeanRatioNumber(entry) == polygonMeanRatioNumber
    @test Smoothing.isAllFixedNodesPolygon(entry) == isFixedPolygon
end

@testset "MinHeapEntry" begin
    polygonIndex = 3
    polygonMeanRatioNumber = 0.4
    isFixedPolygon = false
    initialPenaltySum = 0.06
    entry = Smoothing.MinHeapEntry(polygonIndex, polygonMeanRatioNumber, isFixedPolygon, initialPenaltySum)

    @test Smoothing.getPolygonIndex(entry) == polygonIndex
    @test Smoothing.getPenaltyCorrectedMeanRatioNumber(entry) == polygonMeanRatioNumber + initialPenaltySum
    @test Smoothing.getMeanRatioNumber(entry) == polygonMeanRatioNumber
    @test Smoothing.isAllFixedNodesPolygon(entry) == isFixedPolygon
end

@testset "MinHeapEntry compareOperators" begin
    orderedEntries = [
        Smoothing.MinHeapEntry(14, 0.2, false, 0.0),
        Smoothing.MinHeapEntry(30, 0.18, false, 0.04),
        Smoothing.MinHeapEntry(9, 0.2, false, 0.02),
        Smoothing.MinHeapEntry(11, 0.2, false, 0.02),
        Smoothing.MinHeapEntry(2, 0.22, false, 0.0),
        Smoothing.MinHeapEntry(4, 0.1, false, 0.13),
        Smoothing.MinHeapEntry(1, 0.1, true, 0.0),
        Smoothing.MinHeapEntry(0, 0.2, true, 0.0)
    ]

    for firstIndex in 1:length(orderedEntries)
        for secondIndex in 1:length(orderedEntries)
            firstEntry = orderedEntries[firstIndex]
            secondEntry = orderedEntries[secondIndex]

            @test (firstIndex < secondIndex) == (firstEntry < secondEntry)
            @test (firstIndex == secondIndex) == (firstEntry == secondEntry)
            @test (firstIndex > secondIndex) == (firstEntry > secondEntry)
        end
    end
end

@testset "MinHeapEntry updateMethods_throwIfPolygonIsFixed" begin
    polygonIndex = 1
    polygonMeanRatioNumber = 0.3
    updatedPolygonMeanRatioNumber = 0.37
    entry = Smoothing.MinHeapEntry(polygonIndex, polygonMeanRatioNumber, true, 0.0)

    @test_throws ErrorException Smoothing.updateMeanRatioNumber(entry, updatedPolygonMeanRatioNumber)
    penaltyChange = 0.3
    @test_throws ErrorException Smoothing.updateMeanRatioNumberAndAddToPenaltySum(entry, updatedPolygonMeanRatioNumber, penaltyChange)
end

@testset "MinHeapEntry updateMethods" begin
    polygonIndex = 12
    polygonMeanRatioNumber = 0.33
    entry = Smoothing.MinHeapEntry(polygonIndex, polygonMeanRatioNumber, false, 0.0)

    firstPolygonMeanRatioNumberUpdate = 0.4
    Smoothing.updateMeanRatioNumber(entry, firstPolygonMeanRatioNumberUpdate)
    @test Smoothing.getMeanRatioNumber(entry) == firstPolygonMeanRatioNumberUpdate
    @test Smoothing.getPenaltyCorrectedMeanRatioNumber(entry) == firstPolygonMeanRatioNumberUpdate

    secondPolygonMeanRatioNumberUpdate = 0.6
    penaltyChange = 0.02
    Smoothing.updateMeanRatioNumberAndAddToPenaltySum(entry, secondPolygonMeanRatioNumberUpdate, penaltyChange)
    @test Smoothing.getMeanRatioNumber(entry) == secondPolygonMeanRatioNumberUpdate
    @test Smoothing.getPenaltyCorrectedMeanRatioNumber(entry) == secondPolygonMeanRatioNumberUpdate + penaltyChange

    # Penalty sum is >= 0.0 check.
    secondPenaltyChange = -1.0
    Smoothing.updateMeanRatioNumberAndAddToPenaltySum(entry, firstPolygonMeanRatioNumberUpdate, secondPenaltyChange)
    @test Smoothing.getMeanRatioNumber(entry) == firstPolygonMeanRatioNumberUpdate
    @test Smoothing.getPenaltyCorrectedMeanRatioNumber(entry) == firstPolygonMeanRatioNumberUpdate
end

@testset "MinHeapEntry addToPenaltySum" begin
    polygonIndex = 12
    polygonMeanRatioNumber = 0.33
    entry = Smoothing.MinHeapEntry(polygonIndex, polygonMeanRatioNumber, false, 0.0)

    penaltyChange = 0.4
    Smoothing.addToPenaltySum(entry, penaltyChange)
    @test Smoothing.getPenaltyCorrectedMeanRatioNumber(entry) == polygonMeanRatioNumber + penaltyChange

    # Penalty sum is >= 0.0 check.
    Smoothing.addToPenaltySum(entry, -2.0 * penaltyChange)
    @test Smoothing.getPenaltyCorrectedMeanRatioNumber(entry) == polygonMeanRatioNumber
end

@testset "PolygonQualityMinHeap constructorAndLowestQualityPolygonIndexGetter" begin
    mesh = Testdata.getMixedSampleMesh()
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    try
        Smoothing.isConsistent(minHeap)
        @test true
    catch e
        @test false
    end
    @test Smoothing.getLowestQualityPolygonIndex(minHeap) == 1
    @test Smoothing.isConsistent(minHeap)
end
    
struct TestData
    indexOfPolygonToAdjust::Int
    newPolygonMeanRatioNumber::Float64
    penaltyChange::Float64
    expectedNewWorstPolygonIndex::Int
end

@testset "PolygonQualityMinHeap updateTests" begin
    mesh = Testdata.getMixedSampleMesh()
    testCases = [
        TestData(1, 0.6, 0.0, 1),   # Polygon 1 is still worst after update.
        TestData(1, 0.6, 0.1, 1),   # Polygon 1 is still worst after update.
        TestData(1, 0.71, 0.0, 5),  # Polygon 1 is second worst.
        TestData(1, 1.0, 0.1, 5),   # Polygon 1 is best.
        TestData(6, 0.6, 0.0, 6),   # Make polygon 6 worst.
        TestData(4, 0.6, 0.0, 4),   # Make polygon 4 worst.
    ]

    for testData in testCases
        minHeap = Smoothing.PolygonQualityMinHeap(mesh)
        if testData.penaltyChange == 0.0
            Smoothing.updateMeanRatioNumberIfNotFixedPolygon(minHeap, testData.indexOfPolygonToAdjust, testData.newPolygonMeanRatioNumber)
        else
            Smoothing.updateMeanRatioNumberAndAddToPenaltySum(minHeap, testData.indexOfPolygonToAdjust, testData.newPolygonMeanRatioNumber, testData.penaltyChange)
        end
        newWorstPolygonIndex = Smoothing.getLowestQualityPolygonIndex(minHeap)
        @test testData.expectedNewWorstPolygonIndex == newWorstPolygonIndex
        try
            Smoothing.isConsistent(minHeap)
            @test true
        catch e
            @test false
        end
    end
end
    
@testset "PolygonQualityMinHeap addToPenaltySum" begin
    mesh = Testdata.getMixedSampleMesh()
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    Smoothing.addToPenaltySum(minHeap, 1, 0.1)  # Make polygon 1 second worst.
    @test Smoothing.getLowestQualityPolygonIndex(minHeap) == 5

    Smoothing.addToPenaltySum(minHeap, 5, 0.1)  # Make polygon 1 worst again.
    @test Smoothing.getLowestQualityPolygonIndex(minHeap) == 1

    Smoothing.addToPenaltySum(minHeap, 1, 0.2)  # Polygon 6 becomes worst.
    @test Smoothing.getLowestQualityPolygonIndex(minHeap) == 6

    try
        Smoothing.isConsistent(minHeap)
        @test true
    catch e
        @test false
    end
end

@testset "PolygonQualityMinHeap isAllFixedMesh notAllFixedNodes" begin
    mesh = Testdata.getMixedSampleMesh()
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    @test !Smoothing.isAllFixedMesh(minHeap)
end

@testset "PolygonQualityMinHeap isAllFixedMesh allFixedNodes" begin
    nodes = Testdata.getMixedSampleMeshNodes()
    polygons = Testdata.getMixedSampleMeshPolygons()
    fixedNodeIndices = Testdata.getMixedSampleMeshFixedNodeIndices()
    push!(fixedNodeIndices, 10, 11)  # Indices are 1-based in Julia
    mesh = Mesh.PolygonalMesh(nodes, polygons, fixedNodeIndices)
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    @test Smoothing.isAllFixedMesh(minHeap)
end

@testset "PolygonQualityMinHeap getQMinStar equalsQMin" begin
    mesh = Testdata.getMixedSampleMesh()
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    tolerance = 1.0e-15
    @test abs(Smoothing.getQMinStar(minHeap) - 6.298366572977736e-01) <= tolerance
end

@testset "PolygonQualityMinHeap getQMinStar_differsFromQMin" begin
    nodes = Testdata.getMixedSampleMeshNodes()
    polygons = Testdata.getMixedSampleMeshPolygons()
    fixedNodeIndices = Testdata.getMixedSampleMeshFixedNodeIndices()
    push!(fixedNodeIndices, 11)  # Indices are 1-based in Julia
    mesh = Mesh.PolygonalMesh(nodes, polygons, fixedNodeIndices)
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    tolerance = 1.0e-15
    @test abs(Smoothing.getQMinStar(minHeap) - 7.085662394599952e-01) <= tolerance
end

@testset "PolygonQualityMinHeap getQMinStar_throwIfAllFixedMesh" begin
    nodes =Testdata. getMixedSampleMeshNodes()
    polygons = Testdata.getMixedSampleMeshPolygons()
    fixedNodeIndices = Testdata.getMixedSampleMeshFixedNodeIndices()
    push!(fixedNodeIndices, 10, 11)  # Indices are 1-based in Julia
    mesh = Mesh.PolygonalMesh(nodes, polygons, fixedNodeIndices)
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    @test_throws Exception Smoothing.getQMinStar(minHeap)
end

@testset "PolygonQualityMinHeap containsAnInvalidPolygon" begin
    mesh = Testdata.getMixedSampleMesh()
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    @test !Smoothing.containsAnInvalidPolygon(minHeap)

    mesh = Testdata.getInvalidMixedSampleMesh()
    minHeap = Smoothing.PolygonQualityMinHeap(mesh)
    @test Smoothing.containsAnInvalidPolygon(minHeap)
end
