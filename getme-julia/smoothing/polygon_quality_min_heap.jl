
mutable struct MinHeapEntry
    isFixedPolygon::Bool
    penaltyCorrectedMeanRatioNumber::Float64
    meanRatioNumber::Float64
    qualityPenaltySum::Float64
    polygonIndex::Int

    MinHeapEntry(polygonIndex::Int, initialMeanRatioNumber::Float64, isFixedPolygon::Bool, initialQualityPenaltySum::Float64=0.0) = new(
        isFixedPolygon,
        initialMeanRatioNumber + max(0.0, initialQualityPenaltySum),
        initialMeanRatioNumber,
        initialQualityPenaltySum,
        polygonIndex
    )
end

function getPolygonIndex(entry::MinHeapEntry)
    return entry.polygonIndex
end

function getPenaltyCorrectedMeanRatioNumber(entry::MinHeapEntry)
    return entry.penaltyCorrectedMeanRatioNumber
end

function getMeanRatioNumber(entry::MinHeapEntry)
    return entry.meanRatioNumber
end

function isAllFixedNodesPolygon(entry::MinHeapEntry)
    return entry.isFixedPolygon
end

function Base.:<(a::MinHeapEntry, b::MinHeapEntry)
    if a.isFixedPolygon != b.isFixedPolygon
        return a.isFixedPolygon < b.isFixedPolygon
    elseif a.penaltyCorrectedMeanRatioNumber != b.penaltyCorrectedMeanRatioNumber
        return a.penaltyCorrectedMeanRatioNumber < b.penaltyCorrectedMeanRatioNumber
    elseif a.meanRatioNumber != b.meanRatioNumber
        return a.meanRatioNumber < b.meanRatioNumber
    elseif a.qualityPenaltySum != b.qualityPenaltySum
        return a.qualityPenaltySum < b.qualityPenaltySum
    else
        return a.polygonIndex < b.polygonIndex
    end
end

function Base.:==(a::MinHeapEntry, b::MinHeapEntry)
    return a.isFixedPolygon == b.isFixedPolygon && 
           a.penaltyCorrectedMeanRatioNumber == b.penaltyCorrectedMeanRatioNumber && 
           a.meanRatioNumber == b.meanRatioNumber && 
           a.qualityPenaltySum == b.qualityPenaltySum && 
           a.polygonIndex == b.polygonIndex
end

function updateMeanRatioNumber(entry::MinHeapEntry, newPolygonMeanRatioNumber::Float64)
    if entry.isFixedPolygon
        error("All fixed polygon cannot change quality.")
    end
    entry.meanRatioNumber = newPolygonMeanRatioNumber
    entry.penaltyCorrectedMeanRatioNumber = newPolygonMeanRatioNumber + entry.qualityPenaltySum
end

function updateMeanRatioNumberAndAddToPenaltySum(entry::MinHeapEntry, newPolygonMeanRatioNumber::Float64, penaltyChange::Float64)
    if entry.isFixedPolygon
        error("All fixed polygon cannot change quality.")
    end
    entry.meanRatioNumber = newPolygonMeanRatioNumber
    entry.qualityPenaltySum = max(0.0, entry.qualityPenaltySum + penaltyChange)
    entry.penaltyCorrectedMeanRatioNumber = entry.meanRatioNumber + entry.qualityPenaltySum
end

function addToPenaltySum(entry::MinHeapEntry, penaltyChange::Float64)
    entry.qualityPenaltySum = max(0.0, entry.qualityPenaltySum + penaltyChange)
    entry.penaltyCorrectedMeanRatioNumber = entry.meanRatioNumber + entry.qualityPenaltySum
end


# Min heap for lowest penalty corrected polygon quality lookup. Implementation
# uses a binary tree stored in a vector (cf.
# https://en.wikipedia.org/wiki/Binary_heap). To speed up arbitrary polygon
# quality adjustments, a lookup table for polygon index to heap vector index is
# also kept in sync.

mutable struct PolygonQualityMinHeap
    binaryTree::Vector{MinHeapEntry}
    polygonIndexToBinaryTreeEntryIndex::Vector{Int}

    function PolygonQualityMinHeap(mesh::Mesh.PolygonalMesh)
        polygonIndexToBinaryTreeEntryIndex = fill(typemax(Int), Mesh.getNumberOfPolygons(mesh))
        binaryTree = Vector{MinHeapEntry}(undef,Mesh.getNumberOfPolygons(mesh))

        heap = new(binaryTree, polygonIndexToBinaryTreeEntryIndex)

        meanRatioQualityNumbers = Mesh.computeMeanRatioQualityNumberOfPolygons(Mesh.getPolygons(mesh), Mesh.getNodes(mesh))
        for polygonIndex in 1:(Mesh.getNumberOfPolygons(mesh))
            # TODO chceck how polygonIndex is used
            push!(heap.binaryTree, MinHeapEntry(polygonIndex, meanRatioQualityNumbers[polygonIndex], Mesh.isFixedPolygon(mesh, polygonIndex)))
            heap.polygonIndexToBinaryTreeEntryIndex[polygonIndex] = polygonIndex
            minHeapifyEntryOfPolygon(heap, polygonIndex)
        end
        return heap
    end
end

function getLowestQualityPolygonIndex(heap::PolygonQualityMinHeap)
    return getPolygonIndex(heap.binaryTree[1])
end

function updateMeanRatioNumberIfNotFixedPolygon(heap::PolygonQualityMinHeap, polygonIndex::Int, newPolygonMeanRatioNumber::Float64)
    entry = heap.binaryTree[heap.polygonIndexToBinaryTreeEntryIndex[polygonIndex]]
    if isAllFixedNodesPolygon(entry)
        return
    end
    updateMeanRatioNumber(entry, newPolygonMeanRatioNumber)
    minHeapifyEntryOfPolygon(heap, polygonIndex)
end

function updateMeanRatioNumberAndAddToPenaltySum(heap::PolygonQualityMinHeap, polygonIndex::Int, newPolygonMeanRatioNumber::Float64, penaltyChange::Float64)
    entry = heap.binaryTree[heap.polygonIndexToBinaryTreeEntryIndex[polygonIndex]]
    updateMeanRatioNumberAndAddToPenaltySum(entry, newPolygonMeanRatioNumber, penaltyChange)
    minHeapifyEntryOfPolygon(heap, polygonIndex)
end

function addToPenaltySum(heap::PolygonQualityMinHeap, polygonIndex::Int, penaltyChange::Float64)
    entry = heap.binaryTree[heap.polygonIndexToBinaryTreeEntryIndex[polygonIndex]]
    addToPenaltySum(entry, penaltyChange)
    minHeapifyEntryOfPolygon(heap, polygonIndex)
end

function isConsistent(heap::PolygonQualityMinHeap)::Bool
    doSizesMatch = length(heap.binaryTree) == length(heap.polygonIndexToBinaryTreeEntryIndex)
    return doSizesMatch && isPolygonIndexToBinaryTreeEntryIndexConsistent(heap) && isBinaryTreeConsistent(heap)
end

function isPolygonIndexToBinaryTreeEntryIndexConsistent(heap::PolygonQualityMinHeap)::Bool
    for polygonIndex in 1:length(heap.polygonIndexToBinaryTreeEntryIndex)
        indexInBinaryTree = heap.polygonIndexToBinaryTreeEntryIndex[polygonIndex]
        if polygonIndex != heap.binaryTree[indexInBinaryTree].polygonIndex
            return false
        end
    end
    indexedBinaryTreePolygons = sort(heap.polygonIndexToBinaryTreeEntryIndex) # TODO: check if heap also has to be sorted
    for index in 1:length(heap.polygonIndexToBinaryTreeEntryIndex)
        if index != indexedBinaryTreePolygons[index]
            return false
        end
    end
    return true
end

function isBinaryTreeConsistent(heap::PolygonQualityMinHeap)::Bool
    numberOfEntries = length(heap.binaryTree)
    for entryIndex in 1:numberOfEntries
        leftChildIndex = 2 * entryIndex
        if leftChildIndex <= numberOfEntries && isFirstQualityLower(heap, leftChildIndex, entryIndex)
            return false
        end
        rightChildIndex = leftChildIndex + 1
        if rightChildIndex <= numberOfEntries && isFirstQualityLower(heap, rightChildIndex, entryIndex)
            return false
        end
    end
    return true
end

function isAllFixedMesh(heap::PolygonQualityMinHeap)
    return isAllFixedNodesPolygon(heap.binaryTree[1])
end

function getQMinStar(heap::PolygonQualityMinHeap)::Float64
    if isAllFixedMesh(heap)
        error("QMinStar is not defined for all fixed polygon meshes.")
    end
    qMinStar = Inf
    for entry in heap.binaryTree
        if !isAllFixedNodesPolygon(entry) && entry.meanRatioNumber < qMinStar
            qMinStar = entry.meanRatioNumber
        end
    end
    return qMinStar
end

function containsAnInvalidPolygon(heap::PolygonQualityMinHeap)::Bool
    isInvalidPolygonPredicate = entry -> entry.meanRatioNumber < 0.0
    return any(isInvalidPolygonPredicate, heap.binaryTree)
end

function isFirstQualityLower(heap::PolygonQualityMinHeap, firstEntryIndex::Int, secondEntryIndex::Int)
    return heap.binaryTree[firstEntryIndex] < heap.binaryTree[secondEntryIndex]
end

function swapMinHeapEntriesAndAdjustMapping(heap::PolygonQualityMinHeap, firstEntryIndex::Int, secondEntryIndex::Int)
    heap.binaryTree[firstEntryIndex], heap.binaryTree[secondEntryIndex] = heap.binaryTree[secondEntryIndex], heap.binaryTree[firstEntryIndex]
    heap.polygonIndexToBinaryTreeEntryIndex[heap.binaryTree[firstEntryIndex].polygonIndex] = firstEntryIndex
    heap.polygonIndexToBinaryTreeEntryIndex[heap.binaryTree[secondEntryIndex].polygonIndex] = secondEntryIndex
end

function minHeapifyEntryOfPolygon(heap::PolygonQualityMinHeap, polygonIndex::Int)
    entryIndex = heap.polygonIndexToBinaryTreeEntryIndex[polygonIndex]

    # Up-swaps
    while entryIndex > 1
        parentEntryIndex = div(entryIndex, 2)
        if isFirstQualityLower(heap, parentEntryIndex, entryIndex)
            # Correct order - terminate.
            break
        end
        swapMinHeapEntriesAndAdjustMapping(heap, parentEntryIndex, entryIndex)
        entryIndex = parentEntryIndex
    end

    # Down-swaps
    maxEntryIndex = length(heap.binaryTree)
    while entryIndex <= maxEntryIndex
        leftChildIndex = 2 * entryIndex
        rightChildIndex = leftChildIndex + 1
        doCompareWithLeftChild = leftChildIndex <= maxEntryIndex
        doCompareWithRightChild = rightChildIndex <= maxEntryIndex

        # Start comparing with lowest quality child.
        if doCompareWithLeftChild && doCompareWithRightChild && isFirstQualityLower(heap, rightChildIndex, leftChildIndex)
            leftChildIndex, rightChildIndex = rightChildIndex, leftChildIndex
            doCompareWithLeftChild, doCompareWithRightChild = doCompareWithRightChild, doCompareWithLeftChild
        end

        if doCompareWithLeftChild && isFirstQualityLower(heap, leftChildIndex, entryIndex)
            swapMinHeapEntriesAndAdjustMapping(heap, leftChildIndex, entryIndex)
            entryIndex = leftChildIndex
        elseif doCompareWithRightChild && isFirstQualityLower(heap, rightChildIndex, entryIndex)
            swapMinHeapEntriesAndAdjustMapping(heap, rightChildIndex, entryIndex)
            entryIndex = rightChildIndex
        else
            # Polygon quality is smaller than both children.
            break
        end
    end
end
