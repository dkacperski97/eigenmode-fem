include("../GetmeJulia.jl")

import .Smoothing
import .Mesh

using Test
# using Smoothing

module Testdata
# include("../GetmeJulia.jl")
# import .Smoothing
# import .Mesh

function getMixedSampleMeshNodes()
    return Main.Mathematics.Vector2D[
        Main.Mathematics.Vector2D(0.0, 0.0),   # 0
        Main.Mathematics.Vector2D(5.0, -1.0),  # 1
        Main.Mathematics.Vector2D(7.0, -2.0),  # 2
        Main.Mathematics.Vector2D(9.0, 0.0),   # 3
        Main.Mathematics.Vector2D(9.0, 2.0),   # 4
        Main.Mathematics.Vector2D(9.0, 5.0),   # 5
        Main.Mathematics.Vector2D(6.0, 5.0),   # 6
        Main.Mathematics.Vector2D(3.0, 5.0),   # 7
        Main.Mathematics.Vector2D(0.0, 3.0),   # 8
        Main.Mathematics.Vector2D(6.0, 2.0),   # 9
        Main.Mathematics.Vector2D(3.0, 1.0),   # 10
    ]
end

function getMixedSampleMeshPolygons()
    return Main.Mathematics.Polygon[
        Main.Mathematics.Polygon([1, 2, 11]),       # 0
        Main.Mathematics.Polygon([2, 10, 11]),      # 1
        Main.Mathematics.Polygon([2, 3, 4, 5, 10]), # 2
        Main.Mathematics.Polygon([5, 6, 7, 10]),    # 3
        Main.Mathematics.Polygon([10, 7, 11]),      # 4
        Main.Mathematics.Polygon([7, 8, 9, 11]),    # 5
        Main.Mathematics.Polygon([1, 11, 9]),       # 6
    ]
end

function getMixedSampleMeshFixedNodeIndices()
    return Set([1, 2, 3, 4, 5, 6, 7, 8, 9])
end

function getMixedSampleMesh()
    return Main.Mesh.PolygonalMesh(getMixedSampleMeshNodes(), getMixedSampleMeshPolygons(), getMixedSampleMeshFixedNodeIndices())
end

function getInvalidMixedSampleMesh()
    invalidMesh = getMixedSampleMesh()
    invalidMesh.nodes[10] = Main.Mathematics.Vector2D(17.0, 2.0)
    return invalidMesh
end

end

include("Mesh/mesh_quality_test.jl")

# include("Smoothing/polygon_quality_min_heap_test.jl")
# include("Smoothing/common_algorithms_test.jl")
