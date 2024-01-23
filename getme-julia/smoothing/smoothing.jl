module Smoothing

import ..mathematics.Mathematics
import ..mesh.Mesh

include("common_algorithms.jl")
include("default_configuration.jl")
include("basic_getme_simultaneous_config.jl")
include("basic_laplace_config.jl")
include("getme_sequential_config")
include("getme_simultaneous_config")
include("getme_config")
include("smoothing_result_base.jl")
include("smoothing_result.jl")
include("getme_result.jl")
include("polygon_quality_min_heap.jl")
include("getme_sequential.jl")

include("getme_algorithms.jl")

end