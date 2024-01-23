
# Maximal number of smoothing iterations for algorithms improving all elements
# within one smoothing step.
const maxIterations = 10_000

# Quality based termination threshold for algorithms improving all elements
# within one smoothing step. Terminate if the improvement of the arithmetic
# mean of all element quality numbers of one iteration drops below the
# given threshold.
const qMeanImprovementThreshold = 1.0e-4

# Regularizing transformation set to use for polygons in GETMe variants.
@enum PolygonTransformationSet Generic GETMeBookExamples

# Get default regularizing polygon transformation set vector. Polygons with n
# nodes will be transformed by the n-th entry transformation of this vector.
function getRegularizingPolygonTransformations(
    maxNumberOfPolygonNodes::Int,
    transformationSet::PolygonTransformationSet)
    @assert maxNumberOfPolygonNodes >= 3 "Minimal valid number of polygon nodes is three."
    polygonTransformations = Vector{Mathematics.GeneralizedPolygonTransformation}(undef, maxNumberOfPolygonNodes + 1)
    for numberOfPolygonNodes in 0:maxNumberOfPolygonNodes
        polygonTransformations[numberOfPolygonNodes + 1] = numberOfPolygonNodes
    end
    if transformationSet == PolygonTransformationSet.GETMeBookExamples
        lambda = 0.5
        polygonTransformations[4] = Mathematics.GeneralizedPolygonTransformation(lambda, pi / 4.0)
        if maxNumberOfPolygonNodes >= 4
            polygonTransformations[5] = Mathematics.GeneralizedPolygonTransformation(lambda, pi / 6.0)
        end
    end
    return polygonTransformations
end