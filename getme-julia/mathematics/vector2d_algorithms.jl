function areEqual(first::Vector2D, second::Vector2D, tolerance::Float64)
    return getLengthSquared(first - second) <= tolerance * tolerance
end

function areEqual(first::Vector{Vector2D}, second::Vector{Vector2D}, tolerance::Float64)
    return all([areEqual(f, s, tolerance) for (f, s) in zip(first, second)])
end

function getBoundingBox(points::Vector{Vector2D})
    @assert !isempty(points) "Cannot compute bounding box for empty vector of points."
    xMin, xMax = extrema([p.x for p in points])
    yMin, yMax = extrema([p.y for p in points])
    return BoundingBox(xMin, xMax, yMin, yMax)
end

function getRandomVector(maxVectorLength::Float64)
    @assert maxVectorLength > 0.0 "Radius has to be > 0.0."
    radius = rand() * maxVectorLength
    angle = rand() * 2π
    return Vector2D(radius * cos(angle), radius * sin(angle))
end