struct Vector2D
    x::Float64
    y::Float64
end

getX(v::Vector2D) = v.x

getY(v::Vector2D) = v.y

Base.:+(v1::Vector2D, v2::Vector2D) = Vector2D(v1.x + v2.x, v1.y + v2.y)

Base.:-(v1::Vector2D, v2::Vector2D) = Vector2D(v1.x - v2.x, v1.y - v2.y)

Base.:(==)(v1::Vector2D, v2::Vector2D) = v1.x == v2.x && v1.y == v2.y

Base.:*(factor::Number, v::Vector2D) = Vector2D(factor * v.x, factor * v.y)

Base.:/(v::Vector2D, divisor::Number) = Vector2D(v.x / divisor, v.y / divisor)

getLengthSquared(v::Vector2D) = v.x * v.x + v.y * v.y

getLength(v::Vector2D) = sqrt(v.x * v.x + v.y * v.y)
