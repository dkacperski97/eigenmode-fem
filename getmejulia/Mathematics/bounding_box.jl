struct BoundingBox
    xMin::Float64
    xMax::Float64
    yMin::Float64
    yMax::Float64

    function BoundingBox(xMin::Float64, xMax::Float64, yMin::Float64, yMax::Float64)
        if !(xMin <= xMax)
            throw(DomainError((xMin, xMax), "xMin <= xMax expected."))
        end
        if !(yMin <= yMax)
            throw(DomainError((yMin, yMax), "yMin <= yMax expected."))
        end
        new(xMin, xMax, yMin, yMax)
    end
end

Base.:(==)(b1::BoundingBox, b2::BoundingBox) = b1.xMin == b2.xMin && b1.xMax == b2.xMax && b1.yMin == b2.yMin && b1.yMax == b2.yMax

getXMin(b::BoundingBox) = b.xMin

getXMax(b::BoundingBox) = b.xMax

getYMin(b::BoundingBox) = b.yMin

getYMax(b::BoundingBox) = b.yMax

getXDimension(b::BoundingBox) = b.xMax - b.xMin

getYDimension(b::BoundingBox) = b.yMax - b.yMin