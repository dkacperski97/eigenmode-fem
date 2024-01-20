struct BoundingBox
    xMin::Float64
    xMax::Float64
    yMin::Float64
    yMax::Float64

    function BoundingBox(xMin, xMax, yMin, yMax)
        if !(xMin <= xMax)
            throw(ArgumentError("xMin <= xMax expected."))
        end
        if !(yMin <= yMax)
            throw(ArgumentError("yMin <= yMax expected."))
        end
        new(xMin, xMax, yMin, yMax)
    end
end
