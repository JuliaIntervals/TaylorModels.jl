function calculate_set(U::TaylorN, V::TaylorN, bounds)

    xs = []
    ys = []

    # iterate over 4 sides of initial box

    num_points = 30

    b = bounds[2].lo
    for a in linspace(bounds[1].lo, bounds[2].hi, num_points)
        push!(xs, U([a, b]))
        push!(ys, V([a, b]))
    end

    a = bounds[1].hi
    for b in linspace(bounds[2].lo, bounds[2].hi, num_points)
        push!(xs, U([a, b]))
        push!(ys, V([a, b]))
    end

    b = bounds[2].hi
    for a in linspace(bounds[2].hi, bounds[2].lo, num_points)
        push!(xs, U([a, b]))
        push!(ys, V([a, b]))
    end

    a = bounds[1].lo
    for b in linspace(bounds[2].hi, bounds[2].lo, num_points)
        push!(xs, U([a, b]))
        push!(ys, V([a, b]))
    end

    return mid.(xs), mid.(ys)

end
