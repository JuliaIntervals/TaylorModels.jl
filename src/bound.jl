
doc"""
Compute a rigorous bound for a TaylorModel.
"""
function bound(f::TaylorModel)
    x0, I, p = f.x0, f.I, f.p

    B = zero(I)
    n = degree(p)

    for i = n:-1:0
        B = B * (I - x0) + bound(p[i], f.bounds)
    end

    return B
end

function bound(p::Taylor1, x0, I)

    B = zero(I)
    n = degree(p)

    for i = n:-1:0
        B = B * (I - x0) + p[i]
    end

    return B
end

function bound(p::Taylor1, x0, I, bounds)

    B = zero(I)
    n = degree(p)

    for i = n:-1:0
        B = B * (I - x0) + bound(p[i], bounds)
    end

    return B
end


bound(f::TaylorN, bounds) = evaluate(f, bounds)  # can replace by better polynomial bounder


# bound(f::TaylorModel) = bound(f.p, f.x0, f.I)

bound(x::Interval, bounds) = x
