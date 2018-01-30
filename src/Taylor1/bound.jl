
doc"""
Compute a rigorous bound for a Taylor1Model.
"""
function bound(f::Taylor1Model)
    x0, I, p = f.x0, f.I, f.p

    B = zero(I)
    n = degree(p)

    for i = n:-1:0
        B = B * (I - x0) + bound(p[i])
    end

    return B
end

function bound(p::Taylor1, x0, I)

    B = zero(I)
    n = degree(p)

    for i = n:-1:0
        B = B * (I - x0) + bound(p[i])
    end

    return B
end


bound(f::TaylorN) = evaluate(f)  # can replace by better polynomial bounder


# bound(f::Taylor1Model) = bound(f.p, f.x0, f.I)

bound(x::Interval) = x
