
doc"""
Compute a rigorous bound for a Taylor1Model.
"""
bound(f::Taylor1Model) = bound(f, f.I)

doc"""
Compute a rigorous bound for a `Taylor1Model` `f` evaluated over the interval `I`.
This interval must be a subset of `f.I`.
"""
function bound(f::Taylor1Model, I)
    if !(I ⊆ f.I)
        error("Cannot evaluate `Taylor1Model` on interval $I that is not included in interval of definition, $(f.I)")
    end

    return bound(f.p, f.x0, I) + f.Δ
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
