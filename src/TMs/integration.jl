# integrate

function integrate(a::TM1AbsRem{T}, c0::Interval{T}) where {T}
    order = get_order(a)
    integ_pol = integrate(a.pol, c0)
    # Is this too wide?
    Δ = (a.iI-a.x0)^order
    Δ = (a.iI-a.x0) * ( remainder(a) +  a.pol[order] * Δ )
    return TM1AbsRem( integ_pol, Δ, a.x0, a.iI )
end
integrate(a::TM1AbsRem{T}) where {T} = integrate(a, Interval(zero(T)))
