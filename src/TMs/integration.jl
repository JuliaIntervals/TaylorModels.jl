# integration.jl

# integrate
function integrate(a::TM1AbsRem{T}, c0::Interval{T}) where {T}
    order = get_order(a)
    integ_pol = integrate(a.pol, c0)
    aux = (a.iI-a.x0)
    Δ = aux * remainder(a) +  a.pol[order] * aux^(order+1)
    return TM1AbsRem( integ_pol, Δ, a.x0, a.iI )
end
integrate(a::TM1AbsRem{T}) where {T} = integrate(a, Interval(zero(T)))
