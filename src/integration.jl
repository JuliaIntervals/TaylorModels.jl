# integration.jl

"""
    integrate(a, c0)

Integrates the one-variable Taylor Model (`TaylorModel1`
or `RTaylorModel1`) with respect to the independent variable; `c0` is
the interval representing the integration constant; if omitted
it is considered as the zero interval.
"""
function integrate(a::TaylorModel1{T,S}, c0::T) where {T,S}
    integ_pol = integrate(a.pol, c0)
    δ = a.dom-a.x0
    Δ = bound_integration(a, δ)
    return TaylorModel1( integ_pol, Δ, a.x0, a.dom )
end
function integrate(a::TaylorModel1{T,S}, c0, cc0) where {T<:TaylorN, S}
    integ_pol = integrate(a.pol, c0)
    δ = a.dom - a.x0
    Δ = bound_integration(a, δ, cc0)
    return TaylorModel1(integ_pol, Δ, a.x0, a.dom)
end

integrate(a::TaylorModel1{T,S}) where {T,S} = integrate(a, zero(T))
integrate(a::TaylorModel1{T,S}, δI) where {T<:TaylorN,S} = integrate(a, zero(T), δI)

function integrate(a::RTaylorModel1{T,S}, c0::T) where {T,S}
    integ_pol = integrate(a.pol, c0)
    δ = a.dom - a.x0
    Δ = bound_integration(a, δ)
    return RTaylorModel1( integ_pol, Δ, a.x0, a.dom )
end
function integrate(a::RTaylorModel1{T,S}, c0, cc0) where {T<:TaylorN,S}
    integ_pol = integrate(a.pol, c0)
    δ = a.dom - a.x0
    Δ = bound_integration(a, δ, cc0)
    return RTaylorModel1( integ_pol, Δ, a.x0, a.dom )
end
integrate(a::RTaylorModel1{T,S}) where {T,S} = integrate(a, zero(T))
integrate(a::RTaylorModel1{T,S}, δI) where {T<:TaylorN,S} = integrate(a, zero(T), δI)

function integrate(a::TaylorModel1{TaylorModelN{N,T,S},S},
        c0::TaylorModelN{N,T,S}) where {N,T,S}
    integ_pol = integrate(a.pol, c0)
    δ = a.dom-a.x0

    # Remainder bound after integrating
    Δ = bound_integration(a, δ)
    ΔN = Δ(a[0].dom - a[0].x0)

    return TaylorModel1( integ_pol, ΔN, a.x0, a.dom )
end

"""
    integrate(fT, which)

Integrates a `TaylorModelN` with respect to `which` variable.
The returned `TaylorModelN` corresponds to the Taylor Model
of the definite integral ∫f(x) - ∫f(expansion_point).
"""
function integrate(fT::TaylorModelN, which=1)
    p̂ = integrate(fT.pol, which)
    order = fT.pol.order
    r = TaylorN(p̂.coeffs[1:order+1])
    s = TaylorN(p̂.coeffs[order+2:end])
    Δ = bound_integration(fT, s, which)
    return TaylorModelN(r, Δ, fT.x0, fT.dom)
end
function integrate(fT::TaylorModelN, s::Symbol)
    which = TaylorSeries.lookupvar(s)
    return integrate(fT, which)
end

"""
    bound_integration(xTM1::TaylorModel1{Interval{S},S}, δt::Interval{S})
    bound_integration(xTM1::Vector{TaylorModel1{Interval{S},S}}, δt::Interval{S})

Remainder bound for the integration of a series, given by
``δ * remainder(a) +  a.pol[order] * δ^(order+1) / (order+1)``.
This is tighter that the one used by Berz+Makino, which corresponds to
``Δ = aux * remainder(a) +  a.pol[order] * aux^(order+1)``.

"""
@inline function bound_integration(a::TaylorModel1{T,S}, δ) where {T,S}
    order = get_order(a)
    aux = δ^order / (order+1)
    Δ = δ * (remainder(a) + getcoeff(polynomial(a), order) * aux)
    return Δ
end
@inline function bound_integration(a::TaylorModel1{T, S}, δ, δI) where {T<:TaylorN,S}
    order = get_order(a)
    aux = δ^order / (order+1)
    Δ = δ * (remainder(a) + getcoeff(polynomial(a), order)(δI) * aux)
    return Δ
end

function bound_integration(a::Vector{TaylorModel1{T,S}}, δ) where {T,S}
    order = get_order(a[1])
    aux = δ^order / (order+1)
    Δ = δ .* (remainder.(a) .+ getcoeff.(polynomial.(a), order) .* aux)
    return IntervalBox(Δ)
end
function bound_integration(fT::TaylorModelN, s::TaylorN, which)
    Δ = s(fT.dom - fT.x0) + fT.rem * (fT.dom[which] - fT.x0[which])
    return Δ
end


"""
    bound_integration(xTM1::TaylorModel1{Interval{S},S}, δt::Interval{S})
    bound_integration(xTM1::Vector{TaylorModel1{Interval{S},S}}, δt::Interval{S})

Remainder bound for the integration of a series, given by
``δ * remainder(a) +  a.pol[order] * δ^(order+1) / (order+1)``.
This is tighter that the one used by Berz+Makino, which corresponds to
``Δ = aux * remainder(a) +  a.pol[order] * aux^(order+1)``.

"""
@inline function bound_integration(a::RTaylorModel1{T,S}, δ) where {T,S}
    order = get_order(a)
    Δ = δ * remainder(a)
    Δ = Δ/(order+2) + getcoeff(polynomial(a), order)/(order+1)
    return Δ
end
@inline function bound_integration(a::RTaylorModel1{T,S}, δ, δI) where {T<:TaylorN,S}
    order = get_order(a)
    Δ = δ * remainder(a)
    Δ = Δ/(order+2) + getcoeff(polynomial(a), order)(δI)/(order+1)
    return Δ
end


function picard_lindelof(f!, dxTM1TMN::Vector{TaylorModel1{T,S}},
        xTM1TMN::Vector{TaylorModel1{T,S}}, t, params) where {T,S}

    x_picard = similar(xTM1TMN)
    picard_lindelof!(f!, dxTM1TMN, xTM1TMN, t, x_picard, params)
    return x_picard
end

function picard_lindelof!(f!,
        dxTM1TMN::Vector{TaylorModel1{T,S}},
        xTM1TMN::Vector{TaylorModel1{T,S}},
        x_picard::Vector{TaylorModel1{T,S}}, t, params) where {T,S}

    dof = length(xTM1TMN)
    f!(dxTM1TMN, xTM1TMN, params, t)
    @inbounds for ind = 1:dof
        x_picard[ind] = integrate(dxTM1TMN[ind], xTM1TMN[ind][0])
    end
    return x_picard
end
