# integration.jl

for TM in tupleTMs
    @eval begin
        function integrate(a::$(TM){T,S}) where {T<:TS.NumberNotSeries,S}
            integ_pol = integrate(a.pol)
            Δ = bound_integration(a, centered_dom(a))
            return $(TM)( integ_pol, Δ, expansion_point(a), domain(a) )
        end
        function integrate(a::$(TM){TaylorN{T},S}, cc0) where {T,S}
            integ_pol = integrate(a.pol)
            Δ = bound_integration(a, centered_dom(a), cc0)
            return $(TM)(integ_pol, Δ, expansion_point(a), domain(a))
        end
        integrate(a::$(TM){T,S}, c0) where {T,S} = c0 + integrate(a)
        integrate(a::$(TM){TaylorN{T},S}, c0, δI) where {T,S} = c0 + integrate(a, δI)


        @inline function bound_integration(a::$(TM){T,S}, δ::Interval{S}) where {T<:TS.NumberNotSeries,S}
            order = get_order(a)
            if $TM == TaylorModel1
                aux = pown(δ, order) / interval(order+1)
                Δ = δ * (remainder(a) + interval(getcoeff(polynomial(a), order)) * aux)
            else
                Δ = δ * remainder(a)
                Δ = Δ/interval(order+2) + interval(getcoeff(polynomial(a), order))/interval(order+1)
            end
            return Δ
        end
        @inline function bound_integration(a::$(TM){TaylorN{T}, S}, δ::Interval{S}, δI) where {T,S}
            order = get_order(a)
            if $TM == TaylorModel1
                aux = pown(δ, order) / interval(order+1)
                Δ = δ * (remainder(a) + getcoeff(polynomial(a), order)(δI) * aux)
            else
                Δ = δ * remainder(a)
                Δ = Δ/interval(order+2) + getcoeff(polynomial(a), order)(δI)/interval(order+1)
            end
            return Δ
        end
    end
end


function integrate(a::TaylorModel1{TaylorModelN{N,T,S},S},
        c0::TaylorModelN{N,T,S}) where {N,T,S}
    integ_pol = integrate(a.pol, c0)
    δ = centered_dom(a)

    # Remainder bound after integrating
    Δ = bound_integration(a, δ)
    ΔN = Δ(centered_dom(a[0]))

    return TaylorModel1( integ_pol, ΔN, expansion_point(a), domain(a) )
end

function integrate(fT::TaylorModelN, which=1)
    p̂ = integrate(fT.pol, which)
    order = get_order(fT)
    r = TaylorN(p̂.coeffs[1:order+1])
    s = TaylorN(p̂.coeffs[order+2:end])
    Δ = bound_integration(fT, s, which)
    return TaylorModelN(r, Δ, expansion_point(fT), domain(fT))
end
function integrate(fT::TaylorModelN, s::Symbol)
    which = TaylorSeries.lookupvar(s)
    return integrate(fT, which)
end



@inline function bound_integration(a::Vector{TaylorModel1{T,S}}, δ) where {T,S}
    order = get_order(a[1])
    aux = pown(δ, order) / interval(order+1)
    Δ = δ .* (remainder.(a) .+ interval.(getcoeff.(polynomial.(a), order)) .* aux)
    return Δ
end
@inline function bound_integration(fT::TaylorModelN, s::TaylorN, which)
    Δ = s(centered_dom(fT)) + remainder(fT) * centered_dom(fT)[which]
    return Δ
end


@doc """
    integrate(a::TM{T,S}, c0)
    integrate(a::TM{TaylorN{T},S}, c0, cc0)

Integrates the one-variable Taylor Model (`TaylorModel1` or `RTaylorModel1`) with
respect to the independent variable. `c0` is the integration constant; if omitted
it is taken as zero. When the coefficients of `a` are `TaylorN` variables,
the domain is specified by `cc0::AbstractVector{<:Interval}`.

---

    integrate(fT, which)

Integrates a `fT::TaylorModelN` with respect to `which` variable.
The returned `TaylorModelN` corresponds to the Taylor Model
of the definite integral ∫f(x) - ∫f(expansion_point).
""" integrate


@doc """
    bound_integration(xTM::TaylorModel1{T,S}, δ)
    bound_integration(xTM::Vector{TaylorModel1{T,S}}, δ)

Bound the remainder of the integration of a xTM::TaylorModel1, where δ is the domain used
to bound the integration. The remainder corresponds to
``δ * remainder(a) +  a.pol[order] * δ^(order+1) / (order+1)``.
This is tighter that the one used by Berz+Makino, which corresponds to
``Δ = δ * remainder(a) +  a.pol[order] * δ^(order+1)``.

---
    bound_integration(xTM::RTaylorModel1{T,S}, δ)
    bound_integration(xTM::Vector{RTaylorModel1{T,S}, δ)

Remainder bound for the integration of a xTM::RTaylorModel1, where δ is the domain used
to bound the integration. The remainder corresponds to
``Δ = δ * remainder(a)/(order+2) + getcoeff(polynomial(a), order)/(order+1)``.

""" bound_integration


function picard_lindelof(f!, dxTM1TMN::Vector{TaylorModel1{T,S}},
        xTM1TMN::Vector{TaylorModel1{T,S}}, t, params) where {T,S}
    x_picard = similar(xTM1TMN)
    picard_lindelof!(f!, dxTM1TMN, xTM1TMN, t, x_picard, params)
    return x_picard
end

function picard_lindelof!(f!,
        x_picard::Vector{TaylorModel1{T,S}},
        dxTM1TMN::Vector{TaylorModel1{T,S}},
        xTM1TMN ::Vector{TaylorModel1{T,S}}, t, params) where {T,S}
    f!(dxTM1TMN, xTM1TMN, params, t)
    @inbounds for ind in eachindex(xTM1TMN)
        x_picard[ind] = integrate(dxTM1TMN[ind], xTM1TMN[ind][0])
    end
    return x_picard
end
