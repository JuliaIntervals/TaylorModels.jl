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

    # Remainder bound after integrating.
    Δ = bound_integration(a, δ)

    return TaylorModel1( integ_pol, Δ, a.x0, a.dom )
end
integrate(a::TaylorModel1{T,S}) where {T,S} = integrate(a, zero(T))

function integrate(a::RTaylorModel1{T,S}, c0::T) where {T,S}
    order = get_order(a)
    integ_pol = integrate(a.pol, c0)

    # Remainder bound after integrating...
    Δ = (a.dom-a.x0) * remainder(a)
    Δ = Δ/(order+2) + a.pol[order]/(order+1)

    return RTaylorModel1( integ_pol, Δ, a.x0, a.dom )
end
integrate(a::RTaylorModel1{T,S}) where {T,S} = integrate(a, zero(T))


"""
    bound_integration(xTM1::TaylorModel1{Interval{S},S}, δt::Interval{S})
    bound_integration(xTM1::Vector{TaylorModel1{Interval{S},S}}, δt::Interval{S})

Remainder bound for the integration of a series, given by
``δ * remainder(a) +  a.pol[order] * δ^(order+1) / (order+1)``.
This is tighter that the one used by Berz+Makino, which corresponds to
``Δ = aux * remainder(a) +  a.pol[order] * aux^(order+1)``.

"""
function bound_integration(a::TaylorModel1{T,S}, δ) where {T,S}
    order = get_order(a)
    aux = δ^order / (order+1)
    Δ = δ * (remainder(a) + getcoeff(polynomial(a),order) * aux)
    return Δ
end
function bound_integration(a::Vector{TaylorModel1{T,S}}, δ) where {T,S}
    order = get_order(a[1])
    aux = δ^order / (order+1)
    Δ = δ .* (remainder.(a) .+ getcoeff.(polynomial.(a), order) .* aux)
    return IntervalBox(Δ)
end



"""
    picard_lindelöf(f, tm::T, xm::T, x0::Interval)
    𝒫(f, tm::T, xm::T, x0::Interval)

Returns the application of the Picard-Lindelöf operator
associated to the ODE ``\\dot{x} = f(t,x)``,
with initial condition `x0`. Here, `tm` and `xm` are
(one-variable) Taylor Models (`TaylorModel1` or `RTaylorModel1`).

𝒫 is an abbreviation of this operator, which is obtained
as `\\mscrP<TAB>`.)
"""
picard_lindelöf(f, tm::T, xm::T, x0::Interval) where
    {T<:Union{TaylorModel1, RTaylorModel1}} = integrate(f(tm, xm), x0)

const 𝒫 = picard_lindelöf

"""
    check_existence(f, tm::T, xm::T, x0::Interval, x_test::Interval)

Checks that the range of one iterate of the Picard-Lindelöf
operator is contained in the a-priori interval `x_test` (of
the dependent variable) that bounds the solution of
the ODE defined by `f`. This function returns an interval
of the independent variable where the a-priori solution
is warranted to exist; see [`shrink_for_existance`](@ref).
Here, `tm` and `xm` are Taylor Models (`TaylorModel1` or `RTaylorModel1`)
defined for the independent and dependent variables, and `x0`
is the initial condition.
"""
function check_existence(f, tm::T, xm::T, x0::Interval, x_test::Interval,
        max_steps::Integer=20) where {T<:Union{TaylorModel1, RTaylorModel1}}

    pl = 𝒫(f, tm, xm, x0)
    tt = shrink_for_existance(pl, tm.dom, x_test, max_steps)
    if pl(tt-tm.x0) ⊆ x_test
        return tt
    else
        return emptyinterval(tt)
    end
end

"""
    shrink_for_existance(xm::T, t_interval, x_test, max_steps::Integer=20)

Shrinks the a-priori independent-variable interval `t_interval`
so the range of `xm`, the Taylor Model (`TaylorModel1` or `RTaylorModel1`)
associated with the dependent variable, is contained in the
a-priori interval `x_test`. The method used is some sort
of bisection. If no independent-variable interval is found within
`max_steps`, an empty interval is returned.
"""
function shrink_for_existance(xm::T, t_interval, x_test, max_steps::Integer=20) where
        {T<:Union{TaylorModel1, RTaylorModel1}}

    tt = t_interval
    t1 = t_interval

    nsteps = 0
    xm_pol = xm.pol

    while nsteps < max_steps
        if xm(t1) ⊆ x_test
            (t_interval.hi-t1.hi < 1.0e-3) && return t1
            tt = t1
            mm = (tt.hi+t_interval.hi)*0.5
            t1 = interval(t1.lo, mm)
        else
            t_interval = t1
            if t1 == tt
                mm = (tt.lo+t1.hi)*0.5
                t1 = interval(t1.lo, mm)
                tt = t1
            else
                mm = (tt.hi+t1.hi)*0.5
                t1 = interval(t1.lo, mm)
            end
        end
        nsteps += 1
    end

    xm(t1) ⊆ x_test && return t1
    xm(tt) ⊆ x_test && return tt

    return emptyinterval(tt)
end

"""
    tight_remainder(f, tm::T, xm::T, x0::Interval, max_steps::Integer=20)

Returns a Taylor Model for the dependent variable, with a tighter
remainder, which is obtained after successive iteration
of the Picard-Lindelöf. If the remainder is not tighter
(and identity with the former iterate is not obtained) a
Taylor Model with an empty interval is returned.

"""
function tight_remainder(f, tm::T, xm::T, x0::Interval, max_steps::Integer=20) where
        {T<:Union{TaylorModel1, RTaylorModel1}}

    xOld = deepcopy(xm)
    for ind = 1:max_steps
        xNew = 𝒫(f, tm, xOld, x0)
        if diam(remainder(xNew)) ≥ diam(remainder(xOld))
            xOld == xNew && return xOld
            return T(xOld.pol, emptyinterval(xOld.rem), xOld.x0, xOld.dom)
        end
        xOld = xNew
    end
    return xOld
end
