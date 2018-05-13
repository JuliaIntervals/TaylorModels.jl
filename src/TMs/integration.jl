# integration.jl

# integrate
function integrate(a::TM1AbsRem{T}, c0::Interval{T}) where {T}
    order = get_order(a)
    integ_pol = integrate(a.pol, c0)
    aux = (a.iI-a.x0)

    # Remainder bound after integrating. This is tighter
    # that the one used by Berz+Makino, which corresponds to:
    # Δ = aux * remainder(a) +  a.pol[order] * aux^(order+1)
    Δ = aux * remainder(a) +  a.pol[order] * aux^(order+1) / (order+1)

    return TM1AbsRem( integ_pol, Δ, a.x0, a.iI )
end
integrate(a::TM1AbsRem{T}) where {T} = integrate(a, Interval(zero(T)))

function integrate(a::TM1RelRem{T}, c0::Interval{T}) where {T}
    order = get_order(a)
    integ_pol = integrate(a.pol, c0)
    Δ = (a.iI-a.x0) * remainder(a)

    # Remainder bound after integrating...
    Δ = Δ/(order+2) + a.pol[order]/(order+1)

    return TM1RelRem( integ_pol, Δ, a.x0, a.iI )
end
integrate(a::TM1RelRem{T}) where {T} = integrate(a, Interval(zero(T)))


# \mscrP --> 𝒫
picard_lindelöf(f, tm::T, xm::T, x0::Interval) where
    {T<:Union{TM1AbsRem, TM1RelRem}} = integrate(f(tm, xm), x0)
const 𝒫 = picard_lindelöf

function check_existence(f, tm::T, xm::T, x0::Interval, x_test::Interval,
        max_steps::Int=20) where {T<:Union{TM1AbsRem, TM1RelRem}}

    pl = 𝒫(f, tm, xm, x0)
    tt = shrink_for_existance(pl, tm.iI, x_test, max_steps)
    if pl(tt-tm.x0) ⊆ x_test
        return tt
    else
        return emptyinterval(tt)
    end
end

function shrink_for_existance(xm::T, t_interval, x_test, max_steps::Int=20) where
        {T<:Union{TM1AbsRem, TM1RelRem}}

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

function tight_remainder(f, tm::T, xm::T, x0::Interval, max_steps::Int=20) where
        {T<:Union{TM1AbsRem, TM1RelRem}}

    xOld = deepcopy(xm)
    for ind = 1:max_steps
        xNew = 𝒫(f, tm, xOld, x0)
        if diam(remainder(xNew)) ≥ diam(remainder(xOld))
            xOld == xNew && return xOld
            return T(xOld.pol, emptyinterval(xOld.rem), xOld.x0, xOld.iI)
        end
        xOld = xNew
    end
    return xOld
end
