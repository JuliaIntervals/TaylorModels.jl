# recipe.jl

using RecipesBase

@recipe function g(f::TaylorModel1)
    ffp = fp_rpa(f)
    fT = polynomial(ffp)
    Δ = remainder(ffp)
    ξ0 = ffp.x0

    seriesalpha --> 0.3
    seriestype := :shape

    xs = range(f.dom.lo, stop=f.dom.hi, length=100)
    evals = fT.(xs .- ξ0) .+ Δ

    # make polygon:
    xs = [xs; reverse(xs); xs[1]]
    ys = [inf.(evals); reverse(sup.(evals)); inf(evals[1])]

    xs, ys
end

@recipe function g(f::RTaylorModel1)
    ffp = fp_rpa(f)
    fT = polynomial(ffp)
    Δ = remainder(ffp)
    ξ0 = ffp.x0
    order = get_order(f)+1

    seriesalpha --> 0.5
    seriestype := :shape

    xs = range(f.dom.lo, stop=f.dom.hi, length=100)
    evals = fT.(xs .- ξ0)

    corrs = (xs .- ξ0) .^ order
    Δrel = Δ .* corrs
    evalslo = inf.(evals + Δrel)
    evalshi = sup.(evals + Δrel)

    xs = [xs; reverse(xs); xs[1]]
    ys = [evalslo; reverse(evalshi); evalslo[1]]

    xs, ys
end

@recipe function g(sol::TMSol; vars=(0,1), timediv=1)
    # Some initial checks
    @assert timediv > 0 "`timediv must be 1 or larger"
    @assert length(vars) == length(unique(vars)) == 2 "`vars` must have 2 distinct integer indices"
    @assert 0 ≤ minimum(vars) && maximum(vars) ≤ get_numvars()

    return _plot_intvbox(sol, vars=vars, timediv=timediv)
end

function _plot_intvbox(sol::TMSol; vars=(0,1), timediv=1)
    domT = mince_in_time(sol, var=0, timediv=timediv)

    if 0 ∈ vars
        tup0 = findfirst(0 ∈ vars)
        var1 = vars[3-tup0]
        v1 = _mince_in_time(sol, domT, var1, timediv)
        if tup0 == 1
            return @. IntervalBox(domT, v1)
        else
            return @. IntervalBox(v1, domT)
        end
    end

    var1, var2 = vars
    v1 = _mince_in_time(sol, domT, var1, timediv)
    v2 = _mince_in_time(sol, domT, var2, timediv)
    return @. IntervalBox(v1, v2)
end

function mince_in_time(sol::TMSol; var::Int=0, timediv::Int=0)
    domT = _mince_in_time(sol, Val(true), timediv)
    var == 0 && return domT
    return _mince_in_time(sol, domT, var, timediv)
end

# Mince in time var (var == 0)
function _mince_in_time(sol::TMSol, ::Val{true}, timediv::Int=1)
    # Case timediv == 1
    if timediv == 1
        return @. sol.time + domain(sol.xTM[1,:])
    end

    # Case timediv > 1
    domT = Vector{typeof(domain(sol[1,1]))}(undef, timediv*length(sol))
    i0 = 1
    @inbounds for indT in eachindex(sol)
        i1 = indT*timediv
        domT[i0:i1] .= get_time(sol,indT) .+ mince(domain(sol.xTM[1,indT]), timediv)
        i0 = i1 + 1
    end

    return domT
end

# Mince other var (var > 0)
function _mince_in_time(sol::TMSol, domT::Vector{Interval{T}}, var::Int, timediv::Int=1) where {T}
    N = get_numvars()
    @assert 1 ≤ var ≤ N

    # Case timediv == 1
    if timediv == 1
        return @. getindex(getfield(sol.fp, :v), var)
    end

    # Case timediv > 1
    vv = similar(domT)
    normalized_box = symmetric_box(N, Float64)
    # tTM  = get_time(sol)
    δt = mince(domain(sol[1,1]), timediv)

    i0 = 1
    @inbounds for indT in eachindex(sol)
        i1 = indT*timediv
        δt .= mince(domain(sol[indT,1]), timediv)
        @. vv[i0:i1] = evaluate(evaluate(sol[indT,var], δt), (normalized_box,))
        i0 = i1 + 1
    end

    return vv    
end
