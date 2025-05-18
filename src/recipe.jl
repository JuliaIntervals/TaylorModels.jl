# recipe.jl

using RecipesBase

@recipe function g(f::TaylorModel1)
    ffp = f#fp_rpa(f)
    fT = polynomial(ffp)
    Δ = remainder(ffp)
    ξ0 = expansion_point(ffp)

    seriesalpha --> 0.3
    seriestype := :shape

    xs = range(inf(domain(f)), stop=sup(domain(f)), length=100)
    evals = fT.(xs .- ξ0) .+ Δ

    # make polygon:
    xs = [xs; reverse(xs); xs[1]]
    ys = [inf.(evals); reverse(sup.(evals)); inf(evals[1])]

    return (xs, ys)
end

@recipe function g(f::RTaylorModel1)
    ffp = f#fp_rpa(f)
    fT = polynomial(ffp)
    Δ = remainder(ffp)
    ξ0 = expansion_point(ffp)
    order = get_order(f)+1

    seriesalpha --> 0.5
    seriestype := :shape

    xs = range(inf(domain(f)), stop=sup(domain(f)), length=100)
    evals = fT.(xs .- ξ0)

    corrs = (xs .- ξ0) .^ order
    Δrel = Δ .* corrs
    evalslo = inf.(evals .+ Δrel)
    evalshi = sup.(evals .+ Δrel)

    xs = [xs; reverse(xs); xs[1]]
    ys = [evalslo; reverse(evalshi); evalslo[1]]

    return (xs, ys)
end

@recipe function g(sol::TMSol; vars=(0,1), timediv=1)
    @assert length(vars) == length(unique(vars)) == 2 "`vars` must have 2 distinct integer indices"

    return _plot_intvbox(sol, vars=vars, timediv=timediv)
end

function _plot_intvbox(sol::TMSol; vars=(0,1), timediv=1)
    domT = mince_in_time(sol, var=0, timediv=timediv)

    if 0 ∈ vars
        tup0 = findfirst(0 ∈ vars)
        var1 = vars[3-tup0]
        v1 = _mince_in_time(sol, domT, var1, timediv)
        if tup0 == 1
            return SVector.(domT, v1)
        else
            return SVector.(v1, domT)
        end
    end

    var1, var2 = vars
    v1 = _mince_in_time(sol, domT, var1, timediv)
    v2 = _mince_in_time(sol, domT, var2, timediv)
    return SVector.(v1, v2)
end

"""
    mince_in_time(sol::TMSol; var=0, timediv=1) --> ::Vector{Interval}

For `var=0`, this function divides the time domain of each entry of `sol` in
`timediv` parts (`timediv==1` is the initial domain), and returns the time
intervals where the solution holds. This is useful for plotting or finding
specific events.
For `var ≥ 1`, this function evaluates the flowpipes at the split domain
intervals, which is useful to decrease the overapproximations associated
to the whole time domain.
"""
function mince_in_time(sol::TMSol; var::Int=0, timediv::Int=1)
    @assert timediv > 0 "`timediv must be 1 or larger"
    @assert 0 ≤ var ≤ get_numvars(sol)

    domT = _mince_in_time(sol, Val(true), timediv)
    var == 0 && return domT
    return _mince_in_time(sol, domT, var, timediv)
end

# Mince in time var (var == 0)
function _mince_in_time(sol::TMSol, ::Val{true}, timediv::Int=1)
    # Case timediv == 1
    if timediv == 1
        return expansion_point(sol) .+ domain(sol)
    end

    # Case timediv > 1
    domT = Vector{typeof(domain(sol,1))}(undef, timediv*length(sol))
    i0 = 1
    @inbounds for indT in eachindex(sol)
        i1 = indT*timediv
        domT[i0:i1] .= expansion_point(sol,indT) .+ mince(domain(sol,indT), timediv)
        i0 = i1 + 1
    end

    return domT
end

# Mince other var (var > 0)
function _mince_in_time(sol::TMSol, domT::Vector{Interval{T}}, var::Int,
        timediv::Int=1) where {T}
    N = get_numvars(sol)
    @assert 1 ≤ var ≤ N

    # Case timediv == 1
    if timediv == 1
        return getindex.(flowpipe(sol), var)
    end

    # Case timediv > 1
    vv = similar(domT)
    normalized_box = fill(-1..1, SVector{N})
    δt = mince(domain(sol,1), timediv)

    i0 = 1
    @inbounds for indT in eachindex(sol)
        i1 = indT*timediv
        δt .= mince(domain(sol,indT), timediv)
        @. vv[i0:i1] = evaluate(evaluate(sol[indT,var], δt), (normalized_box,))
        i0 = i1 + 1
    end

    return vv
end
