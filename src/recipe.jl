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

@recipe function g(sol::TMSol; vars=(0,1), ntdiv=1)
    # Some initial checks
    @assert ntdiv > 0 "`ntdiv must be 1 or larger"
    @assert length(vars) == length(unique(vars)) == 2 "`vars` must have 2 distinct integer indices"
    dof = get_numvars()
    @assert 0 ≤ minimum(vars) && maximum(vars) ≤ dof

    return _plot_intvbox(sol, vars=vars, ntdiv=ntdiv)
end

function _plot_intvbox(sol; vars=(0,1), ntdiv=1)
    # Initializations
    normalized_box = symmetric_box(get_numvars(), Float64)
    tTM  = getfield(sol, :time)
    xTMv = getfield(sol, :xTMv)
    plotbox = Vector{IntervalBox{2, Float64}}(undef, ntdiv*(ltime-1))
    v1 = Vector{Interval{Float64}}(undef, length(plotbox))
    v2 = similar(v1)

    # Case in which time is plotted
    if 0 ∈ vars
        tupT = findfirst(0 ∈ vars)
        tup2 = vars[3-tupT]
        ii = 1
        @inbounds for indT in eachindex(tTM)
            domTv = mince(domain(xTMv[tup2,indT]), ntdiv)
            for domT in domTv
                v1[ii] = tTM[indT]+domT
                v2[ii] = evaluate(xTMv[tup2,indT](domT), normalized_box)
                ii += 1
            end
        end

        if vars[1] == 0
            @. plotbox = IntervalBox(v1, v2)
        else
            @. plotbox = IntervalBox(v2, v1)
        end
    else

        # Case in which time is not plotted
        var1 = vars[1]
        var2 = vars[2]
        ii = 1
        @inbounds for indT in eachindex(tTM)
            domTv = mince(domain(xTMv[var1,indT]), ntdiv)
            for domT in domTv
                v1[ii] = evaluate(xTMv[var1,indT](domT), normalized_box)
                v2[ii] = evaluate(xTMv[var2,indT](domT), normalized_box)
                ii += 1
            end
        end
        @. plotbox = IntervalBox(v1, v2)
    end

    return plotbox
end
