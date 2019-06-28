# recipe.jl

using RecipesBase

@recipe function g(f::TaylorModel1)
    ffp = fp_rpa(f)
    fT = polynomial(ffp)
    Δ = remainder(ffp)
    ξ0 = ffp.x0

    alpha --> 0.3
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

    alpha --> 0.5
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
