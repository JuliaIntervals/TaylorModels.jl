
## Taylor coeffcients of different functions
## TODO: Change nth_deriv -> taylor_coeff
nth_deriv(::typeof(exp), i, x) = exp(x) / factorial(i)

function nth_deriv(::typeof(sin), i, x)
    j = i % 4

    (j == 0) && (return sin(x) / factorial(i))
    (j == 1) && (return cos(x) / factorial(i))
    (j == 2) && (return -sin(x) / factorial(i))
    (j == 3) && (return -cos(x) / factorial(i))
end

nth_deriv(::typeof(inv), i, x) = (-1)^i / (x^(i+1))

doc"""
Make a TaylorModel for a given function over a given domain.
"""
function make_Taylor_model(f, n, x0, I::Interval{T}, bounds) where T

    a = zeros(typeof(I), n+1)

    for i in 0:n
        a[i+1] = nth_deriv(f, i, x0)
    end

    p = Taylor1(a)

    Γ = nth_deriv(f, n+1, I)

    if sup(Γ) ≤ 0 || inf(Γ) ≥ 0

        lo = Interval(inf(I))
        hi = Interval(sup(I))

        Δlo = f(lo) - bound(p, x0, lo)
        Δhi = f(hi) - bound(p, x0, hi)
        Δx0 = f(x0) - bound(p, x0, x0)
        Δ = Interval(min(inf(Δlo), inf(Δx0), inf(Δhi)), max(sup(Δlo), sup(Δx0), sup(Δhi)))

    else
        V = (I - x0)^(n+1)
        Δ = V * Γ
    end

    return TaylorModel(n, x0, I, p, Δ, bounds)
end


import Base.zero
zero(::Type{TaylorModel{T}}, n, x0, I::Interval{T}, bounds) where {T<:AbstractFloat} = TaylorModel(n, x0, I, Taylor1{Interval{Float64}}(zeros(n+1)), Interval{T}(0), bounds)


doc"""
Evaluate a polynomial of a TaylorModel.
`b` are the coefficients of the Taylor1nomial.
"""
function poly_eval_of_TM(b, f, I::Interval{T}, x0, n, bounds) where {T<:AbstractFloat}

    #M = TaylorModel{T}(n, x0, I, Taylor1(zeros(n+1)), Interval{T}(0))
    M = zero(TaylorModel{T}, n, x0, I, bounds)

    for i in n:-1:0
        M *= f
        M.p[0] += b[i]  # add constant
    end

    return M
end

doc"""
Calculate the TaylorModel of `(g∘f)` given a function `g` and a TaylorModel `f`.
"""
function TMcomposition(g, f::TaylorModel)
    x0, I, n = f.x0, f.I, f.n

    Bf = bound(f.p, x0, I)
    a, Δf = f.p, f.Δ

    Mg = make_Taylor_model(g, n, a[0], Bf + Δf, f.bounds)

    b, Δg = Mg.p, Mg.Δ

    a[0] = 0  # zero out first element

    M1 = TaylorModel(n, x0, I, a, Δf, f.bounds)
    M = poly_eval_of_TM(b, M1, I, x0, n, f.bounds)

    c, Δ = M.p, M.Δ

    return TaylorModel(n, x0, I, c, Δ+Δg, f.bounds)
end


# allow exp(t) etc. for a TaylorModel t:
for g in (:exp, :sin, :inv)
    @eval $g(f::TaylorModel) = TMcomposition($g, f)
end


function /(f::TaylorModel, g::TaylorModel)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I

    M = f * inv(g)
end
