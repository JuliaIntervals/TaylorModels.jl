
## Taylor coeffcients of different functions
taylor_coeff(::typeof(exp), i, x) = exp(x) / factorial(i)

function taylor_coeff(::typeof(sin), i, x)
    j = i % 4

    (j == 0) && (return sin(x) / factorial(i))
    (j == 1) && (return cos(x) / factorial(i))
    (j == 2) && (return -sin(x) / factorial(i))
    (j == 3) && (return -cos(x) / factorial(i))
end

taylor_coeff(::typeof(inv), i, x) = (-1)^i / (x^(i+1))

doc"""
Make a Taylor1Model for a given function over a given domain.
"""
function Taylor1Model(f, n, x0, I)

    #a = zeros(typeof(I), n+1)
    a = taylor_coeff.(f, 0:n, x0)

    # for i in 0:n
    #     a[i+1] = taylor_coeff(f, i, x0)
    # end

    p = Taylor1(a)

    Γ = taylor_coeff(f, n+1, I)

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

    return Taylor1Model(n, x0, I, p, Δ)
end


import Base.zero
zero(::Type{Taylor1Model{T}}, n, x0, I::Interval{T}) where {T<:AbstractFloat} = Taylor1Model(n, x0, I, Taylor1{Interval{Float64}}(zeros(n+1)), Interval{T}(0))


doc"""
Evaluate a polynomial of a polynomial-like object (Taylor1Model, Taylor1) `f`.
`a` are the coefficients of the polynomial, expressed as a `Taylor1`.
"""
function evaluate_polynomial(a::Taylor1, f)

    M = zero(f)

    # Horner method:
    for i in degree(f):-1:0
        M = (M * f) + a[i]
    end

    return M
end

doc"""
Calculate the Taylor1Model of `(g∘f)` given a function `g` and a Taylor1Model `f`.
"""
function TMcomposition(g, f::Taylor1Model)
    #x0, I, n = f.x0, f.I, f.n
    n = f.n

    f_of_x0 = f[0]

    # calculate interval for g: image of f
    Bf = bound(f)  # Bf + Δf in Joldes thesis
    #a, Δf = f.p, f.Δ

    # location x0 for g is f(x0), which is first coefficient of expansion of f:


    Mg = Taylor1Model(g, n, f_of_x0, Bf)

    b, Δg = Mg.p, Mg.Δ

    #f -= a[0]  #  now is f(x) - f(x0)
    # a[0] = 0  # more efficient


    M1 = f - f_of_x0  # f(x) - f(x0) # Taylor1Model(f, f.p - a[0], Δf)
    M = evaluate_polynomial(b, M1)

    c, Δ = M.p, M.Δ

    return Taylor1Model(f, c, Δ+Δg)
end


# allow exp(t) etc. for a Taylor1Model t:
for g in (:exp, :sin, :inv)
    @eval $g(f::Taylor1Model) = TMcomposition($g, f)
end


function /(f::Taylor1Model, g::Taylor1Model)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I

    M = f * inv(g)
end
