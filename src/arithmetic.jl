

function +(f::TaylorModel, g::TaylorModel)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I

    return TaylorModel(f.n, f.x0, f.I, f.p + g.p, f.Δ + g.Δ, f.bounds)
end

function -(f::TaylorModel, g::TaylorModel)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I

    return TaylorModel(f.n, f.x0, f.I, f.p - g.p, f.Δ - g.Δ, f.bounds)
end

function +(α::Real, f::TaylorModel)
    g = copy(f)

    g.p[0] += (α..α)

    return g
end

-(f::TaylorModel) = (-1 * f)
-(α::Real, f::TaylorModel) = α + (-f)


*(α::Real, f::TaylorModel) = TaylorModel(f.n, f.x0, f.I, Taylor1(α*f.p), (α..α)*f.Δ, f.bounds)


function *(f::TaylorModel, g::TaylorModel)
    @assert f.n == g.n
    @assert f.x0 == g.x0
    @assert f.I == g.I
    @assert f.bounds == g.bounds

    n, x0, I = f.n, f.x0, f.I


    c = zeros(eltype(g.p), 2n+1)

    a = f.p
    b = g.p

    for i in 0:n
        for j in 0:n
            c[i+j+1] += a[i] * b[j]
        end
    end

    d = zeros(eltype(g.p), 2n+1)

    d[n+2:end] = c[n+2:end]

    B = bound(Taylor1(d), x0, I, f.bounds)
    Bf = bound(a, x0, I, f.bounds)
    Bg = bound(b, x0, I, f.bounds)   # assuming g.bounds == f.bounds

    Δ = B + (f.Δ * Bg) + (g.Δ * Bf) + (f.Δ * g.Δ)

    return TaylorModel(n, f.x0, f.I, Taylor1(c[1:n+1]), Δ, f.bounds)

end
