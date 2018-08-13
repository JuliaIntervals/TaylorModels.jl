data(f::Taylor1Model) = (f.n, f.x0, f.I)

function +(f::Taylor1Model, g::Taylor1Model)
    @assert data(f) == data(g)

    return Taylor1Model(f, f.p + g.p, f.Δ + g.Δ)
end

function -(f::Taylor1Model, g::Taylor1Model)
    @assert data(f) == data(g)

    return Taylor1Model(f, f.p - g.p, f.Δ - g.Δ)
end

+(α::Real, f::Taylor1Model) = Taylor1Model(f, f.p + α, f.Δ)
+(f::Taylor1Model, α::Real)= α + f

+(f::Taylor1Model, α::Real)= α + f

-(f::Taylor1Model) = (-1 * f)
-(α::Real, f::Taylor1Model) = α + (-f)
-(f::Taylor1Model, α::Real) = f + (-α)

*(α::Real, f::Taylor1Model) = Taylor1Model(f, α*f.p, α*f.Δ)
*(f::Taylor1Model, α::Real) = α * f

function *(f::Taylor1Model, g::Taylor1Model)
    @assert data(f) == data(g)

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

    B = bound(Taylor1(d), x0, I)
    Bf = bound(a, x0, I)
    Bg = bound(b, x0, I)

    Δ = B + (f.Δ * Bg) + (g.Δ * Bf) + (f.Δ * g.Δ)

    return Taylor1Model(f, Taylor1(c[1:n+1]), Δ)

end
