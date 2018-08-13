
data(f::TaylorNModel) = (f.n, f.x0, f.I)

function +(f::TaylorNModel, g::TaylorNModel)
    @assert data(f) == data(g)

    return TaylorNModel(f.n, f.x0, f.I, f.p + g.p, f.Δ + g.Δ, f.order_bounds
     + g.order_bounds)
end

function -(f::TaylorNModel, g::TaylorNModel)
    @assert data(f) == data(g)

    return TaylorNModel(f.n, f.x0, f.I, f.p - g.p, f.Δ - g.Δ, f.order_bounds
     - g.order_bounds)
end

function +(α::Real, f::TaylorNModel)
    g = copy(f)

    g.p.coeffs[1].coeffs[1] += (α..α)

    return g
end

-(f::TaylorNModel) = (-1 * f)
-(α::Real, f::TaylorNModel) = α + (-f)


*(α::Real, f::TaylorNModel) = TaylorNModel(f.n, f.x0, f.I, α*f.p, (α..α)*f.Δ,
α .* f.order_bounds)

*(α::Interval, f::TaylorNModel) = TaylorNModel(f.n, f.x0, f.I, α*f.p, α*f.Δ,
α .* f.order_bounds)


"""
Truncate a TaylorNModel to order n
"""
function truncate(f::TaylorNModel, n)
    if n >= f.n  # nothing to do
        return f
    end

    coeffs = f.p.coeffs[1:n+1]
    order_bounds = f.order_bounds[1:n+1]

    # add in order_bounds to remainder:
    Δ = f.Δ + sum(f.order_bounds[n+2:end])

    TaylorNModel(n, f.x0, f.I, coeffs, Δ, order_bounds)
end


function *(f::TaylorNModel, g::TaylorNModel)
    @assert f.n == g.n  # not necessary
    @assert f.x0 == g.x0
    @assert f.I == g.I

    n, x0, I = f.n, f.x0, f.I

    # TODO

    product_series = f.p * g.p

    Δ = f.Δ * g.Δ

    for i in 0:n
        for j in (n-i+1):n  # j + k > n
            Δ += f.order_bounds[i] * g.order_bounds[j]
        end
    end

    Δ += (sum(f.order_bounds) * g.Δ)
    Δ += (sum(g.order_bounds) * f.Δ)

    return TaylorNModel(n, x0, I, product_series, Δ)
end
