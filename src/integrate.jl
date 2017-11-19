const ∫ = integrate


doc"""
Integrate a TaylorModel.
`x0` is optional constant to add.
"""
function integrate(f::TaylorModel, x0=0)

    p2 = integrate(f.p)

    high_order_term = f.p[end]
    Δ = integral_bound(f)

    t = TaylorModel(f.n, f.x0, f.I, p2, Δ, f.bounds)
    t.p[0] = x0  # constant term

    return t

end

function integral_bound(f::TaylorModel)
    n = degree(f.p)
    high_order_term = f.p[n]

    coeff = bound(high_order_term, f.bounds)
    power = (f.I - f.x0)^n

    ((coeff * power) + f.Δ) * diam(f.I)
end


function Taylor_step(fs, n, u0, v0, t_interval, bounds)

    u = u0
    v = v0
    u_new = u0   # so exist outside loop
    v_new = v0

    # build up Taylor series by Picard:
    for i in 1:n+1   # how many iterations are required?
        u_new = u0 + ∫(fs[1](u, v))
        v_new = v0 + ∫(fs[2](u, v))

        u, v = u_new, v_new
    end

    # prepare Taylor Model:
    uu = TaylorModel(n, Interval(t_interval.lo), t_interval, u, 0..0, bounds)
    vv = TaylorModel(n, Interval(t_interval.lo), t_interval, v, 0..0, bounds)

    uu_new = ∫( fs[1](uu, vv), u0[0] )
    vv_new = ∫( fs[2](uu, vv), v0[0] )

    uΔ = uu_new.Δ
    vΔ = vv_new.Δ

    # make sure the intervals contain 0:
    m = mag(uΔ)
    uΔ = -m..m

    m = mag(vΔ)
    vΔ = -m..m

    @show uΔ, vΔ

     while ! ((uu_new.Δ ⊆ uu.Δ) && (vv_new.Δ ⊆ vv.Δ))
        uΔ *= 2
        vΔ *= 2

        @show uΔ, vΔ

        uu = TaylorModel(n, Interval(t_interval.lo), t_interval, u, uΔ, bounds)
        vv = TaylorModel(n, Interval(t_interval.lo), t_interval, v, vΔ, bounds)

        uu_new = ∫( fs[1](uu, vv), u0[0] )
        vv_new = ∫( fs[2](uu, vv), v0[0] )

    end


    # only need to bound the integral, not actually carry out
    # the whole integral

    # contract:

    for i in 1:10
       uu, vv = uu_new, vv_new
       uu_new = ∫( fs[1](uu, vv), u0[0] )
       vv_new = ∫( fs[2](uu, vv), v0[0] )
    end

    @show uu_new.Δ, vv_new.Δ

    U = uu_new(t_interval.hi)
    V = vv_new(t_interval.hi)

    return U, V, uu_new, vv_new
end
