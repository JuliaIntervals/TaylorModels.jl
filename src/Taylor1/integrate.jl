const ∫ = integrate


doc"""
Integrate a Taylor1Model.
`x0` is optional constant to add.
"""
function integrate(f::Taylor1Model, x0=0)

    p2 = integrate(f.p)

    Δ = integral_bound(f)

    g = Taylor1Model(f, [0; p2[1:end]], Δ)

    return g

end


function integral_bound(f::Taylor1Model)
    n = f.n
    high_order_term = f.p[n]

    coeff = bound(high_order_term)
    power = (f.I - f.x0)^n

    return ((coeff * power) + f.Δ) * diam(f.I)
end

"""
fs is an array of functions
This works in 2 variables currently.
"""
function Taylor_step(fs, n, u0, v0, t_interval)

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
    uu = Taylor1Model(n, Interval(t_interval.lo), t_interval, u, 0..0)
    vv = Taylor1Model(n, Interval(t_interval.lo), t_interval, v, 0..0)

    uΔ = integral_bound(fs[1](uu, vv))
    vΔ = integral_bound(fs[2](uu, vv))

    # make sure the intervals contain 0:
    m = mag(uΔ)
    uΔ = -m..m

    m = mag(vΔ)
    vΔ = -m..m

    @show uΔ, vΔ

    uu = Taylor1Model(n, Interval(t_interval.lo), t_interval, u, uΔ)
    vv = Taylor1Model(n, Interval(t_interval.lo), t_interval, v, vΔ)

    integ_u_bound = uΔ
    integ_v_bound = vΔ

    while ! ((integ_u_bound ⊆ uu.Δ) && (integ_v_bound ⊆ vv.Δ))
        uΔ *= 2
        vΔ *= 2

        @show uΔ, vΔ

        uu = Taylor1Model(n, Interval(t_interval.lo), t_interval, u, uΔ)
        vv = Taylor1Model(n, Interval(t_interval.lo), t_interval, v, vΔ)

        integ_u_bound = integral_bound(fs[1](uu, vv))
        integ_v_bound = integral_bound(fs[2](uu, vv))

    end


    # only need to bound the integral, not actually carry out
    # the whole integral

    # contract:

    uu_new = ∫( fs[1](uu, vv), u0[0] )
    vv_new = ∫( fs[2](uu, vv), v0[0] )

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
