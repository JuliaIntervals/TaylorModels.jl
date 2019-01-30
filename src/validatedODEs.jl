# Method using TaylorModel1's; necessary to get proper remainder
function TaylorIntegration.jetcoeffs!(eqsdiff!, t::Taylor1{T},
        x::AbstractVector{TaylorModel1{T,S}}, dx::AbstractVector{TaylorModel1{T,S}},
        xaux::AbstractVector{TaylorModel1{T,S}}) where {T,S}

    order = get_order(x[1])
    for ord in 0:order-1
        ordnext = ord+1

        # Set `taux`, auxiliary Taylor1 variable to order `ord`
        @inbounds taux = Taylor1( t.coeffs[1:ordnext] )
        # Set `xaux`, auxiliary vector of TaylorModel1 to order `ord`
        @inbounds for j in eachindex(x)
            xaux[j] = TaylorModel1(
                Taylor1(x[j].pol.coeffs[1:ordnext]), x[j].rem, x[j].x0, x[j].I )
        end

        # Equations of motion
        eqsdiff!(taux, xaux, dx)

        # Recursion relations
        @inbounds for j in eachindex(x)
            x[j].pol[ordnext] = dx[j].pol[ord]/ordnext
        end
    end
    nothing
end


# Stepsize, for Taylor1{TaylorModelN{N,R,T}}
function TaylorIntegration.stepsize(x::Taylor1{TaylorModelN{N,R,T}}, epsilon::T) where {N,R,T<:Real}
    ord = get_order(x)
    h = convert(T, Inf)
    for k in (ord-1, ord)
        @inbounds aux = norm( x[k].pol, Inf)
        # aux == zero(R) && continue
        aux = epsilon / aux
        kinv = one(T)/k
        aux = aux^kinv
        h = min(h, inf(aux))
    end
    return h
end


# # Method for scalar ODEs
# function validate_solution(tT1, xT1Tm, Δtest, t00, t_interval)
#     orderT = get_order(xT1Tm)
#     δt = t_interval-t00
#     aux = (δt^(orderT+1)) / (orderT+1)
#     bx0 = xT1Tm[0].x0
#     δbx0 = xT1Tm[0].iI
#     δI = δbx0 - bx0
#     Δ = zero(bx0)
#     for i = 1:100
#         Δv = Δtest * δt + last_coeff * aux
#         Δ = evaluate(Δv, δI)
#         (Δ ⊂ Δtest) && return Δtest
#         if Δ == Δtest
#             Δtest = interval(prevfloat(Δtest.lo), nextfloat(Δtest.hi))
#             continue
#         end
#         Δtest = hull(Δ, Δtest)
#     end
#     error("Failure to validate existence after 100 steps: $Δ, $Δtest")
# end
# Method for "vectorial" ODEs
function validate_solution(tT1, x, dx, Δtest, δt)
    orderT = get_order(x[1])
    aux = δt^(orderT+1) #/ (orderT+1)
    bx0 = getfield.(getcoeff.(x[1],0), :x0)
    δbx0 = getfield.(getcoeff.(x[1],0), :I)
    δI = δbx0 - bx0
    Δ = zero.(bx0)
    # Δtest = zero.(bx0)
    last_coeff = getcoeff.(dx, orderT)
    vv = Array{eltype(Δtest)}(undef, length(Δtest))
    for i = 1:100
        Δv = Δtest .* δt .+ last_coeff .* aux
        Δ = evaluate(Δv, δI)
        all(Δ .⊂ Δtest) && return Δtest
        # # This is untested
        # if all(Δ .== Δtest)
        #     # for ind in eachindex(Δtest)
        #     #     Δtest[ind] = Interval(prevfloat(Δtest[ind].lo), nextfloat(Δtest[ind].hi))
        #     #     continue
        #     # end
        #     Δtest = IntervalBox( [ Interval(prevfloat(Δ[ind].lo), nextfloat(Δ[ind].hi)) for ind in eachindex(Δ)])
        # end
        for ind in eachindex(vv)
            vv[ind] = Δtest[ind]
            (Δ[ind] ⊂ Δtest[ind]) && continue
            vv[ind] = hull(Δ[ind], Δtest[ind])
            if Δ[ind] == Δtest[ind]
                vv[ind] = Interval(prevfloat(vv[ind].lo), nextfloat(vv[ind].hi))
            end
        end
        Δtest = IntervalBox(vv)
    end
    error("Failure to validate existence after 100 steps: $Δ, $Δtest")
end


function TaylorIntegration.taylorstep!(f!, t::Taylor1{Interval{T}},
        x::Vector{Taylor1{TaylorModelN{N,Interval{T},T}}},
        dx::Vector{Taylor1{TaylorModelN{N,Interval{T},T}}},
        xaux::Vector{Taylor1{TaylorModelN{N,Interval{T},T}}},
        xTM1::Vector{TaylorModel1{Interval{T},T}},
        dxTM1::Vector{TaylorModel1{Interval{T},T}},
        xTM1aux::Vector{TaylorModel1{Interval{T},T}},
        t0::T, t1::T,
        xTMN::Vector{TaylorModelN{N,Interval{T},T}},
        orderT::Int, abstol::T,
        parse_eqs::Bool=true) where {N,T<:Real}

    @assert t1 > t0

    # Some parameters
    bx0  = getfield(getcoeff( x[1],0), :x0)
    δbx0 = getfield(getcoeff(dx[1],0), :I )

    # Compute the Taylor coefficients (non-validated integration);
    # may use specialized method
    TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, x, dx, xaux)

    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(x, abstol)
    δt = min(δt, t1-t0)
    tf = t0 + δt
    time_interval = interval(t0, tf)

    # # Second integration to determine remainder; does not use specialized method
    # TaylorIntegration.jetcoeffs!(f!, time_interval + Taylor1(orderT+1),
    #     xTM1, dxTM1, xTM1aux)
    # # Compute an approx for the reaminder
    # Δtest = bound_integration(xTM1, time_interval)

    # Obtain remainder so solution is properly contained
    # Δ = validate_solution(t, x, dx, Δtest, t0, tf)
    Δ = validate_solution(t, x, dx, zero.(bx0), δt)
    # Validated solution TaylorModel1
    # xTMN .= evaluate( TaylorModel1.(x, Δ, Interval(t0), time_interval), δt )
    xTMN .= evaluate.( x, δt ) .+ Δ#, Interval(t0), time_interval), δt )

    return δt
end

function taylorinteg(f!, q0::IntervalBox{N,T}, δq0::IntervalBox{N,T},
        t0::T, tmax::T, orderQ::Int, orderT::Int, abstol::T;
        maxsteps::Int=500, parse_eqs::Bool=true) where {N, T<:Real}

    # Set orderQ for jet transport
    dof = N
    δq = set_variables("δ", numvars=dof, order=orderQ)
    set_variables("δ", numvars=dof, order=2*orderQ)

    # Allocation for output
    tv = Array{T}(undef, maxsteps+1)
    xv = Array{IntervalBox{N,T}}(undef, maxsteps+1)
    xTMv = Array{TaylorModelN{N,Interval{T},T}}(undef, dof, maxsteps+1)

    # Initialize the vector of Taylor1 expansions
    ti0 = Interval(t0)
    t   = ti0 + Taylor1(orderT)
    x       = Array{Taylor1{TaylorModelN{N,Interval{T},T}}}(undef, dof)
    dx      = Array{Taylor1{TaylorModelN{N,Interval{T},T}}}(undef, dof)
    xaux    = Array{Taylor1{TaylorModelN{N,Interval{T},T}}}(undef, dof)
    xTM1    = Array{TaylorModel1{Interval{T},T}}(undef, dof)
    dxTM1   = Array{TaylorModel1{Interval{T},T}}(undef, dof)
    xTM1aux = Array{TaylorModel1{Interval{T},T}}(undef, dof)
    xTMN    = Array{TaylorModelN{N,Interval{T},T}}(undef, dof)

    # Initial conditions
    δbx0 = q0 + δq0
    zI = zero(Interval{T})
    @inbounds for i in eachindex(x)
        qaux = q0[i] + δq[i]
        x[i] = Taylor1( TaylorModelN(qaux, zI, q0, δbx0), orderT )
        dx[i] = x[i]
        xTM1[i] = TaylorModel1( (qaux)(δq0) + Taylor1(orderT+1), zI, ti0, ti0)
        dxTM1[i] = xTM1[i]
        xTMN[i] = x[i][0]
    end
    @inbounds tv[1] = t0
    @inbounds xv[1] = IntervalBox( evaluate(xTMN, δq0) )
    @inbounds xTMv[:,1] .= xTMN

    # Determine if specialized jetcoeffs! method exists
    parse_eqs = parse_eqs && (length(methods(TaylorIntegration.jetcoeffs!)) > 2)
    if parse_eqs
        try
            TaylorIntegration.jetcoeffs!(Val(f!), t, x, dx)
        catch
            parse_eqs = false
        end
    end

    # Integration
    nsteps = 1
    println(nsteps, "\t", t0, "\t", "\t", diam.(xv[1]))
    while t0 < tmax
        δt = TaylorIntegration.taylorstep!(f!, t, x, dx, xaux,
            xTM1, dxTM1, xTM1aux, t0, tmax, xTMN, orderT, abstol, parse_eqs)
        x0 = IntervalBox( evaluate(xTMN, δq0) )

        # New initial conditions
        t0 += δt
        @inbounds for i in eachindex(x)
            x[i]  = Taylor1( xTMN[i], orderT )
            dx[i] = Taylor1( zero(xTMN[i]), orderT )
        end
        @inbounds t[0] = Interval(t0)
        nsteps += 1
        @inbounds tv[nsteps] = t0
        @inbounds xv[nsteps] = x0
        @inbounds xTMv[:,nsteps] .= xTMN
        # println(nsteps, "\t", t0, "\t", x0)
        println(nsteps, "\t", t0, "\t", diam.(x0))
        if nsteps > maxsteps
            @info("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    return view(tv,1:nsteps), view(xv,1:nsteps)
end
