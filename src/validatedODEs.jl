# Some methods for validated integration of ODEs

# Stepsize, for Taylor1{TaylorModelN{N,R,T}}
function TaylorIntegration.stepsize(x::Taylor1{TaylorModelN{N,R,T}}, epsilon::T) where {N,R,T<:Real}
    ord = get_order(x)
    h = convert(T, Inf)
    for k in (ord-1, ord)
        @inbounds aux = norm( x[k].pol, Inf)
        (0 ∈ aux) && continue
        aux = epsilon / aux
        kinv = one(T)/k
        aux = aux^kinv
        h = min(h, inf(aux))
    end
    return h
end


"""
    remainder_taylorstep(dx, δI, δt)

Returns a remainder for the integration step using Taylor Models, exploiting
repeated iterations of Schauder's fix point theorem (up to 100) of Picard's
integral operator, considering only the remainder. Inputs are: `dx`, the
vector with the RHS of the defining ODEs, `δI` the interval box where the
initial conditions are varied, and `δt` is the integration interval.
"""
function remainder_taylorstep(dx::Vector{Taylor1{TaylorModelN{N,R,T}}},
        δI::IntervalBox{N,T}, δt::T) where {N,R,T}

    orderT = get_order(dx[1])
    aux = δt^orderT / (orderT+1)
    last_coeffTM = getcoeff.(dx, orderT)
    last_coeff_I = evaluate(last_coeffTM, δI) .* aux
    vv = Array{R}(undef, N)
    Δtest = zero.(δI)

    # This mimics the Schauder's fix point theorem (100 iterations)
    for i = 1:100
        # Only the remainders (Picard's integrationto order orderT) are included
        Δ = δt .* (Δtest .+ last_coeff_I )

        # This enlarges Δtest in all directions, and returns
        if Δ == Δtest
            @inbounds for ind in eachindex(vv)
                vv[ind] = Interval(prevfloat(vv[ind].lo), nextfloat(vv[ind].hi))
            end
            Δtest = IntervalBox(vv)
            return Δtest
        end

        # If needed, the tested remainder is enlarged
        @inbounds for ind in eachindex(vv)
            vv[ind] = Δ[ind]
            (Δ[ind] ⊆ Δtest[ind]) && continue
            vv[ind] = hull(Δ[ind], Δtest[ind])
        end
        Δtest = IntervalBox(vv)
    end

    # NOTE: Return after 100 iterations; this should be changed
    # to ensure convergence. Perhaps shrink δt ?
    return Δtest
end


function TaylorIntegration.taylorstep!(f!, t::Taylor1{R},
        x::Vector{Taylor1{TaylorModelN{N,R,T}}},
        dx::Vector{Taylor1{TaylorModelN{N,R,T}}},
        xaux::Vector{Taylor1{TaylorModelN{N,R,T}}},
        t0::T, t1::T,
        orderT::Int, abstol::T,
        parse_eqs::Bool=true) where {N,R,T<:Real}

    @assert t1 > t0

    # Compute the Taylor coefficients (non-validated integration)
    TaylorIntegration.__jetcoeffs!(Val(parse_eqs), f!, t, x, dx, xaux)

    # Compute the step-size of the integration using `abstol`
    δt = TaylorIntegration.stepsize(x, abstol)
    δt = min(δt, t1-t0)

    return δt
end

function taylorinteg(f!, q0::IntervalBox{N,T}, δq0::IntervalBox{N,T},
        t0::T, tmax::T, orderQ::Int, orderT::Int, abstol::T;
        maxsteps::Int=500, parse_eqs::Bool=true) where {N, T<:Real}

    # Set proper parameters for jet transport
    dof = N
    if (get_numvars() != dof || get_order() != 2orderQ)
        set_variables("δ", numvars=dof, order=2*orderQ)
    end

    # Some variables
    R   = Interval{T}
    ti0 = Interval(t0)
    t   = ti0 + Taylor1(orderT)

    # Allocation of vectors
    tv    = Array{T}(undef, maxsteps+1)
    xv    = Array{IntervalBox{N,T}}(undef, maxsteps+1)
    xTMNv = Array{TaylorModelN{N,R,T}}(undef, dof, maxsteps+1)
    x     = Array{Taylor1{TaylorModelN{N,R,T}}}(undef, dof)
    dx    = Array{Taylor1{TaylorModelN{N,R,T}}}(undef, dof)
    xaux  = Array{Taylor1{TaylorModelN{N,R,T}}}(undef, dof)
    xTMN  = Array{TaylorModelN{N,R,T}}(undef, dof)

    # Set initial conditions
    zI = zero(R)
    @inbounds for i in eachindex(x)
        qaux = q0[i] + TaylorN(i, order=orderQ)
        x[i] = Taylor1( TaylorModelN(qaux, zI, q0, q0 + δq0), orderT )
        dx[i] = x[i]
        xTMN[i] = x[i][0]
    end

    # Output vectors
    @inbounds tv[1] = t0
    @inbounds xv[1] = IntervalBox( evaluate(xTMN, δq0) )
    @inbounds xTMNv[:,1] .= xTMN

    # Determine if specialized jetcoeffs! method exists (built by @taylorize)
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
    while t0 < tmax
        # One step integration (non-validated)
        δt = TaylorIntegration.taylorstep!(f!, t, x, dx, xaux,
            t0, tmax, orderT, abstol, parse_eqs)

        # Validate the solution: build a tight remainder (based on Schauder thm)
        # This is to compute dx[:][orderT] (now zero), needed for the remainder
        f!(t, x, dx)
        Δ = remainder_taylorstep(dx, δq0, δt) # remainder of integration step

        # Evaluate the solution (TaylorModelN) at δt including remainder Δ
        xTMN .= evaluate.( x, δt ) .+ Δ
        # Construct IntervalBox for output (wrapping effect?)
        x0 = IntervalBox( evaluate(xTMN, δq0) )

        # New initial conditions and output
        t0 += δt
        @inbounds for i in eachindex(x)
            x[i]  = Taylor1( xTMN[i], orderT )
            dx[i] = Taylor1( zero(xTMN[i]), orderT )
        end
        @inbounds t[0] = Interval(t0)
        nsteps += 1
        @inbounds tv[nsteps] = t0
        @inbounds xv[nsteps] = x0
        @inbounds xTMNv[:,nsteps] .= xTMN
        println(nsteps, "\t", t0, "\t", x0)
        if nsteps > maxsteps
            @info("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end
    end

    return view(tv,1:nsteps), view(xv,1:nsteps), view(transpose(view(xTMNv,:,1:nsteps)),1:nsteps,:)
end
