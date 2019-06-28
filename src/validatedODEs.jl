# Some methods for validated integration of ODEs

"""
    remainder_taylorstep!(f!, t, x, dx, xI, dxI, δI, δt)

Returns a remainder for the integration step for the dependent variables (`x`)
checking that the solution satisfies the criteria for existence and uniqueness.
"""
function remainder_taylorstep!(f!::Function, t::Taylor1{T},
        x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
        xI::Vector{Taylor1{Interval{T}}}, dxI::Vector{Taylor1{Interval{T}}},
        δI::IntervalBox{N,T}, δt::Interval{T}) where {N,T}

    orderT = get_order(dx[1])
    aux = δt^(orderT+1)
    Δx  = IntervalBox( [  xI[i][orderT+1] * aux for i in eachindex(xI)] )
    Δ0  = IntervalBox( [  dx[i][orderT](δI) * aux / (orderT+1) for i in eachindex(x)] )
    Δdx = IntervalBox( [ dxI[i][orderT+1] * aux for i in eachindex(xI)] )
    Δ = Δ0 + Δdx * δt
    Δxold = Δx

    # Checking existence and uniqueness
    iscontractive(Δ, Δx) && return Δx

    # If the check didn't work, compute new remainders. A new Δx is proposed,
    # and the corresponding Δdx is computed
    xxI  = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
    dxxI = Array{Taylor1{TaylorN{Interval{T}}}}(undef, N)
    vv = Array{Interval{T}}(undef, N)
    for its = 1:10
        # Remainder of Picard iteration
        Δ = picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0)

        # Checking existence and uniqueness
        iscontractive(Δ, Δx) && return Δx
        # iscontractive(Δ, Δx) && return _contract_iteration!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δdx, Δ0)

        # Expand Δx in the directions needed
        Δxold = Δx
        if Δ ⊆ Δx
            @inbounds for ind in 1:N
                # Widen the directions where ⊂ does not hold
                vv[ind] = Δx[ind]
                if Δ[ind] == Δx[ind]
                    vv[ind] = widen.(Δ[ind])
                end
            end
            Δx = IntervalBox(vv)
            continue
        end
        Δx = Δ
    end

    # If it didn't work, throw an error
    @format full
    error("""
    Error: it cannot prove existence and unicity of the solution:
        t0 = $(t[0])
        δt = $(δt)
        Δ  = $(Δ)
        Δxo = $(Δxold)
        Δx = $(Δx)
        $(Δ .⊆ Δxold)
    """)
end


"""
    iscontractive(Δ, Δx)

Checks if `Δ .⊂ Δx` is satisfied. If ``Δ ⊆ Δx` is satisfied, it returns
`true` if all cases where `==` holds corresponds to the zero `Interval`.
"""
function iscontractive(Δ::IntervalBox{N,T}, Δx::IntervalBox{N,T}) where{N,T}
    zi = Interval{T}(0, 0)
    @inbounds for ind in 1:N
        Δ[ind] ⊂ Δx[ind] && continue
        Δ[ind] == Δx[ind] == zi && continue
        return false
    end
    return true
end


"""
    picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0)

Return the remainder of Picard operator
"""
function picard_remainder!(f!::Function, t::Taylor1{T},
    x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
    xxI::Vector{Taylor1{TaylorN{Interval{T}}}},
    dxxI::Vector{Taylor1{TaylorN{Interval{T}}}},
    δI::IntervalBox{N,T}, δt::Interval{T},
    Δx::IntervalBox{N,T}, Δ0::IntervalBox{N,T}) where {N,T}

    # Extend `x` and `dx` to have interval coefficients
    zI = zero(Δx[1])
    @inbounds for ind in eachindex(x)
        xxI[ind]  = x[ind] + Δx[ind]
        dxxI[ind] = dx[ind] + zI
    end

    # Compute `dxxI` from the equations of motion
    f!(t, xxI, dxxI)

    # Picard iteration, considering only the bound of `f` and the last coeff of f
    Δdx = IntervalBox( evaluate.( (dxxI - dx)(δt), δI... ) )
    Δ = Δ0 + Δdx * δt
    return Δ
end


# Picard iterations to contract further Δx, once Δ ⊂ Δx holds
# **Currently not used**
function _contract_iteration!(f!::Function, t::Taylor1{T},
        x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}},
        xxI::Vector{Taylor1{TaylorN{Interval{T}}}}, dxxI::Vector{Taylor1{TaylorN{Interval{T}}}},
        δI::IntervalBox{N,T}, δt::Interval{T},
        Δx::IntervalBox{N,T}, Δdx::IntervalBox{N,T}, Δ0::IntervalBox{N,T}) where {N,T}

    # Some abbreviations
    zI = zero(Δx[1])
    Δ = Δ0 + Δdx * δt
    Δxold = Δx

    # Picard contractions
    for its = 1:10
        # Remainder of Picard iteration
        Δ = picard_remainder!(f!, t, x, dx, xxI, dxxI, δI, δt, Δx, Δ0)

        # If contraction doesn't hold, return old bound
        iscontractive(Δ, Δx) || return Δxold

        # Contract estimate
        Δxold = Δx
        Δx = Δ
    end

    return Δxold
end


"""
    absorb_remainder(a::TaylorModelN{N,T,T}) where {N,T}

Returns a TaylorModelN, equivalent to `a`, such that the remainder
is mostly absorbed in the coefficients. The linear shift assumes
that `a` is normalized to the `IntervalBox(-1..1, Val(N))`.

Ref: Xin Chen, Erika Abraham, and Sriram Sankaranarayanan,
"Taylor Model Flowpipe Construction for Non-linear Hybrid
Systems", in Real Time Systems Symposium (RTSS), pp. 183-192 (2012),
IEEE Press.
"""
function absorb_remainder(a::TaylorModelN{N,T,T}) where {N,T}
    Δ = remainder(a)
    orderQ = get_order(a)
    δ = IntervalBox(Interval{T}(-1,1), Val(N))
    aux = diam(Δ)/(2N)
    rem = Interval{T}(0, 0)

    # Linear shift
    lin_shift = mid(Δ) + sum((aux*TaylorN(i, order=orderQ) for i in 1:N))
    bpol = a.pol + lin_shift

    # Compute the new remainder
    aI = a(δ)
    bI = bpol(δ)

    if bI ⊆ aI
        rem = Interval(aI.lo-bI.lo, aI.hi-bI.hi)
    elseif aI ⊆ bI
        rem = Interval(bI.lo-aI.lo, bI.hi-aI.hi)
    else
        r_lo = aI.lo-bI.lo
        r_hi = aI.hi-bI.hi
        if r_lo > 0
            rem = Interval(-r_lo, r_hi)
        else
            rem = Interval( r_lo, -r_hi)
        end
    end

    return TaylorModelN(bpol, rem, a.x0, a.dom)
end


"""
    validated-step!
"""
function validated_step!(f!, t::Taylor1{T},
        x::Vector{Taylor1{TaylorN{T}}}, dx::Vector{Taylor1{TaylorN{T}}}, xaux::Vector{Taylor1{TaylorN{T}}},
        tI::Taylor1{T},
        xI::Vector{Taylor1{Interval{T}}}, dxI::Vector{Taylor1{Interval{T}}}, xauxI::Vector{Taylor1{Interval{T}}},
        t0::T, tmax::T, x0::Vector{TaylorN{T}}, x0I::Vector{Interval{T}},
        xTMN::Vector{TaylorModelN{N,T,T}}, xv::Vector{IntervalBox{N,T}},
        rem::Vector{Interval{T}}, δq_norm::IntervalBox{N,T},
        q0::IntervalBox{N,T}, q0box::IntervalBox{N,T}, nsteps::Int,
        orderT::Int, abstol::T, parse_eqs::Bool,
        check_property::Function=(t, x)->true) where {N,T}
    #
    # One step integration (non-validated)
    δt = TaylorIntegration.taylorstep!(f!, t, x, dx, xaux,
        t0, tmax, x0, orderT, abstol, parse_eqs)
    # One step integration for the initial box
    δtI = TaylorIntegration.taylorstep!(f!, tI, xI, dxI, xauxI,
        t0, tmax, x0I, orderT+1, abstol, parse_eqs)

    # This updates the `dx[:][orderT]` and `dxI[:][orderT+1]`, which are currently zero
    f!(t, x, dx)
    f!(tI, xI, dxI)

    # Test if `check_property` is satisfied; if not, half the integration time.
    # If after 25 checks `check_property` is not satisfied, thow an error.
    nsteps += 1
    issatisfied = false
    # Δ = zero.(δq_norm)
    for nchecks = 1:25
        # Validate the solution: remainder consistent with Schauder thm
        Δ = remainder_taylorstep!(f!, t, x, dx, xI, dxI, δq_norm, Interval(0.0, δt))
        rem .= rem .+ Δ

        # Create TaylorModelN to store remainders and evaluation
        @inbounds begin
            for i in eachindex(x)
                xTMN[i] = fp_rpa( TaylorModelN(x[i](0..δt), rem[i], q0, q0box) )

                # If remainder is still too big, do it again
                j = 0
                while (j < 10) && (mag(rem[i]) > 1.0e-10)
                    j += 1
                    xTMN[i] = absorb_remainder(xTMN[i])
                    rem[i] = remainder(xTMN[i])
                end
            end
            xv[nsteps] = evaluate(xTMN, δq_norm) # IntervalBox

            if !check_property(t0+δt, xv[nsteps])
                δt = δt/2
                continue
            end
        end # @inbounds

        issatisfied = true
        break
    end

    if !issatisfied
        error("""
            `check_property` is not satisfied:
            $t0 $nsteps $δt
            $(xv[nsteps])
            $(check_property(t0+δt, xv[nsteps]))""")
    end

    return δt
end


function validated_integ(f!, qq0::AbstractArray{T,1}, δq0::IntervalBox{N,T},
        t0::T, tmax::T, orderQ::Int, orderT::Int, abstol::T;
        maxsteps::Int=500, parse_eqs::Bool=true,
        check_property::Function=(t, x)->true) where {N, T<:Real}

    # Set proper parameters for jet transport
    @assert N == get_numvars()
    dof = N

    # Some variables
    R   = Interval{T}
    q0 = IntervalBox(qq0)
    t   = t0 + Taylor1(orderT)
    tI  = t0 + Taylor1(orderT+1)
    δq_norm = IntervalBox(Interval{T}(-1, 1), Val(N))
    q0box = q0 .+ δq_norm

    # Allocation of vectors
    # Output
    tv    = Array{T}(undef, maxsteps+1)
    xv    = Array{IntervalBox{N,T}}(undef, maxsteps+1)
    xTM1v = Array{TaylorModel1{TaylorN{T},T}}(undef, dof, maxsteps+1)
    # xTMNv = Array{TaylorModelN{N,T,T}}(undef, dof, maxsteps+1)
    # Internals: jet transport integration
    x     = Array{Taylor1{TaylorN{T}}}(undef, dof)
    dx    = Array{Taylor1{TaylorN{T}}}(undef, dof)
    xaux  = Array{Taylor1{TaylorN{T}}}(undef, dof)
    x0    = Array{TaylorN{T}}(undef, dof)
    xTMN  = Array{TaylorModelN{N,T,T}}(undef, dof)
    # Internals: Taylor1{Interval{T}} integration
    xI    = Array{Taylor1{Interval{T}}}(undef, dof)
    dxI   = Array{Taylor1{Interval{T}}}(undef, dof)
    xauxI = Array{Taylor1{Interval{T}}}(undef, dof)
    x0I   = Array{Interval{T}}(undef, dof)

    # Set initial conditions
    zI = zero(R)
    rem = Array{Interval{T}}(undef, dof)
    @inbounds for i in eachindex(x)
        qaux = normalize_taylor(qq0[i] + TaylorN(i, order=orderQ), δq0, true)
        x[i] = Taylor1( qaux, orderT)
        dx[i] = x[i]
        x0[i] = copy(qaux)
        xTMN[i] = TaylorModelN(x[i][0], zI, q0, q0box)
        #
        xI[i] = Taylor1( q0box[i], orderT+1 )
        dxI[i] = xI[i]
        x0I[i] = qaux(δq_norm)
        rem[i] = zI
        #
        xTM1v[i, 1] = TaylorModel1(deepcopy(x[i]), zI, zI, zI)
        # xTMNv[i, 1] .= xTMN[i]
    end

    # Output vectors
    @inbounds tv[1] = t0
    @inbounds xv[1] = IntervalBox( evaluate(xTMN, δq_norm) )
    # @inbounds xTMNv[:, 1] .= xTMN[:]

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

        # Validated step of the integration
        δt = validated_step!(f!, t, x, dx, xaux, tI, xI, dxI, xauxI,
            t0, tmax, x0, x0I, xTMN, xv, rem, δq_norm,
            q0, q0box, nsteps, orderT, abstol, parse_eqs, check_property)
        # @show(rem)
        # println()

        # New initial conditions and time
        nsteps += 1
        t0 += δt
        @inbounds t[0] = t0
        @inbounds tI[0] = t0
        @inbounds tv[nsteps] = t0
        @inbounds for i in eachindex(x)
            xTM1v[i, nsteps] = TaylorModel1(deepcopy(x[i]), rem[i], 0..0, 0 ..δt)
            # xTMNv[i, nsteps] = xTMN[i]
            aux = x[i](δt)
            x[i]  = Taylor1( aux, orderT )
            dx[i] = Taylor1( zero(aux), orderT )
        end
        # @show(IntervalBox(rem))
        # @inbounds xTMNv[:, nsteps] .= xTMN[:]

        # println(nsteps, "\t", t0, "\t", remainder.(xTMN[:]), "\t", diam(Δ))
        if nsteps > maxsteps
            @info("""
            Maximum number of integration steps reached; exiting.
            """)
            break
        end

    end

    return view(tv,1:nsteps), view(xv,1:nsteps), view(xTM1v, :, 1:nsteps)
end
