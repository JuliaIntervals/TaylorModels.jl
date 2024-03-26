# show.jl

# Is the following really needed?
# function showfull(io::IO, f::Union{TaylorModel1, RTaylorModel1, TaylorModelN})
#     print(io,
#     """Taylor model of degree $(get_order(f)):
#     - x0:  $(expansion_point(f))
#     -  I:  $(domain(f))
#     -  p:  $(polynomial(f))
#     -  Δ:  $(remainder(f))
#     """
#     )
# end
# showfull(x) = showfull(stdout::IO, x)

function show(io::IO, a::Union{TaylorModel1, RTaylorModel1, TaylorModelN})
    if TaylorSeries._show_default[end]
        return Base.show_default(io, a)
    else
        return print(io, pretty_print(a))
    end
end

for T in (:TaylorModel1, :TaylorModelN, :RTaylorModel1)
    @eval function pretty_print(a::$T)
        _bigOnotation = TaylorSeries.bigOnotation[end]
        _bigOnotation && TaylorSeries.displayBigO(false)
        a_rem = remainder(a)
        strout = $T == TaylorModelN ?
            string(a.pol, " + ", a_rem) : string(a.pol, "+ ", a_rem)
        if $T == RTaylorModel1
            _order = get_order(a)
            strout = strout * " t" * TaylorSeries.superscriptify(_order+1)
        end
        _bigOnotation && TaylorSeries.displayBigO(true)
        return strout
    end
end
