# show.jl


function show(io::IO, a::Union{TM1AbsRem, TM1RelRem, TMNAbsRem})
    if TaylorSeries._show_default[end]
        return Base.show_default(io, a)
    else
        return print(io, pretty_print(a))
    end
end

for T in (:TM1AbsRem, :TMNAbsRem, :TM1RelRem)
    @eval function pretty_print(a::$T)
        _bigOnotation = TaylorSeries.bigOnotation[end]
        _bigOnotation && TaylorSeries.displayBigO(false)
        strout = $T == TMNAbsRem ?
            string(a.pol, " + ", a.rem) : string(a.pol, "+ ", a.rem)
        if $T == TM1RelRem
            _order = get_order(a)
            strout = strout * " t" * TaylorSeries.superscriptify(_order+1)
        end
        _bigOnotation && TaylorSeries.displayBigO(true)
        strout
    end
end
