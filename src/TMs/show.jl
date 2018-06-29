# show.jl


function show(io::IO, a::Union{TM1AbsRem, TM1RelRem, TMNAbsRem})
    if TaylorSeries._show_default[end]
        return Base.show_default(io, a)
    else
        return print(io, pretty_print(a))
    end
end

for T in (:TM1AbsRem, :TMNAbsRem)
    @eval function pretty_print(a::$T)
        _bigOnotation = TaylorSeries.bigOnotation[end]
        _bigOnotation && TaylorSeries.displayBigO(false)
        if $T == :TM1AbsRem
            strout = string(a.pol) * "± " * string(a.rem)
        else
            strout = string(a.pol) * " ± " * string(a.rem)
        end
        _bigOnotation && TaylorSeries.displayBigO(true)
        strout
    end
end

function pretty_print(a::TM1RelRem)
    _bigOnotation = TaylorSeries.bigOnotation[end]
    _bigOnotation && TaylorSeries.displayBigO(false)
    _order = get_order(a)
    strout = string(a.pol) * "± " * string(a.rem) * " t" * TaylorSeries.superscriptify(_order+1)
    _bigOnotation && TaylorSeries.displayBigO(true)
    strout
end
