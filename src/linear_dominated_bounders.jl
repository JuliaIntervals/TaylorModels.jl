#=
 source - https://bt.pa.msu.edu/pub/papers/GOM05/GOM05.pdf
=#

using TaylorModels

function univariate_ln(f::Function, dom::Interval, order::Int)
    x0 = Interval(mid(dom))
    x = TaylorModel1(order, x0, dom)
    return f(x - x0)
end

function flip_coordinates_ln(TM::TaylorModel1)
    if getcoeff(linear_polynomial(TM.pol), 1) >= 0
        return TM
    else # change coordinate sign
        coeffs = TM.pol.coeffs
        coeffs_flip = similar(coeffs)
        @inbounds for (i, c) in enumerate(coeffs)
            coeffs_flip[i] = iseven(i) ? -c : c
        end
        pol_flip = Taylor1(coeffs_flip)
        TM_flip = TaylorModel1(pol_flip, TM.rem, -TM.x0, -TM.dom)
        return TM_flip
    end
end

function evaluate_ln(TM::TaylorModel1)
    I1 = linear_polynomial(TM.pol) + constant_term(TM.pol)
    Ih = TM.pol - I1
    bound = evaluate(Ih, TM.dom) + (evaluate(I1, TM.dom)).lo
    return bound
end

function re_expand_dom(d::Interval, ϵ::Number, TM::TaylorModel1)
    dom1 = dom = TM.dom
    a = 2 #value can be discussed
    if d.hi - d.lo <= ϵ
        return dom1, 1
    else
        dom = dom.lo..(dom.lo + (d.hi - d.lo)/(a*TM.pol.coeffs[2].hi))
           #else condition can be discussed
           if (dom.hi - dom.lo) <= (d.hi - d.lo)/(a*TM.pol.coeffs[2].hi)
               return  dom1, 1
           end
    end
    return dom, 0
end

function check_dom(a::Interval, b::Interval, TM::TaylorModel1)
    if a ⊂ b
        return evaluate(TM, -a)
    else
        return evaluate(TM, a)
    end
end

function iterate_ln(f::Function, dom::Interval, order::Int, ϵ::Number)
    Pm =  univariate_ln(f, dom, order)
    Pplus =  flip_coordinates_ln(Pm)
    M = evaluate_ln(Pplus)
    dom, flag = re_expand_dom(M, ϵ, Pplus)
    return dom, flag, Pplus
end
"""
### Examples
If the polynomials have single variable, then this functions exactly calculates 
lower bound.

```julia
julia> using TaylorModels

julia> f(x) = 1 - 5x + x^3/3
f (generic function with 1 method)

julia> I = 2..3
[2, 3]

julia> ord = 4
4

julia> ϵ = 1
1

julia> linear_dominated_bounder(f, I, ord, ϵ)
[-6.91053, -5.29935]

```
"""
function  linear_dominated_bounder(f::Function, dom::Interval, order::Int, ϵ::Number)
    dom1 = dom
    dom, flag, Pplus = iterate_ln(f, dom, order, ϵ)
    for i in range(1, stop = 100) #stop value can be discussed
        dom, flag, Pplus = iterate_ln(f, dom, order, ϵ)
        if flag == 1
            return check_dom(dom, dom1, Pplus)
        end
    end
    return check_dom(dom, dom1, Pplus)
end
