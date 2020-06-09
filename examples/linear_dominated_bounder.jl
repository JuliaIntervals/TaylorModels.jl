# ## Usage example of the Linear Dominated Bounder

using Plots
using BranchAndPrune

## Beale function
fb(x1, x2) = (1.5 - x1 * (1 - x2))^2 + (2.25 - x1 * (1 - x2^2)) * (2.25 - x1 * (1 - x2^2)) +
    (2.625 - x1 * (1 - x2^3)) * (2.625 - x1 * (1 - x2^3))

fb(box) = fb(box[1], box[2])

# In this example we're going to demostrate the usefulness of the 
# linear dominated bounder (ldb) in searching for a minimum.

mutable struct MinimumSearch{N, T} <: AbstractBreadthFirstSearch{IntervalBox{N, T}}
    f::Function
    initial::IntervalBox{N, T}
    tol::Float64
    cutoff::Float64
end

function BranchAndPrune.process(search::MinimumSearch, interval)
    bound = search.f(interval)

    if bound.lo > search.cutoff
        return :discard, interval
    end

    if bound.hi < search.cutoff
        search.cutoff = bound.hi
    end

    if prod(diam.(interval)) < search.tol
        return :store, interval
    else
        return :bisect, interval
    end
end

mutable struct MinimumSearchTM{N, T} <: AbstractBreadthFirstSearch{IntervalBox{N, T}}
    f::Function
    initial::IntervalBox{N, T}
    tol::Float64
    cutoff::Float64
    order::Int64
    ldb_iter::Int64
    ldb_tol::Float64
end

function BranchAndPrune.process(search::MinimumSearchTM, interval)
    X0 = IntervalBox(mid(interval))
    ftm = define_taylor_model(interval, X0, search.f, order=search.order)
    bound = linear_dominated_bounder(ftm, max_iter=search.ldb_iter, ϵ=search.ldb_tol)

    if bound.lo > search.cutoff
        return :discard, interval
    end

    if bound.hi < search.cutoff
        search.cutoff = bound.hi
    end

    if prod(diam.(interval)) < search.tol
        return :store, interval
    else
        return :bisect, interval
    end
end

function BranchAndPrune.bisect(::Union{MinimumSearch, MinimumSearchTM}, box)
    interval_x = box[1]
    interval_y = box[2]
    w_x, w_y = diam.(box)
    mid_x, mid_y = mid.(box)
    if w_x > w_y
        return IntervalBox(interval_x.lo .. mid_x, interval_y),
            IntervalBox(mid_x .. interval_x.hi, interval_y)
    else
        return IntervalBox(interval_x, interval_y.lo .. mid_y),
            IntervalBox(interval_x, mid_y .. interval_y.hi)
    end
end

function run_search(f, interval, tol; bounder=false, order=4, cutoff=1.,
                    ldb_iter=5, ldb_tol=1e-15)
    history = []
    if bounder
        search = MinimumSearchTM(f, interval, tol, cutoff,
                                 order, ldb_iter, ldb_tol)
    else
        search = MinimumSearch(f, interval, tol, cutoff)
    end
    
    local end_tree = nothing
    
    for working_tree in search
        end_tree = working_tree
        push!(history, data(end_tree))
    end
    
    return data(end_tree), history
end

box = IntervalBox(-4.5 .. 4.5, -4.5 .. 4.5)
tol = 1 / 32

d, history = run_search(fb, box, tol)
history = unique([h for hh in history for h in hh])

p1 = plot(history, c=:wheat2, legend=false, ratio=:equal, title="Interval")
plot!(p1, d, c=:blue)
scatter!(p1, [3], [0.5], c=:red, ms=1.5)

d, history = run_search(fb, box, tol, bounder=true, order=6)
history = unique([h for hh in history for h in hh])

p2 = plot(history, c=:wheat2, legend=false, ratio=:equal, title="LDB")
plot!(p2, d, c=:blue)
scatter!(p2, [3], [0.5], c=:red, ms=1.5)

plot(p1, p2, layout=(1, 2), dpi=150)
ylims!(-4.6, 4.6)
