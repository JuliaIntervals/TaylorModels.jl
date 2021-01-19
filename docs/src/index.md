# TaylorModels.jl

This package combines the [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) and
[TaylorSeries.jl](https://github.com/JuliaDiff/TaylorSeries.jl) packages to provide
**Taylor models**, i.e. Taylor polynomials with guaranteed error bounds to approximate functions.

An introduction is available in [this video from JuliaCon 2018](https://www.youtube.com/watch?v=o1h7BUW04NI).

## Introduction to Taylor Models

Taylor models provide a way to rigorously manipulate and evaluate functions using floating-point arithmetic. They have been widely used for validated computing: in global optimization and range bounding, for validated solutions of ODEs, rigorous quadrature, etc.

A Taylor model (TM) of order $n$ for a function $f$ which is supposed to be $n + 1$ times continuously differentiable over an interval $[a,b]$, is a rigorous polynomial approximation of $f$. Specifically, it is a couple $(P, \Delta)$ formed by a polynomial $P$ of degree $n$, and an interval part $\Delta$, such that $f(x) − P(x) \in \Delta$, $\forall x ∈ [a,b]$. Roughly speaking, as their name suggests, the polynomial can be seen as a Taylor expansion of the function at a given point. The interval $\Delta$ (also called interval remainder) provides the validation of the approximation, meaning that it provides an enclosure of all the approximation errors encountered (truncation, roundings).

Here we generate TMs of order 6 and 7 over $I = [-0.5,1.0]$. We can view a TM as a a tube around the actual function.

```julia
using TaylorModels

f(x) = x*(x-1.1)*(x+2)*(x+2.2)*(x+2.5)*(x+3)*sin(1.7*x+0.5)
a =  -0.5 .. 1.0 # Domain
x0 = mid(a)     # Expansion point
tm6 = TaylorModel1(6, interval(x0), a) # Independent variable for Taylor models, order 6
tm7 = TaylorModel1(7, interval(x0), a)  # Independent variable for Taylor models, order 7
# Taylor models corresponding to f(x) of order 6 and 7
ftm6 = f(tm6)
ftm7 = f(tm7)

# Now the plot
using Plots; gr()
plot(range(inf(a), stop=sup(a), length=1000), x->f(x), label="f(x)", lw=2, xaxis="x", yaxis="f(x)")
plot!(ftm6, label="6th order")
plot!(ftm7, label="7th order")
```
![intro_plot](./assets/intro_plot.png)

## Taylor model variables and range bounding

Taylor models can be applied to the problem of range bounding, that is to find an
interval $I \subseteq \mathbb{R}$ such that $f(x) \in I$ on a given domain $x \in D$
(including possible floating point errors, see examples below). There are several ways
to construct a Taylor model. A convenient way is to define a "Taylor model variable",
which is then passed as argument to Julia functions. The following examples should help
to clarify this method. To simplify the presentation we have only considered the univariate case,
but this package can also handle multivariate Taylor models with the `TaylorModelN` type.

Here we construct a Taylor model variable specifying that:

- The truncation order is 3.
- The expansion is around the origin (`interval(0)`).
- The domain is the real interval centered around the origin $[-0.5, 0.5]$.

```@example range_bound_1d
using TaylorModels

t = TaylorModel1(3, interval(0), -0.5..0.5)
```
Here the polynomial part is (the interval) $1$, and the remainder is zero.
We can pass this Taylor model variable to any Julia function, for example:

```@example range_bound_1d
texp = exp(t)
```

This expression is a polynomial of order 3 (in agreement with the truncation order
specified in the construction of `t`), whose coefficients are intervals that are
guaranteed to contain the exact coefficient of the Taylor expansion of the function
$t \mapsto e^t$ in $D : [-0.5, 0.5] \subset \mathbb{R}$. Similarly, we can expand
trigonometric functions such as $t \mapsto \sin(t)$:

```@example range_bound_1d
tsin = sin(t)
```

Common arithmetic operators, such as addition (`+`) and multiplication (`*`) work
with Taylor model variables out-of-the-box:

```@example range_bound_1d
s = texp + tsin
```

```@example range_bound_1d
p = texp * tsin
```

To bound the range of a Taylor model in one variable, use the function `bound_taylor1`:

```@example range_bound_1d
using TaylorModels: bound_taylor1

[bound_taylor1(x) for x in [s, p]]
```

This shows in particular that $0.12499 \leq e^t + \sin(t) \leq 2.12501$ and that
$-0.291667 \leq e^t \sin(t) \leq 0.791667$ for all $t \in D$. Such bounds are in general not tight.
If desired, the common approach to improve the bounds is to evaluate the Taylor model
on a smaller interval, e.g.

```@example range_bound_1d
D = domain(s) # domain -0.5 .. 0.5

E = evaluate(s, D) # original, no mincing
```

```@example range_bound_1d
Dm = mince(D, 8) # split the domain into 8 smaller chunks

Em = evaluate.(s, Dm) # evaluate the Taylor model on each sub-domain

Rm = reduce(hull, Em) # take the convex hull, i.e. the smallest interval that contains them all
```

Here the lower bound has been improved by mincing (or splitting) the domain, and it
may improve by repeating such operation recursively on smaller domains.
In particular, the fact that the lower bound is greater than zero constitutes an
algorithmic proof that $s : t \mapsto e^t + \sin(t)$ is positive on $D$.
Let's visualize the function $s(t)$ and the bounds obtained so far.

```@example range_bound_1d
using Plots

Dt = range(-0.5, 0.5, length=100)

fig = plot(xlab="t", ylab="s(t)", legend=:topleft)
plot!(fig, Dt, t -> exp(t) + sin(t), lab="", c=:black)

# range bounds
plot!(fig, Dt, t -> sup(E), lab="N = 1", c=:blue, style=:dash)
plot!(fig, Dt, t -> inf(E), c=:blue, lab="", style=:dash)

plot!(fig, Dt, t -> sup(Rm), lab="N = 8", c=:red, style=:dash)
plot!(fig, Dt, t -> inf(Rm), c=:red, lab="", style=:dash)

R16 = reduce(hull, evaluate.(s, mince(D, 16)))
plot!(fig, Dt, t -> sup(R16), lab="N = 16", c=:orange, style=:dash)
plot!(fig, Dt, t -> inf(R16), c=:orange, lab="", style=:dash)
```

## Internal representation

Consider again the Taylor model variable from the [Taylor model variables and range bounding](@ref) example.

```@example range_bound_1d
t
```

Such constructor is an alias for

```julia
    TaylorModel1(x0 + Taylor1(eltype(x0), ord), zero(dom), interval(x0), dom)
```

Taylor models in one variable are internally represented using four fields: a Taylor
series (`pol`) in one variable that holds the polynomial approximation of order `ord`;
the interval remainder (`rem`); the expansion point (`x0`), and the interval domain of
interest (`dom`). Getter functions are defined for each of these fields:

```@example range_bound_1d
get_order(t)
```

```@example range_bound_1d
remainder(t)
```

```@example range_bound_1d
polynomial(t)
```

```@example range_bound_1d
domain(t)
```

```@example range_bound_1d
expansion_point(t)
```

Finally, note that the Taylor model type is a parametric structure with two fields,
`T` and `S`. The first field, `T`, refers to the numeric type of the coefficients of
the polynomial, in this case an interval with double precision floating point values
(`Interval{Float64}`). The second field, `S`, refers to the numeric type of the interval
that holds the remainder, expansion point and domain of interest (in this case `Float64`).

```@example range_bound_1d
typeof(t)
```

If we had defined the expansion point using `0.0` instead of `interval(0)`, the
coefficients of (the polynomial part of) this Taylor model variable would be floats instead of intervals.

```@example range_bound_1d
z = TaylorModel1(3, 0.0, -0.5..0.5)
```

```@example range_bound_1d
typeof(z)
```

```@example range_bound_1d
polynomial(z)
```

Using a polynomial with interval coefficients guarantees that all arithmetic
operations involving `t` are conservative, or rigorous, with respect to floating
point arithmetic.

### Authors
- [Luis Benet](http://www.cicc.unam.mx/~benet/), Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders), Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

## Bibliography
- [*Rigorous Polynomial Approximations and Applications*](https://tel.archives-ouvertes.fr/tel-00657843), Mioara Maria Joldes, Ecole normale supérieure de lyon - ENS LYON (2011)

## Acknowledgements ##

Financial support is acknowledged from DGAPA-UNAM PAPIIT grants IN-117117, IG-100616 and IG-100819. DPS acknowledges support through a *Cátedra Marcos Moshinsky* (2018).
