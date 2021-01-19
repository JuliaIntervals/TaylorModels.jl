# TaylorModels.jl

This package combines the `IntervalArithmetic.jl` and `TaylorSeries.jl` packages to provide **Taylor models**, i.e.
Taylor polynomials with guaranteed error bounds to approximate functions.

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

## Rigorous function evaluation

We can use Taylor series and Taylor models to evaluate a function where naive
evaluation using floating-point arithmetic would give incorrect results.
For example, consider the following function:

```math
g : \mathbb{R} \to \mathbb{R},\qquad g(\delta) = \frac{2(1 - δe^{-\delta} - e^{-\delta})}{(1 - e^{-\delta})^2}
```

Using calculus we can prove that $g$ is continuous on its domain and $g(0) = 1$.
However, the $0/0$ limit for $\delta \to 0$ can easily get floating-point computations wrong:

```@example eval_uni
g(δ) = 2(1 - δ*exp(-δ) - exp(-δ)) / (1 - exp(-δ))^2

[g(1/10^x) for x in 7:9]
```

Note that `g(1e-7)` gives approximately `0.9992`, while `g(1e-8)` gives `0.0`.
This is a common problem known as [catastrophic cancellation](https://en.wikipedia.org/wiki/Catastrophic_cancellation)
and it is mainly due to limited accuarcy of double precision floating point arithmetic.

```@example eval_uni
using Plots

fig = plot(ylab="g(δ)", xlab="δ")

dt = range(-10, 10, length=100)
plot!(fig, dt, g, color=:black, lab="", lw=2.0)
```

We can use a Taylor series with rational coefficients to do Taylor series manipulations
and more accurately evaluate $g$ near zero. First we define a Taylor (series) variable
of the desired maximum order:

```@example eval_uni
using TaylorSeries

n = 16  # order
t = Taylor1(Rational{Int}, n) # the coeffs of this Taylor series are rationals
```

Now we define a truncated exponential $u_n = \lfloor e^{-t} \rfloor_n := \sum_{i=0}^{n} \frac{(-t)^i}{i!} + \mathcal{O}(t^{n+1})$:

```@example eval_uni
u = sum((-t)^i / factorial(i) for i in 0:n)
```

Finally, we compute the Taylor series of $g(\delta)$ of the desired order:

```@example eval_uni
g_taylor = 2*(1 - t*u -u) / (1 - u)^2
```

```@example eval_uni
[evaluate(g_taylor, 1/10^x) for x in 7:9]
```

The evaluation for small values of $\delta$ gives a result that is closer to zero, as expected.
We can use the coefficients obtained to make a function that uses the Taylor
series expansion for small $\delta$, and the Julia's built-in function `expm1`
(which accurately computes $e^x - 1$) for larger $\delta$. As a rule of thumb, we take
$\delta = 0.2$ to decide which method to use. To evaluate the Taylor series we use the
`@evalpoly` macro that evaluates a polynomial using Horner's method.

```@example eval_uni
const coeffs = Tuple(Float64.(g_taylor.coeffs))
```

!!! note
    We have chosen to convert the array to a tuple. This way `@evalpoly` unrolls the loops
    of Horner's method at compile time because the number of coefficients of a tuple is statically known.

```@example eval_uni
@evalpoly(1/10^7, coeffs...)
```

```@example eval_uni
function g_poly(δ)
    if abs(δ) < 0.2
        return @evalpoly(δ, coeffs...)
    else
        y = expm1(-δ) # exp(-δ) - 1
        return -2(y + δ*exp(-δ)) / y^2
    end
end

[g_poly(1/10^x) for x in 7:9]
```

Let's now build a [rational function](https://en.wikipedia.org/wiki/Rational_function)
whose numerator and denominator are Taylor models, with guaranteed remainder error bounds.
We take the approach of using a Taylor model variable on the domain
$D = [-0.2, 0.2] \subseteq \mathbb{R}$ of order 16.

```@example eval_uni
using TaylorModels

tm = TaylorModel1(16, interval(0), -0.1..0.1)
```

```@example eval_uni
num(δ) = 2(1 - δ*exp(-δ) - exp(-δ))
denom(δ) = (1 - exp(-δ))^2
```

```@example eval_uni
N = num(tm)
```

```@example eval_uni
remainder(N)
```

```@example eval_uni
D = denom(tm)
```

```@example eval_uni
remainder(D)
```

### Authors
- [Luis Benet](http://www.cicc.unam.mx/~benet/), Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México (UNAM)
- [David P. Sanders](http://sistemas.fciencias.unam.mx/~dsanders), Departamento de Física, Facultad de Ciencias, Universidad Nacional Autónoma de México (UNAM)

## Bibliography
- [*Rigorous Polynomial Approximations and Applications*](https://tel.archives-ouvertes.fr/tel-00657843), Mioara Maria Joldes, Ecole normale supérieure de lyon - ENS LYON (2011)

## Acknowledgements ##

Financial support is acknowledged from DGAPA-UNAM PAPIIT grants IN-117117, IG-100616 and IG-100819. DPS acknowledges support through a *Cátedra Marcos Moshinsky* (2018).
