# ======================
# Univariate functions
# ======================
DAISY_1D["sine"] = Dict("dom" => a,
    "f" => (x -> sin(x)),
    "ref" => Interval(-0.9998434996853951, 2.724232859130124)  # IntervalOptimisation
    )

DAISY_1D["bspline0"] = Dict("dom" => a,
    "f" => (x -> (1 - x) * (1 - x) * (1 - x) / 6.0),
    "ref" => Interval(0.36616666317087393, 27.729165285894855) # MOSEK deg 5
    )

DAISY_1D["bspline1"] = Dict("dom" => a,
    "f" => (x -> (3*x*x*x - 6*x*x + 4) / 6.0),
    "ref" => Interval(-65.14582571514087, 0.5631666653631329) # MOSEK deg 5
    )

DAISY_1D["bspline2"] = Dict("dom" => a,
    "f" => (x -> (-3*x*x*x  + 3*x*x + 3*x + 1) / 6.0),
    "ref" => Interval(0.07407407744904644, 53.60416654726812) # MOSEK deg 5
    )

DAISY_1D["bspline3"] = Dict("dom" => a,
    "f" => (x -> -x*x*x / 6.0),
    "ref" => Interval(0.004500048981347225, 15.18749986989248)  # MOSEK deg 5
    )
