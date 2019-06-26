# ======================
# Univariate functions
# ======================

# (*) can be checked by hand
DAISY_1D["sine"] = Dict("dom" => a,
    "f" => (x -> sin(x)),
    "ref" => Interval(-1.0, 0.977530117665097) # (*)
    )

DAISY_1D["bspline0"] = Dict("dom" => a,
    "f" => (x -> (1 - x) * (1 - x) * (1 - x) / 6.0),
    "ref" => Interval(0.36616666666666675, 27.729166666666668) # (*)
    )

DAISY_1D["bspline1"] = Dict("dom" => a,
    "f" => (x -> (3*x*x*x - 6*x*x + 4) / 6.0),
    "ref" => Interval(-65.14583333333333, 0.5631666666666667) # (*)
    )

DAISY_1D["bspline2"] = Dict("dom" => a,
    "f" => (x -> (-3*x*x*x  + 3*x*x + 3*x + 1) / 6.0),
    "ref" => Interval(0.07516666666666667, 53.604166666666664) # (*)
    )

DAISY_1D["bspline3"] = Dict("dom" => a,
    "f" => (x -> -x*x*x / 6.0),
    "ref" => Interval(0.0045, 15.1875) # (*)
    )
