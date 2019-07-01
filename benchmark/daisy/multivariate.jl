# ======================
# Multivariate functions
# ======================
DAISY["himmilbeau"] = Dict("dom" => a×b,
    "numvars" => 2,
    "vars" => "x1 x2",
    "f" => ((x1, x2) -> (x1*x1 + x2 - 11)*(x1 * x1 + x2 - 11) + (x1 + x2*x2 - 7)*(x1 + x2*x2 - 7)),
    "ref" => Interval(85.46830677734748, 221.7338939301446) # MOSEK deg 6
    )

DAISY["kepler1"] = Dict("dom" => a×b×c×d×e×f,
    "numvars" => 6,
    "vars" => "x1 x2 x3 x4 x5 x6",
    "f" => ((x1, x2, x3, x4, x5, x6) -> x2 * x5 + x3 * x6 - x2 * x3 - x5 * x6 + x1 * (-x1 + x2 + x3 - x4 + x5 + x6)),
    "ref" => Interval(-5.255935494810441, 7.321362422825775)  # MOSEK deg 6
    )

DAISY["kepler2"] = Dict("dom" => a×b×c×d,
    "numvars" => 4,
    "vars" => "x1 x2 x3 x4",
    "f" => ((x1, x2, x3, x4) -> x1 * x4 * (-x1 + x2 + x3 - x4) +
                                x2 * (x1 - x2 + x3 + x4) + x3 * (x1 + x2 - x3 + x4) -x2 * x3 * x4 - x1 * x3 - x1 * x2 - x4),
    "ref" => Interval(-195.36974909125482, 78.3669520644375)  # MOSEK deg 6
    )

DAISY["kepler3"] = Dict("dom" => a×b×c×d×e×f,
    "numvars" => 6,
    "vars" => "x1 x2 x3 x4 x5 x6",
    "f" => ((x1, x2, x3, x4, x5, x6) -> x1 * x4 * (-x1 + x2 + x3 - x4 + x5 + x6) +
                                        x2 * x5 * (x1 - x2 + x3 + x4 - x5 + x6) +x3* x6 * (x1 + x2 - x3 + x4 + x5 - x6) -
                                        x2 * x3 * x4 -x1* x3* x5 - x1 * x2 * x6 - x4 * x5 * x6),
    "ref" => Interval(-309.8484155131222, 17.982082401462407) # MOSEK deg 6
    )

DAISY["Rigidbody1"] = Dict("dom" => a×b×c,
    "numvars" => 3,
    "vars" => "x1 x2 x3",
    "f" => ((x1, x2, x3) -> -x1*x2 - 2*x2*x3 - x1 - x3),
    "ref" => Interval(-20.786552979420335, -0.540012836551535) # MOSEK deg 6
    )

DAISY["Rigidbody2"] = Dict("dom" => a×b×c,
    "numvars" => 3,
    "vars" => "x1 x2 x3",
    "f" => ((x1, x2, x3) -> 2*(x1*x2*x3) + (3*x3*x3) - x2*(x1*x2*x3) + (3*x3*x3) - x2),
    "ref" => Interval(68.81138021006673, 359.98566570476504)  # MOSEK deg 6
    )
