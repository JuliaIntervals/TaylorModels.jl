# ==========
# Arithmetic
# ==========

SUITE["Arithmetic"] = BenchmarkGroup()

# The examples are taken from Althoff, M., Grebenyuk, D., & Kochdumper, N. (2018).
# Implementation of Taylor models in CORA 2018. In Proc. of the 5th International Workshop on
# Applied Verification for Continuous and Hybrid Systems.

t = Taylor1(6)

# taylor models of order 6
b1 = TaylorModel1(0.5 + 1.5t, 0..0, 0.5..0.5, -1..1)
b2 = TaylorModel1(2.5 + 0.5t, 0..0, 0.5..0.5, -1..1)
b3 = TaylorModel1(-5.0 + 1.0t, 0..0, 0.5..0.5, -1..1)
b4 = TaylorModel1(5.0 + 1.0t, 0..0, 0.5..0.5, -1..1)

B1 = [b1, b2]
B2 = [b3, b4]


SUITE["Arithmetic"]["addition"] = @benchmarkable $B1 + $B2
SUITE["Arithmetic"]["matrix multiplication"] = @benchmarkable transpose($B1) * $B2
SUITE["Arithmetic"]["pointwise multiplication"] = @benchmarkable B1 .* B2
SUITE["Arithmetic"]["division by scalar"] = @benchmarkable $B1 / 2
SUITE["Arithmetic"]["pointwise division"] = @benchmarkable $B1 ./ $B2
SUITE["Arithmetic"]["power function"] = @benchmarkable $B1.^3
SUITE["Arithmetic"]["sine function"] = @benchmarkable sin.($B1)
SUITE["Arithmetic"]["combination of functions"] = @benchmarkable sin($B1[1, 1]) + 
                                            $B1[2, 1].^2  - transpose($B1) * $B2
