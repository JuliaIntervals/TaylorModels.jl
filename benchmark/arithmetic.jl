# ==========
# Arithmetic
# ==========

SUITE["Arithmetic-Vector"] = BenchmarkGroup()

SUITE["Arithmetic-Scalar"] = BenchmarkGroup()

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


SUITE["Arithmetic-Vector"]["addition"] = @benchmarkable $B1 + $B2
SUITE["Arithmetic-Vector"]["matrix multiplication"] = @benchmarkable transpose($B1) * $B2
SUITE["Arithmetic-Vector"]["pointwise multiplication"] = @benchmarkable B1 .* B2
SUITE["Arithmetic-Vector"]["division by scalar"] = @benchmarkable $B1 / 2
SUITE["Arithmetic-Vector"]["pointwise division"] = @benchmarkable $B1 ./ $B2
SUITE["Arithmetic-Vector"]["power function"] = @benchmarkable $B1.^3
SUITE["Arithmetic-Vector"]["sine function"] = @benchmarkable sin.($B1)
SUITE["Arithmetic-Vector"]["combination of functions"] = @benchmarkable sin($B1[1, 1]) + $B1[2, 1].^2  - transpose($B1) * $B2

SUITE["Arithmetic-Scalar"]["addition"] = @benchmarkable $b1 + $b2
SUITE["Arithmetic-Scalar"]["multiplication"] = @benchmarkable $b1 * $b2
SUITE["Arithmetic-Scalar"]["division"] = @benchmarkable $b1 / 2
SUITE["Arithmetic-Scalar"]["pointwise division"] = @benchmarkable $b1 / $b2
SUITE["Arithmetic-Scalar"]["power function"] = @benchmarkable $b1^3
SUITE["Arithmetic-Scalar"]["sine function"] = @benchmarkable sin($b1)
SUITE["Arithmetic-Scalar"]["combination of functions"] = @benchmarkable sin($b1[1]) + $b1[2]^2  - $b1 * $b2
