source('baseline/mor.R')
Rcpp::sourceCpp('baseline/mor.cpp')
load("data/testspc.RData")


# Check correct work
cpp_erosion = mor_erosion(testspc[1,], 50)
r_erosion = mor.erosion(testspc[1,,drop=FALSE], 50)
all(cpp_erosion == r_erosion[1,])

cpp_dilation = mor_dilation(testspc[1,], 50)
r_dilation = mor.dilation(testspc[1,,drop=FALSE], 50)
all(cpp_dilation == r_dilation[1,])

cpp_opening = mor_opening(testspc[1,],50)
r_opening = mor.opening(testspc[1,,drop=F],50)
all(cpp_opening == r_opening[1,])

cpp_P = mor_P(testspc[1,],50)
r_P = mor.P(testspc[1,,drop=F],50)
all(cpp_P == r_P[1,])

cpp_morbl = baseline_mor(testspc[1,],50)
r_morbl = baseline.mor(testspc[1,,drop=F],50)
all(cpp_morbl == r_morbl[1,])

cpp_wopt = mor_getwopt(testspc[1,],50)
r_wopt = mor.getwopt(testspc[1,,drop=F],50)
cpp_wopt == (r_wopt-1)

cpp_imor = baseline_imor(testspc[1,],0.0001)
r_imor = baseline.imor(testspc[1,,drop=F],0.0001)
all(cpp_imor == r_imor[1,])


# Benchmarks
library(rbenchmark)
benchmark(
  baseline.mor(testspc,50),
  apply(testspc,1,baseline_mor,w=50),
  baseline.imor(testspc,0.0001),
  apply(testspc,1,baseline_imor,tol=0.0001),
  spc.fit.poly.below(as.hyperSpec(testspc), poly.order = 4),
  replications = 5
)
#                                                        test replications elapsed relative user.self sys.self user.child sys.child
# 4             apply(testspc, 1, baseline_imor, tol = 1e-04)            5   4.645    7.995     4.645    0.000          0         0
# 2                   apply(testspc, 1, baseline_mor, w = 50)            5   0.581    1.000     0.582    0.000          0         0
# 3                             baseline.imor(testspc, 1e-04)            5 171.611  295.372   171.614    0.017          0         0
# 1                                 baseline.mor(testspc, 50)            5  21.206   36.499    21.209    0.000          0         0
# 5 spc.fit.poly.below(as.hyperSpec(testspc), poly.order = 4)            5   0.790    1.360     0.786    0.004          0         0
