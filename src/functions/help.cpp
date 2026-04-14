#include "functions/functions.hpp"
#include <algorithm>
#include <sstream>
#include <vector>

namespace mm::cal {

 static const std::unordered_map<std::string, std::string> kFunctionHelp = {
     // ---- basic_math ----
     {                      "abs",                                                                                  "abs[x]: absolute value. For complex x, returns magnitude |x|."},
     {                     "sign",                                                                                                  "sign[x]: sign of real x. Returns -1, 0, or 1."},
     {                     "sqrt",                                                                    "sqrt[x]: square root. For negative real x, returns complex principal value."},
     {                     "cbrt",                                                                                                                       "cbrt[x]: real cube root."},
     {                    "floor",                                                                                                               "floor[x]: greatest integer <= x."},
     {                     "ceil",                                                                                                                "ceil[x]: smallest integer >= x."},
     {                    "trunc",                                                                                                                "trunc[x]: truncate toward zero."},
     {                      "pow",                                                                                    "pow[a, b]: a raised to the power b. Complex values allowed."},
     {                 "nextpow2",                                                                            "nextpow2[x]: smallest integer n such that 2^n >= x. Requires x > 0."},
     {                     "DtoR",                                                                                                           "DtoR[x]: convert degrees to radians."},
     {                     "DtoG",                                                                                                             "DtoG[x]: convert degrees to grads."},
     {                     "RtoD",                                                                                                           "RtoD[x]: convert radians to degrees."},
     {                     "RtoG",                                                                                                             "RtoG[x]: convert radians to grads."},
     {                     "GtoD",                                                                                                             "GtoD[x]: convert grads to degrees."},
     {                     "GtoR",                                                                                                             "GtoR[x]: convert grads to radians."},

     // ---- exp_log ----
     {                      "exp",                                                                                                              "exp[x]: exponential function e^x."},
     {                    "log10",                                                                                    "log10[x]: base-10 logarithm. Requires x > 0 for real input."},
     {                     "log2",                                                                                      "log2[x]: base-2 logarithm. Requires x > 0 for real input."},
     {                      "log",                          "log[x] or log[base, x]: natural logarithm or logarithm with custom base. Complex principal value is used when needed."},
     {                       "ln",                                                                                                            "ln[x] or ln[base, x]: alias of log."},
     {                    "log1p",                                                                       "log1p[x]: numerically stable log[1 + x]. Requires x > -1 for real input."},
     {                    "expm1",                                                                                                       "expm1[x]: numerically stable exp[x] - 1."},

     // ---- trig ----
     {                      "sin",                                                                                           "sin[x]: sine. Angle input is interpreted in degrees."},
     {                      "cos",                                                                                         "cos[x]: cosine. Angle input is interpreted in degrees."},
     {                      "tan",                                                                                        "tan[x]: tangent. Angle input is interpreted in degrees."},
     {                      "cot",                                                                                      "cot[x]: cotangent. Angle input is interpreted in degrees."},
     {                      "sec",                                                                                         "sec[x]: secant. Angle input is interpreted in degrees."},
     {                      "csc",                                                                                       "csc[x]: cosecant. Angle input is interpreted in degrees."},
     {                     "asin",                                                                                               "asin[x]: inverse sine. Returns angle in degrees."},
     {                     "acos",                                                                                             "acos[x]: inverse cosine. Returns angle in degrees."},
     {                     "atan",                                                                                            "atan[x]: inverse tangent. Returns angle in degrees."},
     {                    "atan2",                                                                                        "atan2[y, x]: quadrant-aware inverse tangent in degrees."},

     // ---- hyperbolic ----
     {                     "sinh",                                                                                                                      "sinh[x]: hyperbolic sine."},
     {                     "cosh",                                                                                                                    "cosh[x]: hyperbolic cosine."},
     {                     "tanh",                                                                                                                   "tanh[x]: hyperbolic tangent."},
     {                    "asinh",                                                                                                             "asinh[x]: inverse hyperbolic sine."},
     {                    "acosh",                                                                                                           "acosh[x]: inverse hyperbolic cosine."},
     {                    "atanh",                                                                                                          "atanh[x]: inverse hyperbolic tangent."},
     {                     "csch",                                                                                                    "csch[x]: hyperbolic cosecant = 1 / sinh[x]."},
     {                     "sech",                                                                                                      "sech[x]: hyperbolic secant = 1 / cosh[x]."},
     {                     "coth",                                                                                             "coth[x]: hyperbolic cotangent = cosh[x] / sinh[x]."},

     // ---- statistics / arithmetic helpers ----
     {                     "fact",                                                                                                      "fact[n]: factorial n! for integer n >= 0."},
     {                      "gcd",                                                                                           "gcd[a, b, ...]: greatest common divisor of integers."},
     {                     "comb",                                                                                                          "comb[n, r]: binomial coefficient nCr."},
     {                     "perm",                                                                                                             "perm[n, r]: permutation count nPr."},
     {                      "lcm",                                                                                             "lcm[a, b, ...]: least common multiple of integers."},
     {                      "sum",                                                                          "sum[x1, x2, ...] or sum[{..}]: sum of values. Complex values allowed."},
     {                     "prod",                                                                    "prod[x1, x2, ...] or prod[{..}]: product of values. Complex values allowed."},
     {                     "mean",                                                                      "mean[x1, x2, ...] or mean[{..}]: arithmetic mean. Complex values allowed."},
     {                      "mod",                                                                                                  "mod[x, y]: modulo using x - y * floor[x / y]."},
     {                  "geomean",                                                                                "geomean[x1, ...]: geometric mean. Real nonnegative values only."},
     {                 "harmmean",                                                                                 "harmmean[x1, ...]: harmonic mean. Zero values are not allowed."},
     {                 "quantile",                                                                        "quantile[p, x1, x2, ...]: linear-interpolated quantile for 0 <= p <= 1."},
     {                      "min",                                                                                                      "min[x1, x2, ...]: minimum of real values."},
     {                      "max",                                                                                                      "max[x1, x2, ...]: maximum of real values."},
     {                    "clamp",                                                                                                       "clamp[x, lo, hi]: clamp x into [lo, hi]."},
     {                    "fract",                                                                                                        "fract[x]: fractional part x - floor[x]."},
     {                    "gamma",                                                                                                                      "gamma[x]: Gamma function."},
     {                   "lgamma",                                                                                               "lgamma[x]: logarithm of absolute Gamma function."},
     {                  "digamma",                                                                                                       "digamma[x]: derivative of log[Gamma[x]]."},
     {                     "zeta",                                                                                                      "zeta[s]: Riemann zeta function for s > 1."},
     {                 "trigamma",                                                                                                         "trigamma[x]: derivative of digamma[x]."},
     {                     "beta",                                                                                                    "beta[a, b]: complete beta function B[a, b]."},
     {                    "ibeta",                                                "ibeta[a, b, x]: regularized incomplete beta function I_x[a, b], with a > 0, b > 0, 0 <= x <= 1."},
     {                   "betaln",                                                                                 "betaln[x, y]:  logarithmic beta function ln(B(x,y)). x, y > 0."},
     {                    "binom", "binom[x, y]: generalized binomial coefficient. If y is an integer, it is evaluated using the series definition via the falling factorial / y!."},
     {              "fallingfact",            "fallingfact[x, n]: descending factorial x(x-1)... Evaluates as a product when n is an integer, and as the Gamma function otherwise."},
     {               "risingfact",   "risingfact[x, n]: factorial of x(x+1)... If n is an integer, it is evaluated as a product; otherwise, it is evaluated as the Gamma function."},
     {                     "mode",                                                                          "mode[x1, x2, ...]: mode of real values. On ties, smallest value wins."},
     {                      "erf",                                                                                                                        "erf[x]: error function."},
     {                     "erfc",                                                                                                         "erfc[x]: complementary error function."},
     {                    "round",                                                                      "round[x] or round[x, n]: round to nearest integer or to n decimal digits."},
     {                      "var",                                                                                                             "var[x1, ...]: population variance."},
     {                     "vars",                                                                                                                "vars[x1, ...]: sample variance."},
     {                   "stddev",                                                                                                "stddev[x1, ...]: population standard deviation."},
     {                  "stddevs",                                                                                                   "stddevs[x1, ...]: sample standard deviation."},
     {                   "median",                                                                                                                       "median[x1, ...]: median."},
     {                      "mad",                                                                                       "mad[x1, ...]: median absolute deviation from the median."},
     {                     "madR",                                                                                          "madR[x1, ...]: mean absolute deviation from the mean."},
     {                     "skew",                                                                                                                       "skew[x1, ...]: skewness."},
     {                    "kurtp",                                                                                                    "kurtp[x1, ...]: population excess kurtosis."},
     {                    "kurts",                                                                                                        "kurts[x1, ...]: sample excess kurtosis."},
     {               "percentile",                                                                                         "percentile[p, x1, ...]: percentile with 0 <= p <= 100."},
     {                      "cov",                                                                                     "cov[x1..xn, y1..yn]: population covariance of two samples."},
     {                     "corr",                                                                                         "corr[x1..xn, y1..yn]: Pearson correlation coefficient."},
     {                      "ave",                                                                                                                       "ave[...]: alias of mean."},
     {                      "rms",                                                                                                                "rms[x1, ...]: root mean square."},
     {                       "cv",                                                                                         "cv[x1, ...]: coefficient of variation = stddev / mean."},
     {                   "stderr",                                                                                                   "stderr[x1, ...]: standard error of the mean."},
     {                   "zscore",                                                                                                            "zscore[x, mu, sigma]: z-score of x."},
     {                      "iqr",                                                                                                     "iqr[x1, ...]: interquartile range Q3 - Q1."},
     {                 "trimmean",                                                                           "trimmean[p, x1, ...]: mean after trimming fraction p from each tail."},
     {                   "winsor",                                                                                           "winsor[p, x1, ...]: winsorized mean with fraction p."},
     {             "corrspearman",                                                                                       "corrspearman[x1..xn, y1..yn]: Spearman rank correlation."},
     {                  "polylog",                                                                                              "polylog[s, z]: polylogarithm Li_s[z] for |z| < 1."},

     // ---- signal ----
     {                     "sinc",                                                    "sinc[x]: sin[x]/x with degree-based angle interpretation and removable singularity handled."},
     {                     "cosc",                                                                    "cosc[x]: [1 - cos[x]] / x style helper with removable singularity handling."},
     {                     "tanc",                                                                          "tanc[x]: tan[x] / x style helper with removable singularity handling."},
     {                    "sinhc",                                                                                     "sinhc[x]: sinh[x] / x with removable singularity handling."},
     {                    "tanhc",                                                                                     "tanhc[x]: tanh[x] / x with removable singularity handling."},
     {                     "expc",                                                                    "expc[x]: [exp[x] - 1] / x style helper with removable singularity handling."},
     {                      "fft",                                                 "fft[v] or fft[M]: FFT of 1D vector or 2D matrix. Falls back to DFT for non-power-of-two sizes."},
     {                     "ifft",                                                                                     "ifft[v] or ifft[M]: inverse FFT of 1D vector or 2D matrix."},
     {                      "dft",                                                                                      "dft[v]: direct discrete Fourier transform of a 1D vector."},
     {                  "hilbert",                                                                                             "hilbert[v]: analytic signal via Hilbert transform."},

     // ---- geometry / vector ----
     {                    "hypot",                                                                                "hypot[x, y]: robust against rounding errors in sqrt[x^2 + y^2]."},
     {                     "norm",                                                                                             "norm[x1, x2, ...]: Euclidean norm of real scalars."},
     {                     "vadd",                                                                                                                   "vadd[a, b]: vector addition."},
     {                     "vsub",                                                                                                                "vsub[a, b]: vector subtraction."},
     {                  "vscalar",                                                                                                  "vscalar[v, s]: multiply vector v by scalar s."},
     {                     "vdot",                                                                                                                "vdot[a, b]: vector dot product."},
     {                   "vcross",                                                                                                         "vcross[a, b]: 3D vector cross product."},
     {                    "vnorm",                                                                                                          "vnorm[v]: Euclidean norm of vector v."},
     {               "vmanhattan",                                                                                          "vmanhattan[a, b]: Manhattan distance between vectors."},
     {               "veuclidean",                                                                                          "veuclidean[a, b]: Euclidean distance between vectors."},
     {                      "dot",                                                                                                "dot[a, b]: dot product of two flat multivalues."},
     {                     "vsum",                                                                                                               "vsum[v]: sum of vector elements."},
     {               "vnormalize",                                                                                             "vnormalize[v]: normalized vector with unit length."},
     {                 "vproject",                                                                                          "vproject[a, b]: projection of vector a onto vector b."},
     {                   "vangle",                                                                                                "vangle[a, b]: angle between vectors in degrees."},
     {                 "vreflect",                                                                                 "vreflect[a, n]: reflection of vector a across normal vector n."},
     {            "vreflect_axis",                                                                              "vreflect_axis[a, n]: reflection of vector a across axis vector n."},
     {                  "vlength",                                                                                                                    "vlength[v]: alias of vnorm."},
     {                "vdistance",                                                                                                          "vdistance[a, b]: alias of veuclidean."},
     {                    "vunit",                                                                                              "vunit[v]: unit vector in the same direction as v."},
     {                   "scalar",                                                                                                                "scalar[v, s]: alias of vscalar."},
     {                     "lerp",                                                                                           "lerp[a, b, t]: linear interpolation between a and b."},
     {                 "distance",                                                                         "distance[x1..xn, y1..yn]: Euclidean distance between two point tuples."},
     {                  "totient",                                                                                                     "totient[n]: Euler's totient function φ[n]."},
     {                "covmatrix",                                                         "covmatrix[x1, ...]: sample variance of one sample. Despite the name, returns a scalar."},
     {               "corrmatrix",                                            "corrmatrix[x1..xn, y1..yn]: Pearson correlation of two samples. Despite the name, returns a scalar."},
     {              "percentrank",                                                                                     "percentrank[v, x1, ...]: percentage of sample values <= v."},
     {                  "winsorR",                                                                                                  "winsorR[p, x1, ...]: winsorized mean variant."},
     {                 "convolve",                                                       "convolve[x1..xn, y1..yn]: convolution-like accumulated sum for two equal-length samples."},

     // ---- matrix / linear algebra ----
     {                   "matrix",                                   "matrix[r1, r2, ...] or matrix[M]: construct a matrix from row vectors or validate an existing 2D multivalue."},
     {                     "madd",                                                                                                                   "madd[A, B]: matrix addition."},
     {                     "mmul",                                                                                                             "mmul[A, B]: matrix multiplication."},
     {               "mtranspose",                                                                                                          "mtranspose[A]: transpose of matrix A."},
     {                   "mtrace",                                                                                                                  "mtrace[A]: trace of matrix A."},
     {                    "mrank",                                                                                                                    "mrank[A]: rank of matrix A."},
     {                     "mdet",                                                                                                       "mdet[A]: determinant of a square matrix."},
     {                 "identity",                                                                                                            "identity[n]: n x n identity matrix."},
     {                    "zeros",                                                                                                  "zeros[rows, cols]: zero matrix of given size."},
     {                     "mget",                                                                                    "mget[A, i, j]: element access at row i, column j [0-based]."},
     {                    "mrows",                                                                                                                      "mrows[A]: number of rows."},
     {                    "mcols",                                                                                                                   "mcols[A]: number of columns."},
     {                    "mdiag",                                                                                                        "mdiag[A]: diagonal entries of matrix A."},
     {                      "mlu",                                                                                                   "mlu[A]: LU decomposition. Returns {P, L, U}."},
     {               "meigenvals",                                                                              "meigenvals[A]: real eigenvalue approximations of a square matrix."},
     {               "meigenvecs",                                                                                       "meigenvecs[A]: eigenvectors for symmetric real matrices."},
     {                 "minverse",                                                     "minverse[A]: matrix inverse using SVD-based inversion. Requires nonsingular square matrix."},
     {                      "mqr",                                                                                                      "mqr[A]: QR decomposition. Returns {Q, R}."},
     {                     "msvd",                                                                                     "msvd[A]: singular value decomposition. Returns {U, S, Vt}."},
     {                    "mcond",                                                                                    "mcond[A]: 2-norm condition number based on singular values."},
     {               "mcondition",                                                                                                                 "mcondition[A]: alias of mcond."},
     {                    "mlsqr",                                                              "mlsqr[A, b]: least-squares solution using QR for overdetermined systems [m >= n]."},
     {                   "mnormf",                                                                                                         "mnormf[A]: Frobenius norm of matrix A."},
     {                    "mnorm",                                                                                                                     "mnorm[A]: alias of mnormf."},
     {          "mmaxeigen_power",                                                                                       "mmaxeigen_power[A]: dominant eigenvalue by power method."},
     {                "mmaxeigen",                                                                                  "mmaxeigen[A]: maximum eigenvalue from eigenvalue computation."},
     {                   "msumsv",                                                                                                             "msumsv[A]: sum of singular values."},
     {          "misse_symmetric",                                                                                       "misse_symmetric[A]: returns 1 if A is symmetric, else 0."},
     {                    "mispd",                                                                  "mispd[A]: returns 1 if A appears positive definite via Cholesky test, else 0."},

     // ---- complex ----
     {                       "re",                                                                                                                              "re[z]: real part."},
     {                     "real",                                                                                                                          "real[z]: alias of re."},
     {                       "im",                                                                                                                         "im[z]: imaginary part."},
     {                     "imag",                                                                                                                          "imag[z]: alias of im."},
     {                      "arg",                                                                                                     "arg[z]: argument / phase angle in degrees."},
     {                     "conj",                                                                                                                    "conj[z]: complex conjugate."},
     {                    "polar",                                                                "polar[r, th]: construct complex value from magnitude r and angle th in degrees."},
     {                     "rect",                                                                                                                   "rect[r, th]: alias of polar."},
     {                      "cis",                                                                                                "cis[th]: cos[th] + i sin[th], angle in degrees."},
     {                     "proj",                                                                                               "proj[z]: Riemann sphere projection of complex z."},
     {                     "unit",                                                                                                       "unit[z]: z normalized to unit magnitude."},
     {                     "csgn",                                                                                                                        "csgn[z]: alias of unit."},
     {                      "mag",                                                                                                                          "mag[z]: alias of abs."},

     // ---- random ----
     {                     "rand",                                                                                            "rand[], rand[hi], rand[lo, hi]: random real number."},
     {                  "randint",                                                                                       "randint[], randint[hi], randint[lo, hi]: random integer."},
     {                    "randn",                                                                                    "randn[], randn[mu], randn[mu, sigma]: normal random number."},
     {                   "choice",                                                                                             "choice[x1, x2, ...]: random choice from arguments."},
     {                      "fma",                                                                                                      "fma[a, b, c]: fused multiply-add a*b + c."},

     // ---- area / volume ----
     {              "area_circle",                                                                                              "area_circle[d]: area of a circle from diameter d."},
     {            "area_triangle",                                                                     "area_triangle[base, height] or area_triangle[a, b, c]: area of a triangle."},
     {           "area_trapezoid",                                                                                                  "area_trapezoid[a, b, h]: area of a trapezoid."},
     {             "area_polygon",                                                                       "area_polygon[x1, y1, x2, y2, ...]: polygon area by the shoelace formula."},
     {             "vol_cylinder",                                                                                  "vol_cylinder[d, h]: cylinder volume from diameter and height."},
     {                 "vol_cone",                                                                                          "vol_cone[d, h]: cone volume from diameter and height."},
     {               "vol_sphere",                                                                                                    "vol_sphere[d]: sphere volume from diameter."},
     {                "vol_prism",                                                                          "vol_prism[x1, y1, ..., h]: prism volume from polygon base and height."},

     // ---- engineering ----
     {                   "stress",                                                                                                           "stress[F, A]: stress = force / area."},
     {                   "strain",                                                                                           "strain[dL, L]: strain = extension / original length."},
     {                    "young",                                                                                      "young[sigma, epsilon]: Young's modulus = stress / strain."},
     {              "moment_rect",                                                                                      "moment_rect[b, h]: second moment of area for a rectangle."},
     {            "moment_circle",                                                                                          "moment_circle[d]: second moment of area for a circle."},
     {          "sectionmod_rect",                                                                                        "sectionmod_rect[b, h]: section modulus for a rectangle."},
     {        "sectionmod_circle",                                                                                            "sectionmod_circle[d]: section modulus for a circle."},
     {         "torsion_J_circle",                                                                                     "torsion_J_circle[d]: polar moment of inertia for a circle."},
     {            "polarZ_circle",                                                                                          "polarZ_circle[d]: polar section modulus for a circle."},
     {              "bolt_stress",                                                                                 "bolt_stress[F, d]: nominal bolt stress from load and diameter."},
     {      "torque_from_preload",                                                                                    "torque_from_preload[F, d, K]: torque estimate from preload."},
     {      "preload_from_torque",                                                                                    "preload_from_torque[T, d, K]: preload estimate from torque."},
     {                 "friction",                                                                                                      "friction[mu, N]: friction force = mu * N."},

     // ---- mold injection ----
     {               "mold_clamp",                                                                                               "mold_clamp[P_avg, A_proj]: clamp force estimate."},
     {          "mold_clamp_safe",                                                                            "mold_clamp_safe[SF, P_avg, A_proj]: clamp force with safety factor."},
     {                "mold_Pinj",                                                                  "mold_Pinj[P_cav, dP_runner, dP_gate, dP_nozzle]: injection pressure estimate."},
     {            "mold_flowrate",                                                                            "mold_flowrate[V, t_fill]: flow rate from part volume and fill time."},
     {       "mold_gate_velocity",                                                                     "mold_gate_velocity[Q, A_gate]: gate velocity from flow rate and gate area."},
     {          "mold_shear_gate",                                                                                            "mold_shear_gate[Q, b, h]: gate shear rate estimate."},
     {        "mold_shear_runner",                                                                                           "mold_shear_runner[Q, D]: runner shear rate estimate."},
     {"mold_pressure_loss_runner",                                                                         "mold_pressure_loss_runner[mu, L, Q, D]: runner pressure loss estimate."},
     {      "mold_eject_friction",                                                                                              "mold_eject_friction[mu, N]: eject friction force."},
     {       "mold_eject_contact",                                                                           "mold_eject_contact[p_contact, A_contact]: contact force on ejection."},
     {         "mold_eject_total",                                                                        "mold_eject_total[mu, p_contact, A_contact]: total eject force estimate."},
     {              "mold_mu_tex",                                                                                          "mold_mu_tex[k, h]: texture-related friction estimate."},
     {          "mold_pin_stress",                                                                                        "mold_pin_stress[F_eject, n, A_pin]: ejector pin stress."},
     {    "mold_plate_deflection",                                                                          "mold_plate_deflection[K, P, a, E, t]: mold plate deflection estimate."},

     // ---- finance ----
     {                   "fin_fv",                                                                                    "fin_fv[rate, nper, pmt]: future value of repeated payments."},
     {                   "fin_pv",                                                                                   "fin_pv[rate, nper, pmt]: present value of repeated payments."},
     {                  "fin_pmt",                                                                                            "fin_pmt[rate, nper, pv]: payment amount per period."},
     {        "fin_total_payment",                                                                       "fin_total_payment[rate, nper, pv]: total amount paid across all periods."},
     {                 "fin_ppmt",                                                                           "fin_ppmt[rate, nper, per, pv]: principal payment for a given period."},
     {                 "fin_ipmt",                                                                            "fin_ipmt[rate, nper, per, pv]: interest payment for a given period."},
     {                  "fin_npv",                                                                                               "fin_npv[rate, cf1, cf2, ...]: net present value."},
     {                  "fin_irr",                                                                                               "fin_irr[cf0, cf1, ...]: internal rate of return."},
     {                 "fin_mirr",                                                                       "fin_mirr[reinvestRate, cf0, cf1, ...]: modified internal rate of return."},
     {                 "fin_nper",                                                                                                    "fin_nper[rate, pmt, pv]: number of periods."},
     {                 "fin_rate",                                                                                      "fin_rate[nper, pmt, pv]: implied rate solved numerically."},
     {              "fin_cumipmt",                                                                   "fin_cumipmt[rate, nper, pv, start, end, type]: cumulative interest payments."},
     {             "fin_cumprinc",                                                                 "fin_cumprinc[rate, nper, pv, start, end, type]: cumulative principal payments."},
     {       "fin_effective_rate",                                             "fin_effective_rate[nominal, npery]: effective annual rate from nominal rate and compounding count."},
     {         "fin_nominal_rate",                                             "fin_nominal_rate[effective, npery]: nominal annual rate from effective rate and compounding count."},
     {                 "fin_cagr",                                                                                          "fin_cagr[start, end, n]: compound annual growth rate."},

     // ---- calculus ----
     {                     "diff",                                        "diff[f, x] or diff[f, x, h]: numerical first derivative by central difference. Pass function name as f."},
     {                    "diff2",                                     "diff2[f, x] or diff2[f, x, h]: numerical second derivative by central difference. Pass function name as f."},
     {                "integrate",                                                          "integrate[f, a, b] or integrate[f, a, b, n]: numerical integral using Simpson's rule."},
     {                  "simpson",                                                   "simpson[f, a, b] or simpson[f, a, b, n]: Simpson-rule numerical integration. n must be even."},
     {                    "trapz",                                                                   "trapz[f, a, b] or trapz[f, a, b, n]: trapezoidal-rule numerical integration."},
     {                "factorint",                                                              "factorint[n]: integer prime factorization. Returns a multivalue of prime factors."},
     {             "primefactors",                                                                 "primefactors[n]: alias of factorint[n]. Returns a multivalue of prime factors."},

     // ---- others / utilities ----
     {                     "cnst",                                                                      "cnst[name]: retrieve named constant by string. e.g. cnst[\"c\"]=299792458"},
     {                      "fib",                                                                                                                      "fib[n]: Fibonacci number."},
     {                  "isprime",                                                                                                    "isprime[n]: primality test. Returns 1 or 0."},
     {                "nextprime",                                                                                              "nextprime[n]: next prime strictly greater than n."},
     {                "prevprime",                                                                                             "prevprime[n]: previous prime strictly less than n."},
     {                       "if",                                                                                                "if[cond, a, b]: returns a if cond != 0, else b."},
     {                  "convert",                                                                   "convert[value, unit] or convert[value, from_unit, to_unit]: unit conversion."},
     {                       "In",                                                                                      "In[n]: re-evaluate the nth input expression from history."},
     {                      "Out",                                                                                               "Out[n]: fetch the nth output value from history."},
     {                     "Prev",                                                                                       "Prev[k]: fetch previous output relative to current base."},
     {                        "%",                                                                               "%: fetch previous output relative to current base. like Prev[k]."},
     {                    "Clear",                                                                                                           "Clear[]: request screen/state clear."},
     {                     "Exit",                                                                                                                     "Exit[]: request REPL exit."},
     {                   "filein",                                                                                              "filein[path]: read file and evaluate its content."},
     {                  "fileout",                                                                                  "fileout[value, path]: serialize value and write it to a file."},
     {                     "clip",                                                                                            "clip[value]: copy formatted value to the clipboard."},
     {                   "silent",                                                                                         "silent[value]: suppress display of the current result."},
     {                  "explain",                                                                                        "explain[value]: show structured explanation of a value."},
 };

 const std::unordered_map<std::string, std::string> &functionHelpMap() { return kFunctionHelp; }

 std::string getFunctionHelp(const std::string &name) {
  auto it = kFunctionHelp.find(name);
  if (it != kFunctionHelp.end()) return it->second;

  return "No help available for function: " + name;
 }

 std::string getFunctionHelpIndex() {
  std::vector<std::string> names;
  names.reserve(kFunctionHelp.size());
  for (const auto &kv : kFunctionHelp) {
   names.push_back(kv.first);
  }
  std::sort(names.begin(), names.end());

  std::ostringstream oss;
  oss << "[Functions explain] (" << names.size() << "):\n";

  std::string currentGroup;
  auto groupOf = [](const std::string &name) -> std::string {
   if (name.rfind("fin_", 0) == 0) return "finance";
   if (name.rfind("mold_", 0) == 0) return "mold";
   if (name.size() >= 2 && name[0] == 'm' && std::isalpha(static_cast<unsigned char>(name[1]))) return "matrix";
   if (name.size() >= 2 && name[0] == 'v' && std::isalpha(static_cast<unsigned char>(name[1]))) return "vector";
   if (name == "fft" || name == "ifft" || name == "dft" || name == "hilbert" || name == "sinc" || name == "cosc" || name == "tanc" || name == "sinhc" || name == "tanhc" || name == "expc") return "signal";
   if (name == "re" || name == "real" || name == "im" || name == "imag" || name == "arg" || name == "conj" || name == "polar" || name == "rect" || name == "cis" || name == "proj" || name == "unit" || name == "csgn" || name == "mag") return "complex";
   return "general";
  };

  for (const auto &name : names) {
   const std::string grp = groupOf(name);
   if (grp != currentGroup) {
    if (!currentGroup.empty()) oss << "\n";
    currentGroup = grp;
    oss << "[" << currentGroup << "]\n";
   }
   oss << "  " << name << "\n";
  }

  oss << "\nUse :help <function> for details.";
  return oss.str();
 }

} // namespace mm::cal