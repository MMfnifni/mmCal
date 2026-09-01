// aggregate・array・elementaryの回帰テスト
#include "aggregate_array_elementary_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

namespace mmcal::tests {
namespace {

std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

error::CalcError evalError(kernel::KernelSession& session, std::string_view source) {
    try { static_cast<void>(session.evaluate(source)); }
    catch (const error::CalcError& e) { return e; }
    throw std::logic_error("Expected CalcError");
}

} // namespace

void runAggregateArrayElementaryTests(TestRunner& tests) {
    kernel::KernelSession session;

    // exact-first aggregate functions.
    tests.expectEqual(eval(session, "sum[]"), std::string{"0"},
        "empty sum uses the additive identity");
    tests.expectEqual(eval(session, "prod[]"), std::string{"1"},
        "empty product uses the multiplicative identity");
    tests.expectEqual(eval(session, "sum[1,2,3]"), std::string{"6"},
        "scalar sum is exact");
    tests.expectEqual(eval(session, "sum[{1,2,3}]"), std::string{"6"},
        "a rank-1 array can be aggregated directly");
    tests.expectEqual(eval(session, "prod[{2,3,4}]"), std::string{"24"},
        "array product is exact");
    tests.expectEqual(eval(session, "mean[1,2,4]"), std::string{"7/3"},
        "mean preserves an exact rational result");
    tests.expectEqual(eval(session, "ave[1,2,4]"), std::string{"7/3"},
        "ave is an alias of mean");
    tests.expectEqual(eval(session, "min[3,1,2]"), std::string{"1"},
        "min compares exact real numbers");
    tests.expectEqual(eval(session, "max[3,1,2]"), std::string{"3"},
        "max compares exact real numbers");
    tests.expectEqual(eval(session, "min[x,3]"), std::string{"min[x, 3]"},
        "min stays symbolic when ordering cannot be proved");

    // Exact sequence generation and explicit element-wise iteration.
    tests.expectEqual(eval(session, "range[5]"), std::string{"{1, 2, 3, 4, 5}"},
        "range[n] constructs an exact inclusive integer sequence");
    tests.expectEqual(eval(session, "range[0,1,1/3]"), std::string{"{0, 1/3, 2/3, 1}"},
        "range supports exact Rational steps without floating drift");
    tests.expectEqual(eval(session, "range[5,1,-2]"), std::string{"{5, 3, 1}"},
        "range supports descending exact sequences with an explicit negative step");
    tests.expectEqual(eval(session, "map[sin,{0,Pi/2,Pi}]"), std::string{"{0, 1, 0}"},
        "map applies a function explicitly to Array scalar leaves");
    static_cast<void>(eval(session, "mapSquare[x]:=x^2+1"));
    tests.expectEqual(eval(session, "map[mapSquare,{{1,2},{3,4}}]"),
        std::string{"{{2, 5}, {10, 17}}"},
        "map accepts user functions while preserving dense Array shape");
    tests.expectEqual(eval(session, "map[mapSquare,{{1,2},{3}}]"),
        std::string{"{{2, 5}, {10}}"},
        "map preserves ragged brace structure while visiting scalar leaves");
    tests.expectEqual(eval(session, "table[i^2,{i,5}]"), std::string{"{1, 4, 9, 16, 25}"},
        "table evaluates a held body over an exact local iterator");
    tests.expectEqual(eval(session, "table[i/2,{i,0,2,1/2}]"),
        std::string{"{0, 1/4, 1/2, 3/4, 1}"},
        "table supports exact Rational iterator ranges");
    tests.expectEqual(eval(session, "table[table[i+j,{j,2}],{i,2}]"),
        std::string{"{{2, 3}, {3, 4}}"},
        "nested table scopes iterator variables independently");
    tests.expectEqual(eval(session, "table[1/(i+j+1),{{i,4},{j,4}}]"),
        std::string{"{{1/3, 1/4, 1/5, 1/6}, {1/4, 1/5, 1/6, 1/7}, {1/5, 1/6, 1/7, 1/8}, {1/6, 1/7, 1/8, 1/9}}"},
        "table accepts a brace value of iterator specifications as nested iteration");
    static_cast<void>(eval(session, "i:=99"));
    tests.expectEqual(eval(session, "table[i,{i,3}]"), std::string{"{1, 2, 3}"},
        "table iterator scope shadows a global definition locally");
    tests.expectEqual(eval(session, "i"), std::string{"99"},
        "table restores the outer definition after iteration");

    tests.expect(evalError(session, "table[log[0],{i,3}]").type() == error::CalcErrorType::Domain,
        "table propagates body errors without swallowing the original diagnostic");
    tests.expectEqual(eval(session, "i"), std::string{"99"},
        "table unwinds its local iterator scope when the body throws");

    // shared array shape layer and exact vector/matrix utilities.
    tests.expectEqual(eval(session, "identity[2]"), std::string{"{{1, 0}, {0, 1}}"},
        "identity constructs an exact matrix");
    tests.expectEqual(eval(session, "zeros[2,3]"), std::string{"{{0, 0, 0}, {0, 0, 0}}"},
        "zeros constructs the requested rectangular matrix");
    tests.expectEqual(eval(session, "mget[{{1,2},{3,4}},0,1]"), std::string{"2"},
        "mget uses zero-based row and column indices");
    tests.expectEqual(eval(session, "trace[{{1,2},{3,4}}]"), std::string{"5"},
        "trace is exact");
    tests.expectEqual(eval(session, "mtrace[{{1,2},{3,4}}]"), std::string{"5"},
        "legacy mtrace aliases trace");
    tests.expectEqual(eval(session, "rows[{{1,2},{3,4}}]"), std::string{"2"},
        "rows reports matrix shape");
    tests.expectEqual(eval(session, "cols[{{1,2},{3,4}}]"), std::string{"2"},
        "cols reports matrix shape");
    tests.expectEqual(eval(session, "diag[{{1,2,3},{4,5,6}}]"), std::string{"{1, 5}"},
        "diag supports rectangular matrices");
    tests.expectEqual(eval(session, "dimensions[{{1,2,3},{4,5,6}}]"),
        std::string{"{2, 3}"}, "dimensions exposes the stored shape");
    tests.expectEqual(eval(session, "arrayRank[{{1,2},{3,4}}]"), std::string{"2"},
        "arrayRank is distinct from matrixRank");
    tests.expectEqual(eval(session, "at[{{1,2},{3,4}},1,0]"), std::string{"3"},
        "at uses zero-based multi-dimensional indices");
    tests.expectEqual(eval(session, "at[{{1,2},{3,4}},1]"), std::string{"{3, 4}"},
        "at supports row-major prefix slicing for higher-level Array results");
    tests.expectEqual(eval(session, "reshape[{1,2,3,4},{2,2}]"),
        std::string{"{{1, 2}, {3, 4}}"}, "reshape preserves row-major element order");
    tests.expectEqual(eval(session, "dimensions[zeros[0,3]]"), std::string{"{0, 3}"},
        "zero-length leading dimensions preserve their trailing shape");
    tests.expectEqual(eval(session, "zeros[0,3]"), std::string{"reshape[{}, {0, 3}]"},
        "formatter preserves zero-length matrix shape through reshape");
    tests.expectEqual(eval(session, "reshape[{}, {0,3}]"),
        std::string{"reshape[{}, {0, 3}]"}, "zero-length reshape is round-trip stable");
    static_cast<void>(eval(session, "arrayPair[x]:={x,x+1}"));
    tests.expectEqual(eval(session, "{arrayPair[1],arrayPair[3]}"),
        std::string{"{{1, 2}, {3, 4}}"},
        "evaluated nested arrays are flattened into one rectangular Array");
    tests.expect(evalError(session, "{identity[1],2}").type() == error::CalcErrorType::Type,
        "evaluated arrays reject mixed scalar and array leaves");
    tests.expect(evalError(session, "{identity[1],identity[2]}").type() == error::CalcErrorType::Type,
        "evaluated arrays reject inconsistent child shapes");
    tests.expectEqual(eval(session, "{{1,2},{3,4}}+{{5,6},{7,8}}"),
        std::string{"{{6, 8}, {10, 12}}"}, "same-shape Array addition is elementwise");
    tests.expectEqual(eval(session, "2*{{1,2},{3,4}}"),
        std::string{"{{2, 4}, {6, 8}}"}, "Array multiplication accepts scalar factors only");
    tests.expect(evalError(session, "{{1,2},{3,4}}*{{5,6},{7,8}}").type()
            == error::CalcErrorType::Type,
        "Array-by-Array multiplication requires explicit dot");
    tests.expectEqual(eval(session, "vadd[{1,2},{3,4}]"), std::string{"{4, 6}"},
        "vector addition is elementwise and exact");
    tests.expect(evalError(session, "vscalar[{1,2},{}]").type() == error::CalcErrorType::Type,
        "vector scaling rejects an array-valued scale");
    tests.expectEqual(eval(session, "vdot[{1,2},{3,4}]"), std::string{"11"},
        "vector dot product is exact");
    tests.expectEqual(eval(session, "vcross[{1,0,0},{0,1,0}]"), std::string{"{0, 0, 1}"},
        "3D cross product is exact");
    tests.expectEqual(eval(session, "vnorm[{3,4}]"), std::string{"5"},
        "real vector norm preserves exact square roots");
    tests.expectEqual(eval(session, "vnormalize[{3,4}]"), std::string{"{3/5, 4/5}"},
        "vector normalization preserves exact rationals");
    tests.expectEqual(eval(session, "vproject[{1,2},{0,1}]"), std::string{"{0, 2}"},
        "vector projection is exact");
    tests.expectEqual(eval(session, "vangle[{1,0},{0,1}]"), std::string{"Pi/2"},
        "vector angle follows the default Radian semantics");
    tests.expectEqual(eval(session, "vreflect[{1,1},{0,1}]"), std::string{"{1, -1}"},
        "reflection about a normal vector is exact");
    tests.expectEqual(eval(session, "vreflect_axis[{1,1},{0,1}]"), std::string{"{-1, 1}"},
        "reflection about an axis vector is exact");
    tests.expectEqual(eval(session, "vsum[{1,2,3}]"), std::string{"6"},
        "vector element sum shares exact aggregation semantics");
    tests.expectEqual(eval(session, "{1,2}/2"), std::string{"{1/2, 1}"},
        "Array division supports scalar scaling without introducing elementwise Array division");
    tests.expect(evalError(session, "2/{1,2}").type() == error::CalcErrorType::Type,
        "scalar-by-Array division remains undefined");
    tests.expectEqual(eval(session, "inner[{1,I},{1,I}]"), std::string{"2"},
        "inner is Hermitian and conjugates its first vector argument");
    tests.expectEqual(eval(session, "outer[{1,2},{3,4,5}]"),
        std::string{"{{3, 4, 5}, {6, 8, 10}}"},
        "outer constructs the bilinear rank-one matrix exactly");
    tests.expectEqual(eval(session, "distance[{1,I},{2,0}]"), std::string{"sqrt[2]"},
        "distance routes vector distance through the Hermitian norm");
    tests.expectEqual(eval(session, "projection[{1,I},{1,0}]"), std::string{"{1, 0}"},
        "projection uses the Hermitian projection convention for complex vectors");
    tests.expectEqual(eval(session, "vproject[{1,I},{1,0}]"), std::string{"{1, 0}"},
        "legacy vproject routes to the canonical Hermitian projection");
    tests.expectEqual(eval(session, "vmanhattan[{1,I},{2,0}]"), std::string{"2"},
        "Manhattan distance extends naturally to complex components through abs");

    tests.expectEqual(eval(session, "grad[x^2+y^2+z^2,{x,y,z}]"),
        std::string{"{2x, 2y, 2z}"},
        "grad computes a Cartesian gradient over explicit coordinate symbols");
    tests.expectEqual(eval(session, "divergence[{x^2,y^2,z^2},{x,y,z}]"),
        std::string{"2x+2y+2z"},
        "divergence contracts field components with matching Cartesian partial derivatives");
    tests.expectEqual(eval(session, "curl[{y*z,z*x,x*y},{x,y,z}]"),
        std::string{"{0, 0, 0}"},
        "curl computes the exact three-dimensional Cartesian rotation");
    tests.expectEqual(eval(session, "laplacian[x^2+y^2+z^2,{x,y,z}]"), std::string{"6"},
        "laplacian sums exact second Cartesian partial derivatives");
    tests.expectEqual(eval(session, "jacobian[{x*y,sin[z]},{x,y,z}]"),
        std::string{"{{y, x, 0}, {0, 0, cos[z]}}"},
        "jacobian returns one derivative row per vector-field component");
    tests.expectEqual(eval(session, "hessian[x^2+x*y+y^2,{x,y}]"),
        std::string{"{{2, 1}, {1, 2}}"},
        "hessian returns the exact square second-derivative matrix");
    tests.expectEqual(eval(session, "divergence[grad[x^2+y^2+z^2,{x,y,z}],{x,y,z}]"),
        std::string{"6"},
        "vector-calculus operators materialize nested held operators without evaluating coordinates");
    tests.expectEqual(eval(session, "curl[grad[x*y*z,{x,y,z}],{x,y,z}]"),
        std::string{"{0, 0, 0}"},
        "curl of an exact gradient reduces to zero");
    tests.expect(evalError(session, "curl[{x,y},{x,y}]").type() == error::CalcErrorType::Domain,
        "curl rejects non-three-dimensional Cartesian fields");
    tests.expect(evalError(session, "grad[x^2,{x,x}]").type() == error::CalcErrorType::Domain,
        "vector-calculus coordinate lists require distinct symbols");

    // stable elementary/cardinal functions remain explicit until certified evaluation.
    tests.expectEqual(eval(session, "expm1[0]"), std::string{"0"},
        "expm1 has an exact zero value");
    tests.expectEqual(eval(session, "log1p[0]"), std::string{"0"},
        "log1p has an exact zero value");
    tests.expectEqual(eval(session, "sinc[0]"), std::string{"1"},
        "sinc fills its removable singularity exactly");
    tests.expectEqual(eval(session, "cosc[0]"), std::string{"0"},
        "cosc fills its removable singularity exactly");
    tests.expectEqual(eval(session, "tanc[0]"), std::string{"1"},
        "tanc fills its removable singularity exactly");
    tests.expectEqual(eval(session, "sinhc[0]"), std::string{"1"},
        "sinhc fills its removable singularity exactly");
    tests.expectEqual(eval(session, "tanhc[0]"), std::string{"1"},
        "tanhc fills its removable singularity exactly");
    tests.expectEqual(eval(session, "expc[0]"), std::string{"1"},
        "expc fills its removable singularity exactly");
    tests.expectEqual(eval(session, "sinc[N[0,5]]"), std::string{"1.0"},
        "sinc certifies a finite-precision interval across its removable singularity");
    tests.expectEqual(eval(session, "cosc[N[0,5]]"), std::string{"0.0"},
        "cosc certifies a finite-precision interval across its removable singularity");
    tests.expectEqual(eval(session, "tanc[N[0,5]]"), std::string{"1.0"},
        "tanc certifies a finite-precision interval across its removable singularity");
    tests.expectEqual(eval(session, "sinhc[N[0,5]]"), std::string{"1.0"},
        "sinhc certifies a finite-precision interval across its removable singularity");
    tests.expectEqual(eval(session, "tanhc[N[0,5]]"), std::string{"1.0"},
        "tanhc certifies a finite-precision interval across its removable singularity");
    tests.expectEqual(eval(session, "expc[N[0,5]]"), std::string{"1.0"},
        "expc certifies a finite-precision interval across its removable singularity");
    tests.expectEqual(eval(session, "sinc[Pi/2]"), std::string{"2/Pi"},
        "sinc reduces exact special angles without approximation");
    tests.expectEqual(eval(session, "tanc[Pi/4]"), std::string{"4/Pi"},
        "tanc preserves an exact special-angle result");
    tests.expect(evalError(session, "tanc[Pi/2]").type() == error::CalcErrorType::Domain,
        "tanc detects an exact tangent pole before numerical evaluation");
    tests.expectEqual(eval(session, "N[sinc[Pi/2],20]"),
        std::string{"0.63661977236758134308"},
        "sinc uses the default Radian argument");
    tests.expectEqual(eval(session, "N[sinc[90 Deg],20]"),
        std::string{"0.63661977236758134308"},
        "sinc explicit Degree input is radian-normalized consistently");
    tests.expectEqual(eval(session, "N[cosc[Pi/3],20]"),
        std::string{"0.47746482927568600731"},
        "cosc is certified without early floating conversion");
    tests.expectEqual(eval(session, "N[tanc[Pi/4],20]"),
        std::string{"1.2732395447351626862"},
        "tanc has certified pole-aware evaluation");
    tests.expectEqual(eval(session, "D[sinc[x],x]"),
        std::string{"cases[(x cos[x]-sin[x])/x^2 if x != 0; 0 if x == 0]"},
        "sinc derivative preserves the removable singularity at zero");
    tests.expectEqual(eval(session, "D[cosc[x],x]"),
        std::string{"cases[(-1+x sin[x]+cos[x])/x^2 if x != 0; 1/2 if x == 0]"},
        "cosc derivative uses its exact one-half limit at zero");
    tests.expectEqual(eval(session, "D[tanc[x],x]"),
        std::string{"cases[(x sec[x]^2-tan[x])/x^2 if x != 0; 0 if x == 0]"},
        "tanc derivative preserves the removable singularity at zero");
    tests.expectEqual(eval(session, "D[sinhc[x],x]"),
        std::string{"cases[(x cosh[x]-sinh[x])/x^2 if x != 0; 0 if x == 0]"},
        "sinhc derivative preserves the removable singularity at zero");
    tests.expectEqual(eval(session, "D[tanhc[x],x]"),
        std::string{"cases[(x sech[x]^2-tanh[x])/x^2 if x != 0; 0 if x == 0]"},
        "tanhc derivative preserves the removable singularity at zero");
    tests.expectEqual(eval(session, "D[expc[x],x]"),
        std::string{"cases[(1+x exp[x]-exp[x])/x^2 if x != 0; 1/2 if x == 0]"},
        "expc derivative uses its exact one-half limit at zero");
    tests.expectEqual(eval(session, "N[expm1[1/10^30],50]"),
        std::string{"0.00000000000000000000000000000100000000000000000000000000000050"},
        "expm1 survives severe cancellation by precision refinement");
    tests.expectEqual(eval(session, "N[log1p[1/10^30],50]"),
        std::string{"0.00000000000000000000000000000099999999999999999999999999999950"},
        "log1p survives severe cancellation by precision refinement");
    tests.expect(evalError(session, "log1p[-1]").type() == error::CalcErrorType::Domain,
        "log1p rejects its exact branch singularity");
    tests.expectEqual(eval(session, "solve[log1p[y]*x+1==0,x]"),
        std::string{"cases[{x == -1/log1p[y]} if log1p[y] != 0; {} if log1p[y] == 0] if y != -1"},
        "Solver receives log1p definedness from MathRegistry");

    // Angle migration and generalized suffix syntax.
    tests.expectEqual(eval(session, "sin[Pi/6]"), std::string{"1/2"},
        "Angle semantics: Radian is the session default");
    tests.expectEqual(eval(session, "sin[Pi/6 Rad]"), std::string{"1/2"},
        "Angle semantics: a suffix can apply to a compound expression");
    tests.expectEqual(eval(session, "sin[30 Deg]"), std::string{"1/2"},
        "Angle semantics: Degree remains explicitly selectable");
    tests.expectEqual(eval(session, "sin[100 Grad]"), std::string{"1"},
        "Angle semantics: Gradian remains explicitly selectable");

    kernel::KernelSession degreeSession;
    degreeSession.setDefaultAngleUnit(mathematics::AngleUnit::Degree);
    tests.expectEqual(eval(degreeSession, "sin[30]"), std::string{"1/2"},
        "Angle semantics: a session can still select Degree as its default");
    kernel::KernelSession gradianSession;
    gradianSession.setDefaultAngleUnit(mathematics::AngleUnit::Gradian);
    tests.expectEqual(eval(gradianSession, "sin[100]"), std::string{"1"},
        "Angle semantics: a session can still select Gradian as its default");
}

} // namespace mmcal::tests
