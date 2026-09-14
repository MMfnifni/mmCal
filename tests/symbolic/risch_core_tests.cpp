// Risch core の exact certificate 回帰テスト
#include "risch_core_tests.hpp"

#include "expression/expr.hpp"
#include "evaluation/builtin_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "symbolic/differential_tower.hpp"
#include "symbolic/risch_core.hpp"
#include "symbolic/risch_differential_equation.hpp"
#include "symbolic/risch_differential_reduction.hpp"
#include "symbolic/risch_expression.hpp"
#include "symbolic/risch_tower_recognizer.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <utility>
#include <vector>

namespace mmcal::tests {
namespace {

using numeric::BigInt;
using numeric::Rational;
using symbolic::RationalPolynomial;

[[nodiscard]] RationalPolynomial polynomial(
    std::initializer_list<std::int64_t> coefficients) {
    std::vector<Rational> exact;
    exact.reserve(coefficients.size());
    for (const std::int64_t coefficient : coefficients)
        exact.emplace_back(BigInt{coefficient});
    return RationalPolynomial{std::move(exact)};
}

[[nodiscard]] RationalPolynomial multiply(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    std::vector<Rational> result(
        lhs.degree() + rhs.degree() + 1, Rational{BigInt{0}});
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            result[i + j] += lhs.coefficient(i) * rhs.coefficient(j);
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial add(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    const std::size_t size = std::max(
        lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<Rational> result(size, Rational{BigInt{0}});
    for (std::size_t i = 0; i < size; ++i)
        result[i] = lhs.coefficient(i) + rhs.coefficient(i);
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] RationalPolynomial scale(
    const RationalPolynomial& value,
    std::int64_t multiplier) {
    std::vector<Rational> result = value.coefficients();
    for (Rational& coefficient : result)
        coefficient *= Rational{BigInt{multiplier}};
    return RationalPolynomial{std::move(result)};
}

[[nodiscard]] bool same(
    const RationalPolynomial& lhs,
    const RationalPolynomial& rhs) {
    return lhs.coefficients() == rhs.coefficients();
}

[[nodiscard]] symbolic::risch::RationalFunction rationalFunction(
    std::initializer_list<std::int64_t> numerator,
    std::initializer_list<std::int64_t> denominator = {1}) {
    return {polynomial(numerator), polynomial(denominator)};
}

void runSharedPolynomialKernelTests(TestRunner& tests) {
    const RationalPolynomial value = polynomial({1, 2, 3});
    tests.expect(
        same(
            symbolic::risch::shiftRationalPolynomialExact(value, 2),
            polynomial({17, 14, 3})),
        "shared Q[x] kernel: exact positive polynomial shift");
    tests.expect(
        same(
            symbolic::risch::shiftRationalPolynomialExact(value, -1),
            polynomial({2, -4, 3})),
        "shared Q[x] kernel: exact negative polynomial shift");
}

[[nodiscard]] std::size_t countCalls(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    evaluation::BuiltinId id) {
    std::size_t count = builtins.isCallTo(expression, id) ? 1 : 0;
    if (expression.isCall())
        for (const expression::Expr& argument : expression.asCall().arguments)
            count += countCalls(argument, builtins, id);
    return count;
}

void runDifferentialTowerTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const expression::Symbol x = symbols.intern("x");
    const expression::Symbol t = symbols.intern("t");
    const expression::Symbol e = symbols.intern("e");
    symbolic::risch::DifferentialTower tower{x, 2};

    const auto primitive = tower.appendPrimitive(
        t, expression::Expr{x},
        expression::Expr{numeric::Number{BigInt{1}}});
    tests.expect(
        primitive && primitive.level == 1 && tower.depth() == 1,
        "Risch tower: appends a primitive extension with a stable level");

    const auto exponential = tower.appendExponential(
        e, expression::Expr{t}, expression::Expr{t});
    tests.expect(
        exponential && exponential.level == 2
            && tower.extensions().back().kind()
                == symbolic::risch::DifferentialExtensionKind::Exponential,
        "Risch tower: appends an exponential extension over the lower field");

    const auto duplicate = tower.appendPrimitive(
        t, expression::Expr{x}, expression::Expr{x});
    tests.expect(
        duplicate.error
            == symbolic::risch::DifferentialTowerError::DuplicateGenerator,
        "Risch tower: rejects duplicate generators");

    symbolic::risch::DifferentialTower selfDependent{x};
    const auto self = selfDependent.appendPrimitive(
        t, expression::Expr{x}, expression::Expr{t});
    tests.expect(
        self.error
            == symbolic::risch::DifferentialTowerError::SelfDependentDefinition,
        "Risch tower: rejects a derivative outside the lower field");
}

void runHermiteTests(TestRunner& tests) {
    const RationalPolynomial quadratic = polynomial({1, 0, 1});
    const RationalPolynomial denominator = multiply(
        multiply(quadratic, quadratic), quadratic);
    const symbolic::risch::RationalFunction input{
        polynomial({1}), denominator};
    const auto reduced = symbolic::risch::hermiteReduceRationalFunction(input);
    tests.expect(
        reduced && reduced.value->exactVerified
            && symbolic::risch::verifyHermiteReduction(input, *reduced.value),
        "Risch Hermite: verifies the derivative identity exactly");
    tests.expect(
        reduced && reduced.value->squareFreePart.denominator.degree() == 2,
        "Risch Hermite: lowers a repeated cubic power to a square-free denominator");

    const RationalPolynomial linear = polynomial({-1, 1});
    const RationalPolynomial mixedDenominator = multiply(
        multiply(linear, linear), denominator);
    const symbolic::risch::RationalFunction mixed{
        polynomial({2, 0, 0, 0, 1}), mixedDenominator};
    const auto mixedReduction = symbolic::risch::hermiteReduceRationalFunction(mixed);
    tests.expect(
        mixedReduction && mixedReduction.value->exactVerified,
        "Risch Hermite: combines CRT blocks without losing its certificate");
}

void runLrtTests(TestRunner& tests) {
    const symbolic::risch::RationalFunction linear{
        polynomial({3}), polynomial({-2, 1})};
    const auto linearLrt = symbolic::risch::lazardRiobooTrager(linear);
    tests.expect(
        linearLrt && linearLrt.value->exactVerified
            && linearLrt.value->logarithmicTerms.size() == 1
            && same(
                linearLrt.value->logarithmicTerms.front().residuePolynomial,
                polynomial({-3, 1})),
        "Risch LRT: handles a constant residue equation for a linear pole");

    // D=(x^2-1)(x^2-4), residues 1 at +/-1 and 2 at +/-2.
    const symbolic::risch::RationalFunction grouped{
        polynomial({0, -12, 0, 6}),
        polynomial({4, 0, -5, 0, 1})};
    const auto lrt = symbolic::risch::lazardRiobooTrager(grouped);
    tests.expect(
        lrt && lrt.value->exactVerified
            && symbolic::risch::verifyLrtResult(grouped, *lrt.value),
        "Risch LRT: verifies resultant and quotient-ring divisibility");
    tests.expect(
        lrt && lrt.value->logarithmicTerms.size() == 1
            && lrt.value->logarithmicTerms.front().poleMultiplicity == 2
            && same(
                lrt.value->logarithmicTerms.front().residuePolynomial,
                polynomial({2, -3, 1})),
        "Risch LRT: groups conjugate equal-multiplicity residues compactly");

    bool generatedFamiliesVerified = true;
    for (std::int64_t degree = 2; degree <= 6; ++degree) {
        std::vector<RationalPolynomial> factors;
        RationalPolynomial generatedDenominator = polynomial({1});
        for (std::int64_t root = 1; root <= degree; ++root) {
            factors.push_back(polynomial({-root, 1}));
            generatedDenominator = multiply(
                generatedDenominator, factors.back());
        }
        RationalPolynomial generatedNumerator;
        for (std::int64_t omitted = 0; omitted < degree; ++omitted) {
            RationalPolynomial cofactor = polynomial({1});
            for (std::int64_t index = 0; index < degree; ++index)
                if (index != omitted)
                    cofactor = multiply(
                        cofactor, factors[static_cast<std::size_t>(index)]);
            generatedNumerator = add(
                generatedNumerator,
                scale(cofactor, 1 + omitted % 3));
        }
        const symbolic::risch::RationalFunction generated{
            generatedNumerator, generatedDenominator};
        const auto generatedLrt =
            symbolic::risch::lazardRiobooTrager(generated);
        generatedFamiliesVerified = generatedFamiliesVerified
            && generatedLrt && generatedLrt.value->exactVerified
            && symbolic::risch::verifyLrtResult(
                generated, *generatedLrt.value);
    }
    tests.expect(
        generatedFamiliesVerified,
        "Risch LRT: certifies generated equal-residue families through degree six");

    // Bronstein の LRT 例。resultant は monic 化後 (z^2+1/4)^3。
    const symbolic::risch::RationalFunction bronstein{
        polynomial({6, 0, -3, 0, 1}),
        polynomial({4, 0, 5, 0, -5, 0, 1})};
    const auto algebraic = symbolic::risch::lazardRiobooTrager(bronstein);
    const Rational quarter{BigInt{1}, BigInt{4}};
    tests.expect(
        algebraic && algebraic.value->logarithmicTerms.size() == 1
            && algebraic.value->logarithmicTerms.front().poleMultiplicity == 3
            && same(
                algebraic.value->logarithmicTerms.front().residuePolynomial,
                RationalPolynomial{{quarter, Rational{BigInt{0}}, Rational{BigInt{1}}}}),
        "Risch LRT: represents algebraic residues without expanding individual roots");

    const symbolic::risch::RischResult integrated =
        symbolic::risch::integrateRationalRisch(bronstein);
    tests.expect(
        integrated.status == symbolic::risch::RischResultStatus::Elementary
            && integrated.exactVerified && integrated.rational.has_value(),
        "Risch result: marks a certified rational integral elementary");

    symbolic::risch::RischOptions limited;
    limited.maximumLrtDegree = 3;
    const symbolic::risch::RischResult stopped =
        symbolic::risch::integrateRationalRisch(bronstein, limited);
    tests.expect(
        stopped.status == symbolic::risch::RischResultStatus::ResourceLimit
            && stopped.failure == symbolic::risch::RischFailure::DegreeLimit,
        "Risch result: reports a resource limit instead of non-elementarity");

    symbols::SymbolTable symbols;
    const evaluation::BuiltinRegistry builtins =
        evaluation::BuiltinRegistry::defaults(symbols);
    const mathematics::MathRegistry mathRegistry =
        mathematics::MathRegistry::defaults(symbols, builtins);
    const mathematics::AngleSemantics angles;
    const expression::Symbol x = symbols.intern("x");
    const auto materialized = symbolic::risch::materializeLrtLogarithms(
        bronstein, *algebraic.value, x, builtins, mathRegistry, angles);
    tests.expect(
        materialized
            && countCalls(*materialized.value, builtins, evaluation::BuiltinId::Log)
                == 2,
        "Risch LRT: materializes one logarithm per algebraic residue, not per pole");
    symbolic::risch::LrtResult uncertified = *algebraic.value;
    uncertified.exactVerified = false;
    const auto rejected = symbolic::risch::materializeLrtLogarithms(
        bronstein, uncertified, x, builtins, mathRegistry, angles);
    tests.expect(
        !rejected
            && rejected.failure == symbolic::risch::RischFailure::CertificateFailed,
        "Risch LRT: refuses to materialize an uncertified residue decomposition");
}

void runTowerRecognizerTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const evaluation::BuiltinRegistry builtins =
        evaluation::BuiltinRegistry::defaults(symbols);
    const mathematics::MathRegistry mathRegistry =
        mathematics::MathRegistry::defaults(symbols, builtins);
    const mathematics::AngleSemantics angles;
    const expression::Symbol x = symbols.intern("x");
    const expression::Expr logX = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Log), {expression::Expr{x}});
    const expression::Expr exponential = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Exp), {
            expression::Expr::call(
                builtins.symbol(evaluation::BuiltinId::Multiply),
                {expression::Expr{x}, logX})});
    const expression::Expr input = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Add), {logX, exponential});
    const auto recognized = symbolic::risch::recognizeDifferentialTower(
        input, x, symbols, builtins, mathRegistry, angles);
    tests.expect(
        recognized.complete() && recognized.tower.depth() == 2
            && recognized.tower.extensions()[0].kind()
                == symbolic::risch::DifferentialExtensionKind::Primitive
            && recognized.tower.extensions()[1].kind()
                == symbolic::risch::DifferentialExtensionKind::Exponential,
        "Risch tower recognizer: builds nested Log/Exp extensions bottom-up");
    tests.expect(
        symbolic::containsSymbol(
            recognized.tower.extensions()[1].differentialCoefficient(),
            recognized.tower.extensions()[0].generator())
            && !symbolic::containsSymbol(
                recognized.tower.extensions()[1].differentialCoefficient(),
                recognized.tower.extensions()[1].generator()),
        "Risch tower recognizer: rewrites logarithmic derivatives into the lower field");

    const expression::Expr duplicated = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Add), {logX, logX});
    const auto reused = symbolic::risch::recognizeDifferentialTower(
        duplicated, x, symbols, builtins, mathRegistry, angles);
    tests.expect(
        reused.complete() && reused.tower.depth() == 1,
        "Risch tower recognizer: reuses structurally equal generators");

    const expression::Expr unsupported = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Log), {
            expression::Expr::call(
                builtins.symbol(evaluation::BuiltinId::Sin),
                {expression::Expr{x}})});
    const auto partial = symbolic::risch::recognizeDifferentialTower(
        unsupported, x, symbols, builtins, mathRegistry, angles);
    tests.expect(
        partial.status
            == symbolic::risch::DifferentialTowerRecognitionStatus::Partial
            && partial.tower.depth() == 0,
        "Risch tower recognizer: leaves non-rational derivatives unsupported");

    symbolic::risch::DifferentialTowerRecognitionOptions oneLevel;
    oneLevel.maximumTowerDepth = 1;
    const auto limited = symbolic::risch::recognizeDifferentialTower(
        exponential, x, symbols, builtins, mathRegistry, angles, oneLevel);
    tests.expect(
        limited.status
            == symbolic::risch::DifferentialTowerRecognitionStatus::ResourceLimit
            && limited.tower.depth() == 1,
        "Risch tower recognizer: reports its depth boundary without discarding lower levels");
}

void runDifferentialReductionTests(TestRunner& tests) {
    using symbolic::risch::DifferentialExtensionKind;
    using symbolic::risch::DifferentialPolynomial;
    const symbolic::risch::DifferentialDerivation primitiveLog{
        DifferentialExtensionKind::Primitive,
        rationalFunction({1}, {0, 1})};
    const DifferentialPolynomial t{{
        rationalFunction({0}), rationalFunction({1})}};
    const auto normal = symbolic::risch::classifyDifferentialPolynomial(
        t, primitiveLog);
    tests.expect(
        normal && normal.value->exactVerified
            && normal.value->classification
                == symbolic::risch::DifferentialPolynomialClass::Normal,
        "Risch differential reduction: classifies Log[x] as a normal generator");

    const auto reduced = symbolic::risch::hermiteReduceNormalDifferentialPower(
        DifferentialPolynomial{{rationalFunction({1})}}, t, 2, primitiveLog);
    tests.expect(
        reduced && reduced.value->steps == 1
            && reduced.value->exactVerified
            && symbolic::risch::verifyNormalDifferentialReduction(
                DifferentialPolynomial{{rationalFunction({1})}},
                t, 2, primitiveLog, *reduced.value),
        "Risch differential reduction: verifies a primitive normal Hermite step");

    const symbolic::risch::DifferentialDerivation exponential{
        DifferentialExtensionKind::Exponential,
        rationalFunction({1})};
    const auto special = symbolic::risch::classifyDifferentialPolynomial(
        t, exponential);
    tests.expect(
        special && special.value->classification
            == symbolic::risch::DifferentialPolynomialClass::Special,
        "Risch differential reduction: classifies an exponential generator as special");

    const DifferentialPolynomial mixedFactor{{
        rationalFunction({0}), rationalFunction({1}), rationalFunction({1})}};
    const auto mixed = symbolic::risch::classifyDifferentialPolynomial(
        mixedFactor, exponential);
    tests.expect(
        mixed && mixed.value->classification
            == symbolic::risch::DifferentialPolynomialClass::Mixed,
        "Risch differential reduction: detects mixed normal/special factors");

    const auto stopped = symbolic::risch::hermiteReduceNormalDifferentialPower(
        DifferentialPolynomial{{rationalFunction({1})}}, t, 2, exponential);
    tests.expect(
        !stopped
            && stopped.failure
                == symbolic::risch::RischFailure::NonNormalDifferentialFactor,
        "Risch differential reduction: does not misclassify special factors as non-elementary");

    symbolic::risch::RischOptions limited;
    limited.maximumDifferentialOperations = 1;
    const auto budgetStop = symbolic::risch::classifyDifferentialPolynomial(
        t, primitiveLog, limited);
    tests.expect(
        !budgetStop
            && budgetStop.failure
                == symbolic::risch::RischFailure::DifferentialOperationLimit,
        "Risch differential reduction: bounds exact coefficient-field work explicitly");
}

void runRationalDifferentialEquationTests(TestRunner& tests) {
    using symbolic::risch::DifferentialPolynomial;
    using symbolic::risch::ExponentialLaurentTerm;
    const auto polynomialRde =
        symbolic::risch::solveRationalDifferentialEquation(
            rationalFunction({1}), rationalFunction({0, 1}));
    tests.expect(
        polynomialRde && polynomialRde.value->exactVerified
            && symbolic::risch::equivalentRationalFunctions(
                polynomialRde.value->solution,
                rationalFunction({-1, 1})),
        "Risch RDE: solves a polynomial hyperexponential coefficient equation exactly");

    const auto repeatedPole =
        symbolic::risch::solveRationalDifferentialEquation(
            rationalFunction({0}), rationalFunction({1}, {0, 0, 1}));
    tests.expect(
        repeatedPole && repeatedPole.value->exactVerified
            && symbolic::risch::equivalentRationalFunctions(
                repeatedPole.value->solution,
                rationalFunction({-1}, {0, 1})),
        "Risch RDE: derives and verifies a repeated-pole normal denominator");

    const auto weaklyNormalized =
        symbolic::risch::solveRationalDifferentialEquation(
            rationalFunction({1}, {0, 1}), rationalFunction({1}));
    const Rational half{BigInt{1}, BigInt{2}};
    tests.expect(
        weaklyNormalized && weaklyNormalized.value->exactVerified
            && weaklyNormalized.value->weakNormalizer.degree() == 1
            && symbolic::risch::equivalentRationalFunctions(
                weaklyNormalized.value->solution,
                symbolic::risch::RationalFunction{
                    RationalPolynomial{{Rational{BigInt{0}}, half}},
                    polynomial({1})}),
        "Risch RDE: removes a positive integral residue before degree bounding");

    const auto noSolution =
        symbolic::risch::solveRationalDifferentialEquation(
            rationalFunction({1}), rationalFunction({1}, {0, 1}));
    tests.expect(
        !noSolution
            && noSolution.failure
                == symbolic::risch::RischFailure::NoRationalSolution,
        "Risch RDE: proves an inconsistent bounded coefficient system has no rational solution");

    const std::vector<symbolic::risch::RationalFunction> coefficients{
        rationalFunction({0}),
        rationalFunction({1}),
        rationalFunction({1}, {0, 1}),
        rationalFunction({2}, {0, 1}),
        rationalFunction({-1}, {0, 1}),
        rationalFunction({0, 1}, {1, 0, 1}),
        rationalFunction({0, 2}, {1, 0, 1}),
        rationalFunction({1}, {0, 0, 1})};
    const std::vector<symbolic::risch::RationalFunction> knownSolutions{
        rationalFunction({1, 1}),
        rationalFunction({1}, {0, 1}),
        rationalFunction({1}, {1, 0, 1}),
        rationalFunction({1, 1}, {1, 0, 1})};
    bool generatedRdesVerified = true;
    for (const auto& coefficient : coefficients) {
        for (const auto& knownSolution : knownSolutions) {
            const auto rightHandSide =
                symbolic::risch::addRationalFunctionsExact(
                    symbolic::risch::differentiateRationalFunctionExact(
                        knownSolution),
                    symbolic::risch::multiplyRationalFunctionsExact(
                        coefficient, knownSolution));
            const auto generated =
                symbolic::risch::solveRationalDifferentialEquation(
                    coefficient, rightHandSide);
            generatedRdesVerified = generatedRdesVerified
                && generated && generated.value->exactVerified
                && symbolic::risch::verifyRationalDifferentialEquation(
                    coefficient, rightHandSide, *generated.value);
        }
    }
    tests.expect(
        generatedRdesVerified,
        "Risch RDE: recovers and certifies a generated rational-solution grid");

    symbolic::risch::RischOptions matrixLimited;
    matrixLimited.maximumRdeMatrixEntries = 1;
    const auto matrixStop =
        symbolic::risch::solveRationalDifferentialEquation(
            rationalFunction({1}), rationalFunction({0, 1}), matrixLimited);
    tests.expect(
        !matrixStop
            && matrixStop.failure
                == symbolic::risch::RischFailure::RdeMatrixSizeLimit,
        "Risch RDE: reports its exact linear-system memory boundary");

    const auto limitedIntegral =
        symbolic::risch::limitedIntegrateRationalFunction(
            rationalFunction({2}, {0, 1}), rationalFunction({1}, {0, 1}));
    tests.expect(
        limitedIntegral && limitedIntegral.value->exactVerified
            && limitedIntegral.value->constantMultiple
                == Rational{BigInt{2}}
            && limitedIntegral.value->rationalPart.numerator.isZero(),
        "Risch primitive case: performs exact limited integration in Q(x)");

    const symbolic::risch::RationalFunction logarithmicDerivative =
        rationalFunction({1}, {0, 1});
    const symbolic::risch::DifferentialPolynomial primitiveInput{{
        rationalFunction({0}), logarithmicDerivative}};
    const auto primitive = symbolic::risch::reducePrimitivePolynomial(
        primitiveInput, logarithmicDerivative);
    tests.expect(
        primitive && primitive.value->exactVerified
            && primitive.value->polynomialPart.degree() == 2
            && primitive.value->lowerFieldRemainder.numerator.isZero(),
        "Risch primitive case: lowers a polynomial by limited leading-coefficient integration");

    using symbolic::risch::DifferentialRationalFunction;
    const DifferentialRationalFunction repeatedPrimitiveDenominator{
        DifferentialPolynomial{{
            rationalFunction({0}), rationalFunction({1})}},
        DifferentialPolynomial{{
            rationalFunction({1}), rationalFunction({2}),
            rationalFunction({1})}}};
    const auto primitiveRational =
        symbolic::risch::reducePrimitiveRationalFunction(
            repeatedPrimitiveDenominator, logarithmicDerivative);
    tests.expect(
        primitiveRational && primitiveRational.value->exactVerified
            && primitiveRational.value->rationalPart.size() == 1
            && primitiveRational.value->logarithmicPart.empty()
            && primitiveRational.value->residualPart.empty()
            && symbolic::risch::verifyPrimitiveRationalReduction(
                repeatedPrimitiveDenominator, logarithmicDerivative,
                *primitiveRational.value),
        "Risch primitive rational case: reduces a repeated polynomial denominator exactly");

    const DifferentialRationalFunction primitiveLogarithmicDerivative{
        DifferentialPolynomial{{
            rationalFunction({1}, {0, 1}),
            rationalFunction({2}, {0, 1})}},
        DifferentialPolynomial{{
            rationalFunction({1}), rationalFunction({1}),
            rationalFunction({1})}}};
    const auto primitiveLogarithm =
        symbolic::risch::reducePrimitiveRationalFunction(
            primitiveLogarithmicDerivative, logarithmicDerivative);
    tests.expect(
        primitiveLogarithm && primitiveLogarithm.value->exactVerified
            && primitiveLogarithm.value->logarithmicPart.size() == 1
            && primitiveLogarithm.value->logarithmicPart.front().coefficient
                == Rational{BigInt{1}}
            && primitiveLogarithm.value->residualPart.empty(),
        "Risch primitive rational case: recognizes an exact logarithmic derivative residue");

    const DifferentialRationalFunction unsupportedPrimitiveResidue{
        DifferentialPolynomial{{rationalFunction({1})}},
        DifferentialPolynomial{{
            rationalFunction({1}), rationalFunction({0}),
            rationalFunction({1})}}};
    const auto retainedPrimitiveResidue =
        symbolic::risch::reducePrimitiveRationalFunction(
            unsupportedPrimitiveResidue, logarithmicDerivative);
    tests.expect(
        retainedPrimitiveResidue
            && retainedPrimitiveResidue.value->exactVerified
            && retainedPrimitiveResidue.value->logarithmicPart.empty()
            && retainedPrimitiveResidue.value->residualPart.size() == 1,
        "Risch primitive rational case: retains an unproved residue instead of claiming non-elementarity");

    const std::vector<ExponentialLaurentTerm> laurentInput{
        {1, rationalFunction({0, 1})},
        {-1, rationalFunction({1})},
        {0, rationalFunction({3}, {1, 1})}};
    const auto exponential =
        symbolic::risch::reduceExponentialLaurentPolynomial(
            laurentInput, rationalFunction({1}));
    tests.expect(
        exponential && exponential.value->exactVerified
            && exponential.value->laurentPart.size() == 2
            && symbolic::risch::verifyExponentialLaurentReduction(
                laurentInput, rationalFunction({1}), *exponential.value),
        "Risch exponential case: reduces positive and negative Laurent powers through exact RDEs");
}

} // namespace

void runRischCoreTests(TestRunner& tests) {
    runSharedPolynomialKernelTests(tests);
    runDifferentialTowerTests(tests);
    runHermiteTests(tests);
    runLrtTests(tests);
    runTowerRecognizerTests(tests);
    runDifferentialReductionTests(tests);
    runRationalDifferentialEquationTests(tests);
}

} // namespace mmcal::tests
