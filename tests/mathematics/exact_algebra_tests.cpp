// exact代数計算の回帰テスト
#include "exact_algebra_tests.hpp"

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "mathematics/exact_algebra.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "symbols/symbol_table.hpp"
#include "symbolic/algebraic_number.hpp"
#include "symbolic/number_field.hpp"
#include "test_framework.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <thread>
#include <vector>

namespace mmcal::tests {
namespace {

expression::Expr integer(std::int64_t value) {
    return expression::Expr{numeric::Number{numeric::BigInt{value}}};
}

expression::Expr sqrtExpr(
    std::int64_t value,
    const evaluation::BuiltinRegistry& builtins) {
    return expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Sqrt),
        {integer(value)});
}

} // namespace

void runExactAlgebraTests(TestRunner& tests) {
    symbols::SymbolTable symbols;
    const evaluation::BuiltinRegistry builtins = evaluation::BuiltinRegistry::defaults(symbols);

    const expression::Expr difference = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Subtract),
        {sqrtExpr(6, builtins), sqrtExpr(2, builtins)});
    const expression::Expr quarter = expression::Expr::call(
        builtins.symbol(evaluation::BuiltinId::Divide),
        {difference, integer(4)});

    const expression::Expr doubled = mathematics::scaleExactExpression(
        numeric::Rational{numeric::BigInt{2}},
        quarter,
        builtins);
    tests.expectEqual(
        formatting::formatExpr(doubled),
        std::string{"(sqrt[6]-sqrt[2])/2"},
        "ExactAlgebra: rational scaling reduces an existing exact denominator");

    const std::vector<expression::Expr> duplicateTerms{quarter, quarter};
    const std::vector<expression::Expr> combined =
        mathematics::combineStructurallyIdenticalTerms(duplicateTerms, builtins);
    tests.expect(combined.size() == 1,
        "ExactAlgebra: structurally identical terms are collected");
    tests.expectEqual(
        formatting::formatExpr(combined.front()),
        std::string{"(sqrt[6]-sqrt[2])/2"},
        "ExactAlgebra: collected identical radical terms keep exact form");

    const expression::Expr sqrt2 = sqrtExpr(2, builtins);
    const expression::Expr sqrt3 = sqrtExpr(3, builtins);
    const std::vector<expression::Expr> distinctTerms{sqrt2, sqrt3};
    const std::vector<expression::Expr> distinct =
        mathematics::combineStructurallyIdenticalTerms(distinctTerms, builtins);
    tests.expect(distinct.size() == 2,
        "ExactAlgebra: mathematically unrelated structures are not guessed equal");

    const std::array<numeric::Rational, 3> sqrt2Polynomial{
        numeric::Rational{numeric::BigInt{-2}},
        numeric::Rational{},
        numeric::Rational{numeric::BigInt{1}}};
    auto sqrt2Root = symbolic::AlgebraicNumber::create(
        sqrt2Polynomial, 2, symbolic::AlgebraicRootDomain::Real);
    tests.expect(sqrt2Root.has_value(),
        "NumberField: creates the selected exact generator embedding");
    if (sqrt2Root) {
        auto field = symbolic::NumberFieldContext::create(*sqrt2Root);
        tests.expect(field && field->degree() == 2,
            "NumberField: context retains the simple-extension degree");
        if (field) {
            std::vector<numeric::Rational> thetaCoordinates(2);
            thetaCoordinates[1] = numeric::Rational{numeric::BigInt{1}};
            auto theta = symbolic::AlgebraicElement::create(
                field, std::move(thetaCoordinates));
            tests.expect(theta.has_value(),
                "NumberField: generator has a power-basis coordinate vector");
            if (theta) {
                const auto square = theta->multiply(*theta);
                tests.expect(square.has_value()
                        && square->coefficients()[0] == numeric::Rational{numeric::BigInt{2}}
                        && square->coefficients()[1].isZero(),
                    "NumberField: multiplication reduces theta^2 modulo its minimal polynomial");

                const auto inverse = theta->divide(*theta);
                tests.expect(inverse.has_value()
                        && inverse->coefficients()[0] == numeric::Rational{numeric::BigInt{1}}
                        && inverse->coefficients()[1].isZero(),
                    "NumberField: exact inverse implements division inside Q(theta)");

                std::vector<numeric::Rational> shiftedCoordinates{
                    numeric::Rational{numeric::BigInt{1}},
                    numeric::Rational{numeric::BigInt{1}}};
                auto shifted = symbolic::AlgebraicElement::create(
                    field, std::move(shiftedCoordinates));
                const auto minimal = shifted ? shifted->minimalPolynomial() : std::nullopt;
                tests.expect(minimal && minimal->size() == 3
                        && (*minimal)[0] == numeric::Rational{numeric::BigInt{-1}}
                        && (*minimal)[1] == numeric::Rational{numeric::BigInt{-2}}
                        && (*minimal)[2] == numeric::Rational{numeric::BigInt{1}},
                    "NumberField: derives the minimal polynomial from persistent field coordinates");

                const auto shiftedInverse = shifted
                    ? field->reciprocal(shifted->coefficients()) : std::nullopt;
                const auto shiftedInverseAgain = shifted
                    ? field->reciprocal(shifted->coefficients()) : std::nullopt;
                const auto shiftedRoundTrip = shiftedInverse
                    ? field->reciprocal(*shiftedInverse) : std::nullopt;
                tests.expect(shifted && shiftedInverse && shiftedInverseAgain
                        && *shiftedInverse == *shiftedInverseAgain
                        && shiftedRoundTrip
                        && *shiftedRoundTrip == std::vector<numeric::Rational>(
                            shifted->coefficients().begin(), shifted->coefficients().end()),
                    "NumberField Stage5-2: reciprocal reuse preserves exact involution");

                tests.expect(theta->exactEquals(*theta) == std::optional<bool>{true},
                    "NumberField Stage3: same-field coordinate equality is exact");
                tests.expect(shifted
                        && theta->exactEquals(*shifted) == std::optional<bool>{false},
                    "NumberField Stage3: distinct same-field coordinates compare unequal exactly");
                tests.expect(theta->exactSign() == symbolic::AlgebraicSign::Positive,
                    "NumberField Stage3: generator sign is certified from the chosen real embedding");

                std::vector<numeric::Rational> negativeCoordinates{
                    numeric::Rational{numeric::BigInt{1}},
                    numeric::Rational{numeric::BigInt{-1}}};
                const auto negative = symbolic::AlgebraicElement::create(
                    field, std::move(negativeCoordinates));
                tests.expect(negative
                        && negative->exactSign() == symbolic::AlgebraicSign::Negative,
                    "NumberField Stage3: polynomial coordinates are signed by exact interval refinement");
            }
        }
    }

    {
        const std::array<numeric::Rational, 4> cubeRoot5Polynomial{
            numeric::Rational{numeric::BigInt{-5}},
            numeric::Rational{},
            numeric::Rational{},
            numeric::Rational{numeric::BigInt{1}}};
        const auto cubeRoot5 = symbolic::AlgebraicNumber::create(
            cubeRoot5Polynomial, 1, symbolic::AlgebraicRootDomain::Real);
        auto field = cubeRoot5 ? symbolic::NumberFieldContext::create(*cubeRoot5) : nullptr;
        std::vector<numeric::Rational> value{
            numeric::Rational{numeric::BigInt{1}},
            numeric::Rational{numeric::BigInt{1}},
            numeric::Rational{numeric::BigInt{1}}};
        std::atomic<bool> reciprocalRaceOk{field != nullptr};
        std::vector<std::thread> workers;
        if (field) {
            for (std::size_t worker = 0; worker < 8; ++worker) {
                workers.emplace_back([field, value, &reciprocalRaceOk] {
                    for (std::size_t iteration = 0; iteration < 32; ++iteration) {
                        const auto inverse = field->reciprocal(value);
                        if (!inverse) {
                            reciprocalRaceOk = false;
                            return;
                        }
                        const auto product = field->multiply(value, *inverse);
                        if (product.empty()
                            || product[0] != numeric::Rational{numeric::BigInt{1}}
                            || !std::all_of(product.begin() + 1, product.end(),
                                [](const numeric::Rational& coefficient) {
                                    return coefficient.isZero();
                                })) {
                            reciprocalRaceOk = false;
                            return;
                        }
                    }
                });
            }
        }
        for (std::thread& worker : workers)
            worker.join();
        tests.expect(reciprocalRaceOk.load(),
            "NumberField Stage5-2: reciprocal cache publish and lookup are thread-safe");

        std::atomic<bool> minimalPolynomialRaceOk{field != nullptr};
        workers.clear();
        if (field) {
            for (std::size_t worker = 0; worker < 8; ++worker) {
                workers.emplace_back([field, value, &minimalPolynomialRaceOk] {
                    const auto element = symbolic::AlgebraicElement::create(field, value);
                    const auto minimal = element ? element->minimalPolynomial() : std::nullopt;
                    const std::array<numeric::Rational, 4> expected{
                        numeric::Rational{numeric::BigInt{-16}},
                        numeric::Rational{numeric::BigInt{-12}},
                        numeric::Rational{numeric::BigInt{-3}},
                        numeric::Rational{numeric::BigInt{1}}};
                    if (!minimal || !std::equal(
                            minimal->begin(), minimal->end(), expected.begin(), expected.end()))
                        minimalPolynomialRaceOk = false;
                });
            }
        }
        for (std::thread& worker : workers)
            worker.join();
        tests.expect(minimalPolynomialRaceOk.load(),
            "NumberField Stage5-3: minimal-polynomial cache miss and publish are thread-safe");
    }

    {
        std::vector<numeric::Rational> degreeSixPolynomial(7);
        degreeSixPolynomial[0] = numeric::Rational{numeric::BigInt{-2}};
        degreeSixPolynomial[6] = numeric::Rational{numeric::BigInt{1}};
        const auto degreeSixRoot = symbolic::AlgebraicNumber::create(
            degreeSixPolynomial, 2, symbolic::AlgebraicRootDomain::Real);
        auto field = degreeSixRoot
            ? symbolic::NumberFieldContext::create(*degreeSixRoot) : nullptr;
        std::vector<numeric::Rational> thetaCubedCoordinates(6);
        thetaCubedCoordinates[3] = numeric::Rational{numeric::BigInt{1}};
        const auto thetaCubed = field
            ? symbolic::AlgebraicElement::create(field, std::move(thetaCubedCoordinates))
            : std::nullopt;
        const auto subfieldMinimal = thetaCubed
            ? thetaCubed->minimalPolynomial() : std::nullopt;
        tests.expect(subfieldMinimal && subfieldMinimal->size() == 3
                && (*subfieldMinimal)[0] == numeric::Rational{numeric::BigInt{-2}}
                && (*subfieldMinimal)[1].isZero()
                && (*subfieldMinimal)[2] == numeric::Rational{numeric::BigInt{1}},
            "NumberField Stage5-3: incremental Krylov elimination finds the first lower-degree relation");

        std::vector<numeric::Rational> rationalCoordinates(6);
        rationalCoordinates[0] = numeric::Rational{
            numeric::BigInt{3}, numeric::BigInt{2}};
        const auto rationalElement = field
            ? symbolic::AlgebraicElement::create(field, std::move(rationalCoordinates))
            : std::nullopt;
        const auto rationalMinimal = rationalElement
            ? rationalElement->minimalPolynomial() : std::nullopt;
        tests.expect(rationalMinimal && rationalMinimal->size() == 2
                && (*rationalMinimal)[0] == numeric::Rational{
                    numeric::BigInt{-3}, numeric::BigInt{2}}
                && (*rationalMinimal)[1] == numeric::Rational{numeric::BigInt{1}},
            "NumberField Stage5-3: incremental Krylov elimination handles degree-one elements");
    }

    if (sqrt2Root) {
        const symbolic::RationalRootInterval positiveInterval{
            numeric::Rational{numeric::BigInt{1}},
            numeric::Rational{numeric::BigInt{2}}};
        const auto directlyIsolated = symbolic::RealAlgebraicNumber::createFromMinimalPolynomialInterval(
            sqrt2Polynomial, positiveInterval);
        tests.expect(directlyIsolated && directlyIsolated->rootIndex() == 2,
            "NumberField Stage5-1: certified real interval identifies the canonical root directly");

        const symbolic::RationalRootInterval ambiguousInterval{
            numeric::Rational{numeric::BigInt{-2}},
            numeric::Rational{numeric::BigInt{2}}};
        tests.expect(!symbolic::RealAlgebraicNumber::createFromMinimalPolynomialInterval(
                sqrt2Polynomial, ambiguousInterval),
            "NumberField Stage5-1: direct real materialization rejects a non-isolating interval");

        const symbolic::AlgebraicNumber rooted = sqrt2Root->withGeneratorField();
        if (const auto* element = rooted.arithmeticElement()) {
            const auto interval = element->refinedRealInterval(96);
            tests.expect(interval
                    && numeric::Rational{numeric::BigInt{1}} < interval->lower
                    && interval->upper < numeric::Rational{numeric::BigInt{2}},
                "NumberField Stage5-1: field coordinates provide a certified real value interval");
            tests.expect(!rooted.exactRationalParts(),
                "NumberField Stage5-1: proven irrational real Roots skip rational probing");
        }
    }

    if (sqrt2Root) {
        const auto duplicateSqrt2 = symbolic::AlgebraicNumber::create(
            sqrt2Polynomial, 2, symbolic::AlgebraicRootDomain::Real);
        const auto firstField = symbolic::NumberFieldContext::create(*sqrt2Root);
        const auto secondField = duplicateSqrt2
            ? symbolic::NumberFieldContext::create(*duplicateSqrt2) : nullptr;
        tests.expect(firstField && secondField && firstField.get() == secondField.get(),
            "NumberField Stage4: independently constructed identical embedded fields are weak-interned");

        const auto negativeSqrt2ForField = symbolic::AlgebraicNumber::create(
            sqrt2Polynomial, 1, symbolic::AlgebraicRootDomain::Real);
        const auto distinctField = negativeSqrt2ForField
            ? symbolic::NumberFieldContext::create(*negativeSqrt2ForField) : nullptr;
        tests.expect(firstField && distinctField && firstField.get() != distinctField.get(),
            "NumberField Stage4: distinct generator embeddings never share a field context");
    }

    if (sqrt2Root) {
        const symbolic::AlgebraicNumber rooted = sqrt2Root->withGeneratorField();
        tests.expect(rooted.arithmeticElement() != nullptr,
            "NumberField Stage2: a proven Root receives a generator-field representation");
        if (rooted.arithmeticElement()) {
            const auto* fieldIdentity = rooted.arithmeticElement()->field().get();
            const auto square = symbolic::AlgebraicNumber::combine(
                rooted, *sqrt2Root, symbolic::AlgebraicBinaryOperation::Multiply);
            tests.expect(square && square->arithmeticElement()
                    && square->arithmeticElement()->field().get() == fieldIdentity,
                "NumberField Stage2: identical Root identity reuses an existing embedded field");

            std::vector<expression::Expr> rootArguments;
            rootArguments.push_back(expression::Expr::rationalArray(
                {sqrt2Polynomial.size()},
                std::vector<numeric::Rational>{sqrt2Polynomial.begin(), sqrt2Polynomial.end()}));
            rootArguments.push_back(integer(2));
            const expression::Expr cachedRoot = expression::Expr::call(
                builtins.symbol(evaluation::BuiltinId::Root),
                rootArguments,
                std::make_shared<const symbolic::AlgebraicNumber>(rooted));
            const expression::Expr rebuilt = expression::Expr::rebuildCall(
                cachedRoot.asCall(), rootArguments);
            tests.expect(rebuilt.asCall().algebraicValue
                    && rebuilt.asCall().algebraicValue->arithmeticElement()
                    && rebuilt.asCall().algebraicValue->arithmeticElement()->field().get()
                        == fieldIdentity,
                "NumberField Stage2: structurally unchanged Call rebuild preserves algebraic cache");

            rootArguments[1] = integer(1);
            const expression::Expr changed = expression::Expr::rebuildCall(
                cachedRoot.asCall(), std::move(rootArguments));
            tests.expect(!changed.asCall().algebraicValue,
                "NumberField Stage2: changed Root arguments invalidate algebraic cache");
        }
    }

    if (sqrt2Root) {
        const auto negativeSqrt2 = symbolic::AlgebraicNumber::create(
            sqrt2Polynomial, 1, symbolic::AlgebraicRootDomain::Real);
        tests.expect(negativeSqrt2
                && sqrt2Root->exactEquals(*negativeSqrt2) == std::optional<bool>{false},
            "Algebraic Stage3: distinct roots of one minimal polynomial are exactly unequal");
        tests.expect(negativeSqrt2
                && sqrt2Root->exactRealCompare(*negativeSqrt2)
                    == symbolic::AlgebraicOrder::Greater,
            "Algebraic Stage3: real root ordering uses certified isolating intervals");
    }

    const auto exactI = symbolic::AlgebraicNumber::fromComplexRational(
        numeric::Rational{}, numeric::Rational{numeric::BigInt{1}});
    const auto exactIParts = exactI ? exactI->exactRationalParts() : std::nullopt;
    tests.expect(exactIParts
            && exactIParts->first.isZero()
            && exactIParts->second == numeric::Rational{numeric::BigInt{1}},
        "NumberField Stage5-1: quadratic complex Rational detection remains exact");

    const std::array<numeric::Rational, 4> cubeRoot3Polynomial{
        numeric::Rational{numeric::BigInt{-3}},
        numeric::Rational{},
        numeric::Rational{},
        numeric::Rational{numeric::BigInt{1}}};
    auto cubeRoot3 = symbolic::AlgebraicNumber::create(
        cubeRoot3Polynomial, 1, symbolic::AlgebraicRootDomain::Real);
    if (sqrt2Root && cubeRoot3) {
        auto sum = symbolic::AlgebraicNumber::combine(
            *sqrt2Root, *cubeRoot3, symbolic::AlgebraicBinaryOperation::Add);
        tests.expect(sum && sum->arithmeticElement(),
            "NumberField: primitive-element arithmetic retains a persistent field element");
        const auto reversedSum = symbolic::AlgebraicNumber::combine(
            *cubeRoot3, *sqrt2Root, symbolic::AlgebraicBinaryOperation::Add);
        tests.expect(sum && reversedSum && sum->hasSameRootIdentity(*reversedSum),
            "NumberField Stage5-1: cached primitive embeddings are reusable with reversed operands");
        if (sum && sum->arithmeticElement()) {
            const auto* fieldIdentity = sum->arithmeticElement()->field().get();

            const auto duplicateSqrt2 = symbolic::AlgebraicNumber::create(
                sqrt2Polynomial, 2, symbolic::AlgebraicRootDomain::Real);
            const auto duplicateCubeRoot3 = symbolic::AlgebraicNumber::create(
                cubeRoot3Polynomial, 1, symbolic::AlgebraicRootDomain::Real);
            const auto duplicateSum = duplicateSqrt2 && duplicateCubeRoot3
                ? symbolic::AlgebraicNumber::combine(
                    *duplicateSqrt2, *duplicateCubeRoot3, symbolic::AlgebraicBinaryOperation::Add)
                : std::nullopt;
            tests.expect(duplicateSum && duplicateSum->arithmeticElement()
                    && duplicateSum->arithmeticElement()->field().get() == fieldIdentity,
                "NumberField Stage4: independent primitive-element construction reuses the interned context");

            auto square = symbolic::AlgebraicNumber::combine(
                *sum, *sum, symbolic::AlgebraicBinaryOperation::Multiply);
            tests.expect(square && square->arithmeticElement()
                    && square->arithmeticElement()->field().get() == fieldIdentity,
                "NumberField: chained same-field arithmetic reuses the existing context");
        }
    }
}

} // namespace mmcal::tests
