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
#include "symbolic/polynomial.hpp"
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

    {
        // (x^4-2^40)(x^4-2^-40): 4根ずつが半径2^-10と2^10へ分かれる。
        // Newton-polygon初期値を使っても，最終結果は従来どおりRoucheで全根を分離する。
        const numeric::BigInt scale = numeric::BigInt{1} << 40;
        std::vector<numeric::Rational> polynomial(9);
        polynomial.front() = numeric::Rational{numeric::BigInt{1}};
        polynomial[4] = -(numeric::Rational{scale}
            + numeric::Rational{numeric::BigInt{1}, scale});
        polynomial.back() = numeric::Rational{numeric::BigInt{1}};
        const auto roots = symbolic::ComplexAlgebraicNumber::isolateAll(polynomial);
        tests.expect(roots && roots->size() == 8
                && std::all_of(roots->begin(), roots->end(), [](const auto& root) {
                    return numeric::Rational{} < root.isolatingDisk().radius;
                }),
            "Complex Root: multi-radius Newton-polygon candidates preserve certified all-root isolation");
    }
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

    {
        const auto q = [](std::int64_t value) {
            return numeric::Rational{numeric::BigInt{value}};
        };
        const auto multiplyPolynomial = [&](
            const symbolic::RationalPolynomial& lhs,
            const symbolic::RationalPolynomial& rhs) {
            std::vector<numeric::Rational> coefficients(
                lhs.degree() + rhs.degree() + 1, q(0));
            for (std::size_t i = 0; i <= lhs.degree(); ++i)
                for (std::size_t j = 0; j <= rhs.degree(); ++j)
                    coefficients[i + j] += lhs.coefficient(i) * rhs.coefficient(j);
            return symbolic::RationalPolynomial{std::move(coefficients)};
        };

        const symbolic::RationalPolynomial first{{q(1), q(0), q(1)}};
        const symbolic::RationalPolynomial second{{q(2), q(1), q(1)}};
        const symbolic::RationalPolynomial third{{q(3), q(-1), q(1)}};
        const symbolic::RationalPolynomial input = multiplyPolynomial(
            multiplyPolynomial(first, second), third);
        symbolic::RationalPolynomialFactorOptions options;
        options.maximumSparseTerms = 0;
        options.maximumBerlekampMatrixEntries = 0;
        options.maximumCombinations = 3;
        const auto partial = symbolic::factorRationalPolynomialOverQ(input, options);

        symbolic::RationalPolynomial reconstructed{{partial.scalar}};
        for (const auto& factor : partial.factors)
            reconstructed = multiplyPolynomial(reconstructed, factor);
        tests.expect(!partial.complete
                && partial.factors.size() == 2
                && partial.factors[0].degree() == 2
                && partial.factors[1].degree() == 4
                && reconstructed.coefficients() == input.coefficients(),
            "Polynomial factor: a recombination budget stop retains certified partial factors and an exact residual");

        options.maximumCombinations = 8;
        const auto complete = symbolic::factorRationalPolynomialOverQ(input, options);
        tests.expect(complete.complete && complete.factors.size() == 3,
            "Polynomial factor: Cantor-Zassenhaus and Hensel complete after sufficient recombination budget");

        options.maximumCombinations = 0;
        options.preferredRecombinationCandidates = 0;
        options.minimumCldModularFactors = 2;
        const auto cld = symbolic::factorRationalPolynomialOverQ(input, options);
        symbolic::RationalPolynomial cldProduct{{cld.scalar}};
        for (const auto& factor : cld.factors)
            cldProduct = multiplyPolynomial(cldProduct, factor);
        tests.expect(cld.complete
                && cld.factors.size() == 3
                && cldProduct.coefficients() == input.coefficients(),
            "Polynomial factor: exact CLD-LLL recombination succeeds with exhaustive subset search disabled");

        options.maximumLllRank = 2;
        const auto limitedCld = symbolic::factorRationalPolynomialOverQ(
            input, options);
        symbolic::RationalPolynomial limitedProduct{{limitedCld.scalar}};
        for (const auto& factor : limitedCld.factors)
            limitedProduct = multiplyPolynomial(limitedProduct, factor);
        tests.expect(!limitedCld.complete
                && limitedProduct.coefficients() == input.coefficients(),
            "Polynomial factor: a CLD lattice budget stop remains an exact partial factorization");

        // p=3は重根を持つbad reduction，p=5では各Q-irreducible quadraticが
        // 二本のlinear modular factorへ分かれる。CLDが4本を2 blockへ戻す必要がある。
        const symbolic::RationalPolynomial cldFirst{{q(1), q(0), q(1)}};
        const symbolic::RationalPolynomial cldSecond{{q(10), q(9), q(1)}};
        const symbolic::RationalPolynomial fourLocalFactors =
            multiplyPolynomial(cldFirst, cldSecond);
        symbolic::RationalPolynomialFactorOptions groupedOptions;
        groupedOptions.maximumSparseTerms = 0;
        groupedOptions.maximumBerlekampMatrixEntries = 0;
        groupedOptions.maximumPrimeTrials = 2;
        groupedOptions.maximumGoodPrimes = 1;
        groupedOptions.maximumCombinations = 4;
        groupedOptions.preferredRecombinationCandidates = 0;
        groupedOptions.minimumCldModularFactors = 2;
        const auto grouped = symbolic::factorRationalPolynomialOverQ(
            fourLocalFactors, groupedOptions);
        tests.expect(grouped.complete
                && grouped.factors.size() == 2
                && grouped.factors[0].degree() == 2
                && grouped.factors[1].degree() == 2,
            "Polynomial factor: CLD-LLL groups several local factors into exact rational blocks");

        groupedOptions.maximumCldLattices = 0;
        const auto withoutCld = symbolic::factorRationalPolynomialOverQ(
            fourLocalFactors, groupedOptions);
        tests.expect(!withoutCld.complete
                && withoutCld.factors.size() == 1
                && withoutCld.factors.front().coefficients()
                    == fourLocalFactors.coefficients(),
            "Polynomial factor: the same subset budget exposes the CLD accelerator boundary without losing exactness");

        std::vector<numeric::Rational> highDegreeCoefficients(68, q(0));
        highDegreeCoefficients[0] = q(-1);
        highDegreeCoefficients[1] = q(-5);
        highDegreeCoefficients[67] = q(1);
        const auto highDegree = symbolic::factorRationalPolynomialOverQ(
            symbolic::RationalPolynomial{std::move(highDegreeCoefficients)});
        tests.expect(highDegree.complete
                && highDegree.factors.size() == 1
                && highDegree.factors.front().degree() == 67,
            "Polynomial factor: multiple finite-field primes certify a high-degree irreducible leaf by degree intersection");
    }
}

} // namespace mmcal::tests
