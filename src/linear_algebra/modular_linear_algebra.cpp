// 31-bit prime field + CRT + rational reconstructionによるexact線形代数backend。
#include "modular_linear_algebra.hpp"

#include "evaluation/evaluation_budget.hpp"
#include "expression/array_utils.hpp"
#include "numeric/integer_algorithms.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra {
namespace {

using numeric::BigInt;
using numeric::Rational;

constexpr std::uint32_t largestPrime31 = 2'147'483'647U;
constexpr std::size_t primePayloadBits = 30;
constexpr std::size_t maximumModularDeterminantPrimes = 4096;
constexpr std::size_t maximumModularSolvePrimes = 8192;

[[nodiscard]] std::size_t saturatedAdd(std::size_t lhs, std::size_t rhs) noexcept {
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs)
        return std::numeric_limits<std::size_t>::max();
    return lhs + rhs;
}

[[nodiscard]] std::size_t saturatedMultiply(std::size_t lhs, std::size_t rhs) noexcept {
    if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs)
        return std::numeric_limits<std::size_t>::max();
    return lhs * rhs;
}

[[nodiscard]] std::size_t ceilLog2(std::size_t value) noexcept {
    if (value <= 1)
        return 0;
    std::size_t result = 0;
    --value;
    while (value != 0) {
        value >>= 1;
        ++result;
    }
    return result;
}

[[nodiscard]] std::size_t maximumCoefficientBits(
    const IntegerMatrixBuffer& matrix,
    std::size_t columns) noexcept {
    columns = std::min(columns, matrix.columns());
    std::size_t result = 0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < columns; ++column)
            result = std::max(result, matrix(row, column).bitLength());
    return result;
}

[[nodiscard]] bool hasUsefulDensity(
    const IntegerMatrixBuffer& matrix,
    std::size_t columns) noexcept {
    columns = std::min(columns, matrix.columns());
    const std::size_t total = saturatedMultiply(matrix.rows(), columns);
    if (total == 0)
        return false;

    std::size_t nonzero = 0;
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < columns; ++column)
            nonzero += matrix(row, column).isZero() ? 0U : 1U;

    // Sparse/diagonal systems are a particularly strong Bareiss case.
    // At least 25% density is required before paying prime-image/CRT overhead.
    return nonzero >= (total + 3) / 4;
}

// |det A| < 2^B を保証する粗いHadamard bound。
// 各要素|a_ij| < 2^h，||row_i||_2 < sqrt(n) 2^h を使い，
// floating pointを介さずceil(log2 n)で上から押さえる。
[[nodiscard]] std::size_t hadamardExponentBound(
    std::size_t order,
    std::size_t coefficientBits) noexcept {
    if (order == 0 || coefficientBits == 0)
        return 0;
    const std::size_t coefficientPart = saturatedMultiply(order, coefficientBits);
    const std::size_t logPart = saturatedMultiply(order, ceilLog2(order));
    const std::size_t normPart = logPart / 2 + (logPart & 1U);
    return saturatedAdd(coefficientPart, normPart);
}

[[nodiscard]] std::size_t primesForBits(std::size_t bits) noexcept {
    if (bits == std::numeric_limits<std::size_t>::max())
        return bits;
    return bits / primePayloadBits + 1;
}

[[nodiscard]] std::size_t matrixElementCount(
    std::size_t rows,
    std::size_t columns) {
    const std::size_t shape[] = {rows, columns};
    return expression::arrayElementCount(shape);
}

[[nodiscard]] std::uint32_t subtractMod(
    std::uint32_t lhs,
    std::uint32_t rhs,
    std::uint32_t modulus) noexcept {
    return lhs >= rhs ? lhs - rhs : static_cast<std::uint32_t>(
        static_cast<std::uint64_t>(lhs) + modulus - rhs);
}

[[nodiscard]] std::uint32_t multiplyMod(
    std::uint32_t lhs,
    std::uint32_t rhs,
    std::uint32_t modulus) noexcept {
    return static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(lhs) * rhs) % modulus);
}

[[nodiscard]] std::uint32_t powerMod(
    std::uint32_t base,
    std::uint32_t exponent,
    std::uint32_t modulus) noexcept {
    std::uint32_t result = 1;
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = multiplyMod(result, base, modulus);
        exponent >>= 1;
        if (exponent != 0)
            base = multiplyMod(base, base, modulus);
    }
    return result;
}

[[nodiscard]] std::uint32_t inverseMod(
    std::uint32_t value,
    std::uint32_t prime) noexcept {
    return powerMod(value, prime - 2, prime);
}

class DescendingPrimeGenerator final {
public:
    [[nodiscard]] std::uint32_t next() {
        while (candidate_ >= 3) {
            candidate_ -= 2;
            if (numeric::isPrimeUint64(candidate_))
                return static_cast<std::uint32_t>(candidate_);
        }
        throw std::overflow_error("31-bit modular prime sequence exhausted");
    }

private:
    std::uint64_t candidate_ = static_cast<std::uint64_t>(largestPrime31) + 2;
};

class ModPrimeMatrix final {
public:
    ModPrimeMatrix(std::size_t rows, std::size_t columns)
        : rows_(rows), columns_(columns) {
        const std::size_t count = matrixElementCount(rows, columns);
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::TemporaryMatrixElement, count);
        values_.resize(count);
    }

    void load(const IntegerMatrixBuffer& source, std::uint32_t prime) {
        prime_ = prime;
        for (std::size_t row = 0; row < rows_; ++row) {
            evaluation::checkEvaluationCancellation();
            for (std::size_t column = 0; column < columns_; ++column)
                (*this)(row, column) = source(row, column).modulo(prime);
        }
    }

    void loadInverseAugmented(
        const IntegerMatrixBuffer& source,
        std::uint32_t prime) {
        if (rows_ != source.rows() || columns_ != source.columns() * 2)
            throw std::logic_error("Modular inverse workspace shape mismatch");
        prime_ = prime;
        const std::size_t n = source.rows();
        for (std::size_t row = 0; row < n; ++row) {
            evaluation::checkEvaluationCancellation();
            for (std::size_t column = 0; column < n; ++column)
                (*this)(row, column) = source(row, column).modulo(prime);
            for (std::size_t column = 0; column < n; ++column)
                (*this)(row, n + column) = row == column ? 1U : 0U;
        }
    }

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] std::uint32_t prime() const noexcept { return prime_; }

    [[nodiscard]] std::uint32_t& operator()(
        std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }

    [[nodiscard]] std::uint32_t operator()(
        std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }

    void swapRows(std::size_t lhs, std::size_t rhs) noexcept {
        if (lhs == rhs)
            return;
        for (std::size_t column = 0; column < columns_; ++column)
            std::swap((*this)(lhs, column), (*this)(rhs, column));
    }

private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::uint32_t prime_ = 0;
    std::vector<std::uint32_t> values_;
};

[[nodiscard]] std::uint32_t determinantModulo(ModPrimeMatrix& matrix) {
    const std::size_t n = matrix.rows();
    const std::uint32_t prime = matrix.prime();
    std::uint32_t determinant = 1;
    bool negative = false;

    for (std::size_t column = 0; column < n; ++column) {
        evaluation::checkEvaluationCancellation();
        std::size_t pivot = column;
        while (pivot < n && matrix(pivot, column) == 0)
            ++pivot;
        if (pivot == n)
            return 0;
        if (pivot != column) {
            matrix.swapRows(pivot, column);
            negative = !negative;
        }

        const std::uint32_t pivotValue = matrix(column, column);
        determinant = multiplyMod(determinant, pivotValue, prime);
        const std::uint32_t inversePivot = inverseMod(pivotValue, prime);
        for (std::size_t row = column + 1; row < n; ++row) {
            const std::uint32_t entry = matrix(row, column);
            if (entry == 0)
                continue;
            const std::uint32_t factor = multiplyMod(entry, inversePivot, prime);
            matrix(row, column) = 0;
            for (std::size_t c = column + 1; c < n; ++c)
                matrix(row, c) = subtractMod(
                    matrix(row, c), multiplyMod(factor, matrix(column, c), prime), prime);
        }
    }

    if (negative && determinant != 0)
        determinant = prime - determinant;
    return determinant;
}

struct ModularSolveImage final {
    bool unique = false;
    bool consistent = true;
    std::vector<std::uint32_t> solution;
};

[[nodiscard]] ModularSolveImage solveModulo(
    ModPrimeMatrix& matrix,
    std::size_t variables) {
    const std::size_t rhsColumns = matrix.columns() - variables;
    const std::uint32_t prime = matrix.prime();
    std::size_t pivotRow = 0;

    for (std::size_t column = 0;
         column < variables && pivotRow < matrix.rows();
         ++column) {
        evaluation::checkEvaluationCancellation();
        std::size_t selected = pivotRow;
        while (selected < matrix.rows() && matrix(selected, column) == 0)
            ++selected;
        if (selected == matrix.rows())
            continue;

        matrix.swapRows(selected, pivotRow);
        const std::uint32_t inversePivot = inverseMod(matrix(pivotRow, column), prime);
        for (std::size_t c = column; c < matrix.columns(); ++c)
            matrix(pivotRow, c) = multiplyMod(matrix(pivotRow, c), inversePivot, prime);

        for (std::size_t row = 0; row < matrix.rows(); ++row) {
            if (row == pivotRow)
                continue;
            const std::uint32_t factor = matrix(row, column);
            if (factor == 0)
                continue;
            matrix(row, column) = 0;
            for (std::size_t c = column + 1; c < matrix.columns(); ++c)
                matrix(row, c) = subtractMod(
                    matrix(row, c), multiplyMod(factor, matrix(pivotRow, c), prime), prime);
        }
        ++pivotRow;
    }

    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        bool coefficientNonZero = false;
        for (std::size_t column = 0; column < variables; ++column)
            if (matrix(row, column) != 0) {
                coefficientNonZero = true;
                break;
            }
        if (coefficientNonZero)
            continue;
        for (std::size_t rhs = 0; rhs < rhsColumns; ++rhs)
            if (matrix(row, variables + rhs) != 0)
                return {false, false, {}};
    }

    if (pivotRow != variables)
        return {false, true, {}};

    std::vector<std::uint32_t> solution(matrixElementCount(variables, rhsColumns));
    for (std::size_t variable = 0; variable < variables; ++variable)
        for (std::size_t rhs = 0; rhs < rhsColumns; ++rhs)
            solution[variable * rhsColumns + rhs] = matrix(variable, variables + rhs);
    return {true, true, std::move(solution)};
}

void appendCrt(
    BigInt& residue,
    BigInt& modulus,
    std::uint32_t image,
    std::uint32_t prime) {
    const std::uint32_t residueMod = residue.modulo(prime);
    const std::uint32_t modulusMod = modulus.modulo(prime);
    const std::uint32_t delta = subtractMod(image, residueMod, prime);
    const std::uint32_t correction = multiplyMod(
        delta, inverseMod(modulusMod, prime), prime);
    residue += modulus * BigInt::fromUnsigned(correction);
    modulus *= BigInt::fromUnsigned(prime);
}

void appendCrtVector(
    std::vector<BigInt>& residues,
    BigInt& modulus,
    const std::vector<std::uint32_t>& images,
    std::uint32_t prime) {
    if (residues.size() != images.size())
        throw std::logic_error("Modular CRT component count mismatch");

    const std::uint32_t modulusMod = modulus.modulo(prime);
    const std::uint32_t inverseModulus = inverseMod(modulusMod, prime);
    const BigInt previousModulus = modulus;
    for (std::size_t index = 0; index < residues.size(); ++index) {
        const std::uint32_t delta = subtractMod(
            images[index], residues[index].modulo(prime), prime);
        const std::uint32_t correction = multiplyMod(delta, inverseModulus, prime);
        residues[index] += previousModulus * BigInt::fromUnsigned(correction);
    }
    modulus = previousModulus * BigInt::fromUnsigned(prime);
}

[[nodiscard]] BigInt centeredRepresentative(BigInt residue, const BigInt& modulus) {
    if (residue + residue > modulus)
        residue -= modulus;
    return residue;
}

[[nodiscard]] std::optional<Rational> rationalReconstruct(
    const BigInt& residue,
    const BigInt& modulus) {
    if (residue.isZero())
        return Rational{};
    if (modulus <= BigInt{2})
        return std::nullopt;

    const BigInt half = modulus / BigInt{2};
    const BigInt bound = numeric::integerSqrt(half).root;
    BigInt r0 = modulus;
    BigInt r1 = residue;
    BigInt t0;
    BigInt t1{1};

    while (r1.abs() > bound) {
        if (r1.isZero())
            return std::nullopt;
        auto division = numeric::divmod(r0, r1);
        BigInt r2 = r0 - division.quotient * r1;
        BigInt t2 = t0 - division.quotient * t1;
        r0 = std::move(r1);
        r1 = std::move(r2);
        t0 = std::move(t1);
        t1 = std::move(t2);
    }

    if (t1.isZero() || t1.abs() > bound)
        return std::nullopt;

    Rational candidate{r1, t1};
    if (candidate.numerator().abs() > bound || candidate.denominator() > bound)
        return std::nullopt;

    const BigInt congruence = residue * candidate.denominator() - candidate.numerator();
    if (!(congruence % modulus).isZero())
        return std::nullopt;
    return candidate;
}

[[nodiscard]] bool verifyIntegerSolution(
    const IntegerMatrixBuffer& augmented,
    std::size_t variables,
    const std::vector<Rational>& solution) {
    const std::size_t rhsColumns = augmented.columns() - variables;
    if (solution.size() != variables * rhsColumns)
        return false;

    for (std::size_t row = 0; row < augmented.rows(); ++row) {
        evaluation::checkEvaluationCancellation();
        for (std::size_t rhs = 0; rhs < rhsColumns; ++rhs) {
            Rational sum;
            for (std::size_t variable = 0; variable < variables; ++variable)
                sum += Rational{augmented(row, variable)}
                    * solution[variable * rhsColumns + rhs];
            if (!(sum == Rational{augmented(row, variables + rhs)}))
                return false;
        }
    }
    return true;
}

[[nodiscard]] bool verifyAdjugate(
    const IntegerMatrixBuffer& matrix,
    const std::vector<BigInt>& adjugate,
    const BigInt& determinant) {
    const std::size_t n = matrix.rows();
    if (n != 0 && (adjugate.size() / n != n || adjugate.size() % n != 0))
        return false;

    for (std::size_t row = 0; row < n; ++row) {
        evaluation::checkEvaluationCancellation();
        for (std::size_t column = 0; column < n; ++column) {
            BigInt sum;
            for (std::size_t k = 0; k < n; ++k)
                sum += matrix(row, k) * adjugate[k * n + column];
            const BigInt expected = row == column ? determinant : BigInt{};
            if (sum != expected)
                return false;
        }
    }
    return true;
}

[[nodiscard]] std::optional<std::vector<Rational>> reconstructSolution(
    const std::vector<BigInt>& residues,
    const BigInt& modulus) {
    std::vector<Rational> solution;
    solution.reserve(residues.size());
    for (const BigInt& residue : residues) {
        const auto value = rationalReconstruct(residue, modulus);
        if (!value)
            return std::nullopt;
        solution.push_back(*value);
    }
    return solution;
}

} // namespace

std::size_t determinantReconstructionBits(const IntegerMatrixBuffer& matrix) noexcept {
    if (matrix.rows() != matrix.columns())
        return std::numeric_limits<std::size_t>::max();
    return hadamardExponentBound(
        matrix.rows(), maximumCoefficientBits(matrix, matrix.columns()));
}

std::size_t solutionReconstructionBits(
    const IntegerMatrixBuffer& augmented,
    std::size_t variables) noexcept {
    if (variables == 0)
        return 0;
    return hadamardExponentBound(
        variables, maximumCoefficientBits(augmented, augmented.columns()));
}

bool preferModularDeterminant(const IntegerMatrixBuffer& matrix) noexcept {
    if (matrix.rows() != matrix.columns())
        return false;
    const std::size_t order = matrix.rows();
    if (order < 4)
        return false;

    const std::size_t bits = determinantReconstructionBits(matrix);
    if (primesForBits(saturatedAdd(bits, 2)) > maximumModularDeterminantPrimes)
        return false;

    if (!hasUsefulDensity(matrix, matrix.columns()))
        return false;

    const std::size_t height = maximumCoefficientBits(matrix, matrix.columns());
    return order >= 48
        || (order >= 32 && height >= 64)
        || (order >= 24 && height >= 192);
}

bool preferModularSolve(
    const IntegerMatrixBuffer& augmented,
    std::size_t variables) noexcept {
    if (variables < 4 || augmented.columns() <= variables)
        return false;

    const std::size_t bits = solutionReconstructionBits(augmented, variables);
    const std::size_t target = saturatedAdd(saturatedMultiply(bits, 2), 3);
    if (primesForBits(target) > maximumModularSolvePrimes)
        return false;

    if (!hasUsefulDensity(augmented, variables))
        return false;

    const std::size_t height = maximumCoefficientBits(augmented, augmented.columns());
    return variables >= 24
        || (variables >= 12 && height >= 96)
        || (variables >= 8 && height >= 256)
        || (variables >= 6 && height >= 512);
}

BigInt modularDeterminant(
    const IntegerMatrixBuffer& matrix,
    ModularLinearAlgebraStats* stats) {
    if (matrix.rows() != matrix.columns())
        throw std::invalid_argument("Modular determinant requires a square matrix");
    if (matrix.rows() == 0)
        return BigInt{1};
    if (matrix.rows() == 1)
        return matrix(0, 0);

    const std::size_t determinantBits = determinantReconstructionBits(matrix);
    const std::size_t targetModulusBits = saturatedAdd(determinantBits, 2);
    if (primesForBits(targetModulusBits) > maximumModularDeterminantPrimes)
        throw std::length_error("Modular determinant reconstruction requires too many primes");

    DescendingPrimeGenerator primes;
    ModPrimeMatrix image{matrix.rows(), matrix.columns()};
    BigInt residue;
    BigInt modulus{1};

    while (modulus.bitLength() < targetModulusBits) {
        evaluation::checkEvaluationCancellation();
        const std::uint32_t prime = primes.next();
        evaluation::recordEvaluationModularPrime();
        if (stats)
            ++stats->primesTried;
        image.load(matrix, prime);
        const std::uint32_t determinant = determinantModulo(image);
        appendCrt(residue, modulus, determinant, prime);
        if (stats) {
            ++stats->primesAccepted;
            stats->reconstructedModulusBits = modulus.bitLength();
        }
    }

    return centeredRepresentative(std::move(residue), modulus);
}

std::optional<std::vector<Rational>> modularSolve(
    const IntegerMatrixBuffer& augmented,
    std::size_t variables,
    ModularLinearAlgebraStats* stats) {
    if (variables > augmented.columns())
        throw std::invalid_argument("Modular solve variable count exceeds matrix columns");
    const std::size_t rhsColumns = augmented.columns() - variables;
    if (rhsColumns == 0)
        throw std::invalid_argument("Modular solve requires at least one right-hand side");
    if (variables == 0)
        return std::vector<Rational>{};

    const std::size_t boundBits = solutionReconstructionBits(augmented, variables);
    const std::size_t targetModulusBits = saturatedAdd(saturatedMultiply(boundBits, 2), 3);
    const std::size_t requiredPrimes = primesForBits(targetModulusBits);
    if (requiredPrimes > maximumModularSolvePrimes)
        return std::nullopt;
    const std::size_t maximumAttempts = std::min(
        maximumModularSolvePrimes, saturatedAdd(requiredPrimes, 64));

    DescendingPrimeGenerator primes;
    ModPrimeMatrix image{augmented.rows(), augmented.columns()};
    const std::size_t componentCount = matrixElementCount(variables, rhsColumns);
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::TemporaryMatrixElement, componentCount);
    std::vector<BigInt> residues(componentCount);
    BigInt modulus{1};
    std::size_t accepted = 0;

    for (std::size_t attempt = 0; attempt < maximumAttempts; ++attempt) {
        evaluation::checkEvaluationCancellation();
        const std::uint32_t prime = primes.next();
        evaluation::recordEvaluationModularPrime();
        if (stats)
            ++stats->primesTried;
        image.load(augmented, prime);
        ModularSolveImage modular = solveModulo(image, variables);
        if (!modular.consistent || !modular.unique)
            continue;

        appendCrtVector(residues, modulus, modular.solution, prime);
        ++accepted;
        if (stats) {
            ++stats->primesAccepted;
            stats->reconstructedModulusBits = modulus.bitLength();
        }

        // 小さい解は1 primeでも復元できる一方，大きいinverse等で毎primeごとに
        // 全componentへextended Euclidを掛けるのは逆効果になる。1,2,4,8,... primeと
        // 最終保証域だけで候補を作り，必ず元の整数系でexact verificationする。
        const bool reconstructionCheckpoint =
            (accepted & (accepted - 1)) == 0
            || modulus.bitLength() >= targetModulusBits;
        if (reconstructionCheckpoint) {
            const auto candidate = reconstructSolution(residues, modulus);
            if (candidate && verifyIntegerSolution(augmented, variables, *candidate))
                return candidate;
        }

        // Cramer/Hadamard由来の一意再構成域まで到達した後もverificationできない場合は，
        // classificationを推測せずBareiss fallbackへ戻す。
        if (modulus.bitLength() >= targetModulusBits)
            break;
    }

    return std::nullopt;
}


std::optional<std::vector<Rational>> modularInverse(
    const IntegerMatrixBuffer& matrix,
    ModularLinearAlgebraStats* stats) {
    if (matrix.rows() != matrix.columns())
        throw std::invalid_argument("Modular inverse requires a square matrix");
    const std::size_t n = matrix.rows();
    if (n == 0)
        return std::vector<Rational>{};
    if (n == 1) {
        if (matrix(0, 0).isZero())
            throw std::domain_error("Matrix is singular");
        return std::vector<Rational>{Rational{BigInt{1}, matrix(0, 0)}};
    }

    const BigInt determinant = preferModularDeterminant(matrix)
        ? modularDeterminant(matrix, stats)
        : bareissDeterminant(matrix);
    if (determinant.isZero())
        throw std::domain_error("Matrix is singular");

    const std::size_t height = maximumCoefficientBits(matrix, matrix.columns());
    const std::size_t adjugateBits = hadamardExponentBound(n - 1, height);
    const std::size_t targetModulusBits = saturatedAdd(adjugateBits, 2);
    const std::size_t requiredPrimes = primesForBits(targetModulusBits);
    if (requiredPrimes > maximumModularSolvePrimes)
        return std::nullopt;

    if (n > std::numeric_limits<std::size_t>::max() / 2)
        throw std::length_error("Modular inverse augmented column count exceeds the size_t range");

    DescendingPrimeGenerator primes;
    ModPrimeMatrix image{n, n * 2};
    const std::size_t componentCount = matrixElementCount(n, n);
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::TemporaryMatrixElement, componentCount);
    std::vector<BigInt> residues(componentCount);
    BigInt modulus{1};
    std::size_t accepted = 0;
    const std::size_t maximumAttempts = std::min(
        maximumModularSolvePrimes, saturatedAdd(requiredPrimes, 64));

    for (std::size_t attempt = 0; attempt < maximumAttempts; ++attempt) {
        evaluation::checkEvaluationCancellation();
        const std::uint32_t prime = primes.next();
        evaluation::recordEvaluationModularPrime();
        if (stats)
            ++stats->primesTried;

        const std::uint32_t determinantImage = determinant.modulo(prime);
        if (determinantImage == 0)
            continue;

        image.loadInverseAugmented(matrix, prime);
        ModularSolveImage inverseImage = solveModulo(image, n);
        if (!inverseImage.unique || !inverseImage.consistent)
            continue;

        for (std::uint32_t& value : inverseImage.solution)
            value = multiplyMod(value, determinantImage, prime);
        appendCrtVector(residues, modulus, inverseImage.solution, prime);
        ++accepted;
        if (stats) {
            ++stats->primesAccepted;
            stats->reconstructedModulusBits = modulus.bitLength();
        }

        // Hadamard boundまで待てば一意復元できるが，小さいadjugateでは過剰になる。
        // 1,2,4,8,... primeごとにcentered candidateを作り，A*adj=det(A)Iを
        // exactに満たした時点で一意なadjugateとして確定できる。
        const bool verificationCheckpoint =
            (accepted & (accepted - 1)) == 0
            || modulus.bitLength() >= targetModulusBits;
        if (!verificationCheckpoint)
            continue;

        std::vector<BigInt> adjugate;
        adjugate.reserve(residues.size());
        for (const BigInt& residue : residues)
            adjugate.push_back(centeredRepresentative(residue, modulus));
        if (!verifyAdjugate(matrix, adjugate, determinant)) {
            if (modulus.bitLength() >= targetModulusBits)
                return std::nullopt;
            continue;
        }

        std::vector<Rational> inverse;
        inverse.reserve(adjugate.size());
        for (BigInt& value : adjugate)
            inverse.emplace_back(std::move(value), determinant);
        return inverse;
    }

    return std::nullopt;
}

} // namespace mmcal::linear_algebra
