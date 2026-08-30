// 特殊函数の保証付き評価
#include "certified_special_functions.hpp"

#include "certification_error.hpp"
#include "certified_atan.hpp"
#include "certified_constants.hpp"
#include "certified_complex_sqrt.hpp"
#include "certified_complex_transcendental.hpp"
#include "certified_exponential.hpp"
#include "certified_elementary_functions.hpp"
#include "certified_logarithm.hpp"
#include "certified_sqrt.hpp"
#include "certified_trigonometry.hpp"
#include "interval_math.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"
#include "numeric/rational_rounding.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <mutex>
#include <optional>
#include <string>
#include <stdexcept>
#include <vector>

namespace mmcal::approximation {
namespace {

using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

inline void consumeCertifiedWork(std::size_t amount = 1) {
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::CertifiedRefinement, amount);
}

[[nodiscard]] std::size_t checkedAdd(
    std::size_t lhs,
    std::size_t rhs,
    const char* message) {
    if (rhs > std::numeric_limits<std::size_t>::max() - lhs)
        throw std::overflow_error(message);
    return lhs + rhs;
}

[[nodiscard]] Rational binaryThreshold(std::size_t bits) {
    BigInt denominator{1};
    denominator <<= bits;
    return Rational{BigInt{1}, std::move(denominator)};
}

[[nodiscard]] Rational absRational(Rational value) {
    return value.numerator().isNegative() ? -value : value;
}


[[nodiscard]] RealInterval exactInterval(
    const Rational& value,
    std::size_t bits) {
    return RealInterval::fromRational(value, bits);
}

[[nodiscard]] RealInterval exactInterval(std::int64_t value, std::size_t bits) {
    return exactInterval(rational(value), bits);
}

[[nodiscard]] bool exactZeroPoint(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

[[nodiscard]] bool exactZeroPoint(const ComplexInterval& value) noexcept {
    return exactZeroPoint(value.real()) && exactZeroPoint(value.imaginary());
}

[[nodiscard]] bool exactNonPositiveIntegerPoint(const ComplexInterval& value) {
    if (!exactZeroPoint(value.imaginary()) || !value.real().isPoint())
        return false;
    const Rational real = value.real().lower().toRational();
    return real.isInteger() && real <= rational(0);
}

/*
旧実装

変更理由：
- 初回Gamma評価だけでB0...B128をすべて生成し，20桁級でも約68 msのcold-startを払っていた。
- 高精度側でStirling項数を増やしたい場合，単純にmaximumを256へ延ばすと生成だけで数百msへ増える。
- Gamma/Stirlingで実際に必要なのはB_2, B_4, ...だけなので，要求された偶数Bernoulli数までexactに遅延生成する方が適している。

元コード：

// Akiyama-Tanigawa法でBernoulli数をexact Rationalとして一度だけ生成する。GammaのStirling剰余評価では偶数添字だけを使う。
[[nodiscard]] const std::vector<Rational>& bernoulliNumbers() {
    static const std::vector<Rational> values = [] {
        constexpr std::size_t maximum = 128;
        std::vector<Rational> a(maximum + 1);
        std::vector<Rational> b(maximum + 1);
        for (std::size_t m = 0; m <= maximum; ++m) {
            a[m] = Rational{BigInt{1}, BigInt::parse(std::to_string(m + 1))};
            for (std::size_t j = m; j >= 1; --j) {
                a[j - 1] = Rational{BigInt::parse(std::to_string(j))} * (a[j - 1] - a[j]);
                if (j == 1)
                    break;
            }
            b[m] = a[0];
        }
        return b;
    }();
    return values;
}

[[nodiscard]] Rational stirlingCoefficient(std::size_t k) {
    const auto& b = bernoulliNumbers();
    const std::size_t n = 2 * k;
    const BigInt denominator = BigInt::parse(std::to_string(n))
        * BigInt::parse(std::to_string(n - 1));
    return b[n] / Rational{denominator};
}
*/

/*
旧実装

変更理由：
- stateful lazy Akiyama-Tanigawaは低～中精度ではeager生成を避けられるが，
  1000 bit級でB_128近傍へ初めて到達した際のexact Rational更新がcold-startを支配した。
- Stirling backendが必要とするB_2...B_128は固定された厳密有理定数であり，
  実行時に再導出する数学的必要はない。
- したがって定数表はexact decimal numerator/denominatorとして保持し，
  実際に参照した値だけBigInt/Rationalへlazy parseする。

class BernoulliEvenCache final {
public:
    [[nodiscard]] Rational get(std::size_t k) {
        std::lock_guard lock{mutex_};
        const std::size_t targetOrder = 2 * k;
        while (a_.size() <= targetOrder)
            appendNext();
        return evenValues_[k];
    }

private:
    void appendNext() {
        const std::size_t m = a_.size();
        a_.push_back(Rational{
            BigInt{1}, BigInt::fromUnsigned(static_cast<std::uint64_t>(m + 1))});
        for (std::size_t j = m; j >= 1; --j) {
            a_[j - 1] = Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(j))}
                * (a_[j - 1] - a_[j]);
            if (j == 1)
                break;
        }
        if ((m & 1U) == 0)
            evenValues_.push_back(a_[0]);
    }

    std::mutex mutex_;
    std::vector<Rational> a_;
    std::vector<Rational> evenValues_;
};

[[nodiscard]] Rational bernoulliEven(std::size_t k) {
    static BernoulliEvenCache cache;
    return cache.get(k);
}
*/

struct BernoulliLiteral final {
    const char* numerator;
    const char* denominator;
};

constexpr std::array<BernoulliLiteral, 64> bernoulliEvenLiterals = {{
    BernoulliLiteral{"1", "6"},
    BernoulliLiteral{"-1", "30"},
    BernoulliLiteral{"1", "42"},
    BernoulliLiteral{"-1", "30"},
    BernoulliLiteral{"5", "66"},
    BernoulliLiteral{"-691", "2730"},
    BernoulliLiteral{"7", "6"},
    BernoulliLiteral{"-3617", "510"},
    BernoulliLiteral{"43867", "798"},
    BernoulliLiteral{"-174611", "330"},
    BernoulliLiteral{"854513", "138"},
    BernoulliLiteral{"-236364091", "2730"},
    BernoulliLiteral{"8553103", "6"},
    BernoulliLiteral{"-23749461029", "870"},
    BernoulliLiteral{"8615841276005", "14322"},
    BernoulliLiteral{"-7709321041217", "510"},
    BernoulliLiteral{"2577687858367", "6"},
    BernoulliLiteral{"-26315271553053477373", "1919190"},
    BernoulliLiteral{"2929993913841559", "6"},
    BernoulliLiteral{"-261082718496449122051", "13530"},
    BernoulliLiteral{"1520097643918070802691", "1806"},
    BernoulliLiteral{"-27833269579301024235023", "690"},
    BernoulliLiteral{"596451111593912163277961", "282"},
    BernoulliLiteral{"-5609403368997817686249127547", "46410"},
    BernoulliLiteral{"495057205241079648212477525", "66"},
    BernoulliLiteral{"-801165718135489957347924991853", "1590"},
    BernoulliLiteral{"29149963634884862421418123812691", "798"},
    BernoulliLiteral{"-2479392929313226753685415739663229", "870"},
    BernoulliLiteral{"84483613348880041862046775994036021", "354"},
    BernoulliLiteral{"-1215233140483755572040304994079820246041491", "56786730"},
    BernoulliLiteral{"12300585434086858541953039857403386151", "6"},
    BernoulliLiteral{"-106783830147866529886385444979142647942017", "510"},
    BernoulliLiteral{"1472600022126335654051619428551932342241899101", "64722"},
    BernoulliLiteral{"-78773130858718728141909149208474606244347001", "30"},
    BernoulliLiteral{"1505381347333367003803076567377857208511438160235", "4686"},
    BernoulliLiteral{"-5827954961669944110438277244641067365282488301844260429", "140100870"},
    BernoulliLiteral{"34152417289221168014330073731472635186688307783087", "6"},
    BernoulliLiteral{"-24655088825935372707687196040585199904365267828865801", "30"},
    BernoulliLiteral{"414846365575400828295179035549542073492199375372400483487", "3318"},
    BernoulliLiteral{"-4603784299479457646935574969019046849794257872751288919656867", "230010"},
    BernoulliLiteral{"1677014149185145836823154509786269900207736027570253414881613", "498"},
    BernoulliLiteral{"-2024576195935290360231131160111731009989917391198090877281083932477", "3404310"},
    BernoulliLiteral{"660714619417678653573847847426261496277830686653388931761996983", "6"},
    BernoulliLiteral{"-1311426488674017507995511424019311843345750275572028644296919890574047", "61410"},
    BernoulliLiteral{"1179057279021082799884123351249215083775254949669647116231545215727922535", "272118"},
    BernoulliLiteral{"-1295585948207537527989427828538576749659341483719435143023316326829946247", "1410"},
    BernoulliLiteral{"1220813806579744469607301679413201203958508415202696621436215105284649447", "6"},
    BernoulliLiteral{"-211600449597266513097597728109824233673043954389060234150638733420050668349987259", "4501770"},
    BernoulliLiteral{"67908260672905495624051117546403605607342195728504487509073961249992947058239", "6"},
    BernoulliLiteral{"-94598037819122125295227433069493721872702841533066936133385696204311395415197247711", "33330"},
    BernoulliLiteral{"3204019410860907078243020782116241775491817197152717450679002501086861530836678158791", "4326"},
    BernoulliLiteral{"-319533631363830011287103352796174274671189606078272738327103470162849568365549721224053", "1590"},
    BernoulliLiteral{"36373903172617414408151820151593427169231298640581690038930816378281879873386202346572901", "642"},
    BernoulliLiteral{"-3469342247847828789552088659323852541399766785760491146870005891371501266319724897592306597338057", "209191710"},
    BernoulliLiteral{"7645992940484742892248134246724347500528752413412307906683593870759797606269585779977930217515", "1518"},
    BernoulliLiteral{"-2650879602155099713352597214685162014443151499192509896451788427680966756514875515366781203552600109", "1671270"},
    BernoulliLiteral{"21737832319369163333310761086652991475721156679090831360806110114933605484234593650904188618562649", "42"},
    BernoulliLiteral{"-309553916571842976912513458033841416869004128064329844245504045721008957524571968271388199595754752259", "1770"},
    BernoulliLiteral{"366963119969713111534947151585585006684606361080699204301059440676414485045806461889371776354517095799", "6"},
    BernoulliLiteral{"-51507486535079109061843996857849983274095170353262675213092869167199297474922985358811329367077682677803282070131", "2328255930"},
    BernoulliLiteral{"49633666079262581912532637475990757438722790311060139770309311793150683214100431329033113678098037968564431", "6"},
    BernoulliLiteral{"-95876775334247128750774903107542444620578830013297336819553512729358593354435944413631943610268472689094609001", "30"},
    BernoulliLiteral{"5556330281949274850616324408918951380525567307126747246796782304333594286400508981287241419934529638692081513802696639", "4357878"},
    BernoulliLiteral{"-267754707742548082886954405585282394779291459592551740629978686063357792734863530145362663093519862048495908453718017", "510"},
}};

class BernoulliEvenCache final {
public:
    [[nodiscard]] Rational get(std::size_t k) {
        if (k == 0 || k > bernoulliEvenLiterals.size())
            throw std::out_of_range("Bernoulli index exceeds the certified Stirling table");

        std::lock_guard lock{mutex_};
        std::optional<Rational>& value = values_[k - 1];
        if (!value) {
            const BernoulliLiteral& literal = bernoulliEvenLiterals[k - 1];
            value.emplace(
                BigInt::parse(literal.numerator),
                BigInt::parse(literal.denominator));
        }
        return *value;
    }

private:
    std::mutex mutex_;
    std::array<std::optional<Rational>, 64> values_;
};

[[nodiscard]] Rational bernoulliEven(std::size_t k) {
    static BernoulliEvenCache cache;
    return cache.get(k);
}

[[nodiscard]] Rational stirlingCoefficient(std::size_t k) {
    const std::size_t n = 2 * k;
    const BigInt denominator = BigInt::fromUnsigned(static_cast<std::uint64_t>(n))
        * BigInt::fromUnsigned(static_cast<std::uint64_t>(n - 1));
    return bernoulliEven(k) / Rational{denominator};
}

struct StirlingPlan final {
    std::size_t shift = 0;
    std::size_t omittedK = 0;
    Rational remainderBound;
};

/*
旧実装

変更理由：
- maximumK=64固定では640 bit級でx=1/3を272段も右へ送る必要があり，recurrence productが支配的になっていた。
- Stirling和は各項でinverseOddを更新する逐次形だったため，Horner形よりinterval multiplicationが多かった。
- exact pointでもrecurrence productをintervalで1因子ずつ掛けており，巨大なshiftで不要なroundingと中間計算が増えていた。

元コード：

[[nodiscard]] StirlingPlan chooseStirlingPlan(
    const Rational& inputLower,
    std::size_t precisionBits) {
    if (inputLower <= rational(0))
        throw std::domain_error("LogGamma requires a positive interval in the Stirling backend");

    const Rational threshold = binaryThreshold(checkedAdd(
        precisionBits, 20, "Gamma precision is too large"));
    constexpr std::size_t maximumK = 64;

    for (std::size_t shift = 0; shift <= 1'000'000; shift += 8) {
        const Rational x = inputLower + Rational{BigInt::parse(std::to_string(shift))};
        if (x < rational(4))
            continue;

        const Rational inverseSquare = rational(1) / (x * x);
        Rational inverseOdd = rational(1) / x;
        for (std::size_t k = 1; k <= maximumK; ++k) {
            const Rational bound = absRational(stirlingCoefficient(k)) * inverseOdd;
            if (bound <= threshold)
                return StirlingPlan{shift, k, bound};
            inverseOdd *= inverseSquare;
        }
    }
    throw std::overflow_error("Gamma precision requires an excessive recurrence shift");
}

[[nodiscard]] RealInterval encloseLogGammaPositive(
    const RealInterval& input,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (input.lower() <= zero)
        throw std::domain_error("LogGamma positive backend requires x > 0");

    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Gamma working precision is too large");
    const StirlingPlan plan = chooseStirlingPlan(
        input.lower().toRational(), workBits);

    const RealInterval shift = exactInterval(
        Rational{BigInt::parse(std::to_string(plan.shift))}, workBits);
    const RealInterval x = add(input.roundedOutward(workBits), shift, workBits);
    const RealInterval logX = encloseLogPositive(x, workBits).interval;
    const RealInterval xMinusHalf = subtract(x, exactInterval(rational(1, 2), workBits), workBits);

    RealInterval result = subtract(
        multiply(xMinusHalf, logX, workBits), x, workBits);

    // 1/2 log(2 Pi)
    const RealInterval twoPi = multiply(
        enclosePi(workBits).interval, exactInterval(2, workBits), workBits);
    const RealInterval halfLogTwoPi = multiply(
        encloseLogPositive(twoPi, workBits).interval,
        exactInterval(rational(1, 2), workBits), workBits);
    result = add(result, halfLogTwoPi, workBits);

    const RealInterval inverseX = divide(exactInterval(1, workBits), x, workBits);
    const RealInterval inverseSquare = multiply(inverseX, inverseX, workBits);
    RealInterval inverseOdd = inverseX;
    for (std::size_t k = 1; k < plan.omittedK; ++k) {
        const RealInterval coefficient = exactInterval(stirlingCoefficient(k), workBits);
        result = add(result, multiply(coefficient, inverseOdd, workBits), workBits);
        inverseOdd = multiply(inverseOdd, inverseSquare, workBits);
    }

    const Rational omittedCoefficient = stirlingCoefficient(plan.omittedK);
    const Rational bound = plan.remainderBound;
    const RealInterval remainder = omittedCoefficient.numerator().isNegative()
        ? RealInterval::fromRationalBounds(-bound, rational(0), workBits)
        : RealInterval::fromRationalBounds(rational(0), bound, workBits);
    result = add(result, remainder, workBits);

    if (plan.shift != 0) {
        RealInterval product = exactInterval(1, workBits);
        for (std::size_t j = 0; j < plan.shift; ++j) {
            const RealInterval factor = add(
                input.roundedOutward(workBits),
                exactInterval(Rational{BigInt::parse(std::to_string(j))}, workBits),
                workBits);
            product = multiply(product, factor, workBits);
        }
        result = subtract(result, encloseLogPositive(product, workBits).interval, workBits);
    }

    return result.roundedOutward(precisionBits);
}
*/


[[nodiscard]] bool stirlingBoundMeetsBinaryThreshold(
    const Rational& inputLower,
    std::size_t k,
    std::size_t shift,
    std::size_t thresholdBits) {
    const Rational x = inputLower
        + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(shift))};
    if (x < rational(4))
        return false;

    const Rational coefficient = absRational(stirlingCoefficient(k));
    const std::uint64_t exponent = static_cast<std::uint64_t>(2 * k - 1);
    const BigInt denominatorPower = numeric::pow(inputLower.denominator(), exponent);
    const BigInt shiftedNumerator = inputLower.numerator()
        + inputLower.denominator()
            * BigInt::fromUnsigned(static_cast<std::uint64_t>(shift));
    const BigInt shiftedPower = numeric::pow(shiftedNumerator, exponent);

    BigInt lhs = coefficient.numerator().abs() * denominatorPower;
    lhs <<= thresholdBits;
    const BigInt rhs = coefficient.denominator() * shiftedPower;
    return lhs <= rhs;
}

[[nodiscard]] StirlingPlan chooseHighPrecisionStirlingPlan(
    const Rational& inputLower,
    std::size_t precisionBits,
    std::size_t maximumK) {
    const std::size_t thresholdBits = checkedAdd(
        precisionBits, 20, "Gamma precision is too large");

    // 高精度ではshiftを8刻みで走査し，各候補でk=1..maximumKをexact Rational評価すると
    // planner自身がStirling本体より高価になり得る。最大kの厳密剰余不等式をBigIntだけで判定し，
    // doubling + binary searchで十分なshiftを直接求める。最終boundもexact Rationalで再構築する。
    std::size_t low = 0;
    std::size_t high = 8;
    while (!stirlingBoundMeetsBinaryThreshold(
        inputLower, maximumK, high, thresholdBits)) {
        consumeCertifiedWork();
        if (high >= 1'000'000 / 2)
            throw std::overflow_error("Gamma precision requires an excessive recurrence shift");
        high *= 2;
    }

    while (low + 1 < high) {
        consumeCertifiedWork();
        const std::size_t middle = low + (high - low) / 2;
        if (stirlingBoundMeetsBinaryThreshold(
            inputLower, maximumK, middle, thresholdBits))
            high = middle;
        else
            low = middle;
    }

    const Rational x = inputLower
        + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(high))};
    const Rational coefficient = absRational(stirlingCoefficient(maximumK));
    const std::uint64_t exponent = static_cast<std::uint64_t>(2 * maximumK - 1);
    const BigInt numerator = coefficient.numerator().abs()
        * numeric::pow(x.denominator(), exponent);
    const BigInt denominator = coefficient.denominator()
        * numeric::pow(x.numerator(), exponent);
    return StirlingPlan{high, maximumK, Rational{numerator, denominator}};
}


[[nodiscard]] std::optional<StirlingPlan> findStirlingPlan(
    const Rational& inputLower,
    const Rational& threshold,
    std::size_t maximumK) {
    for (std::size_t shift = 0; shift <= 1'000'000; shift += 8) {
        consumeCertifiedWork();
        const Rational x = inputLower + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(shift))};
        if (x < rational(4))
            continue;

        const Rational inverseSquare = rational(1) / (x * x);
        Rational inverseOdd = rational(1) / x;
        for (std::size_t k = 1; k <= maximumK; ++k) {
            consumeCertifiedWork();
            const Rational bound = absRational(stirlingCoefficient(k)) * inverseOdd;
            if (bound <= threshold)
                return StirlingPlan{shift, k, bound};
            inverseOdd *= inverseSquare;
        }
    }
    return std::nullopt;
}

[[nodiscard]] StirlingPlan chooseStirlingPlan(
    const Rational& inputLower,
    std::size_t precisionBits) {
    if (inputLower <= rational(0))
        throw std::domain_error("LogGamma requires a positive interval in the Stirling backend");

    struct CacheEntry final {
        Rational inputLower;
        std::size_t precisionBits = 0;
        StirlingPlan plan;
    };
    static thread_local std::vector<CacheEntry> cache;
    for (const CacheEntry& entry : cache) {
        if (entry.precisionBits == precisionBits && entry.inputLower == inputLower)
            return entry.plan;
    }

    const Rational threshold = binaryThreshold(checkedAdd(
        precisionBits, 20, "Gamma precision is too large"));

    constexpr std::size_t highPrecisionPlannerThreshold = 768;
    if (precisionBits > highPrecisionPlannerThreshold) {
        const StirlingPlan plan = chooseHighPrecisionStirlingPlan(
            inputLower, precisionBits, 64);
        constexpr std::size_t maximumCacheEntries = 16;
        if (cache.size() == maximumCacheEntries)
            cache.erase(cache.begin());
        cache.push_back(CacheEntry{inputLower, precisionBits, plan});
        return plan;
    }

    // 低～中精度で最初の小さいx候補からk=64まで総当たりすると，
    // 最終planがk=16程度でも高次係数まで参照・interval化しやすい。precisionに応じて探索上限を絞る。
    const std::size_t maximumK = std::min<std::size_t>(
        64, std::max<std::size_t>(16, (precisionBits + 4) / 5));
    const auto plan = findStirlingPlan(inputLower, threshold, maximumK);
    if (!plan)
        throw std::overflow_error("Gamma precision requires an excessive recurrence shift");

    // planはexact Rational boundだけから決まり再利用しても意味論が変わらない。thread-local bounded cacheでmutexも共有状態も増やさない。
    constexpr std::size_t maximumCacheEntries = 16;
    if (cache.size() == maximumCacheEntries)
        cache.erase(cache.begin());
    cache.push_back(CacheEntry{inputLower, precisionBits, *plan});
    return *plan;
}

/*
旧実装

変更理由：
- exact Rational z=p/q の (z)_n で Rational を各nodeごとに作ると，積のたびにnormalize/GCDが走る。
- gcd(p,q)=1なら各因子 p+qk もqと互いに素なので，分子積とq^nを別々に構築して最後に一度だけRational化できる。
- JohanssonのGamma実装指針でも，rational rising factorialは分子・分母を未約分のままbinary splittingし，最後だけcanonicalizeするのが推奨されている。

[[nodiscard]] Rational balancedRisingProduct(
    const Rational& input,
    std::size_t first,
    std::size_t count) {
    if (count == 0)
        return rational(1);
    if (count == 1)
        return input + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(first))};
    const std::size_t leftCount = count / 2;
    return balancedRisingProduct(input, first, leftCount)
        * balancedRisingProduct(input, first + leftCount, count - leftCount);
}
*/

[[nodiscard]] BigInt balancedArithmeticProgressionProduct(
    const BigInt& numerator,
    const BigInt& denominator,
    std::size_t first,
    std::size_t count) {
    if (count == 0)
        return BigInt{1};
    if (count == 1)
        return numerator + denominator
            * BigInt::fromUnsigned(static_cast<std::uint64_t>(first));

    const std::size_t leftCount = count / 2;
    return balancedArithmeticProgressionProduct(
               numerator, denominator, first, leftCount)
        * balancedArithmeticProgressionProduct(
               numerator, denominator, first + leftCount, count - leftCount);
}

[[nodiscard]] Rational balancedRisingProduct(
    const Rational& input,
    std::size_t first,
    std::size_t count) {
    if (count == 0)
        return rational(1);
    consumeCertifiedWork(count);

    const BigInt numerator = balancedArithmeticProgressionProduct(
        input.numerator(), input.denominator(), first, count);
    const BigInt denominator = numeric::pow(
        input.denominator(), static_cast<std::uint64_t>(count));
    return Rational{numerator, denominator};
}

[[nodiscard]] RealInterval encloseLogGammaPositive(
    const RealInterval& input,
    std::size_t precisionBits,
    const Rational* exactInput = nullptr) {
    const BigFloat zero;
    if (input.lower() <= zero)
        throw std::domain_error("LogGamma positive backend requires x > 0");

    // Gamma(1)=Gamma(2)=1なのでlogGammaはexactに0。Beta(a,b)でa+b=1となる場合にも効く。
    if (exactInput) {
        if (*exactInput == rational(1) || *exactInput == rational(2))
            return exactInterval(0, precisionBits);
    }
    else if (input.isPoint()) {
        const Rational point = input.lower().toRational();
        if (point == rational(1) || point == rational(2))
            return exactInterval(0, precisionBits);
    }

    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Gamma working precision is too large");
    const Rational planInput = exactInput ? *exactInput : input.lower().toRational();
    const StirlingPlan plan = chooseStirlingPlan(planInput, workBits);

    const Rational shiftValue{BigInt::fromUnsigned(static_cast<std::uint64_t>(plan.shift))};
    const RealInterval shift = exactInterval(shiftValue, workBits);
    const RealInterval x = exactInput
        ? exactInterval(*exactInput + shiftValue, workBits)
        : add(input.roundedOutward(workBits), shift, workBits);
    const RealInterval logX = encloseLogPositive(x, workBits).interval;
    const RealInterval xMinusHalf = subtract(x, exactInterval(rational(1, 2), workBits), workBits);

    RealInterval result = subtract(
        multiply(xMinusHalf, logX, workBits), x, workBits);

    // 1/2 log(2 Pi)
    const RealInterval twoPi = multiply(
        enclosePi(workBits).interval, exactInterval(2, workBits), workBits);
    const RealInterval halfLogTwoPi = multiply(
        encloseLogPositive(twoPi, workBits).interval,
        exactInterval(rational(1, 2), workBits), workBits);
    result = add(result, halfLogTwoPi, workBits);

    const RealInterval inverseX = divide(exactInterval(1, workBits), x, workBits);
    if (plan.omittedK > 1) {
        const RealInterval inverseSquare = multiply(inverseX, inverseX, workBits);
        RealInterval series = exactInterval(stirlingCoefficient(plan.omittedK - 1), workBits);
        for (std::size_t k = plan.omittedK - 1; k > 1; --k) {
            consumeCertifiedWork();
            series = add(
                exactInterval(stirlingCoefficient(k - 1), workBits),
                multiply(inverseSquare, series, workBits),
                workBits);
        }
        result = add(result, multiply(inverseX, series, workBits), workBits);
    }

    // 正実軸上のStirling級数の剰余は最初の省略項と同符号で、絶対値はその項を超えない。入力区間ではlower endpointが最大絶対値を与える。
    const Rational omittedCoefficient = stirlingCoefficient(plan.omittedK);
    const Rational bound = plan.remainderBound;
    const RealInterval remainder = omittedCoefficient.numerator().isNegative()
        ? RealInterval::fromRationalBounds(-bound, rational(0), workBits)
        : RealInterval::fromRationalBounds(rational(0), bound, workBits);
    result = add(result, remainder, workBits);

    if (plan.shift != 0) {
        if (exactInput || input.isPoint()) {
            // exact Rational identityが呼出元に残っている場合は，dyadic intervalへ落とした後の
            // lower/upper endpointではなく元の値からrising factorialを構成する。1/3等では
            // RealInterval::fromRationalがpointにならないため，isPoint()だけではこの経路へ入れなかった。
            const Rational product = balancedRisingProduct(
                exactInput ? *exactInput : input.lower().toRational(), 0, plan.shift);
            result = subtract(
                result,
                encloseLogPositive(exactInterval(product, workBits), workBits).interval,
                workBits);
        }
        else {
            // 非point intervalは旧来の包含安全な逐次interval productを維持する。
            RealInterval product = exactInterval(1, workBits);
            for (std::size_t j = 0; j < plan.shift; ++j) {
                consumeCertifiedWork();
                const RealInterval factor = add(
                    input.roundedOutward(workBits),
                    exactInterval(Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(j))}, workBits),
                    workBits);
                product = multiply(product, factor, workBits);
            }
            result = subtract(result, encloseLogPositive(product, workBits).interval, workBits);
        }
    }

    return result.roundedOutward(precisionBits);
}

[[nodiscard]] bool exactNonPositiveIntegerPoint(const RealInterval& input) {
    if (!input.isPoint())
        return false;
    const Rational value = input.lower().toRational();
    return value.isInteger() && value.numerator() <= BigInt{0};
}

[[nodiscard]] RealInterval gammaNegativeByReflection(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (exactNonPositiveIntegerPoint(input))
        throw std::domain_error("gamma is undefined at non-positive integers");

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "Gamma reflection precision is too large");
    const RealInterval oneMinusX = subtract(
        exactInterval(1, workBits), input.roundedOutward(workBits), workBits);
    const RealInterval gammaComplement = encloseExp(
        encloseLogGammaPositive(oneMinusX, workBits), workBits).interval;
    const RealInterval pi = enclosePi(workBits).interval;
    const RealInterval piX = multiply(pi, input.roundedOutward(workBits), workBits);
    const RealInterval sine = encloseSinRadianInterval(piX, workBits).interval;
    if (sine.containsZero())
        throw PrecisionInsufficient{"Gamma reflection cannot yet exclude a non-positive-integer pole"};
    const RealInterval denominator = multiply(sine, gammaComplement, workBits);
    if (denominator.containsZero())
        throw PrecisionInsufficient{"Gamma reflection denominator cannot yet be proven nonzero"};
    return divide(pi, denominator, workBits).roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval gammaNegativeRationalByReflection(
    const Rational& input,
    std::size_t precisionBits) {
    if (input.isInteger() && input.numerator() <= BigInt{0})
        throw std::domain_error("gamma is undefined at non-positive integers");

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "Gamma rational reflection precision is too large");
    const Rational complement = rational(1) - input;
    const RealInterval complementInterval = exactInterval(complement, workBits);
    const RealInterval gammaComplement = encloseExp(
        encloseLogGammaPositive(complementInterval, workBits, &complement),
        workBits).interval;

    // sin(Pi*x)=sin(2*Pi*(x/2))。exact Rational turnsを使い，Pi*xのinterval化と
    // angle-reductionの情報損失を避ける。
    const RealInterval sine = encloseSinTurns(input / rational(2), workBits).interval;
    if (sine.containsZero())
        throw PrecisionInsufficient{"Gamma reflection cannot yet exclude a non-positive-integer pole"};

    const RealInterval pi = enclosePi(workBits).interval;
    const RealInterval denominator = multiply(sine, gammaComplement, workBits);
    if (denominator.containsZero())
        throw PrecisionInsufficient{"Gamma reflection denominator cannot yet be proven nonzero"};
    return divide(pi, denominator, workBits).roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval twoOverSqrtPi(std::size_t bits) {
    const RealInterval rootPi = encloseSqrt(enclosePi(bits).interval, bits).interval;
    return divide(exactInterval(2, bits), rootPi, bits);
}


[[nodiscard]] Rational unsignedRational(std::size_t value) {
    if (value > static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max()))
        throw std::overflow_error("Fresnel series index is too large");
    return Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(value))};
}

[[nodiscard]] RealInterval symmetricError(
    const Rational& bound,
    std::size_t precisionBits) {
    return RealInterval::fromRationalBounds(-bound, bound, precisionBits);
}

[[nodiscard]] Rational intervalAbsUpper(
    const RealInterval& interval,
    std::size_t precisionBits) {
    return absoluteInterval(interval, precisionBits).upper().toRational();
}

[[nodiscard]] RealInterval pointFresnelSeriesPositive(
    const Rational& x,
    bool cosineIntegral,
    std::size_t precisionBits) {
    // FresnelのMaclaurin級数は大きいxでは巨大な中間項が相殺する。
    // 旧+48bit固定guardではx=4程度でも高精度時に包含幅が縮まらないため、
    // x^2に比例したguardを追加して相殺分を明示的に吸収する。
    const Rational x2ForGuard = x * x;
    const BigInt guardQuotient = x2ForGuard.numerator() / x2ForGuard.denominator();
    const auto guardMagnitude = numeric::tryToUint64(guardQuotient);
    const std::size_t cancellationGuard = guardMagnitude
        ? static_cast<std::size_t>(std::min<std::uint64_t>(*guardMagnitude, 100'000ULL)) * 4U
        : 400'000U;
    const std::size_t workBits = checkedAdd(
        precisionBits,
        checkedAdd(64, cancellationGuard, "Fresnel cancellation guard is too large"),
        "Fresnel working precision is too large");
    const RealInterval pi = enclosePi(workBits).interval;
    const Rational piUpper = pi.upper().toRational();
    const Rational x2 = x * x;
    const Rational x4 = x2 * x2;
    const Rational commonUpper = piUpper * piUpper * x4 / rational(4);
    const RealInterval common = divide(
        multiply(multiply(pi, pi, workBits), exactInterval(x4, workBits), workBits),
        exactInterval(4, workBits), workBits);

    RealInterval term = cosineIntegral
        ? exactInterval(x, workBits)
        : divide(
            multiply(pi, exactInterval(x * x2, workBits), workBits),
            exactInterval(6, workBits), workBits);
    RealInterval sum = term;
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 12, "Fresnel target precision is too large"));

    constexpr std::size_t maximumTerms = 1'000'000;
    for (std::size_t n = 0; n < maximumTerms; ++n) {
        consumeCertifiedWork();
        const std::size_t numeratorIndex = cosineIntegral ? 4 * n + 1 : 4 * n + 3;
        const std::size_t d0 = cosineIntegral ? 2 * n + 1 : 2 * n + 2;
        const std::size_t d1 = cosineIntegral ? 2 * n + 2 : 2 * n + 3;
        const std::size_t d2 = cosineIntegral ? 4 * n + 5 : 4 * n + 7;
        const Rational ratioUpper = commonUpper * unsignedRational(numeratorIndex)
            / (unsignedRational(d0) * unsignedRational(d1) * unsignedRational(d2));

        // 現項より後の比が1未満に入れば以後は単調減少する。
        // 最初の未加算項を等比級数で上から押さえ、Taylor剰余を明示的に区間へ足す。
        if (ratioUpper < rational(1)) {
            const Rational nextBound = intervalAbsUpper(term, workBits) * ratioUpper;
            const Rational tailBound = nextBound / (rational(1) - ratioUpper);
            if (tailBound <= target) {
                sum = add(sum, symmetricError(tailBound, workBits), workBits);
                return sum.roundedOutward(precisionBits);
            }
        }

        const RealInterval ratio = divide(
            multiply(common, exactInterval(unsignedRational(numeratorIndex), workBits), workBits),
            exactInterval(unsignedRational(d0) * unsignedRational(d1) * unsignedRational(d2), workBits),
            workBits);
        term = negate(multiply(term, ratio, workBits));
        sum = add(sum, term, workBits);
    }
    throw std::overflow_error("Fresnel series requires too many terms");
}

struct FresnelPair final {
    RealInterval c;
    RealInterval s;
};

[[nodiscard]] FresnelPair pointFresnelAsymptoticPositive(
    const Rational& x,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 64, "Fresnel asymptotic precision is too large");
    const RealInterval pi = enclosePi(workBits).interval;
    const RealInterval xInterval = exactInterval(x, workBits);
    const RealInterval x2 = exactInterval(x * x, workBits);
    const RealInterval phase = divide(
        multiply(pi, x2, workBits), exactInterval(2, workBits), workBits);
    const RealInterval sine = encloseSinRadianInterval(phase, workBits).interval;
    const RealInterval cosine = encloseCosRadianInterval(phase, workBits).interval;

    RealInterval amplitude = divide(
        exactInterval(1, workBits),
        multiply(pi, xInterval, workBits), workBits);
    RealInterval tailReal = exactInterval(0, workBits);
    RealInterval tailImag = exactInterval(0, workBits);
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 14, "Fresnel asymptotic target precision is too large"));

    constexpr std::size_t maximumTerms = 4096;
    Rational remainderBound;
    for (std::size_t m = 0; m < maximumTerms; ++m) {
        consumeCertifiedWork();
        RealInterval realPart = exactInterval(0, workBits);
        RealInterval imagPart = exactInterval(0, workBits);
        switch (m & 3U) {
        case 0: // i A_m
            realPart = negate(multiply(sine, amplitude, workBits));
            imagPart = multiply(cosine, amplitude, workBits);
            break;
        case 1: // +A_m
            realPart = multiply(cosine, amplitude, workBits);
            imagPart = multiply(sine, amplitude, workBits);
            break;
        case 2: // -i A_m
            realPart = multiply(sine, amplitude, workBits);
            imagPart = negate(multiply(cosine, amplitude, workBits));
            break;
        case 3: // -A_m
            realPart = negate(multiply(cosine, amplitude, workBits));
            imagPart = negate(multiply(sine, amplitude, workBits));
            break;
        }
        tailReal = add(tailReal, realPart, workBits);
        tailImag = add(tailImag, imagPart, workBits);

        // m+1項まで展開した部分積分公式の剰余は、最後に加えたA_m以下。
        // x>=4ではA_mが必要精度まで減少する範囲で打ち切るため、発散域へ進まない。
        remainderBound = intervalAbsUpper(amplitude, workBits);
        if (remainderBound <= target) {
            const RealInterval error = symmetricError(remainderBound, workBits);
            const RealInterval half = exactInterval(rational(1, 2), workBits);
            return FresnelPair{
                add(subtract(half, tailReal, workBits), error, workBits).roundedOutward(precisionBits),
                add(subtract(half, tailImag, workBits), error, workBits).roundedOutward(precisionBits)};
        }

        const Rational odd = unsignedRational(2 * m + 1);
        const RealInterval scale = divide(
            exactInterval(odd, workBits),
            multiply(pi, x2, workBits), workBits);
        const RealInterval nextAmplitude = multiply(amplitude, scale, workBits);
        if (intervalAbsUpper(nextAmplitude, workBits) >= remainderBound)
            break;
        amplitude = nextAmplitude;
    }
    throw PrecisionInsufficient{"Fresnel asymptotic expansion did not reach the requested precision"};
}

[[nodiscard]] RealInterval pointFresnel(
    Rational x,
    bool cosineIntegral,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactInterval(0, precisionBits);
    if (x.numerator().isNegative())
        return negate(pointFresnel(-x, cosineIntegral, precisionBits));

    if (x < rational(8))
        return pointFresnelSeriesPositive(x, cosineIntegral, precisionBits);

    // 漸近級数は固定xで任意精度まで収束する級数ではない。
    // 要求精度に届かない場合は正則なMaclaurin級数へ戻し、速度のために保証を捨てない。
    try {
        const FresnelPair pair = pointFresnelAsymptoticPositive(x, precisionBits);
        return cosineIntegral ? pair.c : pair.s;
    }
    catch (const PrecisionInsufficient&) {
        return pointFresnelSeriesPositive(x, cosineIntegral, precisionBits);
    }
}

[[nodiscard]] RealInterval encloseFresnelRealImpl(
    const RealInterval& input,
    bool cosineIntegral,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Fresnel interval precision is too large");
    const Rational lower = input.lower().toRational();
    const Rational upper = input.upper().toRational();
    RealInterval value = pointFresnel(lower, cosineIntegral, workBits);

    // |C'(x)|=|cos(pi x^2/2)|<=1, |S'(x)|<=1。
    // 入力が丸め区間でもlower endpointからの距離だけ膨らませれば真値を必ず包含できる。
    const Rational width = upper - lower;
    if (!width.isZero())
        value = add(value, symmetricError(width, workBits), workBits);
    return value.roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval pointErfSeriesPositive(
    const Rational& x,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "erf working precision is too large");
    const Rational threshold = binaryThreshold(checkedAdd(
        workBits, 12, "erf precision is too large"));
    const Rational x2 = x * x;

    Rational term = x;
    Rational sum = term;
    Rational tailBound;
    std::uint64_t n = 0;
    for (;;) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        const std::uint64_t next = n + 1;
        const Rational ratio = -x2
            * Rational{BigInt::parse(std::to_string(2 * n + 1))}
            / Rational{
                BigInt::parse(std::to_string(next))
                * BigInt::parse(std::to_string(2 * n + 3))};
        const Rational nextTerm = term * ratio;

        // 次項以降の絶対比は x^2/(n+2) より小さい。
        const Rational q = x2 / Rational{BigInt::parse(std::to_string(n + 2))};
        if (q < rational(1)) {
            tailBound = absRational(nextTerm) / (rational(1) - q);
            if (tailBound <= threshold)
                break;
        }

        sum += nextTerm;
        term = nextTerm;
        n = next;
        if (n > 1'000'000)
            throw std::overflow_error("erf series did not converge within the iteration limit");
    }

    const RealInterval sumWithTail = RealInterval::fromRationalBounds(
        sum - tailBound, sum + tailBound, workBits);
    return multiply(sumWithTail, twoOverSqrtPi(workBits), workBits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] std::optional<RealInterval> pointErfcAsymptoticPositive(
    const Rational& x,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "erfc working precision is too large");
    const Rational threshold = binaryThreshold(checkedAdd(
        workBits, 12, "erfc precision is too large"));
    const Rational x2 = x * x;

    const RealInterval xInterval = exactInterval(x, workBits);
    const RealInterval exponential = encloseExp(
        exactInterval(-x2, workBits), workBits).interval;
    const RealInterval rootPi = encloseSqrt(enclosePi(workBits).interval, workBits).interval;
    const RealInterval prefactor = divide(
        exponential, multiply(xInterval, rootPi, workBits), workBits);

    Rational term{BigInt{1}};
    Rational sum = term;
    Rational previousMagnitude = absRational(term);
    for (std::uint64_t n = 0; n < 1'000'000; ++n) {
        consumeCertifiedWork();
        const Rational nextTerm = -term
            * Rational{BigInt::parse(std::to_string(2 * n + 1))}
            / (rational(2) * x2);
        const Rational nextMagnitude = absRational(nextTerm);
        const RealInterval error = multiply(
            prefactor, exactInterval(nextMagnitude, workBits), workBits);
        if (error.upper().toRational() <= threshold) {
            const RealInterval series = RealInterval::fromRationalBounds(
                sum - nextMagnitude, sum + nextMagnitude, workBits);
            return multiply(prefactor, series, workBits).roundedOutward(precisionBits);
        }

        // 漸近級数は最小項を越えると発散する。要求精度へ届く前に項が増加へ転じた場合はMaclaurin側へfallbackする。
        if (nextMagnitude >= previousMagnitude)
            return std::nullopt;

        sum += nextTerm;
        term = nextTerm;
        previousMagnitude = nextMagnitude;
    }
    return std::nullopt;
}

[[nodiscard]] RealInterval pointErf(
    const Rational& x,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactInterval(0, precisionBits);
    if (x.numerator().isNegative())
        return negate(pointErf(-x, precisionBits));

    if (x < rational(4))
        return pointErfSeriesPositive(x, precisionBits);

    if (const auto erfc = pointErfcAsymptoticPositive(x, precisionBits))
        return subtract(exactInterval(1, precisionBits), *erfc, precisionBits);
    return pointErfSeriesPositive(x, precisionBits);
}

[[nodiscard]] RealInterval pointErfc(
    const Rational& x,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactInterval(1, precisionBits);
    if (x.numerator().isNegative())
        return subtract(exactInterval(2, precisionBits), pointErfc(-x, precisionBits), precisionBits);
    if (x < rational(4))
        return subtract(exactInterval(1, precisionBits), pointErfSeriesPositive(x, precisionBits), precisionBits);
    if (const auto asymptotic = pointErfcAsymptoticPositive(x, precisionBits))
        return *asymptotic;
    return subtract(exactInterval(1, precisionBits), pointErfSeriesPositive(x, precisionBits), precisionBits);
}

enum class LambertProductOrder {
    BelowTarget,
    AboveTarget
};

// f(w)=w exp(w) を point w で保証評価し，exact Rational target との大小を証明する。
// enclosureがtargetと重なる場合だけguard precisionを増やし，所定回数で分離できなければ
// 推測せずPrecisionInsufficientとして上位のadaptive precisionへ返す。
[[nodiscard]] LambertProductOrder compareLambertProduct(
    const Rational& w,
    const Rational& target,
    std::size_t precisionBits) {
    std::size_t workBits = checkedAdd(
        precisionBits, 24, "Lambert W comparison precision is too large");
    for (std::size_t attempt = 0; attempt < 12; ++attempt) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        const RealInterval wInterval = exactInterval(w, workBits);
        const RealInterval exponential = encloseExp(wInterval, workBits).interval;
        const RealInterval product = multiply(wInterval, exponential, workBits);
        if (product.upper().toRational() < target)
            return LambertProductOrder::BelowTarget;
        if (product.lower().toRational() > target)
            return LambertProductOrder::AboveTarget;
        workBits = checkedAdd(
            workBits, std::max<std::size_t>(32, workBits / 2),
            "Lambert W comparison precision is too large");
    }
    throw PrecisionInsufficient{"Lambert W product comparison requires more precision"};
}

[[nodiscard]] RealInterval pointLambertWReal(
    const Rational& target,
    int branch,
    std::size_t precisionBits) {
    if (branch == 0 && target.isZero())
        return exactInterval(0, precisionBits);

    Rational lower;
    Rational upper;
    if (branch == 0) {
        if (target > rational(0)) {
            lower = rational(0);
            // W_0(x) <= x for x>=0 because exp(W_0(x))>=1.
            upper = target;
        }
        else {
            lower = rational(-1);
            upper = rational(0);
        }
    }
    else if (branch == -1) {
        upper = rational(-1);
        lower = rational(-2);
        // W_-1(x) -> -infinity as x -> 0-.  Find a certified left bracket by
        // doubling its magnitude until f(lower)>target on the decreasing branch.
        for (std::size_t i = 0; i < 64; ++i) {
            consumeCertifiedWork();
            if (compareLambertProduct(lower, target, precisionBits)
                == LambertProductOrder::AboveTarget)
                break;
            lower *= rational(2);
            if (i == 63)
                throw PrecisionInsufficient{"Lambert W lower branch bracket is too wide"};
        }
    }
    else {
        throw std::domain_error("Certified real Lambert W supports only branches 0 and -1");
    }

    const Rational targetWidth = binaryThreshold(checkedAdd(
        precisionBits, 8, "Lambert W target precision is too large"));
    const std::size_t maximumIterations = checkedAdd(
        precisionBits, 4096, "Lambert W iteration budget is too large");

    for (std::size_t iteration = 0; upper - lower > targetWidth; ++iteration) {
        consumeCertifiedWork();
        if (iteration >= maximumIterations)
            throw PrecisionInsufficient{"Lambert W bisection did not converge"};
        const Rational midpoint = (lower + upper) / rational(2);
        const LambertProductOrder order = compareLambertProduct(
            midpoint, target, precisionBits);

        if (branch == 0) {
            if (order == LambertProductOrder::BelowTarget)
                lower = midpoint;
            else
                upper = midpoint;
        }
        else {
            // f is strictly decreasing on (-infinity,-1].
            if (order == LambertProductOrder::AboveTarget)
                lower = midpoint;
            else
                upper = midpoint;
        }
    }

    return RealInterval::fromRationalBounds(lower, upper, precisionBits);
}

[[nodiscard]] RealInterval lambertBranchPoint(std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Lambert W branch-point precision is too large");
    // -1/e = -exp(-1).  exp backendの保証区間をそのまま反転してbranch domain判定に使う。
    return negate(encloseExp(exactInterval(-1, workBits), workBits).interval);
}

} // namespace


[[nodiscard]] std::optional<std::uint64_t> ceilAbsToUint64(const Rational& value) {
    const BigInt numerator = value.numerator().abs();
    const BigInt& denominator = value.denominator();
    if (numerator.isZero())
        return 0;
    const BigInt quotient = (numerator + denominator - BigInt{1}) / denominator;
    return numeric::tryToUint64(quotient);
}

[[nodiscard]] bool nonPositiveInteger(const Rational& value) {
    return value.isInteger() && !value.numerator().isPositive();
}

[[nodiscard]] RealInterval pointHypergeometric1F1(
    const Rational& a,
    const Rational& b,
    const Rational& z,
    std::size_t precisionBits) {
    if (nonPositiveInteger(b))
        throw std::domain_error("hypergeometric1F1 has a pole at a non-positive integer b");
    if (z.isZero())
        return exactInterval(1, precisionBits);

    const auto absA = ceilAbsToUint64(a);
    const auto absB = ceilAbsToUint64(b);
    const auto absZ = ceilAbsToUint64(z);
    if (!absA || !absB || !absZ)
        throw CertifiedBackendUnsupported{"hypergeometric1F1 argument is too large for the series backend"};

    // large |z|ではexact Rationalの項・部分和が巨大化し，級数項数より
    // Rational正規化が支配する。項は保証付きdyadic intervalで保持し，
    // zに比例するguardで値の指数成長・符号交代時の相殺を吸収する。
    const bool intervalFastPath = absRational(z) >= rational(16);

    // N>=2|a|,2|b|なら将来項比は 3|z|/(N+1) で上から押さえられる。
    // この上界が1未満へ入れば，残差は幾何級数として保証できる。
    // 固定の|z|境界は設けない。必要項数の保守上界がterm budgetへ収まる限り
    // series backendを試し，それを超える入力だけresource limitとして退く。
    constexpr std::uint64_t maximumTerms = 250000;

    const Rational tolerance = binaryThreshold(checkedAdd(
        precisionBits, 16, "hypergeometric1F1 precision is too large"));

    if (intervalFastPath) {
        if (*absA > maximumTerms / 2 || *absB > maximumTerms / 2
            || *absZ > maximumTerms / 3)
            throw CertifiedBackendUnsupported{
                "hypergeometric1F1 requires too many certified series terms"};
        const std::size_t magnitudeGuard = checkedAdd(
            static_cast<std::size_t>(*absZ) * 2U, 64,
            "hypergeometric1F1 interval-series guard is too large");
        const std::size_t workBits = checkedAdd(
            precisionBits, magnitudeGuard,
            "hypergeometric1F1 interval-series precision is too large");
        RealInterval term = exactInterval(1, workBits);
        RealInterval sum = term;

        for (std::uint64_t n = 0; n < maximumTerms; ++n) {
            consumeCertifiedWork();
            const Rational numeratorFactor = a + Rational{BigInt::fromUnsigned(n)};
            const Rational denominatorFactor = b + Rational{BigInt::fromUnsigned(n)};
            const Rational ratio = numeratorFactor * z
                / (denominatorFactor * Rational{BigInt::fromUnsigned(n + 1)});
            const RealInterval next = multiply(term, exactInterval(ratio, workBits), workBits);

            // N>=2|a|,2|b|なら，j>=Nで
            // |a+j|<=3j/2, |b+j|>=j/2 より将来項比は
            // 3|z|/(N+1)以下。符号・複素位相に依存しないmajorantなので，
            // large negative zの相殺もguard付きintervalで安全に扱える。
            const std::uint64_t Nvalue = n + 1;
            if (Nvalue >= 2U * *absA && Nvalue >= 2U * *absB
                && Nvalue >= 3U * *absZ && ((Nvalue - 3U * *absZ) & 15U) == 0) {
                const Rational N{BigInt::fromUnsigned(Nvalue)};
                const Rational Q = rational(3) * absRational(z) / (N + rational(1));
                if (Q < rational(1)) {
                    const Rational tailBound = absoluteInterval(next, workBits).upper().toRational()
                        / (rational(1) - Q);
                    if (tailBound <= tolerance)
                        return add(sum, symmetricError(tailBound, workBits), workBits)
                            .roundedOutward(precisionBits);
                }
            }

            term = next;
            sum = add(sum, term, workBits);
        }
        throw CertifiedBackendUnsupported{
            "hypergeometric1F1 interval series did not converge within the term limit"};
    }

    const std::uint64_t ratioStart = std::max({
        2U * *absA + 2U,
        2U * *absB + 2U,
        6U * *absZ + 2U});
    if (ratioStart > maximumTerms)
        throw CertifiedBackendUnsupported{"hypergeometric1F1 requires too many series terms"};

    Rational term{BigInt{1}};
    Rational sum{BigInt{1}};
    for (std::uint64_t n = 0; n < maximumTerms; ++n) {
        consumeCertifiedWork();
        const Rational numeratorFactor = a + Rational{BigInt::fromUnsigned(n)};
        const Rational denominatorFactor = b + Rational{BigInt::fromUnsigned(n)};
        if (denominatorFactor.isZero())
            throw std::domain_error("hypergeometric1F1 denominator parameter reaches a pole");

        const Rational next = term * numeratorFactor * z
            / (denominatorFactor * Rational{BigInt::fromUnsigned(n + 1)});

        if (n >= ratioStart) {
            const Rational tailBound = rational(2) * absRational(next);
            if (tailBound <= tolerance)
                return RealInterval::fromRationalBounds(
                    sum - tailBound, sum + tailBound, precisionBits);
        }

        term = next;
        sum += term;
        if (term.isZero())
            return exactInterval(sum, precisionBits);
    }

    throw CertifiedBackendUnsupported{"hypergeometric1F1 series did not converge within the term limit"};
}


[[nodiscard]] RealInterval pointHypergeometric2F1(
    const Rational& a,
    const Rational& b,
    const Rational& c,
    const Rational& z,
    std::size_t precisionBits) {
    if (nonPositiveInteger(c))
        throw std::domain_error("hypergeometric2F1 has a pole at a non-positive integer c");
    if (z.isZero())
        return exactInterval(1, precisionBits);

    const Rational absZ = absRational(z);
    if (absZ >= rational(1))
        throw CertifiedBackendUnsupported{
            "hypergeometric2F1 Gauss series currently requires |z| < 1"};

    const auto absA = ceilAbsToUint64(a);
    const auto absB = ceilAbsToUint64(b);
    const auto absC = ceilAbsToUint64(c);
    if (!absA || !absB || !absC)
        throw CertifiedBackendUnsupported{"hypergeometric2F1 parameter is too large for the series backend"};

    // exact RationalのGauss級数は反復ごとに分子・分母が肥大化し，
    // 小さい分母パラメータの近傍では数値結果だけが必要でも計算量を支配する。
    // tail bound用のパラメータはexactのまま保ち，級数和だけ外向き丸め区間で累積する。
    const std::size_t workBits = checkedAdd(
        precisionBits, 64, "hypergeometric2F1 working precision is too large");
    constexpr std::uint64_t maximumTerms = 250000;
    const Rational tolerance = binaryThreshold(checkedAdd(
        precisionBits, 20, "hypergeometric2F1 precision is too large"));
    RealInterval term = exactInterval(1, workBits);
    RealInterval sum = term;
    const RealInterval zInterval = exactInterval(z, workBits);

    const auto intervalAbsUpper = [](const RealInterval& value) {
        const Rational lower = absRational(value.lower().toRational());
        const Rational upper = absRational(value.upper().toRational());
        return lower > upper ? lower : upper;
    };

    for (std::uint64_t k = 0; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const Rational index{BigInt::fromUnsigned(k)};
        const Rational denominatorFactor = c + index;
        if (denominatorFactor.isZero())
            throw std::domain_error("hypergeometric2F1 denominator parameter reaches a pole");

        RealInterval next = multiply(
            term, exactInterval(a + index, workBits), workBits);
        next = multiply(next, exactInterval(b + index, workBits), workBits);
        next = multiply(next, zInterval, workBits);
        next = divide(next, exactInterval(denominatorFactor, workBits), workBits);
        next = divide(next,
            exactInterval(Rational{BigInt::fromUnsigned(k + 1)}, workBits), workBits);

        // j>=Nなら |a+j|<=j+A, |b+j|<=j+B, |c+j|>=j-C。
        // N>Cを満たす領域では、将来全ての項比を一つのq<1で押さえられる。
        const std::uint64_t n = k + 1;
        if (n > *absC + 1) {
            const Rational N{BigInt::fromUnsigned(n)};
            const Rational q = absZ
                * (Rational{BigInt{1}} + Rational{BigInt::fromUnsigned(*absA)} / N)
                * (Rational{BigInt{1}} + Rational{BigInt::fromUnsigned(*absB)} / N)
                / (Rational{BigInt{1}} - Rational{BigInt::fromUnsigned(*absC)} / N);
            if (q < rational(1)) {
                const Rational tail = intervalAbsUpper(next) / (rational(1) - q);
                if (tail <= tolerance) {
                    sum = add(sum, symmetricError(tail, workBits), workBits);
                    return sum.roundedOutward(precisionBits);
                }
            }
        }

        term = std::move(next);
        sum = add(sum, term, workBits);
    }
    throw CertifiedBackendUnsupported{"hypergeometric2F1 series did not converge within the term limit"};
}

[[nodiscard]] RealInterval intervalTimesRational(
    const RealInterval& interval,
    const Rational& value,
    std::size_t bits) {
    return multiply(interval, exactInterval(value, bits), bits);
}

[[nodiscard]] RealInterval evenSinePowerIntegral(
    std::size_t k,
    const RealInterval& cosine,
    RealInterval previous,
    RealInterval& sineOdd,
    const RealInterval& sineSquared,
    std::size_t bits) {
    if (k == 0)
        return previous;
    const Rational denominator{BigInt::fromUnsigned(2 * k)};
    const Rational recurrence{BigInt::fromUnsigned(2 * k - 1), BigInt::fromUnsigned(2 * k)};
    const RealInterval boundary = divide(
        multiply(sineOdd, cosine, bits), exactInterval(denominator, bits), bits);
    RealInterval result = subtract(
        intervalTimesRational(previous, recurrence, bits), boundary, bits);
    sineOdd = multiply(sineOdd, sineSquared, bits);
    return result;
}

enum class EllipticSeriesKind { F, E, Pi };


enum class CarlsonPower { Half, ThreeHalves };

[[nodiscard]] RealInterval carlsonEqualArgumentBounds(
    std::initializer_list<RealInterval> arguments,
    CarlsonPower power,
    std::size_t bits) {
    if (arguments.size() == 0)
        throw std::invalid_argument("Carlson residual requires at least one argument");

    Rational minimum = arguments.begin()->lower().toRational();
    Rational maximum = arguments.begin()->upper().toRational();
    for (const RealInterval& argument : arguments) {
        minimum = std::min(minimum, argument.lower().toRational());
        maximum = std::max(maximum, argument.upper().toRational());
    }
    if (minimum <= rational(0))
        throw CertifiedBackendUnsupported{
            "Carlson equal-argument residual requires positive arguments"};

    const auto valueAt = [&](const Rational& t) {
        const RealInterval point = exactInterval(t, bits);
        const RealInterval root = encloseSqrt(point, bits).interval;
        if (power == CarlsonPower::Half)
            return divide(exactInterval(1, bits), root, bits);
        return divide(exactInterval(1, bits), multiply(point, root, bits), bits);
    };

    const RealInterval atMaximum = valueAt(maximum);
    const RealInterval atMinimum = valueAt(minimum);
    return RealInterval::fromRationalBounds(
        atMaximum.lower().toRational(), atMinimum.upper().toRational(), bits);
}

[[nodiscard]] RealInterval carlsonRCPositive(
    const RealInterval& x,
    const RealInterval& y,
    std::size_t bits) {
    const Rational xLower = x.lower().toRational();
    const Rational xUpper = x.upper().toRational();
    const Rational yLower = y.lower().toRational();
    const Rational yUpper = y.upper().toRational();
    if (xLower <= rational(0) || yLower <= rational(0))
        throw CertifiedBackendUnsupported{"Carlson RC requires positive real arguments"};

    const auto nearEqualSeries = [&](const RealInterval& u, bool alternating) {
        const Rational uUpper = u.upper().toRational();
        if (u.lower().toRational() < rational(0) || uUpper > rational(1, 16))
            throw CertifiedBackendUnsupported{"Carlson RC near-equal series precondition failed"};

        const RealInterval one = exactInterval(1, bits);
        RealInterval sum = one;
        RealInterval power = one;
        const Rational target = binaryThreshold(bits > 8 ? bits - 8 : bits);
        Rational tailBound = rational(1);

        for (std::size_t k = 1; k < 4096; ++k) {
            consumeCertifiedWork();
            power = multiply(power, u, bits);
            const RealInterval term = divide(
                power, exactInterval(Rational{BigInt::fromUnsigned(2 * k + 1)}, bits), bits);
            if (alternating && (k & 1U))
                sum = subtract(sum, term, bits);
            else
                sum = add(sum, term, bits);

            const RealInterval nextPower = multiply(power, u, bits);
            const RealInterval nextTerm = divide(
                nextPower,
                exactInterval(Rational{BigInt::fromUnsigned(2 * k + 3)}, bits), bits);
            tailBound = nextTerm.upper().toRational();
            if (!alternating)
                tailBound /= rational(1) - uUpper;

            const Rational scale = std::max(rational(1), intervalAbsUpper(sum, bits));
            if (tailBound <= scale * target) {
                const RealInterval tail = RealInterval::fromRationalBounds(
                    alternating ? -tailBound : rational(0), tailBound, bits);
                const RealInterval quotient = add(sum, tail, bits);
                const RealInterval rootX = encloseSqrt(x, bits).interval;
                return divide(quotient, rootX, bits);
            }
        }
        throw CertifiedBackendUnsupported{"Carlson RC near-equal series did not converge"};
    };

    // DLMF 19.2.18--19.2.19。x≈yではdeltaを分母にも再利用する初等函数表示が
    // interval dependencyを増幅するため，u=(|x-y|/x)の正則級数へ切り替える。
    if (xUpper < yLower) {
        const RealInterval delta = subtract(y, x, bits);
        const RealInterval u = divide(delta, x, bits);
        if (u.upper().toRational() <= rational(1, 16))
            return nearEqualSeries(u, true);

        const RealInterval rootDelta = encloseSqrt(delta, bits).interval;
        const RealInterval rootRatio = encloseSqrt(u, bits).interval;
        return divide(encloseAtan(rootRatio, bits).interval, rootDelta, bits);
    }
    if (yUpper < xLower) {
        const RealInterval delta = subtract(x, y, bits);
        const RealInterval u = divide(delta, x, bits);
        if (u.upper().toRational() <= rational(1, 16))
            return nearEqualSeries(u, false);

        const RealInterval rootDelta = encloseSqrt(delta, bits).interval;
        const RealInterval rootX = encloseSqrt(x, bits).interval;
        const RealInterval rootY = encloseSqrt(y, bits).interval;
        const RealInterval numerator = add(rootX, rootDelta, bits);
        const RealInterval logarithm = encloseLogPositive(
            divide(numerator, rootY, bits), bits).interval;
        return divide(logarithm, rootDelta, bits);
    }

    return carlsonEqualArgumentBounds({x, y, y}, CarlsonPower::Half, bits);
}

[[nodiscard]] std::size_t carlsonDuplicationIterations(std::size_t precisionBits) {
    // 単調残差boundだけでも引数差は1回のduplicationで約1/4になる。
    // p bitに対してp/2回＋guardを上限とし，固定parameter thresholdへ逃げない。
    return checkedAdd(precisionBits / 2, 12, "Carlson iteration count is too large");
}

[[nodiscard]] RealInterval carlsonRFPositive(
    RealInterval x,
    RealInterval y,
    RealInterval z,
    std::size_t precisionBits) {
    const std::size_t bits = checkedAdd(
        precisionBits, 40, "Carlson RF working precision is too large");
    x = x.roundedOutward(bits);
    y = y.roundedOutward(bits);
    z = z.roundedOutward(bits);

    const auto nonnegative = [](const RealInterval& value) {
        return value.lower().toRational() >= rational(0);
    };
    if (!nonnegative(x) || !nonnegative(y) || !nonnegative(z))
        throw CertifiedBackendUnsupported{"Carlson RF real backend requires nonnegative arguments"};

    const RealInterval four = exactInterval(4, bits);
    const std::size_t maximumIterations = carlsonDuplicationIterations(precisionBits);
    for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration) {
        consumeCertifiedWork();
        const RealInterval sx = encloseSqrt(x, bits).interval;
        const RealInterval sy = encloseSqrt(y, bits).interval;
        const RealInterval sz = encloseSqrt(z, bits).interval;
        const RealInterval lambda = add(
            add(multiply(sx, sy, bits), multiply(sy, sz, bits), bits),
            multiply(sz, sx, bits), bits);
        x = divide(add(x, lambda, bits), four, bits);
        y = divide(add(y, lambda, bits), four, bits);
        z = divide(add(z, lambda, bits), four, bits);

        try {
            const RealInterval bound = carlsonEqualArgumentBounds(
                {x, y, z}, CarlsonPower::Half, bits);
            const Rational width = bound.upper().toRational() - bound.lower().toRational();
            const Rational scale = std::max(rational(1), intervalAbsUpper(bound, bits));
            if (width <= scale * binaryThreshold(checkedAdd(
                    precisionBits, 16, "Carlson RF target precision is too large")))
                return bound.roundedOutward(precisionBits);
        }
        catch (const CertifiedBackendUnsupported&) {
            // complete caseのx=0等は最初のduplicationで正領域へ移るため続行する。
        }
    }
    return carlsonEqualArgumentBounds(
        {x, y, z}, CarlsonPower::Half, bits).roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval carlsonRDPositive(
    RealInterval x,
    RealInterval y,
    RealInterval z,
    std::size_t precisionBits) {
    const std::size_t bits = checkedAdd(
        precisionBits, 40, "Carlson RD working precision is too large");
    x = x.roundedOutward(bits);
    y = y.roundedOutward(bits);
    z = z.roundedOutward(bits);
    if (x.lower().toRational() < rational(0)
        || y.lower().toRational() < rational(0)
        || z.lower().toRational() <= rational(0))
        throw CertifiedBackendUnsupported{"Carlson RD real backend requires x,y>=0 and z>0"};

    const RealInterval four = exactInterval(4, bits);
    RealInterval accumulated = exactInterval(0, bits);
    Rational weight = rational(1);
    const std::size_t maximumIterations = carlsonDuplicationIterations(precisionBits);
    for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration) {
        consumeCertifiedWork();
        const RealInterval sx = encloseSqrt(x, bits).interval;
        const RealInterval sy = encloseSqrt(y, bits).interval;
        const RealInterval sz = encloseSqrt(z, bits).interval;
        const RealInterval lambda = add(
            add(multiply(sx, sy, bits), multiply(sy, sz, bits), bits),
            multiply(sz, sx, bits), bits);

        // DLMF 19.26.20をscaled duplicationへ直すと
        // RD(x,y,z)=3/(sqrt(z)(z+lambda)) + RD(x',y',z')/4。
        const RealInterval correction = divide(
            exactInterval(3, bits),
            multiply(sz, add(z, lambda, bits), bits), bits);
        accumulated = add(
            accumulated,
            multiply(correction, exactInterval(weight, bits), bits), bits);
        weight *= rational(1, 4);

        x = divide(add(x, lambda, bits), four, bits);
        y = divide(add(y, lambda, bits), four, bits);
        z = divide(add(z, lambda, bits), four, bits);

        const RealInterval residual = carlsonEqualArgumentBounds(
            {x, y, z}, CarlsonPower::ThreeHalves, bits);
        const RealInterval result = add(
            accumulated,
            multiply(residual, exactInterval(weight, bits), bits), bits);
        const Rational width = result.upper().toRational() - result.lower().toRational();
        const Rational scale = std::max(rational(1), intervalAbsUpper(result, bits));
        if (width <= scale * binaryThreshold(checkedAdd(
                precisionBits, 16, "Carlson RD target precision is too large")))
            return result.roundedOutward(precisionBits);
    }

    const RealInterval residual = carlsonEqualArgumentBounds(
        {x, y, z}, CarlsonPower::ThreeHalves, bits);
    return add(accumulated,
        multiply(residual, exactInterval(weight, bits), bits), bits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval carlsonRJPositive(
    RealInterval x,
    RealInterval y,
    RealInterval z,
    RealInterval p,
    std::size_t precisionBits) {
    const std::size_t bits = checkedAdd(
        precisionBits, 48, "Carlson RJ working precision is too large");
    x = x.roundedOutward(bits);
    y = y.roundedOutward(bits);
    z = z.roundedOutward(bits);
    p = p.roundedOutward(bits);
    if (x.lower().toRational() < rational(0)
        || y.lower().toRational() < rational(0)
        || z.lower().toRational() < rational(0)
        || p.lower().toRational() <= rational(0))
        throw CertifiedBackendUnsupported{
            "Carlson RJ real backend requires x,y,z>=0 and p>0"};

    const RealInterval four = exactInterval(4, bits);
    RealInterval accumulated = exactInterval(0, bits);
    Rational weight = rational(1);
    const std::size_t maximumIterations = carlsonDuplicationIterations(precisionBits);
    for (std::size_t iteration = 0; iteration < maximumIterations; ++iteration) {
        consumeCertifiedWork();
        const RealInterval sx = encloseSqrt(x, bits).interval;
        const RealInterval sy = encloseSqrt(y, bits).interval;
        const RealInterval sz = encloseSqrt(z, bits).interval;
        const RealInterval sp = encloseSqrt(p, bits).interval;
        const RealInterval lambda = add(
            add(multiply(sx, sy, bits), multiply(sy, sz, bits), bits),
            multiply(sz, sx, bits), bits);

        // DLMF 19.26.22--23。RC補正を各段で厳密包含し，残りのRJだけ1/4へ縮小する。
        const RealInterval rootSum = add(add(sx, sy, bits), sz, bits);
        const RealInterval alpha = add(
            multiply(p, rootSum, bits),
            multiply(multiply(sx, sy, bits), sz, bits), bits);
        const RealInterval beta = multiply(sp, add(p, lambda, bits), bits);
        const RealInterval rc = carlsonRCPositive(
            squareInterval(alpha, bits), squareInterval(beta, bits), bits);
        const RealInterval correction = multiply(exactInterval(3, bits), rc, bits);
        accumulated = add(
            accumulated,
            multiply(correction, exactInterval(weight, bits), bits), bits);
        weight *= rational(1, 4);

        x = divide(add(x, lambda, bits), four, bits);
        y = divide(add(y, lambda, bits), four, bits);
        z = divide(add(z, lambda, bits), four, bits);
        p = divide(add(p, lambda, bits), four, bits);

        const RealInterval residual = carlsonEqualArgumentBounds(
            {x, y, z, p}, CarlsonPower::ThreeHalves, bits);
        const RealInterval result = add(
            accumulated,
            multiply(residual, exactInterval(weight, bits), bits), bits);
        const Rational width = result.upper().toRational() - result.lower().toRational();
        const Rational scale = std::max(rational(1), intervalAbsUpper(result, bits));
        if (width <= scale * binaryThreshold(checkedAdd(
                precisionBits, 16, "Carlson RJ target precision is too large")))
            return result.roundedOutward(precisionBits);
    }

    const RealInterval residual = carlsonEqualArgumentBounds(
        {x, y, z, p}, CarlsonPower::ThreeHalves, bits);
    return add(accumulated,
        multiply(residual, exactInterval(weight, bits), bits), bits)
        .roundedOutward(precisionBits);
}

struct ReducedEllipticAmplitude final {
    BigInt periods;
    RealInterval reduced;
};

[[nodiscard]] ReducedEllipticAmplitude reduceEllipticAmplitude(
    const RealInterval& phi,
    std::size_t bits) {
    const RealInterval pi = enclosePi(bits).interval;
    const RealInterval ratio = divide(phi, pi, bits);
    const RealInterval shifted = add(ratio, exactInterval(rational(1, 2), bits), bits);
    const BigInt lower = numeric::floorToInteger(shifted.lower().toRational());
    const BigInt upper = numeric::floorToInteger(shifted.upper().toRational());
    if (lower != upper)
        throw PrecisionInsufficient{
            "Elliptic amplitude InformationEnclosure cannot determine the Pi-period reduction",
            PrecisionInsufficientKind::InputInformation};

    const RealInterval period = multiply(
        exactInterval(Rational{lower}, bits), pi, bits);
    return ReducedEllipticAmplitude{lower, subtract(phi, period, bits)};
}

[[nodiscard]] RealInterval reducedEllipticCarlsonWithReduction(
    EllipticSeriesKind kind,
    const RealInterval& n,
    const ReducedEllipticAmplitude& reduction,
    const RealInterval& m,
    std::size_t precisionBits) {
    const std::size_t bits = checkedAdd(
        precisionBits, 72, "elliptic Carlson working precision is too large");
    const std::size_t carlsonPrecision = checkedAdd(
        precisionBits, 24, "elliptic Carlson precision is too large");
    const RealInterval mWork = m.roundedOutward(bits);
    const RealInterval nWork = n.roundedOutward(bits);

    const bool hasPeriods = !reduction.periods.isZero();
    const Rational mLower = mWork.lower().toRational();
    const Rational mUpper = mWork.upper().toRational();
    const Rational nLower = nWork.lower().toRational();
    const Rational nUpper = nWork.upper().toRational();

    // m=1ではEだけが実軸上で有限に退化する。還元区間[-Pi/2,Pi/2]では
    // sqrt(1-sin(phi)^2)=cos(phi)>=0なのでE(phi|1)=sin(phi)，
    // Pi周期ごとの増分は2である。F/Piは端点にbranch/poleを持つため別扱いしない。
    const bool exactMOne = mWork.isPoint() && mLower == rational(1);
    if (kind == EllipticSeriesKind::E && exactMOne) {
        const RealInterval sine = encloseSinRadianInterval(reduction.reduced, bits).interval;
        if (!hasPeriods)
            return sine.roundedOutward(precisionBits);
        const Rational periodIncrement{reduction.periods * BigInt{2}};
        return add(sine, exactInterval(periodIncrement, bits), bits)
            .roundedOutward(precisionBits);
    }

    if (hasPeriods) {
        if (mLower < rational(1) && mUpper >= rational(1))
            throw PrecisionInsufficient{
                "Elliptic parameter InformationEnclosure cannot exclude a branch point in period reduction",
                PrecisionInsufficientKind::InputInformation};
        if (kind == EllipticSeriesKind::Pi
            && nLower < rational(1) && nUpper >= rational(1))
            throw PrecisionInsufficient{
                "Elliptic Pi InformationEnclosure cannot exclude a pole in period reduction",
                PrecisionInsufficientKind::InputInformation};
        if (mUpper >= rational(1)
            || (kind == EllipticSeriesKind::Pi && nUpper >= rational(1)))
            throw CertifiedBackendUnsupported{
                "real Carlson elliptic backend cannot cross a period containing a branch point or pole"};
    }

    const RealInterval sine = encloseSinRadianInterval(reduction.reduced, bits).interval;
    const RealInterval cosine = encloseCosRadianInterval(reduction.reduced, bits).interval;
    const RealInterval sineSquared = squareInterval(sine, bits);
    const RealInterval cosineSquared = squareInterval(cosine, bits);
    const RealInterval one = exactInterval(1, bits);
    const RealInterval y = subtract(one, multiply(mWork, sineSquared, bits), bits);
    if (y.lower().toRational() < rational(0)) {
        if (y.upper().toRational() >= rational(0))
            throw PrecisionInsufficient{
                "Elliptic InformationEnclosure cannot determine the real principal branch",
                PrecisionInsufficientKind::InputInformation};
        throw CertifiedBackendUnsupported{
            "Certified complex elliptic continuation is not implemented"};
    }

    const RealInterval rf = carlsonRFPositive(cosineSquared, y, one, carlsonPrecision);
    RealInterval reduced = multiply(sine, rf, bits);
    if (kind == EllipticSeriesKind::E) {
        const RealInterval rd = carlsonRDPositive(cosineSquared, y, one, carlsonPrecision);
        const RealInterval sineCubed = multiply(sine, sineSquared, bits);
        reduced = subtract(reduced,
            divide(multiply(multiply(mWork, sineCubed, bits), rd, bits),
                exactInterval(3, bits), bits), bits);
    }
    else if (kind == EllipticSeriesKind::Pi) {
        const RealInterval p = subtract(one, multiply(nWork, sineSquared, bits), bits);
        if (p.lower().toRational() <= rational(0)) {
            if (p.upper().toRational() > rational(0))
                throw PrecisionInsufficient{
                    "Elliptic Pi InformationEnclosure cannot exclude a pole",
                    PrecisionInsufficientKind::InputInformation};
            throw CertifiedBackendUnsupported{
                "real ellipticPi Carlson backend does not cross principal-value poles"};
        }
        const RealInterval rj = carlsonRJPositive(cosineSquared, y, one, p, carlsonPrecision);
        const RealInterval sineCubed = multiply(sine, sineSquared, bits);
        reduced = add(reduced,
            divide(multiply(multiply(nWork, sineCubed, bits), rj, bits),
                exactInterval(3, bits), bits), bits);
    }

    if (!hasPeriods)
        return reduced.roundedOutward(precisionBits);

    const RealInterval completeY = subtract(one, mWork, bits);
    if (completeY.lower().toRational() <= rational(0))
        throw CertifiedBackendUnsupported{
            "real complete elliptic Carlson backend requires m < 1"};
    const RealInterval completeRF = carlsonRFPositive(
        exactInterval(0, bits), completeY, one, carlsonPrecision);
    RealInterval complete = completeRF;
    if (kind == EllipticSeriesKind::E) {
        const RealInterval completeRD = carlsonRDPositive(
            exactInterval(0, bits), completeY, one, carlsonPrecision);
        complete = subtract(completeRF,
            divide(multiply(mWork, completeRD, bits), exactInterval(3, bits), bits), bits);
    }
    else if (kind == EllipticSeriesKind::Pi) {
        const RealInterval completeP = subtract(one, nWork, bits);
        if (completeP.lower().toRational() <= rational(0))
            throw CertifiedBackendUnsupported{
                "real complete ellipticPi Carlson backend requires n < 1"};
        const RealInterval completeRJ = carlsonRJPositive(
            exactInterval(0, bits), completeY, one, completeP, carlsonPrecision);
        complete = add(completeRF,
            divide(multiply(nWork, completeRJ, bits), exactInterval(3, bits), bits), bits);
    }

    const Rational periodMultiplier{reduction.periods * BigInt{2}};
    return add(reduced,
        multiply(complete, exactInterval(periodMultiplier, bits), bits), bits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval reducedEllipticCarlson(
    EllipticSeriesKind kind,
    const RealInterval& n,
    const RealInterval& phi,
    const RealInterval& m,
    std::size_t precisionBits) {
    const std::size_t bits = checkedAdd(
        precisionBits, 72, "elliptic Carlson working precision is too large");
    const ReducedEllipticAmplitude reduction = reduceEllipticAmplitude(
        phi.roundedOutward(bits), bits);
    return reducedEllipticCarlsonWithReduction(kind, n, reduction, m, precisionBits);
}

[[nodiscard]] ReducedEllipticAmplitude exactPiMultipleReduction(
    const Rational& piCoefficient,
    std::size_t bits) {
    const BigInt periods = numeric::floorToInteger(piCoefficient + rational(1, 2));
    const Rational reducedCoefficient = piCoefficient - Rational{periods};
    return ReducedEllipticAmplitude{periods, multiply(
        exactInterval(reducedCoefficient, bits), enclosePi(bits).interval, bits)};
}

[[nodiscard]] RealInterval pointEllipticSeries(
    EllipticSeriesKind kind,
    const Rational& n,
    const Rational& phi,
    const Rational& m,
    std::size_t precisionBits) {
    if (phi.isZero())
        return exactInterval(0, precisionBits);
    if (phi.numerator().isNegative())
        return negate(pointEllipticSeries(kind, n, -phi, m, precisionBits));

    const Rational absM = absRational(m);
    const Rational absN = absRational(n);

    // 係数・m^k・tail majorantをexact Rationalで反復すると、|m|,|n|→0.9で
    // 分子分母が項数に比例して肥大化する。級数値は元からinterval評価なので、
    // これらも固定work precisionの外向きintervalへ落として包含保証だけを維持する。
    const std::size_t bits = checkedAdd(
        precisionBits, 48, "elliptic working precision is too large");
    const auto sinBox = encloseSinRadian(phi, bits).interval;
    const auto cosBox = encloseCosRadian(phi, bits).interval;
    const RealInterval sinSquared = multiply(sinBox, sinBox, bits);
    RealInterval sineOdd = sinBox;
    RealInterval integral = exactInterval(phi, bits); // I_0(phi)=phi
    RealInterval sum = integral;

    const RealInterval mBox = exactInterval(m, bits);
    const RealInterval nBox = exactInterval(n, bits);
    RealInterval c = exactInterval(1, bits); // (1/2)_k/k!
    RealInterval e = exactInterval(1, bits); // coefficients of sqrt(1-x)
    RealInterval mPower = exactInterval(1, bits);
    RealInterval q = exactInterval(1, bits); // Pi combined coefficient

    const Rational r = kind == EllipticSeriesKind::Pi
        ? (absM < absN ? absN : absM) : absM;
    const Rational tolerance = binaryThreshold(checkedAdd(
        precisionBits, 18, "elliptic precision is too large"));

    // |phi|<=Pi/2では I_k(phi)=Integral[sin(t)^(2k),{0,phi}]
    // <= |phi| sin(|phi|)^(2k)。従来のI_k<=|phi|だけより遥かに鋭く、
    // m≈0.9でも小振幅なら実効収束率を m*sin(phi)^2 まで下げられる。
    Rational tailRatio = r;
    const Rational halfPiLower = divide(
        enclosePi(bits).interval, exactInterval(2, bits), bits).lower().toRational();
    if (phi <= halfPiLower) {
        const Rational sineSquaredBound = intervalAbsUpper(sinSquared, bits);
        const Rational candidate = r * sineSquaredBound;
        if (candidate < tailRatio)
            tailRatio = candidate;
    }
    // seriesをparameter値そのものではなく，実際のtail収束率で選別する。
    // 小振幅では|m|や|n|が1を越えてもr*sin(phi)^2が十分小さければ
    // 積分路上の級数は一様収束する。逆に収束率が遅い場合はCarlsonへ渡す。
    if (tailRatio > rational(9, 10))
        throw CertifiedBackendUnsupported{
            "elliptic series convergence is too slow for the fast path"};
    const RealInterval rBox = exactInterval(tailRatio, bits);
    const RealInterval gap = exactInterval(rational(1) - tailRatio, bits);
    const RealInterval phiMagnitude = exactInterval(absRational(phi), bits);
    const RealInterval tailScale = kind == EllipticSeriesKind::Pi
        ? divide(phiMagnitude, multiply(gap, gap, bits), bits)
        : divide(phiMagnitude, gap, bits);

    RealInterval rPower = rBox;
    constexpr std::size_t maximumTerms = 200000;
    constexpr std::size_t certificateStride = 8;
    for (std::size_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        integral = evenSinePowerIntegral(
            k, cosBox, integral, sineOdd, sinSquared, bits);
        mPower = multiply(mPower, mBox, bits);
        c = intervalTimesRational(
            c, Rational{BigInt::fromUnsigned(2 * k - 1), BigInt::fromUnsigned(2 * k)}, bits);

        RealInterval coefficient = exactInterval(0, bits);
        if (kind == EllipticSeriesKind::F) {
            coefficient = multiply(c, mPower, bits);
        }
        else if (kind == EllipticSeriesKind::E) {
            e = intervalTimesRational(
                e, Rational{BigInt{static_cast<std::int64_t>(2 * k) - 3},
                    BigInt::fromUnsigned(2 * k)}, bits);
            coefficient = multiply(e, mPower, bits);
        }
        else {
            q = add(multiply(nBox, q, bits), multiply(c, mPower, bits), bits);
            coefficient = q;
        }
        sum = add(sum, multiply(integral, coefficient, bits), bits);

        rPower = multiply(rPower, rBox, bits);
        if ((k % certificateStride) != 0)
            continue;

        RealInterval tail = multiply(tailScale, rPower, bits);
        if (kind == EllipticSeriesKind::Pi) {
            const RealInterval linear = subtract(
                exactInterval(Rational{BigInt::fromUnsigned(k + 2)}, bits),
                multiply(
                    exactInterval(Rational{BigInt::fromUnsigned(k + 1)}, bits),
                    rBox, bits),
                bits);
            tail = multiply(tail, linear, bits);
        }

        const Rational tailUpper = tail.upper().toRational();
        if (tailUpper <= tolerance)
            return add(sum, symmetricError(tailUpper, bits), bits).roundedOutward(precisionBits);
    }
    throw CertifiedBackendUnsupported{"elliptic series did not converge within the term limit"};
}

[[nodiscard]] RealInterval intervalEllipticSeries(
    EllipticSeriesKind kind,
    const RealInterval& n,
    const RealInterval& phi,
    const RealInterval& m,
    std::size_t precisionBits) {
    const std::size_t bits = checkedAdd(
        precisionBits, 56, "elliptic interval working precision is too large");
    const RealInterval nWork = n.roundedOutward(bits);
    const RealInterval phiWork = phi.roundedOutward(bits);
    const RealInterval mWork = m.roundedOutward(bits);

    const RealInterval absMInterval = absoluteInterval(mWork, bits);
    const RealInterval absNInterval = absoluteInterval(nWork, bits);
    const Rational absM = absMInterval.upper().toRational();
    const Rational absN = absNInterval.upper().toRational();
    const RealInterval sine = encloseSinRadianInterval(phiWork, bits).interval;
    const RealInterval cosine = encloseCosRadianInterval(phiWork, bits).interval;
    const RealInterval sineSquared = multiply(sine, sine, bits);
    RealInterval sineOdd = sine;
    RealInterval integral = phiWork;
    RealInterval sum = integral;

    RealInterval c = exactInterval(1, bits);
    RealInterval e = exactInterval(1, bits);
    RealInterval mPower = exactInterval(1, bits);
    RealInterval q = exactInterval(1, bits);
    const Rational r = kind == EllipticSeriesKind::Pi
        ? (absM < absN ? absN : absM) : absM;
    const Rational tolerance = binaryThreshold(checkedAdd(
        precisionBits, 18, "elliptic interval precision is too large"));

    Rational tailRatio = r;
    const Rational phiMagnitudeBound = intervalAbsUpper(phiWork, bits);
    const Rational halfPiLower = divide(
        enclosePi(bits).interval, exactInterval(2, bits), bits).lower().toRational();
    if (phiMagnitudeBound <= halfPiLower) {
        const Rational sineSquaredBound = intervalAbsUpper(sineSquared, bits);
        const Rational candidate = r * sineSquaredBound;
        if (candidate < tailRatio)
            tailRatio = candidate;
    }
    if (tailRatio > rational(9, 10))
        throw CertifiedBackendUnsupported{
            "elliptic interval series convergence is too slow for the fast path"};
    const RealInterval rBox = exactInterval(tailRatio, bits);
    const RealInterval gap = exactInterval(rational(1) - tailRatio, bits);
    const RealInterval phiMagnitude = exactInterval(phiMagnitudeBound, bits);
    const RealInterval tailScale = kind == EllipticSeriesKind::Pi
        ? divide(phiMagnitude, multiply(gap, gap, bits), bits)
        : divide(phiMagnitude, gap, bits);

    RealInterval rPower = rBox;
    constexpr std::size_t maximumTerms = 200000;
    constexpr std::size_t certificateStride = 8;
    for (std::size_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        integral = evenSinePowerIntegral(
            k, cosine, integral, sineOdd, sineSquared, bits);
        mPower = multiply(mPower, mWork, bits);
        c = intervalTimesRational(
            c, Rational{BigInt::fromUnsigned(2 * k - 1), BigInt::fromUnsigned(2 * k)}, bits);

        RealInterval coefficient = exactInterval(0, bits);
        if (kind == EllipticSeriesKind::F) {
            coefficient = multiply(c, mPower, bits);
        }
        else if (kind == EllipticSeriesKind::E) {
            e = intervalTimesRational(
                e, Rational{BigInt{static_cast<std::int64_t>(2 * k) - 3},
                    BigInt::fromUnsigned(2 * k)}, bits);
            coefficient = multiply(e, mPower, bits);
        }
        else {
            q = add(
                multiply(nWork, q, bits),
                multiply(c, mPower, bits),
                bits);
            coefficient = q;
        }
        sum = add(sum, multiply(integral, coefficient, bits), bits);

        rPower = multiply(rPower, rBox, bits);
        if ((k % certificateStride) != 0)
            continue;

        RealInterval tail = multiply(tailScale, rPower, bits);
        if (kind == EllipticSeriesKind::Pi) {
            const RealInterval linear = subtract(
                exactInterval(Rational{BigInt::fromUnsigned(k + 2)}, bits),
                multiply(
                    exactInterval(Rational{BigInt::fromUnsigned(k + 1)}, bits),
                    rBox, bits),
                bits);
            tail = multiply(tail, linear, bits);
        }

        const Rational tailUpper = tail.upper().toRational();
        if (tailUpper <= tolerance)
            return add(sum, symmetricError(tailUpper, bits), bits).roundedOutward(precisionBits);
    }
    throw CertifiedBackendUnsupported{
        "elliptic interval series did not converge within the term limit"};
}


[[nodiscard]] Rational eulerGammaRemainderBound(std::uint64_t n) {
    constexpr std::size_t omittedK = 64;
    const Rational coefficient = absRational(bernoulliEven(omittedK))
        / Rational{BigInt::fromUnsigned(2 * omittedK)};
    const BigInt nPower = numeric::pow(BigInt::fromUnsigned(n), 2 * omittedK);
    return coefficient / Rational{nPower};
}

[[nodiscard]] RealInterval encloseEulerGamma(std::size_t precisionBits) {
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "EulerGamma target precision is too large"));

    std::uint64_t n = 64;
    while (eulerGammaRemainderBound(n) > target) {
        consumeCertifiedWork();
        if (n >= (1ULL << 20))
            throw CertifiedBackendUnsupported{"EulerGamma backend does not support this precision yet"};
        n *= 2;
    }

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "EulerGamma working precision is too large");
    RealInterval harmonic = exactInterval(0, workBits);
    for (std::uint64_t k = 1; k <= n; ++k) {
        consumeCertifiedWork();
        harmonic = add(harmonic,
            exactInterval(Rational{BigInt{1}, BigInt::fromUnsigned(k)}, workBits), workBits);
    }

    const Rational nR{BigInt::fromUnsigned(n)};
    RealInterval result = subtract(
        harmonic,
        encloseLogPositive(exactInterval(nR, workBits), workBits).interval,
        workBits);
    result = subtract(result,
        exactInterval(Rational{BigInt{1}, BigInt::fromUnsigned(2 * n)}, workBits),
        workBits);

    const Rational inverseSquare = Rational{BigInt{1}}
        / (nR * nR);
    Rational inverseEven = inverseSquare;
    for (std::size_t k = 1; k < 64; ++k) {
        consumeCertifiedWork();
        const Rational coefficient = bernoulliEven(k)
            / Rational{BigInt::fromUnsigned(2 * k)};
        result = add(result,
            exactInterval(coefficient * inverseEven, workBits), workBits);
        inverseEven *= inverseSquare;
    }

    // Euler-Maclaurin剰余は次のBernoulli項の絶対値以下で押さえる。
    // 符号へ依存せず対称区間を足し、保証を優先する。
    result = add(result,
        symmetricError(eulerGammaRemainderBound(n), workBits), workBits);
    return result.roundedOutward(precisionBits);
}


struct SineCosineIntegralPair final {
    RealInterval si;
    RealInterval ci;
};

[[nodiscard]] RealInterval signedAsymptoticRemainder(
    const Rational& bound,
    bool positive,
    std::size_t bits) {
    if (bound < rational(0))
        throw std::invalid_argument("asymptotic remainder bound must be nonnegative");
    return positive
        ? RealInterval::fromRationalBounds(rational(0), bound, bits)
        : RealInterval::fromRationalBounds(-bound, rational(0), bits);
}

[[nodiscard]] std::optional<SineCosineIntegralPair>
pointSineCosineIntegralAsymptoticPositive(
    const Rational& x,
    std::size_t precisionBits) {
    if (x <= rational(0))
        throw std::invalid_argument("Si/Ci asymptotic backend requires x > 0");

    // DLMF 6.12.3--6.12.8。正実軸ではf,gの剰余は最初の未使用項以下で，
    // その項と同符号になる。したがって最適打切り前だけを使えば，
    // 発散漸近級数でも片側剰余を含む保証区間を直接構成できる。
    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "Si/Ci asymptotic precision is too large");
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 28, "Si/Ci asymptotic target precision is too large"));
    const Rational x2 = x * x;

    RealInterval fTerm = exactInterval(rational(1) / x, workBits);
    RealInterval gTerm = exactInterval(rational(1) / x2, workBits);
    RealInterval f = fTerm;
    RealInterval g = gTerm;
    Rational previousF = intervalAbsUpper(fTerm, workBits);
    Rational previousG = intervalAbsUpper(gTerm, workBits);

    constexpr std::size_t maximumTerms = 4096;
    for (std::size_t m = 0; m < maximumTerms; ++m) {
        consumeCertifiedWork();
        const Rational fRatio = -unsignedRational(2 * m + 1)
            * unsignedRational(2 * m + 2) / x2;
        const Rational gRatio = -unsignedRational(2 * m + 2)
            * unsignedRational(2 * m + 3) / x2;
        const RealInterval nextF = multiply(fTerm, exactInterval(fRatio, workBits), workBits);
        const RealInterval nextG = multiply(gTerm, exactInterval(gRatio, workBits), workBits);
        const Rational fBound = intervalAbsUpper(nextF, workBits);
        const Rational gBound = intervalAbsUpper(nextG, workBits);

        if (fBound <= target && gBound <= target) {
            const bool remainderPositive = ((m + 1) & 1U) == 0;
            f = add(f, signedAsymptoticRemainder(
                fBound, remainderPositive, workBits), workBits);
            g = add(g, signedAsymptoticRemainder(
                gBound, remainderPositive, workBits), workBits);

            const RealInterval xInterval = exactInterval(x, workBits);
            const RealInterval sine = encloseSinRadianInterval(xInterval, workBits).interval;
            const RealInterval cosine = encloseCosRadianInterval(xInterval, workBits).interval;
            const RealInterval halfPi = divide(
                enclosePi(workBits).interval, exactInterval(2, workBits), workBits);
            const RealInterval si = subtract(
                subtract(halfPi, multiply(f, cosine, workBits), workBits),
                multiply(g, sine, workBits), workBits);
            const RealInterval ci = subtract(
                multiply(f, sine, workBits),
                multiply(g, cosine, workBits), workBits);
            return SineCosineIntegralPair{
                si.roundedOutward(precisionBits),
                ci.roundedOutward(precisionBits)};
        }

        // 最小項を越えてから続けても漸近級数は改善しない。Taylorへfallbackする。
        if (fBound >= previousF || gBound >= previousG)
            return std::nullopt;

        f = add(f, nextF, workBits);
        g = add(g, nextG, workBits);
        fTerm = nextF;
        gTerm = nextG;
        previousF = fBound;
        previousG = gBound;
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<RealInterval> pointPositiveExponentialIntegralEiAsymptotic(
    const Rational& x,
    std::size_t precisionBits) {
    if (x <= rational(0))
        throw std::invalid_argument("positive Ei asymptotic backend requires x > 0");

    // DLMF 6.12.2。n項で打ち切った剰余は次項の(1+chi(n+1))倍以下。
    // bracketはEi(x)/(exp(x)/x)~1なので，ここを相対精度相当まで囲ってから
    // prefactorを掛ければ巨大なEi(x)でもsignificant-digits契約を保てる。
    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "positive Ei asymptotic precision is too large");
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 24, "positive Ei asymptotic target precision is too large"));
    const RealInterval xInterval = exactInterval(x, workBits);
    const RealInterval prefactor = divide(
        encloseExp(xInterval, workBits).interval, xInterval, workBits);

    RealInterval term = exactInterval(1, workBits);
    RealInterval sum = term;
    Rational previousMagnitude = rational(1);

    // chi(2)=2, chi(3)=3*pi/4。またchi(t+2)=chi(t)*(t+2)/(t+1)。
    RealInterval chiEven = exactInterval(2, workBits);
    RealInterval chiOdd = multiply(
        enclosePi(workBits).interval, exactInterval(rational(3, 4), workBits), workBits);

    constexpr std::size_t maximumTerms = 4096;
    for (std::size_t m = 0; m < maximumTerms; ++m) {
        consumeCertifiedWork();
        const Rational ratio = unsignedRational(m + 1) / x;
        const RealInterval next = multiply(term, exactInterval(ratio, workBits), workBits);
        const Rational nextMagnitude = intervalAbsUpper(next, workBits);
        RealInterval chi = (m & 1U) == 0 ? chiEven : chiOdd; // chi(m+2)
        const Rational factorUpper = add(
            exactInterval(1, workBits), chi, workBits).upper().toRational();
        const Rational remainder = nextMagnitude * factorUpper;
        if (remainder <= target) {
            const RealInterval bracket = add(
                sum, symmetricError(remainder, workBits), workBits);
            return multiply(prefactor, bracket, workBits).roundedOutward(precisionBits);
        }

        if (nextMagnitude >= previousMagnitude)
            return std::nullopt;

        sum = add(sum, next, workBits);
        term = next;
        previousMagnitude = nextMagnitude;

        if ((m & 1U) == 0) {
            // 次に必要なeven chiはchi(m+4)。
            chiEven = multiply(
                chiEven,
                exactInterval(unsignedRational(m + 4) / unsignedRational(m + 3), workBits),
                workBits);
        } else {
            chiOdd = multiply(
                chiOdd,
                exactInterval(unsignedRational(m + 4) / unsignedRational(m + 3), workBits),
                workBits);
        }
    }
    return std::nullopt;
}

[[nodiscard]] RealInterval pointNegativeExponentialIntegralEiAsymptotic(
    const Rational& positiveX,
    std::size_t precisionBits) {
    // x>0について E1(x)=exp(-x) integral_0^inf exp(-t)/(x+t) dt。
    // 1/(1+t/x)を有限幾何展開すると，n項後の剰余は符号が既知で
    // 絶対値 <= n!/x^(n+1) となる。したがって Ei(-x)=-E1(x) を
    // significant-digitsを失わず直接囲える。
    if (positiveX <= rational(0))
        throw std::invalid_argument("negative Ei asymptotic requires x > 0");

    const std::size_t workBits = checkedAdd(
        precisionBits, 48, "negative Ei asymptotic precision is too large");
    const Rational relativeTarget = binaryThreshold(checkedAdd(
        precisionBits, 20, "negative Ei asymptotic target precision is too large"));
    const Rational bracketTarget = relativeTarget / positiveX;

    RealInterval term = exactInterval(rational(1) / positiveX, workBits);
    RealInterval sum = term;
    constexpr std::uint64_t maximumTerms = 4096;
    for (std::uint64_t k = 0; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const Rational ratio = -Rational{BigInt::fromUnsigned(k + 1)} / positiveX;
        const RealInterval next = multiply(term, exactInterval(ratio, workBits), workBits);
        const Rational remainder = absoluteInterval(next, workBits).upper().toRational();
        if (remainder <= bracketTarget) {
            const RealInterval bracket = add(sum, symmetricError(remainder, workBits), workBits);
            const RealInterval exponential = encloseExp(
                exactInterval(-positiveX, workBits), workBits).interval;
            return negate(multiply(exponential, bracket, workBits))
                .roundedOutward(precisionBits);
        }
        term = next;
        sum = add(sum, term, workBits);

        // 漸近級数はk~xを越えると項が再増大する。そこまでに要求精度へ
        // 到達しない場合はTaylor側へ無理に戻さずbackend limitとして退く。
        if (Rational{BigInt::fromUnsigned(k + 2)} >= positiveX)
            break;
    }
    throw CertifiedBackendUnsupported{
        "negative Ei asymptotic expansion could not certify the requested precision"};
}

[[nodiscard]] RealInterval pointExponentialIntegralEi(
    Rational x,
    std::size_t precisionBits) {
    if (x.isZero())
        throw std::domain_error("Ei is undefined at zero");
    const Rational absX = absRational(x);
    if (x > rational(0) && x >= rational(16)) {
        if (const auto asymptotic = pointPositiveExponentialIntegralEiAsymptotic(x, precisionBits))
            return *asymptotic;
    }
    if (x < rational(0) && absX >= rational(16)) {
        try {
            return pointNegativeExponentialIntegralEiAsymptotic(absX, precisionBits);
        } catch (const CertifiedBackendUnsupported&) {
            // 固定precisionの漸近級数が最適打切り前に要求幅へ届かない場合は，
            // 収束Taylor級数へ戻す。これはbackend未対応ではなく算法選択の問題である。
        }
    }

    constexpr std::uint64_t maximumTerms = 200'000;
    const auto magnitude = ceilAbsToUint64(absX);
    if (!magnitude || *magnitude >= maximumTerms)
        throw CertifiedBackendUnsupported{
            "Ei requires too many certified series terms"};

    // EiのTaylor項は|x|付近まで増大し得る。exact Rationalで保持すると
    // 大引数で分子・分母が肥大化するため，値・tailともoutward intervalで運ぶ。
    const std::size_t cancellationBits = static_cast<std::size_t>(*magnitude)
        * (x < rational(0) ? 4U : 2U);
    const std::size_t workBits = checkedAdd(
        checkedAdd(precisionBits, 56, "Ei working precision is too large"),
        cancellationBits, "Ei cancellation precision is too large");
    RealInterval sum = add(
        encloseEulerGamma(workBits),
        encloseLogPositive(exactInterval(absX, workBits), workBits).interval,
        workBits);

    RealInterval term = exactInterval(x, workBits); // k=1: x/(1*1!)
    Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "Ei target precision is too large"));
    if (x < rational(0)) {
        // E1(a)=exp(-a) E[1/(a+T)] >= exp(-a)/(a+1), T~Exp(1)。
        // このlower boundでtailの絶対誤差をscaleし，微小なEi(-a)でも
        // requested significant digitsを0.0へ潰さない。
        const RealInterval exponential = encloseExp(exactInterval(x, workBits), workBits).interval;
        const Rational valueLower = exponential.lower().toRational()
            / (absX + rational(1));
        target *= valueLower;
    }
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        sum = add(sum, term, workBits);

        const BigInt kBig = BigInt::fromUnsigned(k);
        const BigInt nextBig = BigInt::fromUnsigned(k + 1);
        const Rational ratio = x * Rational{kBig, nextBig * nextBig};
        const RealInterval next = multiply(term, exactInterval(ratio, workBits), workBits);

        // |t_(j+1)/t_j|=|x|j/(j+1)^2 <= |x|/(j+1)。
        // j>=k以降の一様上界qが1未満ならtailを幾何級数で保証できる。
        const Rational q = absX / Rational{nextBig};
        if (q < rational(1) && (k & 7U) == 0) {
            const Rational nextBound = absoluteInterval(next, workBits).upper().toRational();
            const Rational tail = nextBound / (rational(1) - q);
            if (tail <= target)
                return add(sum, symmetricError(tail, workBits), workBits)
                    .roundedOutward(precisionBits);
        }
        term = next;
    }
    throw CertifiedBackendUnsupported{"Ei series did not converge within the term limit"};
}

[[nodiscard]] RealInterval pointSineIntegralSi(
    Rational x,
    std::size_t precisionBits) {
    if (x.isZero())
        return exactInterval(0, precisionBits);
    if (x.numerator().isNegative())
        return negate(pointSineIntegralSi(-x, precisionBits));

    if (x >= rational(32)) {
        if (const auto asymptotic = pointSineCosineIntegralAsymptoticPositive(x, precisionBits))
            return asymptotic->si;
    }

    constexpr std::uint64_t maximumTerms = 200'000;
    const auto magnitude = ceilAbsToUint64(x);
    if (!magnitude || *magnitude > maximumTerms * 2U - 4U)
        throw CertifiedBackendUnsupported{"Si requires too many certified series terms"};
    const std::size_t cancellationBits = static_cast<std::size_t>(*magnitude) * 2U;
    const std::size_t workBits = checkedAdd(
        checkedAdd(precisionBits, 48, "Si working precision is too large"),
        cancellationBits, "Si cancellation precision is too large");
    const Rational x2 = x * x;
    RealInterval term = exactInterval(x, workBits);
    RealInterval sum = term;
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "Si target precision is too large"));

    for (std::uint64_t k = 0; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const std::uint64_t a = 2 * k + 1;
        const std::uint64_t b = 2 * k + 2;
        const std::uint64_t c = 2 * k + 3;
        const Rational ratio = -x2 * Rational{BigInt::fromUnsigned(a)}
            / Rational{BigInt::fromUnsigned(c) * BigInt::fromUnsigned(c) * BigInt::fromUnsigned(b)};
        const RealInterval next = multiply(term, exactInterval(ratio, workBits), workBits);
        const Rational q = x2
            / Rational{BigInt::fromUnsigned(b) * BigInt::fromUnsigned(c)};
        if (q < rational(1) && (k & 7U) == 0) {
            const Rational nextBound = absoluteInterval(next, workBits).upper().toRational();
            const Rational tail = nextBound / (rational(1) - q);
            if (tail <= target)
                return add(sum, symmetricError(tail, workBits), workBits)
                    .roundedOutward(precisionBits);
        }
        term = next;
        sum = add(sum, term, workBits);
    }
    throw CertifiedBackendUnsupported{"Si series did not converge within the term limit"};
}

[[nodiscard]] RealInterval pointCosineIntegralCiPositive(
    const Rational& x,
    std::size_t precisionBits) {
    if (x <= rational(0))
        throw std::domain_error("Ci real certified backend requires x > 0");

    if (x >= rational(32)) {
        if (const auto asymptotic = pointSineCosineIntegralAsymptoticPositive(x, precisionBits))
            return asymptotic->ci;
    }

    constexpr std::uint64_t maximumTerms = 200'000;
    const auto magnitude = ceilAbsToUint64(x);
    if (!magnitude || *magnitude > maximumTerms * 2U - 4U)
        throw CertifiedBackendUnsupported{"Ci requires too many certified series terms"};
    const std::size_t cancellationBits = static_cast<std::size_t>(*magnitude) * 2U;
    const std::size_t workBits = checkedAdd(
        checkedAdd(precisionBits, 56, "Ci working precision is too large"),
        cancellationBits, "Ci cancellation precision is too large");
    RealInterval sum = add(
        encloseEulerGamma(workBits),
        encloseLogPositive(exactInterval(x, workBits), workBits).interval,
        workBits);

    const Rational x2 = x * x;
    RealInterval term = exactInterval(-(x2 / rational(4)), workBits); // k=1
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "Ci target precision is too large"));
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        sum = add(sum, term, workBits);
        const std::uint64_t a = 2 * k;
        const std::uint64_t b = 2 * k + 1;
        const std::uint64_t c = 2 * k + 2;
        const Rational ratio = -x2 * Rational{BigInt::fromUnsigned(a)}
            / Rational{BigInt::fromUnsigned(c) * BigInt::fromUnsigned(c) * BigInt::fromUnsigned(b)};
        const RealInterval next = multiply(term, exactInterval(ratio, workBits), workBits);
        const Rational q = x2
            / Rational{BigInt::fromUnsigned(b) * BigInt::fromUnsigned(c)};
        if (q < rational(1) && (k & 7U) == 0) {
            const Rational nextBound = absoluteInterval(next, workBits).upper().toRational();
            const Rational tail = nextBound / (rational(1) - q);
            if (tail <= target)
                return add(sum, symmetricError(tail, workBits), workBits)
                    .roundedOutward(precisionBits);
        }
        term = next;
    }
    throw CertifiedBackendUnsupported{"Ci series did not converge within the term limit"};
}

[[nodiscard]] RealInterval pointZetaGreaterThanOne(
    const Rational& s,
    std::size_t precisionBits);

[[nodiscard]] Rational polylogRationalPower(
    Rational base,
    std::size_t exponent) {
    Rational result{BigInt{1}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result *= base;
        exponent >>= 1U;
        if (exponent != 0)
            base *= base;
    }
    return result;
}

[[nodiscard]] std::optional<RealInterval> pointPolylogPositiveIntegerNearOne(
    std::uint64_t order,
    const RealInterval& z,
    std::size_t precisionBits) {
    // DLMF 25.12.12のs->n (n=2,3,...) 極限。
    // mu=log(z), 0<z<1ではbranch ambiguityがなく，
    // Li_n(e^mu)=sum_{k!=n-1} zeta(n-k) mu^k/k!
    //   + mu^(n-1)/(n-1)! (H_{n-1}-log(-mu))
    // と書ける。ここではdirect seriesが遅くなるz≈1だけをfast pathにする。
    if (order < 3 || order > 12
        || z.lower().toRational() <= rational(0)
        || z.upper().toRational() >= rational(1))
        return std::nullopt;

    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "polylog near-one working precision is too large");
    const RealInterval mu = encloseLogPositive(z.roundedOutward(workBits), workBits).interval;
    if (mu.upper() >= BigFloat{})
        return std::nullopt;

    const Rational muAbs = absoluteInterval(mu, workBits).upper().toRational();
    // これは能力境界ではなく，direct seriesとの内部dispatch条件である。
    // 実測上|mu|<=1/20ではdirect seriesの反復costを上回らず，
    // Bernoulli tailも急速に減少する。より離れた点と高orderは既存seriesへ戻す。
    if (muAbs > rational(1, 20))
        return std::nullopt;

    const RealInterval pi = enclosePi(workBits).interval;
    const Rational piLower = pi.lower().toRational();
    const Rational fourPiSquared = rational(4) * piLower * piLower;
    const Rational tailRatio = muAbs * muAbs / fourPiSquared;
    if (tailRatio >= rational(1))
        return std::nullopt;

    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "polylog near-one target precision is too large"));

    RealInterval result = exactInterval(0, workBits);
    RealInterval powerOverFactorial = exactInterval(1, workBits); // mu^k/k!

    // k=0,...,n-2: zeta(n-k)は通常のs>1 certified backendで独立に評価する。
    for (std::uint64_t k = 0; k + 1 < order; ++k) {
        consumeCertifiedWork();
        if (k != 0) {
            powerOverFactorial = multiply(
                powerOverFactorial,
                multiply(
                    mu,
                    exactInterval(Rational{BigInt{1}, BigInt::fromUnsigned(k)}, workBits),
                    workBits),
                workBits);
        }
        const Rational zetaArgument{BigInt::fromUnsigned(order - k)};
        result = add(
            result,
            multiply(
                pointZetaGreaterThanOne(zetaArgument, workBits),
                powerOverFactorial,
                workBits),
            workBits);
    }

    // k=n-1のzeta(1) poleはGamma項との極限でlog(-mu)へ相殺される。
    const std::uint64_t singularK = order - 1;
    powerOverFactorial = multiply(
        powerOverFactorial,
        multiply(
            mu,
            exactInterval(
                Rational{BigInt{1}, BigInt::fromUnsigned(singularK)}, workBits),
            workBits),
        workBits);

    Rational harmonic;
    for (std::uint64_t k = 1; k < order; ++k)
        harmonic += Rational{BigInt{1}, BigInt::fromUnsigned(k)};
    const RealInterval logMinusMu = encloseLogPositive(negate(mu), workBits).interval;
    result = add(
        result,
        multiply(
            powerOverFactorial,
            subtract(exactInterval(harmonic, workBits), logMinusMu, workBits),
            workBits),
        workBits);

    // k=n: zeta(0)=-1/2。
    powerOverFactorial = multiply(
        powerOverFactorial,
        multiply(
            mu,
            exactInterval(Rational{BigInt{1}, BigInt::fromUnsigned(order)}, workBits),
            workBits),
        workBits);
    result = add(
        result,
        multiply(powerOverFactorial, exactInterval(rational(-1, 2), workBits), workBits),
        workBits);

    std::uint64_t currentK = order;
    for (std::size_t r = 1; r <= bernoulliEvenLiterals.size(); ++r) {
        consumeCertifiedWork();
        const std::uint64_t k = order + static_cast<std::uint64_t>(2 * r - 1);
        while (currentK < k) {
            ++currentK;
            powerOverFactorial = multiply(
                powerOverFactorial,
                multiply(
                    mu,
                    exactInterval(
                        Rational{BigInt{1}, BigInt::fromUnsigned(currentK)}, workBits),
                    workBits),
                workBits);
        }

        // zeta(1-2r)=-B_2r/(2r)。偶数の負整数zetaは0なので明示的に飛ばす。
        const Rational coefficient = -bernoulliEven(r)
            / Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * r))};
        result = add(
            result,
            multiply(powerOverFactorial, exactInterval(coefficient, workBits), workBits),
            workBits);

        // |B_2r| = 2(2r)! zeta(2r)/(2pi)^(2r), zeta(2r)<2 を使う。
        // 次の非零項majorant M_{r+1}以降は比 <= (|mu|/(2pi))^2 なので
        // M_{r+1}/(1-q) で残差全体を包含できる。
        const std::size_t nextR = r + 1;
        const std::uint64_t nextK = order
            + static_cast<std::uint64_t>(2 * nextR - 1);
        const BigInt numeratorFactorial = numeric::factorial(
            static_cast<std::uint64_t>(2 * nextR - 1));
        const BigInt denominatorFactorial = numeric::factorial(nextK);
        const Rational majorant = Rational{BigInt{4} * numeratorFactorial, denominatorFactorial}
            * polylogRationalPower(muAbs, static_cast<std::size_t>(nextK))
            / polylogRationalPower(rational(2) * piLower, 2 * nextR);
        const Rational tail = majorant / (rational(1) - tailRatio);
        if (tail <= target)
            return add(result, symmetricError(tail, workBits), workBits)
                .roundedOutward(precisionBits);
    }

    // Bernoulli表内で証明できなければ，能力を落とさず既存direct seriesへfallbackする。
    return std::nullopt;
}

[[nodiscard]] RealInterval pointPolylogPositiveOrder(
    std::uint64_t order,
    Rational z,
    std::size_t precisionBits) {
    if (order == 0)
        throw std::invalid_argument("polylog certified series requires positive order");
    const Rational absZ = absRational(z);
    if (z.isZero())
        return exactInterval(0, precisionBits);

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "polylog working precision is too large");

    if (order >= 3) {
        if (const auto nearOne = pointPolylogPositiveIntegerNearOne(
                order, exactInterval(z, precisionBits), precisionBits))
            return *nearOne;
    }

    // DLMF 25.12.3: 負実軸では w=z/(z-1) が (0,1) に入り，branch cutを跨がない。
    // |z|がある程度大きい場合だけ変換して，小引数での不要な相殺は避ける。
    if (order == 2 && z < rational(-1, 3)) {
        const Rational transformed = z / (z - rational(1));
        const RealInterval logOneMinusZ = encloseLogPositive(
            exactInterval(rational(1) - z, workBits), workBits).interval;
        RealInterval result = multiply(logOneMinusZ, logOneMinusZ, workBits);
        result = multiply(result, exactInterval(rational(-1, 2), workBits), workBits);
        result = subtract(
            result, pointPolylogPositiveOrder(2, transformed, workBits), workBits);
        return result.roundedOutward(precisionBits);
    }

    if (absZ >= rational(1))
        throw CertifiedBackendUnsupported{
            "polylog certified power series requires |z| < 1"};

    // DLMF 25.12.6: 0<z<1 では Li_2(z) を 1-z 側へ反射できる。
    // unit circle近傍を何万項も直接足す代わりに，小さい1-zの級数へ移す。
    if (order == 2 && z > rational(1, 2) && z < rational(1)) {
        const Rational complement = rational(1) - z;
        const RealInterval pi = enclosePi(workBits).interval;
        RealInterval result = multiply(pi, pi, workBits);
        result = multiply(result, exactInterval(rational(1, 6), workBits), workBits);
        const RealInterval logZ = encloseLogPositive(exactInterval(z, workBits), workBits).interval;
        const RealInterval logComplement = encloseLogPositive(
            exactInterval(complement, workBits), workBits).interval;
        result = subtract(result, multiply(logZ, logComplement, workBits), workBits);
        result = subtract(
            result, pointPolylogPositiveOrder(2, complement, workBits), workBits);
        return result.roundedOutward(precisionBits);
    }
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "polylog target precision is too large"));
    RealInterval term = exactInterval(z, workBits); // k=1
    RealInterval sum = term;

    // z^kをexact Rationalとして保持するとunit circle近傍で分子・分母が急成長する。
    // Li_sの項比 z*(k/(k+1))^s だけを小さいRationalとして区間へ掛ける。
    constexpr std::uint64_t maximumTerms = 1'000'000;
    for (std::uint64_t k = 1; k + 1 < maximumTerms; ++k) {
        consumeCertifiedWork();
        const BigInt kPower = numeric::pow(BigInt::fromUnsigned(k), order);
        const BigInt nextPower = numeric::pow(BigInt::fromUnsigned(k + 1), order);
        const Rational ratio = z * Rational{kPower, nextPower};
        const RealInterval next = multiply(term, exactInterval(ratio, workBits), workBits);

        // 後続の項比の絶対値は常に|z|以下なので幾何級数でtailをmajorizeできる。
        const Rational nextBound = absoluteInterval(next, workBits).upper().toRational();
        const Rational tail = nextBound / (rational(1) - absZ);
        sum = add(sum, next, workBits);
        if (tail <= target)
            return add(sum, symmetricError(tail, workBits), workBits)
                .roundedOutward(precisionBits);
        term = next;
    }
    throw CertifiedBackendUnsupported{"polylog series did not converge within the term limit"};
}

[[nodiscard]] RealInterval positiveIntegerBasePower(
    std::uint64_t base,
    const Rational& exponent,
    std::size_t precisionBits) {
    if (base == 1)
        return exactInterval(1, precisionBits);
    if (exponent.isInteger()) {
        const BigInt& e = exponent.numerator();
        const auto magnitude = numeric::tryToUint64(e.abs());
        if (magnitude && *magnitude <= 100000) {
            const BigInt powered = numeric::pow(BigInt::fromUnsigned(base), *magnitude);
            const Rational exact = e.isNegative()
                ? Rational{BigInt{1}, powered}
                : Rational{powered};
            return exactInterval(exact, precisionBits);
        }
    }

    // 半整数指数はLog/Expを経由せず，exact整数冪とcertified sqrtへ分解する。
    if (exponent.denominator() == BigInt{2}) {
        const BigInt& numerator = exponent.numerator();
        const auto magnitude = numeric::tryToUint64(numerator.abs());
        if (magnitude && (*magnitude & 1U) != 0 && *magnitude <= 100000) {
            const BigInt integerPower = numeric::pow(
                BigInt::fromUnsigned(base), (*magnitude - 1) / 2);
            const RealInterval root = encloseSqrt(
                exactInterval(Rational{BigInt::fromUnsigned(base)}, precisionBits),
                precisionBits).interval;
            const RealInterval value = multiply(
                exactInterval(Rational{integerPower}, precisionBits), root, precisionBits);
            if (numerator.isNegative())
                return divide(exactInterval(1, precisionBits), value, precisionBits);
            return value;
        }
    }

    const RealInterval logarithm = encloseLogPositive(
        exactInterval(Rational{BigInt::fromUnsigned(base)}, precisionBits),
        precisionBits).interval;
    const RealInterval scaled = multiply(
        logarithm, exactInterval(exponent, precisionBits), precisionBits);
    return encloseExp(scaled, precisionBits).interval;
}

[[nodiscard]] std::uint64_t smallestFactor(std::uint64_t value) noexcept {
    if ((value & 1U) == 0)
        return value == 2 ? value : 2;
    for (std::uint64_t factor = 3; factor <= value / factor; factor += 2) {
        if (value % factor == 0)
            return factor;
    }
    return value;
}

[[nodiscard]] RealInterval pointZetaGreaterThanOne(
    const Rational& s,
    std::size_t precisionBits) {
    if (s == rational(1))
        throw std::domain_error("zeta has a pole at s = 1");
    if (s < rational(1))
        throw CertifiedBackendUnsupported{"zeta certified real backend currently requires s > 1"};
    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "zeta working precision is too large");
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "zeta target precision is too large"));

    std::uint64_t chosenN = 0;
    std::size_t chosenK = 0;
    Rational chosenRising;
    RealInterval chosenPower = exactInterval(0, workBits);
    for (std::uint64_t n = 8; n <= 4096 && chosenN == 0; n *= 2) {
        consumeCertifiedWork();
        const BigInt nBig = BigInt::fromUnsigned(n);
        const RealInterval inverseNSquared = exactInterval(
            Rational{BigInt{1}, nBig * nBig}, workBits);
        Rational rising = s * (s + rational(1)) * (s + rational(2));
        BigInt factorial = numeric::factorial(4);
        RealInterval power = positiveIntegerBasePower(n, -(s + rational(3)), workBits);

        for (std::size_t k = 2; k <= 64; ++k) {
            consumeCertifiedWork();
            const Rational coefficient = absRational(bernoulliEven(k))
                * rising / Rational{factorial};
            const Rational bound = coefficient * power.upper().toRational();
            if (bound <= target) {
                chosenN = n;
                chosenK = k;
                chosenRising = rising;
                chosenPower = power;
                break;
            }
            if (k == 64)
                break;

            const std::uint64_t first = static_cast<std::uint64_t>(2 * k - 1);
            rising *= s + Rational{BigInt::fromUnsigned(first)};
            rising *= s + Rational{BigInt::fromUnsigned(first + 1)};
            factorial *= BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k + 1));
            factorial *= BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k + 2));
            power = multiply(power, inverseNSquared, workBits);
        }
    }
    if (chosenN == 0)
        throw CertifiedBackendUnsupported{"zeta Euler-Maclaurin budget is insufficient for the requested precision"};

    RealInterval result = exactInterval(0, workBits);
    std::vector<std::optional<RealInterval>> dirichletPowers(chosenN);
    dirichletPowers[1] = exactInterval(1, workBits);
    for (std::uint64_t n = 1; n < chosenN; ++n) {
        consumeCertifiedWork();
        RealInterval term = exactInterval(1, workBits);
        if (n > 1) {
            const std::uint64_t factor = smallestFactor(n);
            if (factor == n) {
                term = positiveIntegerBasePower(n, -s, workBits);
            } else {
                term = multiply(
                    *dirichletPowers[factor],
                    *dirichletPowers[n / factor],
                    workBits);
            }
            dirichletPowers[n] = term;
        }
        result = add(result, term, workBits);
    }

    const RealInterval integralTail = divide(
        positiveIntegerBasePower(chosenN, rational(1) - s, workBits),
        exactInterval(s - rational(1), workBits), workBits);
    result = add(result, integralTail, workBits);
    result = add(result, multiply(
        exactInterval(rational(1, 2), workBits),
        positiveIntegerBasePower(chosenN, -s, workBits), workBits), workBits);

    const BigInt chosenNBig = BigInt::fromUnsigned(chosenN);
    const RealInterval inverseNSquared = exactInterval(
        Rational{BigInt{1}, chosenNBig * chosenNBig}, workBits);
    Rational rising = s;
    BigInt factorial = numeric::factorial(2);
    RealInterval power = positiveIntegerBasePower(chosenN, -(s + rational(1)), workBits);
    for (std::size_t k = 1; k < chosenK; ++k) {
        consumeCertifiedWork();
        const Rational coefficient = bernoulliEven(k) * rising / Rational{factorial};
        result = add(result, multiply(
            exactInterval(coefficient, workBits), power, workBits), workBits);

        const std::uint64_t first = static_cast<std::uint64_t>(2 * k - 1);
        rising *= s + Rational{BigInt::fromUnsigned(first)};
        rising *= s + Rational{BigInt::fromUnsigned(first + 1)};
        factorial *= BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k + 1));
        factorial *= BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k + 2));
        power = multiply(power, inverseNSquared, workBits);
    }

    const BigInt omittedFactorial = numeric::factorial(
        static_cast<std::uint64_t>(2 * chosenK));
    const Rational omittedCoefficient = absRational(bernoulliEven(chosenK))
        * chosenRising / Rational{omittedFactorial};
    const Rational remainderBound = omittedCoefficient * chosenPower.upper().toRational();
    result = add(result, symmetricError(remainderBound, workBits), workBits);
    return result.roundedOutward(precisionBits);
}

struct PsiPlan final {
    std::size_t shift = 0;
    std::size_t omittedK = 0;
    Rational remainderBound;
};

[[nodiscard]] PsiPlan choosePsiPlan(
    const Rational& input,
    std::size_t precisionBits,
    bool trigamma) {
    if (input <= rational(0))
        throw std::domain_error("psi certified backend requires a positive argument");
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "psi target precision is too large"));
    for (std::size_t shift = 0; shift <= 1'000'000; shift += 4) {
        consumeCertifiedWork();
        const Rational x = input + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(shift))};
        if (x < rational(4))
            continue;
        const Rational xSquared = x * x;
        Rational power = trigamma ? xSquared * x : xSquared;
        for (std::size_t k = 1; k <= 64; ++k) {
            consumeCertifiedWork();
            Rational bound = absRational(bernoulliEven(k));
            if (!trigamma)
                bound /= Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k))};
            bound /= power;
            if (bound <= target)
                return PsiPlan{shift, k, bound};
            power *= xSquared;
        }
    }
    throw CertifiedBackendUnsupported{"psi asymptotic budget is insufficient for the requested precision"};
}

[[nodiscard]] RealInterval pointDigammaPositive(
    const Rational& input,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "digamma working precision is too large");
    const PsiPlan plan = choosePsiPlan(input, workBits, false);
    const Rational x = input + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(plan.shift))};
    RealInterval result = encloseLogPositive(exactInterval(x, workBits), workBits).interval;
    result = subtract(result, exactInterval(rational(1, 2) / x, workBits), workBits);
    const Rational xSquared = x * x;
    Rational power = xSquared;
    for (std::size_t k = 1; k < plan.omittedK; ++k) {
        consumeCertifiedWork();
        Rational term = bernoulliEven(k)
            / Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k))};
        term /= power;
        result = subtract(result, exactInterval(term, workBits), workBits);
        power *= xSquared;
    }

    const Rational omitted = -bernoulliEven(plan.omittedK)
        / Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * plan.omittedK))}
        / power;
    const RealInterval remainder = omitted.numerator().isNegative()
        ? RealInterval::fromRationalBounds(-plan.remainderBound, rational(0), workBits)
        : RealInterval::fromRationalBounds(rational(0), plan.remainderBound, workBits);
    result = add(result, remainder, workBits);

    for (std::size_t j = 0; j < plan.shift; ++j) {
        consumeCertifiedWork();
        const Rational divisor = input + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(j))};
        result = subtract(result, exactInterval(rational(1) / divisor, workBits), workBits);
    }
    return result.roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval pointTrigammaPositive(
    const Rational& input,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "trigamma working precision is too large");
    const PsiPlan plan = choosePsiPlan(input, workBits, true);
    const Rational x = input + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(plan.shift))};
    RealInterval result = exactInterval(rational(1) / x, workBits);
    const Rational xSquared = x * x;
    result = add(result, exactInterval(rational(1, 2) / xSquared, workBits), workBits);
    Rational power = xSquared * x;
    for (std::size_t k = 1; k < plan.omittedK; ++k) {
        consumeCertifiedWork();
        const Rational term = bernoulliEven(k) / power;
        result = add(result, exactInterval(term, workBits), workBits);
        power *= xSquared;
    }
    result = add(result, symmetricError(plan.remainderBound, workBits), workBits);
    for (std::size_t j = 0; j < plan.shift; ++j) {
        consumeCertifiedWork();
        const Rational divisor = input + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(j))};
        result = add(result, exactInterval(rational(1) / (divisor * divisor), workBits), workBits);
    }
    return result.roundedOutward(precisionBits);
}

namespace {

[[nodiscard]] ComplexInterval exactComplex(
    const Rational& real,
    const Rational& imaginary,
    std::size_t precisionBits) {
    return ComplexInterval{
        exactInterval(real, precisionBits),
        exactInterval(imaginary, precisionBits)};
}

[[nodiscard]] ComplexInterval exactComplex(
    std::int64_t real,
    std::size_t precisionBits) {
    return exactComplex(rational(real), rational(0), precisionBits);
}


[[nodiscard]] Rational complexAbsUpper(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const RealInterval magnitudeSquared = add(
        squareInterval(value.real(), precisionBits),
        squareInterval(value.imaginary(), precisionBits),
        precisionBits);
    return encloseSqrt(magnitudeSquared, precisionBits).interval.upper().toRational();
}

[[nodiscard]] Rational intervalAbsLower(const RealInterval& value) {
    if (value.containsZero())
        return rational(0);
    const Rational lower = absRational(value.lower().toRational());
    const Rational upper = absRational(value.upper().toRational());
    return lower < upper ? lower : upper;
}

[[nodiscard]] Rational complexAbsLower(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const Rational realLower = intervalAbsLower(value.real());
    const Rational imaginaryLower = intervalAbsLower(value.imaginary());
    const Rational squared = realLower * realLower + imaginaryLower * imaginaryLower;
    return encloseSqrt(exactInterval(squared, precisionBits), precisionBits)
        .interval.lower().toRational();
}

[[nodiscard]] ComplexInterval inflateComplex(
    const ComplexInterval& value,
    const Rational& radius,
    std::size_t precisionBits) {
    if (radius.isZero())
        return value.roundedOutward(precisionBits);
    const RealInterval error = symmetricError(radius, precisionBits);
    return ComplexInterval{
        add(value.real(), error, precisionBits),
        add(value.imaginary(), error, precisionBits)};
}

[[nodiscard]] ComplexInterval multiplyByRational(
    const ComplexInterval& value,
    const Rational& factor,
    std::size_t precisionBits) {
    return multiply(
        value,
        ComplexInterval::fromReal(exactInterval(factor, precisionBits)),
        precisionBits);
}

[[nodiscard]] ComplexInterval addReal(
    const ComplexInterval& value,
    const RealInterval& real,
    std::size_t precisionBits) {
    return add(value, ComplexInterval::fromReal(real), precisionBits);
}

[[nodiscard]] Rational rationalPower(
    Rational base,
    std::size_t exponent) {
    Rational result{BigInt{1}};
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result *= base;
        exponent >>= 1U;
        if (exponent != 0)
            base *= base;
    }
    return result;
}

[[nodiscard]] RealInterval eulerGammaInterval(std::size_t precisionBits) {
    return encloseEulerGamma(precisionBits);
}

[[nodiscard]] bool crossesNegativeRealCut(const ComplexInterval& value) {
    const BigFloat zero;
    return value.real().lower() < zero && value.imaginary().containsZero();
}

[[nodiscard]] bool intervalContains(
    const RealInterval& outer,
    const RealInterval& inner) {
    return outer.lower() <= inner.lower() && outer.upper() >= inner.upper();
}

[[nodiscard]] bool intervalContains(
    const ComplexInterval& outer,
    const ComplexInterval& inner) {
    return intervalContains(outer.real(), inner.real())
        && intervalContains(outer.imaginary(), inner.imaginary());
}

// branch point -1/e の局所座標 u=W+1 では
//
//   q = e z + 1 = (u-1)exp(u)+1 = u^2 A(u)/2,
//   A(u) = sum_{m>=0} 2(m+1) u^m/(m+2)!.
//
// A(0)=1 なので p=sqrt(2q) を固定し，
//
//   u = +/- p / sqrt(A(u))
//
// と書けば branch point でも縮小写像が退化しない。sqrt(A) は1近傍だけを通るため，
// q自体のsqrt branch cutを反復中に跨ぐ問題も避けられる。
[[nodiscard]] ComplexInterval encloseLambertBranchFactor(
    const ComplexInterval& u,
    std::size_t precisionBits) {
    const Rational radius = complexAbsUpper(u, precisionBits);
    if (radius >= rational(3, 4))
        throw CertifiedBackendUnsupported{
            "Lambert W branch-point factor requires |W+1| < 3/4"};

    ComplexInterval term = exactComplex(1, precisionBits);
    ComplexInterval sum = term;
    const Rational ratioBound = rational(2, 3) * radius;
    const Rational target = binaryThreshold(
        precisionBits > 12 ? precisionBits - 12 : std::size_t{1});

    constexpr std::uint64_t maximumTerms = 4096;
    for (std::uint64_t m = 0; m < maximumTerms; ++m) {
        consumeCertifiedWork();
        // a_{m+1}/a_m = (m+2)/((m+1)(m+3)).
        const Rational ratio{
            BigInt::fromUnsigned(m + 2),
            BigInt::fromUnsigned(m + 1) * BigInt::fromUnsigned(m + 3)};
        term = multiplyByRational(
            multiply(term, u, precisionBits), ratio, precisionBits);
        sum = add(sum, term, precisionBits);
        const Rational termBound = complexAbsUpper(term, precisionBits);
        const Rational tail = ratioBound.isZero()
            ? rational(0)
            : termBound * ratioBound / (rational(1) - ratioBound);
        if (tail <= target)
            return inflateComplex(sum, tail, precisionBits);
    }
    throw CertifiedBackendUnsupported{
        "Lambert W branch-point factor series did not converge"};
}

[[nodiscard]] ComplexInterval encloseLambertBranchFactorDerivative(
    const ComplexInterval& u,
    std::size_t precisionBits) {
    const Rational radius = complexAbsUpper(u, precisionBits);
    if (radius >= rational(2, 3))
        throw CertifiedBackendUnsupported{
            "Lambert W branch-point derivative requires |W+1| < 2/3"};

    // A'(u)=sum b_j u^j, b_0=2/3,
    // b_{j+1}/b_j=(j+3)/((j+1)(j+4)).
    ComplexInterval term = exactComplex(rational(2, 3), rational(0), precisionBits);
    ComplexInterval sum = term;
    const Rational ratioBound = rational(3, 4) * radius;
    const Rational target = binaryThreshold(
        precisionBits > 12 ? precisionBits - 12 : std::size_t{1});

    constexpr std::uint64_t maximumTerms = 4096;
    for (std::uint64_t j = 0; j < maximumTerms; ++j) {
        consumeCertifiedWork();
        const Rational ratio{
            BigInt::fromUnsigned(j + 3),
            BigInt::fromUnsigned(j + 1) * BigInt::fromUnsigned(j + 4)};
        term = multiplyByRational(
            multiply(term, u, precisionBits), ratio, precisionBits);
        sum = add(sum, term, precisionBits);
        const Rational termBound = complexAbsUpper(term, precisionBits);
        const Rational tail = ratioBound.isZero()
            ? rational(0)
            : termBound * ratioBound / (rational(1) - ratioBound);
        if (tail <= target)
            return inflateComplex(sum, tail, precisionBits);
    }
    throw CertifiedBackendUnsupported{
        "Lambert W branch-point derivative series did not converge"};
}

[[nodiscard]] ComplexInterval encloseInverseSqrtNearOne(
    const ComplexInterval& value,
    std::size_t precisionBits) {
    const ComplexInterval one = exactComplex(1, precisionBits);
    const ComplexInterval delta = subtract(value, one, precisionBits);
    const Rational radius = complexAbsUpper(delta, precisionBits);
    if (radius >= rational(1, 2))
        throw CertifiedBackendUnsupported{
            "Lambert W branch-point inverse sqrt is outside its local series region"};

    ComplexInterval term = one;
    ComplexInterval sum = one;
    const Rational target = binaryThreshold(
        precisionBits > 12 ? precisionBits - 12 : std::size_t{1});

    constexpr std::uint64_t maximumTerms = 4096;
    for (std::uint64_t n = 0; n < maximumTerms; ++n) {
        consumeCertifiedWork();
        // (1+x)^(-1/2): a_{n+1}/a_n=-(2n+1)/(2n+2).
        const Rational ratio{
            -BigInt::fromUnsigned(2 * n + 1),
            BigInt::fromUnsigned(2 * n + 2)};
        term = multiplyByRational(
            multiply(term, delta, precisionBits), ratio, precisionBits);
        sum = add(sum, term, precisionBits);
        // |a_{k+1}/a_k|<1 なので以後の比をradiusで一様に抑えられる。
        // current term自体をoutward intervalで評価済みなので，別の巨大exact
        // Rational majorantを育てず，そのノルムから残りを直接majorizeする。
        const Rational termBound = complexAbsUpper(term, precisionBits);
        const Rational tail = radius.isZero()
            ? rational(0)
            : termBound * radius / (rational(1) - radius);
        if (tail <= target)
            return inflateComplex(sum, tail, precisionBits);
    }
    throw CertifiedBackendUnsupported{
        "Lambert W branch-point inverse sqrt series did not converge"};
}

[[nodiscard]] ComplexInterval mapLambertBranchPoint(
    const ComplexInterval& p,
    const ComplexInterval& u,
    bool negativeLocalBranch,
    std::size_t precisionBits) {
    const ComplexInterval factor = encloseLambertBranchFactor(u, precisionBits);
    if (factor.containsZero())
        throw CertifiedBackendUnsupported{
            "Lambert W branch-point factor may contain zero"};
    const ComplexInterval inverseRoot = encloseInverseSqrtNearOne(
        factor, precisionBits);
    ComplexInterval mapped = multiply(p, inverseRoot, precisionBits);
    if (negativeLocalBranch)
        mapped = negate(mapped);
    return mapped;
}

[[nodiscard]] ComplexInterval pointLambertWBranchPointContractionFromQ(
    const ComplexInterval& qInput,
    bool negativeLocalBranch,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 72, "Lambert W branch-point precision is too large");
    const ComplexInterval q = qInput.roundedOutward(workBits);
    const Rational qAbs = complexAbsUpper(q, workBits);
    // これはcapability boundaryではなく，局所backendを選ぶための安全なfast-path条件。
    // 外側は従来のexp/log contractionへfallbackする。
    if (qAbs >= rational(1, 16))
        throw CertifiedBackendUnsupported{
            "Lambert W input is outside the branch-point contraction neighborhood"};

    const ComplexInterval p = enclosePrincipalComplexSqrt(
        multiplyByRational(q, rational(2), workBits), workBits);
    const auto midpoint = [&](const ComplexInterval& value) {
        const Rational real = (value.real().lower().toRational()
            + value.real().upper().toRational()) * rational(1, 2);
        const Rational imaginary = (value.imaginary().lower().toRational()
            + value.imaginary().upper().toRational()) * rational(1, 2);
        return exactComplex(real, imaginary, workBits);
    };
    const ComplexInterval centerP = midpoint(p);
    ComplexInterval signedP = negativeLocalBranch ? negate(centerP) : centerP;

    // DLMF 4.13.9_1 の最初の項でseedし，証明前の候補だけNewtonで高精度化する。
    const ComplexInterval p2 = multiply(signedP, signedP, workBits);
    ComplexInterval power = p2;
    ComplexInterval candidate = subtract(
        signedP, multiplyByRational(power, rational(1, 3), workBits), workBits);
    const std::array<Rational, 10> puiseuxCoefficients{
        rational(11, 72), rational(-43, 540), rational(769, 17280),
        rational(-221, 8505), Rational{BigInt{680863}, BigInt{43545600}},
        rational(-1963, 204120), Rational{BigInt{226287557}, BigInt{37623398400LL}},
        Rational{-BigInt{5776369}, BigInt{1515591000}},
        Rational{BigInt{169709463197LL}, BigInt{69528040243200LL}},
        Rational{-BigInt{1118511313}, BigInt{709296588000LL}}};
    for (const Rational& coefficient : puiseuxCoefficients) {
        power = multiply(power, signedP, workBits);
        candidate = add(candidate,
            multiplyByRational(power, coefficient, workBits), workBits);
    }
    const ComplexInterval centerQ = midpoint(q);
    for (std::size_t iteration = 0; iteration < 4; ++iteration) {
        consumeCertifiedWork();
        candidate = midpoint(candidate);
        const ComplexInterval factor = encloseLambertBranchFactor(candidate, workBits);
        const ComplexInterval derivative = encloseLambertBranchFactorDerivative(candidate, workBits);
        const ComplexInterval u2 = multiply(candidate, candidate, workBits);
        const ComplexInterval residual = subtract(
            multiplyByRational(multiply(u2, factor, workBits), rational(1, 2), workBits),
            centerQ, workBits);
        const ComplexInterval slope = add(
            multiply(candidate, factor, workBits),
            multiplyByRational(multiply(u2, derivative, workBits), rational(1, 2), workBits),
            workBits);
        if (slope.containsZero())
            break;
        candidate = subtract(candidate, divide(residual, slope, workBits), workBits);
    }
    candidate = midpoint(candidate);

    const Rational targetWidth = binaryThreshold(checkedAdd(
        precisionBits, 8, "Lambert W branch-point target precision is too large"));
    Rational radius = complexAbsUpper(candidate, workBits) * targetWidth;
    if (radius.isZero())
        radius = targetWidth;

    for (std::size_t radiusAttempt = 0; radiusAttempt < 512 && radius < rational(1, 2);
         ++radiusAttempt, radius *= rational(2)) {
        ComplexInterval box = inflateComplex(candidate, radius, workBits);
        if (complexAbsUpper(box, workBits) >= rational(1, 2))
            continue;

        ComplexInterval mapped = exactComplex(0, workBits);
        try {
            mapped = mapLambertBranchPoint(
                p, box, negativeLocalBranch, workBits);
        } catch (const CertifiedBackendUnsupported&) {
            continue;
        }
        if (!intervalContains(box, mapped))
            continue;

        const ComplexInterval factor = encloseLambertBranchFactor(box, workBits);
        const Rational factorLower = complexAbsLower(factor, workBits);
        if (factorLower.isZero())
            continue;
        const ComplexInterval derivative = encloseLambertBranchFactorDerivative(box, workBits);
        const Rational contraction = complexAbsUpper(mapped, workBits)
            * complexAbsUpper(derivative, workBits)
            / (rational(2) * factorLower);
        if (contraction >= rational(1))
            continue;

        // T_Q(B) subset B と一様縮小率を証明した後は，Qが表す各入力に対する
        // 一意な局所根を保ったまま像へ縮められる。
        box = mapped;
        for (std::size_t iteration = 0; iteration < 4; ++iteration) {
            consumeCertifiedWork();
            const ComplexInterval next = mapLambertBranchPoint(
                p, box, negativeLocalBranch, workBits);
            if (!intervalContains(box, next))
                break;
            box = next;
            const auto componentConverged = [&](const RealInterval& component) {
                const Rational width = component.upper().toRational()
                    - component.lower().toRational();
                if (width.isZero())
                    return true;

                // branch pointへ極端に近いと u=W+1 の非零成分は10^-80等まで
                // 小さくなる。W全体の絶対幅だけで止めると，その微小成分の
                // significant digitsが不足し，外側adaptive Nが局所kernelを何度も
                // 呼び直す。符号を証明できる非零成分は自身の大きさに対する
                // relative widthで止め，0を含む成分だけabsolute targetを使う。
                if (component.containsZero())
                    return width <= targetWidth;
                const Rational scale = intervalAbsUpper(component, workBits);
                return width <= scale * targetWidth;
            };
            if (componentConverged(box.real())
                && componentConverged(box.imaginary()))
                break;
        }
        return add(exactComplex(-1, workBits), box, workBits)
            .roundedOutward(precisionBits);
    }

    throw CertifiedBackendUnsupported{
        "Lambert W branch-point root could not be certified by the local contraction backend"};
}

[[nodiscard]] ComplexInterval pointLambertWBranchPointContraction(
    const ComplexInterval& z,
    bool negativeLocalBranch,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 72, "Lambert W branch-point precision is too large");
    const ComplexInterval zWork = z.roundedOutward(workBits);
    const RealInterval e = encloseExp(exactInterval(1, workBits), workBits).interval;
    const ComplexInterval q = add(
        ComplexInterval{
            multiply(zWork.real(), e, workBits),
            multiply(zWork.imaginary(), e, workBits)},
        exactComplex(1, workBits), workBits);
    return pointLambertWBranchPointContractionFromQ(
        q, negativeLocalBranch, precisionBits);
}

[[nodiscard]] ComplexInterval pointLambertWPrincipalSeries(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 48, "complex Lambert W series precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    if (q.isZero())
        return exactComplex(0, precisionBits);
    // |a_{n+1}/a_n|=((n+1)/n)^(n-1)<e<3 なので |z|<1/3 なら
    // W_0(z)=sum (-n)^(n-1) z^n/n! のtailを単純な幾何級数で保証できる。
    if (rational(3) * q >= rational(1))
        throw CertifiedBackendUnsupported{
            "principal Lambert W series requires |z| < 1/3"};

    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "complex Lambert W series target precision is too large"));
    ComplexInterval term = z.roundedOutward(workBits);
    ComplexInterval sum = term;
    constexpr std::uint64_t maximumTerms = 200000;
    for (std::uint64_t n = 1; n < maximumTerms; ++n) {
        consumeCertifiedWork();
        const BigInt nBig = BigInt::fromUnsigned(n);
        const BigInt nextBig = BigInt::fromUnsigned(n + 1);
        const BigInt numerator = numeric::pow(nextBig, n - 1);
        const BigInt denominator = numeric::pow(nBig, n - 1);
        const Rational coefficientRatio = -Rational{numerator, denominator};
        const ComplexInterval next = multiplyByRational(
            multiply(term, z, workBits), coefficientRatio, workBits);
        const Rational tail = complexAbsUpper(next, workBits)
            / (rational(1) - rational(3) * q);
        if (tail <= target)
            return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
        term = next;
        sum = add(sum, term, workBits);
    }
    throw CertifiedBackendUnsupported{
        "principal Lambert W series did not converge within the term limit"};
}

[[nodiscard]] ComplexInterval pointLambertWLogContraction(
    const ComplexInterval& z,
    const BigInt& branch,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 64, "complex Lambert W contraction precision is too large");
    if (z.containsZero())
        throw PrecisionInsufficient{"Lambert W branch input may contain zero"};
    if (crossesNegativeRealCut(z) && !exactZeroPoint(z.imaginary()))
        throw CertifiedBackendUnsupported{
            "complex Lambert W logarithmic branch backend does not cross the negative-real cut"};

    ComplexInterval logarithm = enclosePrincipalComplexLog(
        z.roundedOutward(workBits), workBits).interval;
    const RealInterval twoPi = multiply(
        exactInterval(2, workBits), enclosePi(workBits).interval, workBits);
    const RealInterval branchOffset = multiply(
        twoPi, exactInterval(Rational{branch}, workBits), workBits);
    logarithm = add(logarithm, ComplexInterval{
        exactInterval(0, workBits), branchOffset}, workBits);
    if (logarithm.containsZero())
        throw CertifiedBackendUnsupported{
            "Lambert W logarithmic seed is too close to zero"};

    ComplexInterval candidate = subtract(
        logarithm, enclosePrincipalComplexLog(logarithm, workBits).interval, workBits);
    const ComplexInterval one = exactComplex(1, workBits);
    for (std::size_t iteration = 0; iteration < 10; ++iteration) {
        consumeCertifiedWork();
        const Rational centerReal = (candidate.real().lower().toRational()
            + candidate.real().upper().toRational()) * rational(1, 2);
        const Rational centerImaginary = (candidate.imaginary().lower().toRational()
            + candidate.imaginary().upper().toRational()) * rational(1, 2);
        candidate = exactComplex(centerReal, centerImaginary, workBits);

        const ComplexInterval exponential = encloseComplexExp(candidate, workBits).interval;
        const ComplexInterval numerator = subtract(
            multiply(candidate, exponential, workBits), z.roundedOutward(workBits), workBits);
        const ComplexInterval denominator = multiply(
            exponential, add(candidate, one, workBits), workBits);
        if (denominator.containsZero())
            break;
        candidate = subtract(candidate, divide(numerator, denominator, workBits), workBits);
    }

    const Rational targetWidth = binaryThreshold(checkedAdd(
        precisionBits, 8, "complex Lambert W contraction target precision is too large"));

    // Rectangle反復は各段でcertified Logを再評価するため高precisionでは高価になる。
    // 先に中心diskのBanach証明を試し，残差と縮小率だけで閉じれば反復を省略する。
    const auto tryDiskCertificate = [&]() -> std::optional<ComplexInterval> {
        const Rational centerReal = (candidate.real().lower().toRational()
            + candidate.real().upper().toRational()) * rational(1, 2);
        const Rational centerImaginary = (candidate.imaginary().lower().toRational()
            + candidate.imaginary().upper().toRational()) * rational(1, 2);
        const ComplexInterval center = exactComplex(
            centerReal, centerImaginary, workBits);
        if (center.containsZero() || crossesNegativeRealCut(center))
            return std::nullopt;

        const ComplexInterval centerMapped = subtract(
            logarithm, enclosePrincipalComplexLog(center, workBits).interval, workBits);
        const Rational residual = complexAbsUpper(
            subtract(centerMapped, center, workBits), workBits);
        const Rational centerMagnitude = complexAbsLower(center, workBits);
        for (const Rational& radius : {
                rational(1, 1024), rational(1, 512), rational(1, 256), rational(1, 128),
                rational(1, 64), rational(1, 32), rational(1, 16), rational(1, 8),
                rational(1, 4)}) {
            if (centerMagnitude <= rational(1) + radius)
                continue;
            const ComplexInterval box = inflateComplex(center, radius, workBits);
            if (box.containsZero() || crossesNegativeRealCut(box))
                continue;
            const Rational contraction = rational(1) / (centerMagnitude - radius);
            if (contraction >= rational(1)
                || residual + contraction * radius > radius)
                continue;
            const Rational certifiedRadius = residual / (rational(1) - contraction);
            // ここでは存在・一意性が証明できたdiskを返す。要求幅への到達判定は
            // 外側のCertifiedEvaluatorが行い，不足時だけ高precisionで再評価する。
            return inflateComplex(center, certifiedRadius, workBits)
                .roundedOutward(precisionBits);
        }
        return std::nullopt;
    };
    if (const auto certified = tryDiskCertificate())
        return *certified;

    // T(w)=Log_k(z)-Log(w).  T(D)⊂D かつ inf|w|>1 ならTはD上の縮小写像で，
    // asymptotic branch seedに対応するW_kの根がD内に一意に存在する。
    for (const Rational& radius : {
            rational(1, 64), rational(1, 32), rational(1, 16), rational(1, 8),
            rational(1, 4), rational(1, 2), rational(1)}) {
        ComplexInterval box = inflateComplex(candidate, radius, workBits);
        if (box.containsZero() || crossesNegativeRealCut(box)
            || complexAbsLower(box, workBits) <= rational(1))
            continue;
        ComplexInterval mapped = subtract(
            logarithm, enclosePrincipalComplexLog(box, workBits).interval, workBits);
        if (!intervalContains(box, mapped))
            continue;

        box = mapped;
        for (std::size_t iteration = 0; iteration < 256; ++iteration) {
            consumeCertifiedWork();
            if (box.containsZero() || crossesNegativeRealCut(box))
                break;
            const ComplexInterval next = subtract(
                logarithm, enclosePrincipalComplexLog(box, workBits).interval, workBits);
            if (!intervalContains(box, next))
                break;
            box = next;
            const Rational realWidth = box.real().upper().toRational()
                - box.real().lower().toRational();
            const Rational imaginaryWidth = box.imaginary().upper().toRational()
                - box.imaginary().lower().toRational();
            if (realWidth <= targetWidth && imaginaryWidth <= targetWidth)
                break;
        }
        return box.roundedOutward(precisionBits);
    }

    throw CertifiedBackendUnsupported{
        "complex Lambert W branch could not be certified by the logarithmic contraction backend"};
}

[[nodiscard]] ComplexInterval pointLambertWPrincipalExpContraction(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 64, "principal Lambert W contraction precision is too large");
    const ComplexInterval zWork = z.roundedOutward(workBits);

    // Newton seedは入力区間の中心値だけから作る。中心値は候補探索専用であり，
    // 証明は後段の入力区間全体に対する縮小写像評価だけで行う。
    const auto midpoint = [&](const ComplexInterval& value) {
        const Rational real = (value.real().lower().toRational()
            + value.real().upper().toRational()) * rational(1, 2);
        const Rational imaginary = (value.imaginary().lower().toRational()
            + value.imaginary().upper().toRational()) * rational(1, 2);
        return exactComplex(real, imaginary, workBits);
    };
    const ComplexInterval centerZ = midpoint(zWork);
    ComplexInterval candidate = centerZ;
    const ComplexInterval one = exactComplex(1, workBits);
    for (std::size_t iteration = 0; iteration < 10; ++iteration) {
        consumeCertifiedWork();
        candidate = midpoint(candidate);
        const ComplexInterval fixedPointValue = multiply(
            centerZ, encloseComplexExp(negate(candidate), workBits).interval, workBits);
        const ComplexInterval denominator = add(candidate, one, workBits);
        if (denominator.containsZero())
            break;
        candidate = subtract(candidate, divide(
            subtract(candidate, fixedPointValue, workBits), denominator, workBits), workBits);
    }
    candidate = midpoint(candidate);

    // T(w)=z exp(-w)。中心cの周囲のdisk |w-c|<=R 上では
    // |T'(w)|<=|z| exp(-Re(c)+R)=q。d=|T(c)-c| として d+qR<=R,
    // q<1 を証明できればBanachの縮小写像定理で一意な固定点を得る。
    // さらに |c|+R<1 ならその固定点は単位円板内にあり，w exp(w) がそこで単葉なためW_0である。
    const Rational candidateAbs = complexAbsUpper(candidate, workBits);
    const Rational zAbs = complexAbsUpper(zWork, workBits);
    const ComplexInterval residual = subtract(
        multiply(zWork, encloseComplexExp(negate(candidate), workBits).interval, workBits),
        candidate,
        workBits);
    const Rational residualAbs = complexAbsUpper(residual, workBits);
    const Rational candidateReal = candidate.real().lower().toRational();

    for (const Rational& radius : {
            rational(1, 1024), rational(1, 512), rational(1, 256), rational(1, 128),
            rational(1, 64), rational(1, 32), rational(1, 16), rational(1, 8),
            rational(1, 4), rational(1, 2)}) {
        if (candidateAbs + radius >= rational(1))
            continue;
        const RealInterval exponential = encloseExp(
            exactInterval(-candidateReal + radius, workBits), workBits).interval;
        const Rational contraction = zAbs * exponential.upper().toRational();
        if (contraction >= rational(1))
            continue;
        if (residualAbs + contraction * radius > radius)
            continue;

        const Rational certifiedRadius = residualAbs / (rational(1) - contraction);
        return inflateComplex(candidate, certifiedRadius, workBits)
            .roundedOutward(precisionBits);
    }
    throw CertifiedBackendUnsupported{
        "principal Lambert W could not be certified by the exponential contraction backend"};
}

[[nodiscard]] ComplexInterval pointLambertWComplex(
    const ComplexInterval& z,
    const BigInt& branch,
    std::size_t precisionBits) {
    // -1/e近傍は通常のNewton/log固定点がW'の発散で悪条件になる。
    // W0はprincipal sqrt側，W-1は負実軸の上側（およびcut上の規約値），
    // W1は下側からのみ局所平方根branchへ接続する。
    const BigFloat zero;
    const bool principalLocal = branch.isZero();
    const bool minusOneLocal = branch == BigInt{-1}
        && z.imaginary().lower() >= zero;
    const bool plusOneLocal = branch == BigInt{1}
        && z.imaginary().upper() < zero;
    if (principalLocal || minusOneLocal || plusOneLocal) {
        try {
            return pointLambertWBranchPointContraction(
                z, minusOneLocal || plusOneLocal, precisionBits);
        } catch (const CertifiedBackendUnsupported&) {
            // 局所近傍外または包含が閉じない場合は既存backendへ続ける。
        }
    }

    if (branch.isZero()) {
        const Rational magnitude = complexAbsUpper(z, checkedAdd(
            precisionBits, 24, "Lambert W boundary precision is too large"));
        if (magnitude < rational(1, 4))
            return pointLambertWPrincipalSeries(z, precisionBits);
        try {
            return pointLambertWPrincipalExpContraction(z, precisionBits);
        } catch (const CertifiedBackendUnsupported&) {
            // branch point近傍ではexp contractionが弱くなる。Maclaurinの保証域内なら
            // 級数へfallbackし，それ以外の|W_0|>=1側はlog contractionへ続ける。
            if (magnitude < rational(1, 3))
                return pointLambertWPrincipalSeries(z, precisionBits);
        }
    }
    return pointLambertWLogContraction(z, branch, precisionBits);
}

[[nodiscard]] ComplexInterval pointErfComplex(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "complex erf working precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    const Rational qSquared = q * q;
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "complex erf target precision is too large"));

    ComplexInterval term = z.roundedOutward(workBits);
    ComplexInterval sum = term;
    Rational majorant = q;
    constexpr std::uint64_t maximumTerms = 200000;
    for (std::uint64_t n = 0; n < maximumTerms; ++n) {
        consumeCertifiedWork();
        const Rational ratio = qSquared
            * Rational{BigInt::fromUnsigned(2 * n + 1)}
            / Rational{
                BigInt::fromUnsigned(n + 1)
                * BigInt::fromUnsigned(2 * n + 3)};
        const Rational futureRatio = qSquared
            / Rational{BigInt::fromUnsigned(n + 1)};
        const Rational nextMajorant = majorant * ratio;
        if (futureRatio < rational(1)) {
            const Rational tail = nextMajorant / (rational(1) - futureRatio);
            if (tail <= target) {
                const ComplexInterval enclosed = inflateComplex(sum, tail, workBits);
                const RealInterval scale = divide(
                    exactInterval(2, workBits),
                    encloseSqrt(enclosePi(workBits).interval, workBits).interval,
                    workBits);
                return multiply(
                    enclosed, ComplexInterval::fromReal(scale), workBits)
                    .roundedOutward(precisionBits);
            }
        }

        const Rational signedRatio = -Rational{BigInt::fromUnsigned(2 * n + 1)}
            / Rational{
                BigInt::fromUnsigned(n + 1)
                * BigInt::fromUnsigned(2 * n + 3)};
        term = multiplyByRational(
            multiply(term, multiply(z, z, workBits), workBits),
            signedRatio, workBits);
        sum = add(sum, term, workBits);
        majorant = nextMajorant;
    }
    throw CertifiedBackendUnsupported{"complex erf series did not converge within the term limit"};
}

[[nodiscard]] std::optional<ComplexInterval> pointE1ComplexAsymptotic(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (input.containsZero())
        return std::nullopt;

    // DLMF 6.12.1。右半平面では剰余は最初の未使用項以下である。
    // 左半平面でもnegative real cutを避け，虚部が0から分離できれば
    // csc(|arg z|)=|z|/|Im z| を用いて同じ項からrigorous boundを作れる。
    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "complex E1 asymptotic working precision is too large");
    const ComplexInterval z = input.roundedOutward(workBits);
    const BigFloat zero;

    Rational remainderFactor{BigInt{1}};
    if (z.real().lower() < zero) {
        const Rational imaginaryLower = intervalAbsLower(z.imaginary());
        if (imaginaryLower.isZero())
            return std::nullopt;
        remainderFactor = complexAbsUpper(z, workBits) / imaginaryLower;
        if (remainderFactor < rational(1))
            remainderFactor = rational(1);
    }

    // 漸近級数を使う価値がある大きさだけを対象にする。これは能力境界ではなく
    // dispatchであり，届かなければ呼出し側の収束級数へ戻る。
    if (complexAbsLower(z, workBits) < rational(16))
        return std::nullopt;

    ComplexInterval term = divide(
        encloseComplexExp(negate(z), workBits).interval,
        z, workBits);
    ComplexInterval sum = term;
    Rational previousTerm = complexAbsUpper(term, workBits);
    const Rational relativeTarget = binaryThreshold(checkedAdd(
        precisionBits, 18, "complex E1 asymptotic target precision is too large"));

    constexpr std::uint64_t maximumTerms = 4096;
    for (std::uint64_t k = 0; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const Rational coefficient = -Rational{BigInt::fromUnsigned(k + 1)};
        const ComplexInterval next = divide(
            multiplyByRational(term, coefficient, workBits), z, workBits);
        const Rational nextMagnitude = complexAbsUpper(next, workBits);
        const Rational remainder = nextMagnitude * remainderFactor;
        const Rational sumLower = complexAbsLower(sum, workBits);

        // 真値の大きさは少なくとも |partial|-remainder。
        // remainder <= relativeTarget*|partial|/4 としておけば，
        // whole-complex significant precisionに十分な余裕を持つ。
        if (!sumLower.isZero()
            && remainder * rational(4) <= relativeTarget * sumLower) {
            return inflateComplex(sum, remainder, workBits)
                .roundedOutward(precisionBits);
        }

        // Poincare級数は最適打切りを越えると再び増大する。
        // そこまでに証明できなければ「このprecisionではこのbackendを使わない」。
        if (k != 0 && nextMagnitude >= previousTerm)
            return std::nullopt;

        term = next;
        sum = add(sum, term, workBits);
        previousTerm = nextMagnitude;
    }
    return std::nullopt;
}

[[nodiscard]] Rational intervalAbsUpper(const RealInterval& value) {
    const Rational lower = absRational(value.lower().toRational());
    const Rational upper = absRational(value.upper().toRational());
    return lower > upper ? lower : upper;
}

[[nodiscard]] std::optional<ComplexInterval> pointEiComplexAcrossPositiveRealAxis(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (input.real().lower() <= zero || !input.imaginary().containsZero())
        return std::nullopt;

    // 正実軸はprincipal Eiのbranch cutではない。imaginary InformationEnclosureが
    // 0を跨ぐ場合にE1(-z)へ写すと，補助函数側のcutを人工的に跨いでしまう。
    // そこで実軸上のEiを基準に，縦線分上の
    //   Ei'(z)=exp(z)/z
    // を積分して入力長方形全体を直接包含する。
    const std::size_t workBits = checkedAdd(
        precisionBits, 48, "positive-axis complex Ei working precision is too large");
    const ComplexInterval z = input.roundedOutward(workBits);
    const RealInterval base = encloseExponentialIntegralEiReal(z.real(), workBits);
    const Rational yRadius = intervalAbsUpper(z.imaginary());
    if (yRadius.isZero())
        return ComplexInterval::fromReal(base).roundedOutward(precisionBits);

    const RealInterval exponential = encloseExp(z.real(), workBits).interval;
    const Rational denominatorLower = z.real().lower().toRational();
    if (denominatorLower <= rational(0))
        return std::nullopt;
    const Rational derivativeUpper = exponential.upper().toRational() / denominatorLower;
    const Rational radius = yRadius * derivativeUpper;
    const RealInterval error = symmetricError(radius, workBits);
    return ComplexInterval{
        add(base, error, workBits),
        error}.roundedOutward(precisionBits);
}

[[nodiscard]] std::optional<ComplexInterval> pointEiComplexAsymptotic(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "complex Ei asymptotic working precision is too large");
    const ComplexInterval workZ = z.roundedOutward(workBits);

    if (const auto acrossAxis = pointEiComplexAcrossPositiveRealAxis(workZ, precisionBits))
        return acrossAxis;

    // Ei(z) = -E1(-z) + Log(z) - Log(-z)
    // は現在のprincipal Log規約と整合する。E1側のStokes/cut補正を
    // 手書きの +/- i Pi にせず，二つのprincipal Logの差へ集約する。
    const ComplexInterval minusZ = negate(workZ);
    const auto e1 = pointE1ComplexAsymptotic(minusZ, workBits);
    if (!e1)
        return std::nullopt;

    const ComplexInterval logZ = enclosePrincipalComplexLog(workZ, workBits).interval;
    const ComplexInterval logMinusZ = enclosePrincipalComplexLog(minusZ, workBits).interval;
    ComplexInterval result = subtract(negate(*e1), logMinusZ, workBits);
    result = add(result, logZ, workBits);
    return result.roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval pointChiPositive(
    const Rational& y,
    std::size_t precisionBits) {
    if (y <= rational(0))
        throw std::invalid_argument("Chi positive backend requires y > 0");
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Chi working precision is too large");
    const RealInterval positive = pointExponentialIntegralEi(y, workBits);
    const RealInterval negative = pointExponentialIntegralEi(-y, workBits);
    return multiply(
        add(positive, negative, workBits),
        exactInterval(rational(1, 2), workBits), workBits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] std::optional<ComplexInterval> pointCiComplexPureImaginary(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (!z.real().isPoint() || !z.real().lower().isZero())
        return std::nullopt;
    if (z.imaginary().containsZero())
        return std::nullopt;

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "pure-imaginary Ci working precision is too large");
    const bool positiveImaginary = z.imaginary().lower() > zero;
    const Rational yLower = positiveImaginary
        ? z.imaginary().lower().toRational()
        : -z.imaginary().upper().toRational();
    const Rational yUpper = positiveImaginary
        ? z.imaginary().upper().toRational()
        : -z.imaginary().lower().toRational();
    if (yLower <= rational(0) || yUpper < yLower)
        return std::nullopt;

    // Chi'(y)=sinh(y)/y>0 (y>0) なのでinterval端点だけで全入力を包含できる。
    const RealInterval lowerChi = pointChiPositive(yLower, workBits);
    const RealInterval upperChi = pointChiPositive(yUpper, workBits);
    const RealInterval chi = hull(lowerChi, upperChi).roundedOutward(workBits);

    RealInterval halfPi = multiply(
        enclosePi(workBits).interval,
        exactInterval(rational(1, 2), workBits), workBits);
    if (!positiveImaginary)
        halfPi = negate(halfPi);
    return ComplexInterval{chi, halfPi}.roundedOutward(precisionBits);
}

[[nodiscard]] std::optional<ComplexInterval> pointCiComplexAsymptoticRightHalfPlane(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (z.real().lower() <= zero)
        return std::nullopt;

    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "complex Ci asymptotic working precision is too large");
    const ComplexInterval workZ = z.roundedOutward(workBits);
    const ComplexInterval iz{
        negate(workZ.imaginary()), workZ.real()};
    const ComplexInterval minusIz{
        workZ.imaginary(), negate(workZ.real())};

    // DLMF 6.5.6: Ci(z)=-1/2(E1(i z)+E1(-i z)), |arg z|<Pi/2。
    const auto positive = pointE1ComplexAsymptotic(iz, workBits);
    const auto negative = pointE1ComplexAsymptotic(minusIz, workBits);
    if (!positive || !negative)
        return std::nullopt;
    return multiplyByRational(
        add(*positive, *negative, workBits), rational(-1, 2), workBits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] std::optional<ComplexInterval> pointCiComplexAsymptotic(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const std::size_t boundaryBits = checkedAdd(
        precisionBits, 24, "complex Ci asymptotic boundary precision is too large");
    if (complexAbsLower(z, boundaryBits) < rational(16))
        return std::nullopt;

    if (const auto pureImaginary = pointCiComplexPureImaginary(z, precisionBits))
        return pureImaginary;

    const BigFloat zero;
    if (z.real().lower() > zero)
        return pointCiComplexAsymptoticRightHalfPlane(z, precisionBits);

    if (z.real().upper() < zero) {
        const std::size_t workBits = checkedAdd(
            precisionBits, 32, "complex Ci continuation precision is too large");
        const ComplexInterval workZ = z.roundedOutward(workBits);
        const ComplexInterval minusZ = negate(workZ);
        const auto reflected = pointCiComplexAsymptoticRightHalfPlane(minusZ, workBits);
        if (!reflected)
            return std::nullopt;

        // Cinはevenなので DLMF 6.2.13 から
        // Ci(z)=Ci(-z)+Log(z)-Log(-z)。principal branch補正をLogへ集約する。
        const ComplexInterval logZ = enclosePrincipalComplexLog(workZ, workBits).interval;
        const ComplexInterval logMinusZ = enclosePrincipalComplexLog(minusZ, workBits).interval;
        ComplexInterval result = add(*reflected, logZ, workBits);
        result = subtract(result, logMinusZ, workBits);
        return result.roundedOutward(precisionBits);
    }

    return std::nullopt;
}

[[nodiscard]] ComplexInterval pointEiComplex(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (z.containsZero())
        throw std::domain_error("Ei is singular at zero");
    if (crossesNegativeRealCut(z))
        throw PrecisionInsufficient{"complex Ei input crosses the principal branch cut"};

    const std::size_t boundaryBits = checkedAdd(
        precisionBits, 24, "complex Ei boundary precision is too large");
    const Rational boundaryUpper = complexAbsUpper(z, boundaryBits);
    if (complexAbsLower(z, boundaryBits) >= rational(16)) {
        if (const auto asymptotic = pointEiComplexAsymptotic(z, precisionBits))
            return *asymptotic;
    }
    const auto magnitude = ceilAbsToUint64(boundaryUpper);
    if (!magnitude || *magnitude > std::numeric_limits<std::size_t>::max() / 2)
        throw CertifiedBackendUnsupported{"complex Ei argument is too large for cancellation planning"};
    const std::size_t workBits = checkedAdd(
        checkedAdd(precisionBits, 40, "complex Ei working precision is too large"),
        static_cast<std::size_t>(*magnitude) * 2,
        "complex Ei cancellation precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "complex Ei target precision is too large"));
    ComplexInterval sum = addReal(
        enclosePrincipalComplexLog(z, workBits).interval,
        eulerGammaInterval(workBits), workBits);
    ComplexInterval term = z.roundedOutward(workBits); // k=1
    sum = add(sum, term, workBits);

    // tail証明で必要なのは項絶対値の厳密な値ではなく上界だけである。
    // qを含むRationalを毎項乗算すると分母が指数的に肥大化するため，
    // majorantは固定working precisionのoutward intervalとして更新する。
    RealInterval majorant = exactInterval(q, workBits);

    constexpr std::uint64_t maximumTerms = 200000;
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const BigInt nextIndex = BigInt::fromUnsigned(k + 1);
        const Rational ratio = q * Rational{BigInt::fromUnsigned(k)}
            / Rational{nextIndex * nextIndex};
        const RealInterval nextMajorant = multiply(
            majorant, exactInterval(ratio, workBits), workBits);
        const Rational futureRatio = q / Rational{nextIndex};
        if (futureRatio < rational(1) && (k & 7U) == 0) {
            const RealInterval tailInterval = divide(
                nextMajorant,
                exactInterval(rational(1) - futureRatio, workBits),
                workBits);
            const Rational tail = tailInterval.upper().toRational();
            if (tail <= target)
                return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
        }
        term = multiplyByRational(
            multiply(term, z, workBits),
            Rational{BigInt::fromUnsigned(k), nextIndex * nextIndex},
            workBits);
        sum = add(sum, term, workBits);
        majorant = nextMajorant;
    }
    throw CertifiedBackendUnsupported{"complex Ei series did not converge within the term limit"};
}

[[nodiscard]] ComplexInterval pointSiComplex(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 48, "complex Si working precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    const Rational qSquared = q * q;
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "complex Si target precision is too large"));
    ComplexInterval term = z.roundedOutward(workBits);
    ComplexInterval sum = term;
    Rational majorant = q;

    constexpr std::uint64_t maximumTerms = 200000;
    for (std::uint64_t n = 0; n < maximumTerms; ++n) {
        consumeCertifiedWork();
        const BigInt a = BigInt::fromUnsigned(2 * n + 1);
        const BigInt b = BigInt::fromUnsigned(2 * n + 2);
        const BigInt c = BigInt::fromUnsigned(2 * n + 3);
        const Rational ratio = qSquared * Rational{a}
            / Rational{b * c * c};
        const Rational nextMajorant = majorant * ratio;
        const Rational futureRatio = qSquared / Rational{b * c};
        if (futureRatio < rational(1)) {
            const Rational tail = nextMajorant / (rational(1) - futureRatio);
            if (tail <= target)
                return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
        }
        term = multiplyByRational(
            multiply(term, multiply(z, z, workBits), workBits),
            -Rational{a, b * c * c}, workBits);
        sum = add(sum, term, workBits);
        majorant = nextMajorant;
    }
    throw CertifiedBackendUnsupported{"complex Si series did not converge within the term limit"};
}

[[nodiscard]] ComplexInterval pointCiComplex(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (z.containsZero())
        throw std::domain_error("Ci is singular at zero");
    if (crossesNegativeRealCut(z))
        throw PrecisionInsufficient{"complex Ci input crosses the principal branch cut"};

    const std::size_t boundaryBits = checkedAdd(
        precisionBits, 24, "complex Ci boundary precision is too large");
    const Rational boundaryUpper = complexAbsUpper(z, boundaryBits);
    if (complexAbsLower(z, boundaryBits) >= rational(16)) {
        if (const auto asymptotic = pointCiComplexAsymptotic(z, precisionBits))
            return *asymptotic;
    }
    const auto magnitude = ceilAbsToUint64(boundaryUpper);
    if (!magnitude || *magnitude > std::numeric_limits<std::size_t>::max() / 2)
        throw CertifiedBackendUnsupported{"complex Ci argument is too large for cancellation planning"};
    const std::size_t workBits = checkedAdd(
        checkedAdd(precisionBits, 56, "complex Ci working precision is too large"),
        static_cast<std::size_t>(*magnitude) * 2,
        "complex Ci cancellation precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    const Rational qSquared = q * q;
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "complex Ci target precision is too large"));
    ComplexInterval sum = addReal(
        enclosePrincipalComplexLog(z, workBits).interval,
        eulerGammaInterval(workBits), workBits);
    ComplexInterval zSquared = multiply(z, z, workBits);
    ComplexInterval term = multiplyByRational(zSquared, rational(-1, 4), workBits);
    sum = add(sum, term, workBits);
    RealInterval majorant = exactInterval(qSquared / rational(4), workBits);

    constexpr std::uint64_t maximumTerms = 200000;
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const BigInt a = BigInt::fromUnsigned(2 * k);
        const BigInt b = BigInt::fromUnsigned(2 * k + 1);
        const BigInt c = BigInt::fromUnsigned(2 * k + 2);
        const Rational ratio = qSquared * Rational{a}
            / Rational{b * c * c};
        const RealInterval nextMajorant = multiply(
            majorant, exactInterval(ratio, workBits), workBits);
        const Rational futureRatio = qSquared / Rational{b * c};
        if (futureRatio < rational(1) && (k & 7U) == 0) {
            const RealInterval tailInterval = divide(
                nextMajorant,
                exactInterval(rational(1) - futureRatio, workBits),
                workBits);
            const Rational tail = tailInterval.upper().toRational();
            if (tail <= target)
                return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
        }
        term = multiplyByRational(
            multiply(term, zSquared, workBits),
            -Rational{a, b * c * c}, workBits);
        sum = add(sum, term, workBits);
        majorant = nextMajorant;
    }
    throw CertifiedBackendUnsupported{"complex Ci series did not converge within the term limit"};
}

struct ComplexFresnelPair final {
    ComplexInterval c;
    ComplexInterval s;
};

[[nodiscard]] ComplexInterval rotateComplexQuarterTurns(
    const ComplexInterval& value,
    unsigned turns) {
    switch (turns & 3U) {
    case 0:
        return value;
    case 1:
        return ComplexInterval{negate(value.imaginary()), value.real()};
    case 2:
        return negate(value);
    default:
        return ComplexInterval{value.imaginary(), negate(value.real())};
    }
}

[[nodiscard]] std::optional<std::pair<ComplexInterval, unsigned>>
mapFresnelToAsymptoticWedge(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    // C(i z)=i C(z), S(i z)=-i S(z) を使い，入力を正実軸まわりの
    // |tan(arg z)|<=1/4 のwedgeへ90度単位で回転する。このwedgeは
    // |arg z|<Pi/8に含まれるため，DLMF 7.12.6/7の剰余を最初の未使用項で押さえられる。
    for (unsigned turns = 0; turns < 4; ++turns) {
        const ComplexInterval rotated = rotateComplexQuarterTurns(input, turns);
        const Rational realLower = rotated.real().lower().toRational();
        if (realLower <= rational(0))
            continue;
        const Rational imaginaryUpper = intervalAbsUpper(rotated.imaginary(), precisionBits);
        if (imaginaryUpper * rational(4) <= realLower)
            return std::pair<ComplexInterval, unsigned>{rotated, turns};
    }
    return std::nullopt;
}

[[nodiscard]] ComplexFresnelPair pointFresnelComplexAsymptoticWedge(
    const ComplexInterval& z,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 80, "complex Fresnel asymptotic precision is too large");
    const ComplexInterval zWork = z.roundedOutward(workBits);
    const ComplexInterval one = exactComplex(1, workBits);
    const ComplexInterval half = exactComplex(rational(1, 2), rational(0), workBits);
    const RealInterval pi = enclosePi(workBits).interval;
    const RealInterval piSquared = multiply(pi, pi, workBits);
    const ComplexInterval zSquared = multiply(zWork, zWork, workBits);
    const ComplexInterval zCubed = multiply(zSquared, zWork, workBits);
    const ComplexInterval zFourth = multiply(zSquared, zSquared, workBits);
    const ComplexInterval piZ = multiply(ComplexInterval::fromReal(pi), zWork, workBits);
    const ComplexInterval piSquaredZCubed = multiply(
        ComplexInterval::fromReal(piSquared), zCubed, workBits);
    const ComplexInterval inversePiSquaredZFourth = divide(
        one,
        multiply(ComplexInterval::fromReal(piSquared), zFourth, workBits),
        workBits);

    ComplexInterval fTerm = divide(one, piZ, workBits);
    ComplexInterval gTerm = divide(one, piSquaredZCubed, workBits);
    ComplexInterval f = fTerm;
    ComplexInterval g = gTerm;

    const ComplexInterval phase = multiply(
        ComplexInterval::fromReal(divide(pi, exactInterval(2, workBits), workBits)),
        zSquared,
        workBits);
    const ComplexInterval sine = encloseComplexSinRadian(phase, workBits);
    const ComplexInterval cosine = encloseComplexCosRadian(phase, workBits);
    const Rational sineMagnitude = complexAbsUpper(sine, workBits);
    const Rational cosineMagnitude = complexAbsUpper(cosine, workBits);
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "complex Fresnel asymptotic target precision is too large"));

    Rational previousFTerm = complexAbsUpper(fTerm, workBits);
    Rational previousGTerm = complexAbsUpper(gTerm, workBits);
    constexpr std::size_t maximumTerms = 4096;
    for (std::size_t m = 0; m < maximumTerms; ++m) {
        consumeCertifiedWork();
        const Rational fFactor = rational(-1)
            * unsignedRational(4 * m + 1) * unsignedRational(4 * m + 3);
        const Rational gFactor = rational(-1)
            * unsignedRational(4 * m + 3) * unsignedRational(4 * m + 5);
        const ComplexInterval nextF = multiply(
            multiplyByRational(fTerm, fFactor, workBits),
            inversePiSquaredZFourth,
            workBits);
        const ComplexInterval nextG = multiply(
            multiplyByRational(gTerm, gFactor, workBits),
            inversePiSquaredZFourth,
            workBits);
        const Rational fTail = complexAbsUpper(nextF, workBits);
        const Rational gTail = complexAbsUpper(nextG, workBits);

        // |arg z|<Pi/8ではf,gの剰余は各々最初の未使用項以下。
        // 最終C/Sへの増幅まで含めた上界で要求精度を判定する。
        const Rational cError = fTail * sineMagnitude + gTail * cosineMagnitude;
        const Rational sError = fTail * cosineMagnitude + gTail * sineMagnitude;
        if (cError <= target && sError <= target) {
            const ComplexInterval enclosedF = inflateComplex(f, fTail, workBits);
            const ComplexInterval enclosedG = inflateComplex(g, gTail, workBits);
            return ComplexFresnelPair{
                add(
                    subtract(
                        half,
                        multiply(enclosedG, cosine, workBits),
                        workBits),
                    multiply(enclosedF, sine, workBits),
                    workBits).roundedOutward(precisionBits),
                subtract(
                    subtract(
                        half,
                        multiply(enclosedF, cosine, workBits),
                        workBits),
                    multiply(enclosedG, sine, workBits),
                    workBits).roundedOutward(precisionBits)};
        }

        // 漸近級数が最小項を過ぎたら，これ以上の展開は改善しない。
        // entireなMaclaurin経路へ戻し，保証を捨てて無理に続けない。
        if (fTail >= previousFTerm && gTail >= previousGTerm)
            throw PrecisionInsufficient{
                "complex Fresnel asymptotic expansion did not reach the requested precision"};
        f = add(f, nextF, workBits);
        g = add(g, nextG, workBits);
        fTerm = nextF;
        gTerm = nextG;
        previousFTerm = fTail;
        previousGTerm = gTail;
    }
    throw PrecisionInsufficient{
        "complex Fresnel asymptotic expansion exceeded the term limit"};
}

[[nodiscard]] ComplexInterval pointFresnelComplex(
    const ComplexInterval& z,
    std::size_t precisionBits,
    bool sineVariant) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "complex Fresnel working precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    // 複素Fresnel級数はentireであり，固定の|z|境界は数学的には不要である。
    // 大引数かつ軸近傍では保証付きf/g漸近展開を優先し，それ以外は
    // term上限と共通EvaluationBudgetで制御したMaclaurin級数へ戻す。
    if (q >= rational(8)) {
        if (const auto mapped = mapFresnelToAsymptoticWedge(z, workBits)) {
            try {
                const ComplexFresnelPair pair = pointFresnelComplexAsymptoticWedge(
                    mapped->first, precisionBits);
                // mapped = i^k z。C(i^k z)=i^k C(z)，
                // S(i^k z)=(-i)^k S(z)より元のorientationへ戻す。
                return sineVariant
                    ? rotateComplexQuarterTurns(pair.s, mapped->second)
                    : rotateComplexQuarterTurns(pair.c, (4U - mapped->second) & 3U);
            } catch (const PrecisionInsufficient&) {
                // 要求精度へ届かない小～中引数ではentireな級数へfallbackする。
            }
        }
    }
    const Rational qSquared = q * q;
    const Rational qFourth = qSquared * qSquared;
    const RealInterval pi = enclosePi(workBits).interval;
    const RealInterval piSquaredQuarter = divide(
        multiply(pi, pi, workBits), exactInterval(4, workBits), workBits);
    const Rational piUpper = pi.upper().toRational();
    const Rational baseMajorant = rationalPower(piUpper / rational(2), 2) * qFourth;
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "complex Fresnel target precision is too large"));
    const ComplexInterval zSquared = multiply(z, z, workBits);
    const ComplexInterval zFourth = multiply(zSquared, zSquared, workBits);
    const ComplexInterval common = multiply(
        zFourth, ComplexInterval::fromReal(piSquaredQuarter), workBits);

    ComplexInterval term = z.roundedOutward(workBits);
    // tail majorantをexact Rationalで反復乗算すると，|z|→8で分母が指数的に肥大化する。
    // majorant自体も上界であればよいので，固定precisionの正のdyadic intervalとして伝播する。
    RealInterval majorant = exactInterval(q, workBits);
    if (sineVariant) {
        term = multiply(term, zSquared, workBits);
        term = multiply(term, ComplexInterval::fromReal(pi), workBits);
        term = multiplyByRational(term, rational(1, 6), workBits);
        majorant = exactInterval(piUpper * q * qSquared / rational(6), workBits);
    }
    ComplexInterval sum = term;

    constexpr std::uint64_t maximumTerms = 200000;
    for (std::uint64_t n = 0; n < maximumTerms; ++n) {
        consumeCertifiedWork();
        BigInt numerator;
        BigInt denominator;
        Rational futureRatio;
        if (!sineVariant) {
            numerator = BigInt::fromUnsigned(4 * n + 1);
            denominator = BigInt::fromUnsigned(2 * n + 1)
                * BigInt::fromUnsigned(2 * n + 2)
                * BigInt::fromUnsigned(4 * n + 5);
            futureRatio = baseMajorant
                / Rational{BigInt::fromUnsigned(2 * n + 3)
                    * BigInt::fromUnsigned(2 * n + 4)};
        } else {
            numerator = BigInt::fromUnsigned(4 * n + 3);
            denominator = BigInt::fromUnsigned(2 * n + 2)
                * BigInt::fromUnsigned(2 * n + 3)
                * BigInt::fromUnsigned(4 * n + 7);
            futureRatio = baseMajorant
                / Rational{BigInt::fromUnsigned(2 * n + 4)
                    * BigInt::fromUnsigned(2 * n + 5)};
        }

        ComplexInterval next = multiply(term, common, workBits);
        next = multiplyByRational(next, -Rational{numerator, denominator}, workBits);
        const RealInterval nextMajorant = multiply(
            majorant,
            exactInterval(baseMajorant * Rational{numerator, denominator}, workBits),
            workBits);
        if (futureRatio < rational(1) && (n & 7U) == 0) {
            const RealInterval tailInterval = divide(
                nextMajorant, exactInterval(rational(1) - futureRatio, workBits), workBits);
            const Rational tail = tailInterval.upper().toRational();
            if (tail <= target)
                return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
        }
        term = std::move(next);
        sum = add(sum, term, workBits);
        majorant = nextMajorant;
    }
    throw CertifiedBackendUnsupported{"complex Fresnel series did not converge within the term limit"};
}

[[nodiscard]] ComplexInterval pointHypergeometric1F1Complex(
    const ComplexInterval& a,
    const ComplexInterval& b,
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (exactNonPositiveIntegerPoint(b))
        throw std::domain_error(
            "hypergeometric1F1 has a pole at a non-positive integer b");

    const std::size_t baseBits = checkedAdd(
        precisionBits, 56, "complex hypergeometric1F1 working precision is too large");
    const Rational initialQ = complexAbsUpper(z, baseBits);
    const Rational initialAbsA = complexAbsUpper(a, baseBits);
    const Rational initialAbsB = complexAbsUpper(b, baseBits);
    const auto absZCeil = ceilAbsToUint64(initialQ);
    const auto absACeil = ceilAbsToUint64(initialAbsA);
    const auto absBCeil = ceilAbsToUint64(initialAbsB);

    // real backendと同じく固定の|z|境界は置かない。
    // N>=2|a|,2|b|,3|z|で得る将来項比majorantがterm budget内に
    // 到達できるかだけを事前に判定し，巨大入力の無駄なguard確保を防ぐ。
    constexpr std::uint64_t maximumTerms = 250000;
    if (!absZCeil || !absACeil || !absBCeil
        || *absACeil > maximumTerms / 2 || *absBCeil > maximumTerms / 2
        || *absZCeil > maximumTerms / 3)
        throw CertifiedBackendUnsupported{
            "complex hypergeometric1F1 requires too many certified series terms"};

    const std::size_t workBits = checkedAdd(
        baseBits, static_cast<std::size_t>(*absZCeil) * 2U,
        "complex hypergeometric1F1 cancellation guard is too large");
    const Rational q = complexAbsUpper(z, workBits);
    const Rational absA = complexAbsUpper(a, workBits);
    const Rational absB = complexAbsUpper(b, workBits);
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 22, "complex hypergeometric1F1 target precision is too large"));

    ComplexInterval term = exactComplex(1, workBits);
    ComplexInterval sum = term;
    for (std::uint64_t k = 0; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const ComplexInterval index = exactComplex(static_cast<std::int64_t>(k), workBits);
        const ComplexInterval denominatorFactor = add(b, index, workBits);
        if (denominatorFactor.containsZero()) {
            if (exactZeroPoint(denominatorFactor))
                throw std::domain_error(
                    "hypergeometric1F1 denominator parameter reaches a pole");
            throw PrecisionInsufficient{
                "hypergeometric1F1 denominator parameter interval may contain a pole"};
        }
        ComplexInterval next = multiply(term, add(a, index, workBits), workBits);
        next = multiply(next, z, workBits);
        next = divide(next, denominatorFactor, workBits);
        next = divide(next,
            exactComplex(static_cast<std::int64_t>(k + 1), workBits), workBits);

        const Rational N{BigInt::fromUnsigned(k + 1)};
        if (N >= rational(2) * absA && N >= rational(2) * absB) {
            const Rational Q = rational(3) * q / (N + rational(1));
            if (Q < rational(1)) {
                const Rational tail = complexAbsUpper(next, workBits)
                    / (rational(1) - Q);
                if (tail <= target)
                    return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
            }
        }
        term = std::move(next);
        sum = add(sum, term, workBits);
    }
    throw CertifiedBackendUnsupported{"complex hypergeometric1F1 series did not converge within the term limit"};
}

[[nodiscard]] ComplexInterval positiveIntegerComplexPower(
    std::uint64_t base,
    const ComplexInterval& exponent,
    std::size_t precisionBits) {
    if (base == 1)
        return exactComplex(1, precisionBits);
    if (exactZeroPoint(exponent.imaginary()) && exponent.real().isPoint()) {
        const Rational realExponent = exponent.real().lower().toRational();
        return ComplexInterval::fromReal(
            positiveIntegerBasePower(base, realExponent, precisionBits));
    }

    // real partが整数/半整数なら振幅をexact power/sqrt側で作り，
    // complex Exp内部の実指数評価を省く。位相だけLog+sin/cosでcertifyする。
    if (exponent.real().isPoint() && exponent.imaginary().isPoint()) {
        const Rational realExponent = exponent.real().lower().toRational();
        if (realExponent.denominator() == BigInt{1}
            || realExponent.denominator() == BigInt{2}) {
            const RealInterval logarithm = encloseLogPositive(
                exactInterval(Rational{BigInt::fromUnsigned(base)}, precisionBits),
                precisionBits).interval;
            const RealInterval phase = multiply(
                logarithm, exponent.imaginary(), precisionBits);
            const RealInterval magnitude = positiveIntegerBasePower(
                base, realExponent, precisionBits);
            const RealInterval cosine = encloseCosRadianInterval(
                phase, precisionBits).interval;
            const RealInterval sine = encloseSinRadianInterval(
                phase, precisionBits).interval;
            return ComplexInterval{
                multiply(magnitude, cosine, precisionBits),
                multiply(magnitude, sine, precisionBits)};
        }
    }

    const RealInterval logarithm = encloseLogPositive(
        exactInterval(Rational{BigInt::fromUnsigned(base)}, precisionBits),
        precisionBits).interval;
    return encloseComplexExp(
        multiply(exponent, ComplexInterval::fromReal(logarithm), precisionBits),
        precisionBits).interval;
}

[[nodiscard]] ComplexInterval pointGammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 88, "complex Gamma working precision is too large");
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 24, "complex Gamma target precision is too large"));

    // Stirling剰余を64 Bernoulli項以内で保証できる最小の右半平面へ移す。
    // 固定shiftでは高精度Nが到達不能な剰余上界を何度も再試行するため，
    // requested precisionから理論上必要な下限を先に選ぶ。
    std::optional<Rational> targetReal;
    for (const std::uint64_t candidate : {40ULL, 80ULL, 160ULL, 320ULL,
             640ULL, 1280ULL, 2560ULL, 4096ULL}) {
        const Rational lower{BigInt::fromUnsigned(candidate)};
        bool sufficient = false;
        for (std::size_t n = 1; n <= 64; ++n) {
            consumeCertifiedWork();
            const Rational coefficient = absRational(bernoulliEven(n))
                / Rational{
                    BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * n))
                    * BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * n - 1))};
            BigInt secPower{1};
            secPower <<= n;
            const Rational bound = coefficient * Rational{secPower}
                / rationalPower(lower, 2 * n - 1);
            if (bound <= target) {
                sufficient = true;
                break;
            }
        }
        if (sufficient) {
            targetReal = lower;
            break;
        }
    }
    if (!targetReal)
        throw CertifiedBackendUnsupported{
            "Certified complex Gamma precision exceeds the current Stirling budget"};

    const Rational lowerReal = input.real().lower().toRational();
    std::uint64_t shift = 0;
    if (lowerReal < *targetReal) {
        const Rational needed = *targetReal - lowerReal;
        BigInt quotient = needed.numerator() / needed.denominator();
        if (needed.numerator() % needed.denominator() != BigInt{0})
            quotient += BigInt{1};
        const auto count = numeric::tryToUint64(quotient);
        if (!count || *count > 4096)
            throw CertifiedBackendUnsupported{"complex Gamma requires too large a recurrence shift"};
        shift = *count;
    }

    ComplexInterval shifted = input.roundedOutward(workBits);
    if (shift != 0)
        shifted = add(shifted, exactComplex(static_cast<std::int64_t>(shift), workBits), workBits);
    const Rational positiveLower = shifted.real().lower().toRational();
    if (positiveLower <= rational(0))
        throw CertifiedBackendUnsupported{"complex Gamma could not move the argument into the right half-plane"};

    // Stirling remainderに必要なのはRe(z)の正の下界だけである。
    // finite/exact non-dyadic入力をBigFloatへ外向き変換した下端をそのまま
    // Rational冪乗すると，precisionに比例して巨大な2進分母を(2n-1)乗してしまう。
    // floor(lower)はより粗いが厳密な下界なので，remainder証明には十分である。
    const BigInt lowerInteger = positiveLower.numerator() / positiveLower.denominator();
    if (lowerInteger <= BigInt{0})
        throw CertifiedBackendUnsupported{"complex Gamma Stirling lower bound is not positive"};
    const Rational stirlingLower{lowerInteger};

    const ComplexInterval logarithm = enclosePrincipalComplexLog(shifted, workBits).interval;
    const ComplexInterval shiftedMinusHalf = subtract(
        shifted, exactComplex(rational(1, 2), rational(0), workBits), workBits);
    ComplexInterval logGamma = subtract(
        multiply(shiftedMinusHalf, logarithm, workBits), shifted, workBits);
    const RealInterval twoPi = multiply(
        enclosePi(workBits).interval, exactInterval(2, workBits), workBits);
    const RealInterval halfLogTwoPi = multiply(
        encloseLogPositive(twoPi, workBits).interval,
        exactInterval(rational(1, 2), workBits), workBits);
    logGamma = addReal(logGamma, halfLogTwoPi, workBits);

    std::size_t omittedN = 0;
    Rational remainderBound;
    for (std::size_t n = 1; n <= 64; ++n) {
        consumeCertifiedWork();
        const Rational coefficient = absRational(bernoulliEven(n))
            / Rational{
                BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * n))
                * BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * n - 1))};
        BigInt secPower{1};
        secPower <<= n; // Re(z)>0 => sec(arg(z)/2)^(2n) <= 2^n.
        const Rational denominator = rationalPower(stirlingLower, 2 * n - 1);
        const Rational bound = coefficient * Rational{secPower} / denominator;
        if (bound <= target) {
            omittedN = n;
            remainderBound = bound;
            break;
        }
    }
    if (omittedN == 0)
        throw CertifiedBackendUnsupported{"complex Gamma Stirling remainder did not reach the target precision"};

    const ComplexInterval inverse = divide(exactComplex(1, workBits), shifted, workBits);
    const ComplexInterval inverseSquared = multiply(inverse, inverse, workBits);
    ComplexInterval inversePower = inverse;
    for (std::size_t k = 1; k < omittedN; ++k) {
        consumeCertifiedWork();
        const Rational coefficient = bernoulliEven(k)
            / Rational{
                BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k))
                * BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k - 1))};
        logGamma = add(
            logGamma, multiplyByRational(inversePower, coefficient, workBits), workBits);
        inversePower = multiply(inversePower, inverseSquared, workBits);
    }
    logGamma = inflateComplex(logGamma, remainderBound, workBits);
    ComplexInterval gamma = encloseComplexExp(logGamma, workBits).interval;

    if (shift != 0) {
        ComplexInterval product = exactComplex(1, workBits);
        for (std::uint64_t j = 0; j < shift; ++j) {
            consumeCertifiedWork();
            const ComplexInterval factor = add(
                input.roundedOutward(workBits),
                exactComplex(static_cast<std::int64_t>(j), workBits), workBits);
            if (factor.containsZero())
                throw std::domain_error("gamma is undefined at a non-positive integer");
            product = multiply(product, factor, workBits);
        }
        gamma = divide(gamma, product, workBits);
    }
    return gamma.roundedOutward(precisionBits);
}

struct ComplexPsiPlan final {
    std::uint64_t shift = 0;
    std::size_t omittedK = 0;
    Rational remainderBound;
};

[[nodiscard]] ComplexPsiPlan chooseComplexPsiPlan(
    const ComplexInterval& input,
    std::size_t precisionBits,
    bool trigamma) {
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 24, "complex psi target precision is too large"));
    const Rational inputLower = input.real().lower().toRational();

    // Re(z)>0 へ移した後に評価する。digammaの複素漸近剰余は
    // sec(arg(z)/2)^(2n+1) 倍で抑えられる。Re(z)>0なら
    // sec(arg/2)<sqrt(2)なので，2^(n+1)を安全な有理上界として使う。
    // trigammaはEuler-Maclaurin remainderを実軸方向に積分し，
    // |R_n| <= |B_2n| / Re(z)^(2n+1) を使う。
    for (const std::uint64_t targetReal : {20ULL, 40ULL, 80ULL, 160ULL,
             320ULL, 640ULL, 1280ULL, 2560ULL, 4096ULL}) {
        const Rational targetLower{BigInt::fromUnsigned(targetReal)};
        std::uint64_t shift = 0;
        if (inputLower < targetLower) {
            const Rational needed = targetLower - inputLower;
            BigInt quotient = needed.numerator() / needed.denominator();
            if (needed.numerator() % needed.denominator() != BigInt{0})
                quotient += BigInt{1};
            const auto converted = numeric::tryToUint64(quotient);
            if (!converted || *converted > 4096)
                continue;
            shift = *converted;
        }

        const Rational lower = inputLower
            + Rational{BigInt::fromUnsigned(shift)};
        if (lower <= rational(0))
            continue;

        for (std::size_t k = 1; k <= 64; ++k) {
            consumeCertifiedWork();
            Rational bound = absRational(bernoulliEven(k));
            if (trigamma) {
                bound /= rationalPower(lower, 2 * k + 1);
            } else {
                bound /= Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k))};
                BigInt sectorFactor{1};
                sectorFactor <<= k + 1; // > 2^(k+1/2)
                bound *= Rational{sectorFactor};
                bound /= rationalPower(lower, 2 * k);
            }
            if (bound <= target)
                return ComplexPsiPlan{shift, k, bound};
        }
    }
    throw CertifiedBackendUnsupported{
        trigamma
            ? "Certified complex trigamma precision exceeds the current Euler-Maclaurin budget"
            : "Certified complex digamma precision exceeds the current asymptotic budget"};
}

[[nodiscard]] bool isExactComplexZero(const ComplexInterval& value) {
    return value.real().isPoint() && value.real().lower().isZero()
        && value.imaginary().isPoint() && value.imaginary().lower().isZero();
}

[[nodiscard]] ComplexInterval pointDigammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 72, "complex digamma working precision is too large");
    const ComplexPsiPlan plan = chooseComplexPsiPlan(input, workBits, false);

    ComplexInterval shifted = input.roundedOutward(workBits);
    if (plan.shift != 0)
        shifted = add(
            shifted, exactComplex(static_cast<std::int64_t>(plan.shift), workBits), workBits);
    if (shifted.real().lower().toRational() <= rational(0))
        throw CertifiedBackendUnsupported{
            "complex digamma could not move the argument into the right half-plane"};

    ComplexInterval result = enclosePrincipalComplexLog(shifted, workBits).interval;
    const ComplexInterval inverse = divide(exactComplex(1, workBits), shifted, workBits);
    result = subtract(result, multiplyByRational(inverse, rational(1, 2), workBits), workBits);
    const ComplexInterval inverseSquared = multiply(inverse, inverse, workBits);
    ComplexInterval inversePower = inverseSquared;
    for (std::size_t k = 1; k < plan.omittedK; ++k) {
        consumeCertifiedWork();
        const Rational coefficient = bernoulliEven(k)
            / Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k))};
        result = subtract(
            result, multiplyByRational(inversePower, coefficient, workBits), workBits);
        inversePower = multiply(inversePower, inverseSquared, workBits);
    }
    result = inflateComplex(result, plan.remainderBound, workBits);

    for (std::uint64_t j = 0; j < plan.shift; ++j) {
        consumeCertifiedWork();
        const ComplexInterval divisor = add(
            input.roundedOutward(workBits),
            exactComplex(static_cast<std::int64_t>(j), workBits), workBits);
        if (isExactComplexZero(divisor))
            throw std::domain_error("digamma is undefined at a non-positive integer");
        if (divisor.containsZero())
            throw PrecisionInsufficient{"complex digamma recurrence interval may contain a pole"};
        result = subtract(
            result, divide(exactComplex(1, workBits), divisor, workBits), workBits);
    }
    return result.roundedOutward(precisionBits);
}

[[nodiscard]] ComplexInterval pointTrigammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 72, "complex trigamma working precision is too large");
    const ComplexPsiPlan plan = chooseComplexPsiPlan(input, workBits, true);

    ComplexInterval shifted = input.roundedOutward(workBits);
    if (plan.shift != 0)
        shifted = add(
            shifted, exactComplex(static_cast<std::int64_t>(plan.shift), workBits), workBits);
    if (shifted.real().lower().toRational() <= rational(0))
        throw CertifiedBackendUnsupported{
            "complex trigamma could not move the argument into the right half-plane"};

    const ComplexInterval inverse = divide(exactComplex(1, workBits), shifted, workBits);
    const ComplexInterval inverseSquared = multiply(inverse, inverse, workBits);
    ComplexInterval result = add(
        inverse, multiplyByRational(inverseSquared, rational(1, 2), workBits), workBits);
    ComplexInterval inversePower = multiply(inverseSquared, inverse, workBits);
    for (std::size_t k = 1; k < plan.omittedK; ++k) {
        consumeCertifiedWork();
        result = add(
            result, multiplyByRational(inversePower, bernoulliEven(k), workBits), workBits);
        inversePower = multiply(inversePower, inverseSquared, workBits);
    }
    result = inflateComplex(result, plan.remainderBound, workBits);

    for (std::uint64_t j = 0; j < plan.shift; ++j) {
        consumeCertifiedWork();
        const ComplexInterval divisor = add(
            input.roundedOutward(workBits),
            exactComplex(static_cast<std::int64_t>(j), workBits), workBits);
        if (isExactComplexZero(divisor))
            throw std::domain_error("trigamma is undefined at a non-positive integer");
        if (divisor.containsZero())
            throw PrecisionInsufficient{"complex trigamma recurrence interval may contain a pole"};
        const ComplexInterval inverseDivisor = divide(
            exactComplex(1, workBits), divisor, workBits);
        result = add(
            result, multiply(inverseDivisor, inverseDivisor, workBits), workBits);
    }
    return result.roundedOutward(precisionBits);
}

[[nodiscard]] bool excludesRealRayFrom(
    const ComplexInterval& value,
    const Rational& start) {
    const Rational imaginaryLower = value.imaginary().lower().toRational();
    const Rational imaginaryUpper = value.imaginary().upper().toRational();
    return imaginaryLower > rational(0)
        || imaginaryUpper < rational(0)
        || value.real().upper().toRational() < start;
}

[[nodiscard]] ComplexInterval complexPiSquaredOverSix(std::size_t precisionBits) {
    const RealInterval pi = enclosePi(precisionBits).interval;
    return ComplexInterval::fromReal(multiply(
        multiply(pi, pi, precisionBits),
        exactInterval(rational(1, 6), precisionBits), precisionBits));
}

[[nodiscard]] ComplexInterval pointPolylogComplex(
    std::uint64_t order,
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (order == 0)
        throw std::invalid_argument("polylog certified complex series requires positive order");
    const std::size_t workBits = checkedAdd(
        precisionBits, 48, "complex polylog working precision is too large");

    // finite-precision real inputでも整数orderのmu=log(z)近傍展開を区間のまま使う。
    // midpointへ戻さず，入力区間幅をそのまま結果へ伝播する。
    if (order >= 3 && z.isProvablyReal()) {
        if (const auto nearOne = pointPolylogPositiveIntegerNearOne(
                order, z.real(), precisionBits))
            return ComplexInterval::fromReal(*nearOne);
    }

    // DLMF 25.12.3。負実軸の大きめの引数は w=z/(z-1) へ写すと
    // unit disk内のより小さい正実数になる。区間全体が負実軸上の場合だけ使う。
    if (order == 2 && z.isProvablyReal()
        && z.real().upper().toRational() < rational(0)) {
        const ComplexInterval denominator = subtract(
            z.roundedOutward(workBits), exactComplex(1, workBits), workBits);
        const ComplexInterval transformed = divide(
            z.roundedOutward(workBits), denominator, workBits);
        const Rational sourceMagnitude = complexAbsUpper(z, workBits);
        const Rational transformedMagnitude = complexAbsUpper(transformed, workBits);
        if (transformedMagnitude < rational(1)
            && transformedMagnitude * rational(4) < sourceMagnitude * rational(3)) {
            const RealInterval oneMinusZ = subtract(
                exactInterval(1, workBits), z.real().roundedOutward(workBits), workBits);
            const RealInterval logarithm = encloseLogPositive(oneMinusZ, workBits).interval;
            RealInterval result = multiply(logarithm, logarithm, workBits);
            result = multiply(result, exactInterval(rational(-1, 2), workBits), workBits);
            const ComplexInterval reflected = pointPolylogComplex(2, transformed, workBits);
            result = subtract(result, reflected.real(), workBits);
            return ComplexInterval::fromReal(result.roundedOutward(precisionBits));
        }
    }

    // DLMF 25.12.6を実軸上のcertified intervalへそのまま適用する。
    // finite-precision入力でも区間全体が(1/2,1)に入る場合だけ使用するため，
    // hidden midpointやbranch-side情報を発明しない。
    if (order == 2 && z.isProvablyReal()
        && z.real().lower().toRational() > rational(1, 2)
        && z.real().upper().toRational() < rational(1)) {
        const RealInterval oneMinus = subtract(
            exactInterval(1, workBits), z.real().roundedOutward(workBits), workBits);
        const ComplexInterval reflected = pointPolylogComplex(
            2, ComplexInterval::fromReal(oneMinus), workBits);
        const RealInterval pi = enclosePi(workBits).interval;
        RealInterval result = multiply(pi, pi, workBits);
        result = multiply(result, exactInterval(rational(1, 6), workBits), workBits);
        const RealInterval logZ = encloseLogPositive(z.real(), workBits).interval;
        const RealInterval logComplement = encloseLogPositive(oneMinus, workBits).interval;
        result = subtract(result, multiply(logZ, logComplement, workBits), workBits);
        result = subtract(result, reflected.real(), workBits);
        return ComplexInterval::fromReal(result.roundedOutward(precisionBits));
    }

    const ComplexInterval workZ = z.roundedOutward(workBits);
    const Rational q = complexAbsUpper(workZ, workBits);

    if (order == 2 && !workZ.containsZero()) {
        const ComplexInterval zMinusOne = subtract(
            workZ, exactComplex(1, workBits), workBits);

        // DLMF 25.12.3。z/(z-1) が直接unit disk内のより小さい引数へ
        // 移る場合は一段のconnection formulaを優先する。これにより例えばz=Iも
        // unit-circle級数を直接扱わずcertifyできる。
        if (excludesRealRayFrom(workZ, rational(1)) && !zMinusOne.containsZero()) {
            const ComplexInterval transformed = divide(workZ, zMinusOne, workBits);
            const Rational transformedMagnitude = complexAbsUpper(transformed, workBits);
            if (transformedMagnitude < rational(1)
                && transformedMagnitude * rational(4) < q * rational(3)) {
                const ComplexInterval oneMinusZ = subtract(
                    exactComplex(1, workBits), workZ, workBits);
                const ComplexInterval logarithm = enclosePrincipalComplexLog(
                    oneMinusZ, workBits).interval;
                ComplexInterval result = multiplyByRational(
                    multiply(logarithm, logarithm, workBits),
                    rational(-1, 2), workBits);
                result = subtract(
                    result, pointPolylogComplex(2, transformed, workBits), workBits);
                return result.roundedOutward(precisionBits);
            }
        }

        // DLMF 25.12.3 + 25.12.4。
        // u=(z-1)/z が元のzより十分小さく，両connection formulaのcutを
        // 区間全体が避けることを証明できる場合だけ二段変換する。
        if (excludesRealRayFrom(workZ, rational(1))) {
            const ComplexInterval transformed = divide(zMinusOne, workZ, workBits); // 1/w
            const Rational transformedMagnitude = complexAbsUpper(transformed, workBits);
            if (transformedMagnitude < rational(1)
                && transformedMagnitude * rational(4) < q * rational(3)) {
                const ComplexInterval w = divide(workZ, zMinusOne, workBits);
                if (excludesRealRayFrom(w, rational(0))) {
                    const ComplexInterval oneMinusZ = subtract(
                        exactComplex(1, workBits), workZ, workBits);
                    const ComplexInterval logOneMinusZ = enclosePrincipalComplexLog(
                        oneMinusZ, workBits).interval;
                    const ComplexInterval logMinusW = enclosePrincipalComplexLog(
                        negate(w), workBits).interval;
                    ComplexInterval result = complexPiSquaredOverSix(workBits);
                    result = subtract(result, multiplyByRational(
                        multiply(logOneMinusZ, logOneMinusZ, workBits),
                        rational(1, 2), workBits), workBits);
                    result = add(result, multiplyByRational(
                        multiply(logMinusW, logMinusW, workBits),
                        rational(1, 2), workBits), workBits);
                    result = add(
                        result, pointPolylogComplex(2, transformed, workBits), workBits);
                    return result.roundedOutward(precisionBits);
                }
            }
        }

        // DLMF 25.12.4。unit disk外で正実軸cutを避ける場合は1/zへ反転する。
        if (complexAbsLower(workZ, workBits) > rational(1)
            && excludesRealRayFrom(workZ, rational(0))) {
            const ComplexInterval inverse = divide(
                exactComplex(1, workBits), workZ, workBits);
            if (complexAbsUpper(inverse, workBits) < rational(1)) {
                const ComplexInterval logMinusZ = enclosePrincipalComplexLog(
                    negate(workZ), workBits).interval;
                ComplexInterval result = negate(pointPolylogComplex(2, inverse, workBits));
                result = subtract(result, complexPiSquaredOverSix(workBits), workBits);
                result = subtract(result, multiplyByRational(
                    multiply(logMinusZ, logMinusZ, workBits),
                    rational(1, 2), workBits), workBits);
                return result.roundedOutward(precisionBits);
            }
        }
    }

    if (q >= rational(1)) {
        if (complexAbsLower(z, workBits) >= rational(1))
            throw CertifiedBackendUnsupported{
                "polylog certified power series requires |z| < 1"};
        throw PrecisionInsufficient{
            "complex polylog magnitude interval straddles the unit-disk series boundary"};
    }
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "complex polylog target precision is too large"));
    ComplexInterval term = z.roundedOutward(workBits); // k=1
    ComplexInterval sum = term;
    constexpr std::uint64_t maximumTerms = 1'000'000;
    for (std::uint64_t k = 1; k + 1 < maximumTerms; ++k) {
        consumeCertifiedWork();
        const BigInt kPower = numeric::pow(BigInt::fromUnsigned(k), order);
        const BigInt nextPower = numeric::pow(BigInt::fromUnsigned(k + 1), order);
        ComplexInterval next = multiply(term, z, workBits);
        next = multiplyByRational(next, Rational{kPower, nextPower}, workBits);

        // |term_{j+1}/term_j|<=|z|より，次項から先を幾何級数でmajorizeする。
        const Rational nextBound = complexAbsUpper(next, workBits);
        const Rational tail = nextBound / (rational(1) - q);
        sum = add(sum, next, workBits);
        if (tail <= target)
            return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
        term = std::move(next);
    }
    throw CertifiedBackendUnsupported{"complex polylog series did not converge within the term limit"};
}

[[nodiscard]] ComplexInterval pointHypergeometric2F1Complex(
    const ComplexInterval& a,
    const ComplexInterval& b,
    const ComplexInterval& c,
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (exactNonPositiveIntegerPoint(c))
        throw std::domain_error(
            "hypergeometric2F1 has a pole at a non-positive integer c");

    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "complex hypergeometric2F1 working precision is too large");
    const Rational q = complexAbsUpper(z, workBits);

    if (q < rational(1)) {
        const Rational absA = complexAbsUpper(a, workBits);
        const Rational absB = complexAbsUpper(b, workBits);
        const Rational absC = complexAbsUpper(c, workBits);
        const Rational target = binaryThreshold(checkedAdd(
            precisionBits, 22, "complex hypergeometric2F1 target precision is too large"));

        ComplexInterval term = exactComplex(1, workBits);
        ComplexInterval sum = term;
        constexpr std::uint64_t maximumTerms = 250000;
        for (std::uint64_t k = 0; k < maximumTerms; ++k) {
            consumeCertifiedWork();
            const ComplexInterval index = exactComplex(static_cast<std::int64_t>(k), workBits);
            const ComplexInterval denominatorFactor = add(c, index, workBits);
            if (denominatorFactor.containsZero()) {
                if (exactZeroPoint(denominatorFactor))
                    throw std::domain_error(
                        "hypergeometric2F1 denominator parameter reaches a pole");
                throw PrecisionInsufficient{
                    "hypergeometric2F1 denominator parameter interval may contain a pole"};
            }
            ComplexInterval next = multiply(term, add(a, index, workBits), workBits);
            next = multiply(next, add(b, index, workBits), workBits);
            next = multiply(next, z, workBits);
            next = divide(next, denominatorFactor, workBits);
            next = divide(next,
                exactComplex(static_cast<std::int64_t>(k + 1), workBits), workBits);

            const std::uint64_t n = k + 1;
            const Rational N{BigInt::fromUnsigned(n)};
            if (N > absC) {
                const Rational Q = q
                    * (rational(1) + absA / N)
                    * (rational(1) + absB / N)
                    / (rational(1) - absC / N);
                if (Q < rational(1)) {
                    const Rational tail = complexAbsUpper(next, workBits)
                        / (rational(1) - Q);
                    if (tail <= target)
                        return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
                }
            }
            term = std::move(next);
            sum = add(sum, term, workBits);
        }
        throw CertifiedBackendUnsupported{"complex hypergeometric2F1 series did not converge within the term limit"};
    }

    const Rational qLower = complexAbsLower(z, workBits);
    if (qLower < rational(1))
        throw PrecisionInsufficient{
            "complex hypergeometric2F1 magnitude interval straddles the |z| = 1 convergence boundary"};
    if (qLower == rational(1)) {
        if (q > rational(1))
            throw PrecisionInsufficient{
                "complex hypergeometric2F1 magnitude interval straddles the |z| = 1 continuation boundary"};

        // z=1 かつ Re(c-a-b)>0 ならGauss summationでunit-circle境界を直接閉じる。
        // その他のunit-circle点はparameter依存の境界収束を別途扱う必要がある。
        const bool exactZOne = exactZeroPoint(z.imaginary()) && z.real().isPoint()
            && z.real().lower().toRational() == rational(1);
        if (exactZOne) {
            const ComplexInterval cMinusA = subtract(c, a, workBits);
            const ComplexInterval cMinusB = subtract(c, b, workBits);
            const ComplexInterval delta = subtract(cMinusA, b, workBits);
            if (delta.real().lower().toRational() > rational(0)) {
                const auto gammaValue = [&](const ComplexInterval& argument) {
                    if (argument.isProvablyReal())
                        return ComplexInterval::fromReal(
                            encloseGammaReal(argument.real(), workBits));
                    return pointGammaComplex(argument, workBits);
                };
                try {
                    const ComplexInterval numerator = multiply(
                        gammaValue(c), gammaValue(delta), workBits);
                    const ComplexInterval denominator = multiply(
                        gammaValue(cMinusA), gammaValue(cMinusB), workBits);
                    return divide(numerator, denominator, workBits)
                        .roundedOutward(precisionBits);
                } catch (const std::domain_error&) {
                    throw CertifiedBackendUnsupported{
                        "hypergeometric2F1 Gauss summation is degenerate for these parameters"};
                }
            }
        }
        throw CertifiedBackendUnsupported{
            "complex hypergeometric2F1 evaluation on |z| = 1 requires a supported boundary formula"};
    }

    // |z|>1 ではprincipal-valueの1/z connection formulaへ送る。
    // a-bが整数の退化公式はdigammaを含む別展開になるため，Gamma poleを検出したら
    // このbackendでは扱わず未評価へ戻す。
    if (z.containsZero())
        throw CertifiedBackendUnsupported{"hypergeometric2F1 continuation requires z != 0"};
    const ComplexInterval reciprocalZ = divide(exactComplex(1, workBits), z, workBits);
    if (complexAbsUpper(reciprocalZ, workBits) >= rational(1))
        throw CertifiedBackendUnsupported{
            "Certified complex hypergeometric2F1 currently supports |z|<1 or a provable |z|>1 connection"};

    const ComplexInterval minusZ = negate(z);
    // z>1 のexact real pointは principal continuation の規約値を持つ。
    // -z が負実軸上のpointである場合は principal Log の +Pi 側をそのまま使う。
    // 幅を持つ区間がcutを跨ぐ場合だけ保守的に未対応へ戻す。
    if (crossesNegativeRealCut(minusZ) && !exactZeroPoint(minusZ.imaginary()))
        throw CertifiedBackendUnsupported{
            "hypergeometric2F1 1/z connection crosses the principal (-z) branch cut"};
    const ComplexInterval logMinusZ = enclosePrincipalComplexLog(minusZ, workBits).interval;
    const auto principalPower = [&](const ComplexInterval& exponent) {
        return encloseComplexExp(multiply(exponent, logMinusZ, workBits), workBits).interval;
    };

    try {
        // continuation係数にはGammaが7個現れるが，a/bやその差がprovably realなら
        // complex Stirling backendへ送る必要はない。実Gammaのreflection/positive
        // backendを使い，結果だけComplexIntervalへ持ち上げる。
        const auto gammaValue = [&](const ComplexInterval& argument) {
            if (argument.isProvablyReal())
                return ComplexInterval::fromReal(
                    encloseGammaReal(argument.real(), workBits));
            return pointGammaComplex(argument, workBits);
        };

        const ComplexInterval gammaC = gammaValue(c);
        const ComplexInterval gammaBMinusA = gammaValue(subtract(b, a, workBits));
        const ComplexInterval gammaAMinusB = gammaValue(subtract(a, b, workBits));
        const ComplexInterval gammaB = gammaValue(b);
        const ComplexInterval gammaA = gammaValue(a);
        const ComplexInterval gammaCMinusA = gammaValue(subtract(c, a, workBits));
        const ComplexInterval gammaCMinusB = gammaValue(subtract(c, b, workBits));

        const ComplexInterval one = exactComplex(1, workBits);
        const ComplexInterval firstHyper = pointHypergeometric2F1Complex(
            a,
            add(subtract(one, c, workBits), a, workBits),
            add(subtract(one, b, workBits), a, workBits),
            reciprocalZ, workBits);
        const ComplexInterval secondHyper = pointHypergeometric2F1Complex(
            b,
            add(subtract(one, c, workBits), b, workBits),
            add(subtract(one, a, workBits), b, workBits),
            reciprocalZ, workBits);

        ComplexInterval firstCoefficient = multiply(gammaC, gammaBMinusA, workBits);
        firstCoefficient = divide(firstCoefficient,
            multiply(gammaB, gammaCMinusA, workBits), workBits);
        ComplexInterval secondCoefficient = multiply(gammaC, gammaAMinusB, workBits);
        secondCoefficient = divide(secondCoefficient,
            multiply(gammaA, gammaCMinusB, workBits), workBits);

        const ComplexInterval firstPower = principalPower(negate(a));
        const ComplexInterval secondPower = principalPower(negate(b));
        const ComplexInterval first = multiply(
            multiply(firstCoefficient, firstPower, workBits), firstHyper, workBits);
        const ComplexInterval second = multiply(
            multiply(secondCoefficient, secondPower, workBits), secondHyper, workBits);
        return add(first, second, workBits).roundedOutward(precisionBits);
    } catch (const std::domain_error&) {
        throw CertifiedBackendUnsupported{
            "hypergeometric2F1 1/z connection is degenerate for these parameters"};
    }
}

[[nodiscard]] ComplexInterval pointZetaEulerMaclaurinComplex(
    const ComplexInterval& s,
    std::size_t precisionBits) {
    const Rational sigma = s.real().lower().toRational();
    const std::size_t workBits = checkedAdd(
        precisionBits, 72, "complex zeta working precision is too large");
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 24, "complex zeta target precision is too large"));
    const Rational absS = complexAbsUpper(s, workBits);

    // Bernoulli remainderのplannerはexact Rationalで(2Pi)^(-2k)を育てない。
    // 外向きinterval recurrenceなら高precisionでも証明幅だけを保持できる。
    const RealInterval oneReal = exactInterval(1, workBits);
    const RealInterval twoPi = multiply(
        exactInterval(2, workBits), enclosePi(workBits).interval, workBits);
    const RealInterval inverseTwoPiSquared = divide(
        oneReal, squareInterval(twoPi, workBits), workBits);

    std::uint64_t chosenN = 0;
    std::size_t chosenK = 0;
    Rational remainderBound;
    for (std::uint64_t N : {8ULL, 16ULL, 32ULL, 64ULL, 128ULL}) {
        const BigInt nBig = BigInt::fromUnsigned(N);
        const RealInterval inverseNSquared = exactInterval(
            Rational{BigInt{1}, nBig * nBig}, workBits);
        RealInterval risingBound = exactInterval(1, workBits);
        RealInterval nPower = positiveIntegerBasePower(N, rational(1) - sigma, workBits);
        RealInterval inverseFourierPower = exactInterval(1, workBits);

        for (std::size_t k = 1; k <= 48; ++k) {
            consumeCertifiedWork();
            const std::uint64_t first = static_cast<std::uint64_t>(2 * k - 2);
            risingBound = multiply(risingBound,
                exactInterval(absS + Rational{BigInt::fromUnsigned(first)}, workBits),
                workBits);
            risingBound = multiply(risingBound,
                exactInterval(absS + Rational{BigInt::fromUnsigned(first + 1)}, workBits),
                workBits);
            nPower = multiply(nPower, inverseNSquared, workBits);
            inverseFourierPower = multiply(
                inverseFourierPower, inverseTwoPiSquared, workBits);
            if (k < 2)
                continue;

            const Rational remainderExponent = sigma + Rational{BigInt::fromUnsigned(
                static_cast<std::uint64_t>(2 * k - 1))};
            // Euler-Maclaurin remainder integral requires Re(s)+2k-1>0.
            // Critical stripではkを増やせば満たせるため，Re(s)>1という人工境界は不要。
            if (remainderExponent <= rational(0))
                continue;

            RealInterval bound = multiply(risingBound, nPower, workBits);
            bound = multiply(bound, inverseFourierPower, workBits);
            bound = multiply(exactInterval(4, workBits), bound, workBits);
            bound = divide(bound, exactInterval(remainderExponent, workBits), workBits);
            const Rational upperBound = bound.upper().toRational();
            if (upperBound <= target) {
                chosenN = N;
                chosenK = k;
                remainderBound = upperBound;
                break;
            }
        }
        if (chosenN != 0)
            break;
    }
    if (chosenN == 0)
        throw CertifiedBackendUnsupported{"complex zeta Euler-Maclaurin budget is insufficient"};

    ComplexInterval result = exactComplex(0, workBits);
    const ComplexInterval negativeS = negate(s);
    std::vector<std::optional<ComplexInterval>> dirichletPowers(chosenN);
    dirichletPowers[1] = exactComplex(1, workBits);
    for (std::uint64_t n = 1; n < chosenN; ++n) {
        consumeCertifiedWork();
        ComplexInterval term = exactComplex(1, workBits);
        if (n > 1) {
            const std::uint64_t factor = smallestFactor(n);
            if (factor == n) {
                term = positiveIntegerComplexPower(n, negativeS, workBits);
            } else {
                term = multiply(
                    *dirichletPowers[factor],
                    *dirichletPowers[n / factor],
                    workBits);
            }
            dirichletPowers[n] = term;
        }
        result = add(result, term, workBits);
    }

    const ComplexInterval oneMinusS = subtract(exactComplex(1, workBits), s, workBits);
    const ComplexInterval integralNumerator = positiveIntegerComplexPower(
        chosenN, oneMinusS, workBits);
    const ComplexInterval sMinusOne = subtract(s, exactComplex(1, workBits), workBits);
    result = add(result, divide(integralNumerator, sMinusOne, workBits), workBits);
    result = add(result, multiplyByRational(
        positiveIntegerComplexPower(chosenN, negativeS, workBits), rational(1, 2), workBits), workBits);

    const BigInt chosenNBig = BigInt::fromUnsigned(chosenN);
    const ComplexInterval inverseNSquared = ComplexInterval::fromReal(exactInterval(
        Rational{BigInt{1}, chosenNBig * chosenNBig}, workBits));
    ComplexInterval rising = s;
    ComplexInterval power = positiveIntegerComplexPower(
        chosenN, negate(add(s, exactComplex(1, workBits), workBits)), workBits);
    BigInt factorial = numeric::factorial(2);
    for (std::size_t k = 1; k < chosenK; ++k) {
        consumeCertifiedWork();
        const Rational coefficient = bernoulliEven(k) / Rational{factorial};
        ComplexInterval correction = multiply(rising, power, workBits);
        correction = multiplyByRational(correction, coefficient, workBits);
        result = add(result, correction, workBits);

        const std::uint64_t first = static_cast<std::uint64_t>(2 * k - 1);
        rising = multiply(rising,
            add(s, exactComplex(static_cast<std::int64_t>(first), workBits), workBits),
            workBits);
        rising = multiply(rising,
            add(s, exactComplex(static_cast<std::int64_t>(first + 1), workBits), workBits),
            workBits);
        power = multiply(power, inverseNSquared, workBits);
        factorial *= BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k + 1));
        factorial *= BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k + 2));
    }
    return inflateComplex(result, remainderBound, workBits).roundedOutward(precisionBits);
}

[[nodiscard]] ComplexInterval pointZetaComplex(
    const ComplexInterval& s,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 80, "zeta continuation working precision is too large");
    const ComplexInterval one = exactComplex(1, workBits);
    const ComplexInterval sWork = s.roundedOutward(workBits);
    const ComplexInterval sMinusOne = subtract(sWork, one, workBits);
    if (sMinusOne.containsZero()) {
        if (sWork.real().isPoint() && sWork.imaginary().isPoint()
            && sWork.real().lower().toRational() == rational(1)
            && sWork.imaginary().lower().isZero())
            throw std::domain_error("zeta has a pole at s = 1");
        throw PrecisionInsufficient{
            "zeta input interval may contain the pole at s = 1",
            PrecisionInsufficientKind::InputInformation};
    }

    const Rational sigmaLower = sWork.real().lower().toRational();
    const Rational sigmaUpper = sWork.real().upper().toRational();

    // Re(s)>=0ではEuler-Maclaurinを直接使う。DLMF 25.11.7の適用条件
    // Re(s)>-2nをremainder planner内で検査するため，critical stripも同じbackendで閉じる。
    if (sigmaLower >= rational(0))
        return pointZetaEulerMaclaurinComplex(sWork, precisionBits);

    // 区間全体が左半平面ならfunctional equationでRe(1-s)>1へ写す。
    // ζ(s)=2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s).
    if (sigmaUpper < rational(0)) {
        const ComplexInterval reflected = subtract(one, sWork, workBits);
        const ComplexInterval reflectedZeta = pointZetaEulerMaclaurinComplex(
            reflected, workBits);
        const ComplexInterval twoPower = positiveIntegerComplexPower(2, sWork, workBits);
        const RealInterval pi = enclosePi(workBits).interval;
        const RealInterval logPi = encloseLogPositive(pi, workBits).interval;
        const ComplexInterval piPower = encloseComplexExp(
            multiply(sMinusOne, ComplexInterval::fromReal(logPi), workBits),
            workBits).interval;
        const ComplexInterval halfPiS = multiply(
            sWork, ComplexInterval::fromReal(multiply(
                pi, exactInterval(rational(1, 2), workBits), workBits)), workBits);
        const ComplexInterval sine = encloseComplexSinRadian(halfPiS, workBits);
        const ComplexInterval gamma = pointGammaComplex(reflected, workBits);
        ComplexInterval result = multiply(twoPower, piPower, workBits);
        result = multiply(result, sine, workBits);
        result = multiply(result, gamma, workBits);
        result = multiply(result, reflectedZeta, workBits);
        return result.roundedOutward(precisionBits);
    }

    // Re(s)=0を跨ぐ狭い入力はfunctional equationと直接式を混在させず，
    // Euler-Maclaurin側で一括包含する。64 Bernoulli項の範囲で証明できない場合だけbounded-work扱い。
    return pointZetaEulerMaclaurinComplex(sWork, precisionBits);
}

} // namespace

/*
旧実装

変更理由：
- xがexact pointでもpublic interval wrapperからlower/upperを別々に呼び，同一計算を2回実行していた。
- lower/upper endpointで共通なBeta(a,b) normalizationまで各回で再計算していた。
- 代表例ibeta[1/3,2/3,1/4]では2F1よりBeta/Gamma normalizationが支配的だったため，
  endpoint計算の共有だけでexactnessを変えず大きく短縮できる。

元コード：

[[nodiscard]] RealInterval pointIncompleteBetaRegularized(
    const Rational& a,
    const Rational& b,
    const Rational& x,
    std::size_t precisionBits) {
    if (a <= rational(0) || b <= rational(0))
        throw std::domain_error("ibeta certified backend requires a > 0 and b > 0");
    if (x < rational(0) || x > rational(1))
        throw std::domain_error("ibeta certified backend requires x in [0,1]");
    if (x.isZero())
        return exactInterval(0, precisionBits);
    if (x == rational(1))
        return exactInterval(1, precisionBits);

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "ibeta working precision is too large");
    if (x > rational(1, 2)) {
        const RealInterval complement = pointIncompleteBetaRegularized(
            b, a, rational(1) - x, workBits);
        return subtract(exactInterval(1, workBits), complement, workBits)
            .roundedOutward(precisionBits);
    }

    const RealInterval logX = encloseLogPositive(exactInterval(x, workBits), workBits).interval;
    const RealInterval xPower = encloseExp(
        multiply(logX, exactInterval(a, workBits), workBits), workBits).interval;
    const RealInterval hyper = encloseHypergeometric2F1Real(
        a, rational(1) - b, a + rational(1), x, workBits);
    const RealInterval numerator = divide(
        multiply(xPower, hyper, workBits), exactInterval(a, workBits), workBits);
    const RealInterval beta = encloseBetaPositive(
        exactInterval(a, workBits), exactInterval(b, workBits), workBits);
    return divide(numerator, beta, workBits).roundedOutward(precisionBits);
}
*/

[[nodiscard]] RealInterval pointIncompleteBetaRegularized(
    const Rational& a,
    const Rational& b,
    const Rational& x,
    const RealInterval& betaNormalization,
    std::size_t precisionBits) {
    if (a <= rational(0) || b <= rational(0))
        throw std::domain_error("ibeta certified backend requires a > 0 and b > 0");
    if (x < rational(0) || x > rational(1))
        throw std::domain_error("ibeta certified backend requires x in [0,1]");
    if (x.isZero())
        return exactInterval(0, precisionBits);
    if (x == rational(1))
        return exactInterval(1, precisionBits);

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "ibeta working precision is too large");
    if (x > rational(1, 2)) {
        // B(a,b)=B(b,a)なのでcomplement側でも同じnormalizationを再利用できる。
        const RealInterval complement = pointIncompleteBetaRegularized(
            b, a, rational(1) - x, betaNormalization, workBits);
        return subtract(exactInterval(1, workBits), complement, workBits)
            .roundedOutward(precisionBits);
    }

    const RealInterval logX = encloseLogPositive(exactInterval(x, workBits), workBits).interval;
    const RealInterval xPower = encloseExp(
        multiply(logX, exactInterval(a, workBits), workBits), workBits).interval;
    const RealInterval hyper = encloseHypergeometric2F1Real(
        a, rational(1) - b, a + rational(1), x, workBits);
    const RealInterval numerator = divide(
        multiply(xPower, hyper, workBits), exactInterval(a, workBits), workBits);
    return divide(
        numerator, betaNormalization.roundedOutward(workBits), workBits)
        .roundedOutward(precisionBits);
}


RealInterval encloseGammaRational(
    const Rational& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Gamma precision must be at least one bit");

    if (input > rational(0)) {
        const RealInterval interval = exactInterval(input, checkedAdd(
            precisionBits, 64, "Gamma rational input precision is too large"));
        return encloseExp(
            encloseLogGammaPositive(interval, precisionBits, &input),
            precisionBits).interval;
    }

    return gammaNegativeRationalByReflection(input, precisionBits);
}

RealInterval encloseLogGammaRational(
    const Rational& input,
    std::size_t precisionBits) {
    if (input > rational(0)) {
        const RealInterval interval = exactInterval(input, checkedAdd(
            precisionBits, 64, "LogGamma rational input precision is too large"));
        return encloseLogGammaPositive(interval, precisionBits, &input);
    }

    const RealInterval gamma = encloseGammaRational(input, checkedAdd(
        precisionBits, 24, "lgamma rational working precision is too large"));
    const RealInterval magnitude = absoluteInterval(gamma, checkedAdd(
        precisionBits, 16, "lgamma rational absolute-value precision is too large"));
    if (magnitude.containsZero())
        throw PrecisionInsufficient{"lgamma could not yet prove Gamma away from zero"};
    return encloseLogPositive(magnitude, precisionBits).interval;
}

RealInterval encloseLambertWReal(
    const RealInterval& input,
    int branch,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Lambert W precision must be at least one bit");
    if (branch != 0 && branch != -1)
        throw std::domain_error("Certified real Lambert W supports only branches 0 and -1");

    const std::size_t domainBits = checkedAdd(
        precisionBits, 32, "Lambert W domain precision is too large");
    const RealInterval branchPoint = lambertBranchPoint(domainBits);
    const Rational inputLower = input.lower().toRational();
    const Rational inputUpper = input.upper().toRational();
    const Rational branchPointLower = branchPoint.lower().toRational();
    const Rational branchPointUpper = branchPoint.upper().toRational();

    if (inputUpper < branchPointLower)
        throw CertifiedBackendUnsupported{
            "Certified complex Lambert W is not implemented for z < -1/e"};
    if (inputLower < branchPointUpper)
        throw PrecisionInsufficient{"Lambert W input is too close to the branch point -1/e"};

    if (branch == -1) {
        if (inputLower > rational(0))
            throw CertifiedBackendUnsupported{
                "Certified complex Lambert W branch -1 is not implemented for z > 0"};
        if (inputLower == rational(0) && inputUpper == rational(0))
            throw std::domain_error("Lambert W branch -1 is singular at z = 0");
        if (inputUpper >= rational(0))
            throw PrecisionInsufficient{"Lambert W branch -1 input is too close to zero"};
    }

    const std::size_t workBits = checkedAdd(
        precisionBits, 16, "Lambert W working precision is too large");
    if (branch == 0) {
        const RealInterval lowerValue = pointLambertWReal(inputLower, 0, workBits);
        const RealInterval upperValue = pointLambertWReal(inputUpper, 0, workBits);
        return RealInterval{lowerValue.lower(), upperValue.upper()}.roundedOutward(precisionBits);
    }

    // W_-1は実区間上で単調減少。入力上端が出力下端に対応する。
    const RealInterval lowerValue = pointLambertWReal(inputUpper, -1, workBits);
    const RealInterval upperValue = pointLambertWReal(inputLower, -1, workBits);
    return RealInterval{lowerValue.lower(), upperValue.upper()}.roundedOutward(precisionBits);
}

ComplexInterval encloseLambertWBranchPointOffset(
    const ComplexInterval& offset,
    const BigInt& branch,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Lambert W precision must be at least one bit");

    const BigFloat zero;
    const bool principalLocal = branch.isZero();
    const bool minusOneLocal = branch == BigInt{-1}
        && offset.imaginary().lower() >= zero;
    const bool plusOneLocal = branch == BigInt{1}
        && offset.imaginary().upper() <= zero
        && offset.imaginary().lower() < zero;
    if (!principalLocal && !minusOneLocal && !plusOneLocal)
        throw CertifiedBackendUnsupported{
            "Lambert W branch is not connected to the requested branch-point side"};

    const std::size_t workBits = checkedAdd(
        precisionBits, 72, "Lambert W branch-point precision is too large");
    const RealInterval e = encloseExp(exactInterval(1, workBits), workBits).interval;

    // 実係数函数の共役対称性を使い，branch point下側のprincipal branchとW_1は
    // 上側の既証明kernelへ写してから共役で戻す。極小な負虚部を高precisionへ
    // 丸めてから反転すると，負値側のoutward roundingだけに不要なcostが出るため，
    // 共役はprecision拡張より先に行う。W_1(z)=conj(W_-1(conj(z))) が局所対応。
    const bool reflectFromLowerHalfPlane = offset.imaginary().upper() <= BigFloat{}
        && offset.imaginary().lower() < BigFloat{}
        && (principalLocal || plusOneLocal);
    const ComplexInterval canonicalOffset = reflectFromLowerHalfPlane
        ? ComplexInterval{offset.real(), negate(offset.imaginary())}
        : offset;
    const ComplexInterval offsetWork = canonicalOffset.roundedOutward(workBits);

    const ComplexInterval q{
        multiply(offsetWork.real(), e, workBits),
        multiply(offsetWork.imaginary(), e, workBits)};
    ComplexInterval result = pointLambertWBranchPointContractionFromQ(
        q, minusOneLocal || plusOneLocal, precisionBits);
    if (reflectFromLowerHalfPlane)
        result = ComplexInterval{result.real(), negate(result.imaginary())};
    return result;
}

ComplexInterval encloseLambertWComplex(
    const ComplexInterval& input,
    const BigInt& branch,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Lambert W precision must be at least one bit");
    if (!branch.isZero() && input.containsZero()) {
        if (input.real().isPoint() && input.imaginary().isPoint()
            && input.real().lower().isZero() && input.imaginary().lower().isZero())
            throw std::domain_error("non-principal Lambert W branches are singular at z = 0");
        throw PrecisionInsufficient{
            "Lambert W input interval may contain the logarithmic branch point z = 0",
            PrecisionInsufficientKind::InputInformation};
    }
    return pointLambertWComplex(input, branch, precisionBits);
}

RealInterval encloseGammaReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Gamma precision must be at least one bit");
    if (exactNonPositiveIntegerPoint(input))
        throw std::domain_error("gamma is undefined at non-positive integers");

    if (input.isPoint()) {
        const Rational point = input.lower().toRational();
        if (point == rational(1) || point == rational(2))
            return exactInterval(1, precisionBits);
    }

    const BigFloat zero;
    if (input.lower() > zero)
        return encloseExp(encloseLogGammaPositive(input, precisionBits), precisionBits).interval;
    if (input.upper() < zero)
        return gammaNegativeByReflection(input, precisionBits);

    throw PrecisionInsufficient{"Gamma interval straddles zero or a possible pole"};
}

RealInterval encloseLogGammaReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (exactNonPositiveIntegerPoint(input))
        throw std::domain_error("lgamma is undefined at non-positive integers");

    const BigFloat zero;
    if (input.lower() > zero)
        return encloseLogGammaPositive(input, precisionBits);

    const RealInterval gamma = encloseGammaReal(input, checkedAdd(
        precisionBits, 24, "lgamma working precision is too large"));
    const RealInterval magnitude = absoluteInterval(gamma, checkedAdd(
        precisionBits, 16, "lgamma absolute-value precision is too large"));
    if (magnitude.containsZero())
        throw PrecisionInsufficient{"lgamma could not yet prove Gamma away from zero"};
    return encloseLogPositive(magnitude, precisionBits).interval;
}

RealInterval encloseErfReal(
    const RealInterval& input,
    std::size_t precisionBits) { // erfは実軸上で単調増加。
    const RealInterval lower = pointErf(input.lower().toRational(), precisionBits);
    const RealInterval upper = pointErf(input.upper().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

RealInterval encloseErfcReal(
    const RealInterval& input,
    std::size_t precisionBits) { // erfcは実軸上で単調減少。
    const RealInterval lower = pointErfc(input.upper().toRational(), precisionBits);
    const RealInterval upper = pointErfc(input.lower().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

ComplexInterval encloseErfComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("erf precision must be at least one bit");
    return pointErfComplex(input, precisionBits);
}

ComplexInterval encloseErfcComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 24, "complex erfc working precision is too large");
    return subtract(
        exactComplex(1, workBits), pointErfComplex(input, workBits), workBits)
        .roundedOutward(precisionBits);
}

ComplexInterval encloseGammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Gamma precision must be at least one bit");
    return pointGammaComplex(input, precisionBits);
}

[[nodiscard]] RealInterval encloseBetaLogPositiveRational(
    const Rational& a,
    const Rational& b,
    std::size_t precisionBits) {
    if (a <= rational(0) || b <= rational(0))
        throw std::domain_error("betaln certified evaluation currently requires a > 0 and b > 0");
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Beta working precision is too large");
    const Rational sum = a + b;
    const RealInterval aInterval = exactInterval(a, workBits);
    const RealInterval bInterval = exactInterval(b, workBits);
    const RealInterval sumInterval = exactInterval(sum, workBits);
    return subtract(
        add(
            encloseLogGammaPositive(aInterval, workBits, &a),
            encloseLogGammaPositive(bInterval, workBits, &b),
            workBits),
        encloseLogGammaPositive(sumInterval, workBits, &sum),
        workBits)
        .roundedOutward(precisionBits);
}

[[nodiscard]] RealInterval encloseBetaPositiveRational(
    const Rational& a,
    const Rational& b,
    std::size_t precisionBits) {
    return encloseExp(
        encloseBetaLogPositiveRational(a, b, precisionBits),
        precisionBits).interval;
}

RealInterval encloseBetaLogPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (a.lower() <= zero || b.lower() <= zero)
        throw std::domain_error("betaln certified evaluation currently requires a > 0 and b > 0");
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "Beta working precision is too large");
    const RealInterval sum = add(a, b, workBits);
    return subtract(
        add(encloseLogGammaPositive(a, workBits), encloseLogGammaPositive(b, workBits), workBits),
        encloseLogGammaPositive(sum, workBits), workBits)
        .roundedOutward(precisionBits);
}


RealInterval encloseFresnelCReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Fresnel precision must be at least one bit");
    return encloseFresnelRealImpl(input, true, precisionBits);
}

RealInterval encloseFresnelSReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Fresnel precision must be at least one bit");
    return encloseFresnelRealImpl(input, false, precisionBits);
}


ComplexInterval encloseFresnelCComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Fresnel C precision must be at least one bit");
    return pointFresnelComplex(input, precisionBits, false);
}

ComplexInterval encloseFresnelSComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Fresnel S precision must be at least one bit");
    return pointFresnelComplex(input, precisionBits, true);
}

RealInterval encloseHypergeometric1F1Real(
    const Rational& a,
    const Rational& b,
    const Rational& z,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("hypergeometric1F1 precision must be at least one bit");
    return pointHypergeometric1F1(a, b, z, precisionBits);
}

ComplexInterval encloseHypergeometric1F1Complex(
    const ComplexInterval& a,
    const ComplexInterval& b,
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("hypergeometric1F1 precision must be at least one bit");
    return pointHypergeometric1F1Complex(a, b, z, precisionBits);
}

RealInterval encloseHypergeometric2F1Real(
    const Rational& a,
    const Rational& b,
    const Rational& c,
    const Rational& z,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("hypergeometric2F1 precision must be at least one bit");
    return pointHypergeometric2F1(a, b, c, z, precisionBits);
}

ComplexInterval encloseHypergeometric2F1Complex(
    const ComplexInterval& a,
    const ComplexInterval& b,
    const ComplexInterval& c,
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("hypergeometric2F1 precision must be at least one bit");
    return pointHypergeometric2F1Complex(a, b, c, z, precisionBits);
}

RealInterval encloseEllipticFReal(
    const Rational& phi, const Rational& m, std::size_t precisionBits) {
    try {
        return pointEllipticSeries(EllipticSeriesKind::F, rational(0), phi, m, precisionBits);
    }
    catch (const CertifiedBackendUnsupported&) {
        return reducedEllipticCarlson(
            EllipticSeriesKind::F, exactInterval(0, precisionBits),
            exactInterval(phi, precisionBits), exactInterval(m, precisionBits), precisionBits);
    }
}

RealInterval encloseEllipticEReal(
    const Rational& phi, const Rational& m, std::size_t precisionBits) {
    try {
        return pointEllipticSeries(EllipticSeriesKind::E, rational(0), phi, m, precisionBits);
    }
    catch (const CertifiedBackendUnsupported&) {
        return reducedEllipticCarlson(
            EllipticSeriesKind::E, exactInterval(0, precisionBits),
            exactInterval(phi, precisionBits), exactInterval(m, precisionBits), precisionBits);
    }
}

RealInterval encloseEllipticPiReal(
    const Rational& n, const Rational& phi, const Rational& m, std::size_t precisionBits) {
    try {
        return pointEllipticSeries(EllipticSeriesKind::Pi, n, phi, m, precisionBits);
    }
    catch (const CertifiedBackendUnsupported&) {
        return reducedEllipticCarlson(
            EllipticSeriesKind::Pi, exactInterval(n, precisionBits),
            exactInterval(phi, precisionBits), exactInterval(m, precisionBits), precisionBits);
    }
}

RealInterval encloseEllipticFReal(
    const RealInterval& phi, const RealInterval& m, std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticF precision must be at least one bit");
    try {
        return intervalEllipticSeries(
            EllipticSeriesKind::F, exactInterval(0, precisionBits), phi, m, precisionBits);
    }
    catch (const CertifiedBackendUnsupported&) {
        return reducedEllipticCarlson(
            EllipticSeriesKind::F, exactInterval(0, precisionBits), phi, m, precisionBits);
    }
}

RealInterval encloseEllipticEReal(
    const RealInterval& phi, const RealInterval& m, std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticE precision must be at least one bit");
    try {
        return intervalEllipticSeries(
            EllipticSeriesKind::E, exactInterval(0, precisionBits), phi, m, precisionBits);
    }
    catch (const CertifiedBackendUnsupported&) {
        return reducedEllipticCarlson(
            EllipticSeriesKind::E, exactInterval(0, precisionBits), phi, m, precisionBits);
    }
}

RealInterval encloseEllipticPiReal(
    const RealInterval& n, const RealInterval& phi, const RealInterval& m,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticPi precision must be at least one bit");
    try {
        return intervalEllipticSeries(
            EllipticSeriesKind::Pi, n, phi, m, precisionBits);
    }
    catch (const CertifiedBackendUnsupported&) {
        return reducedEllipticCarlson(
            EllipticSeriesKind::Pi, n, phi, m, precisionBits);
    }
}

RealInterval encloseEllipticFRealPiMultiple(
    const Rational& piCoefficient,
    const RealInterval& m,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticF precision must be at least one bit");
    const std::size_t bits = checkedAdd(
        precisionBits, 72, "elliptic Carlson working precision is too large");
    return reducedEllipticCarlsonWithReduction(
        EllipticSeriesKind::F, exactInterval(0, bits),
        exactPiMultipleReduction(piCoefficient, bits), m, precisionBits);
}

RealInterval encloseEllipticERealPiMultiple(
    const Rational& piCoefficient,
    const RealInterval& m,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticE precision must be at least one bit");
    const std::size_t bits = checkedAdd(
        precisionBits, 72, "elliptic Carlson working precision is too large");
    return reducedEllipticCarlsonWithReduction(
        EllipticSeriesKind::E, exactInterval(0, bits),
        exactPiMultipleReduction(piCoefficient, bits), m, precisionBits);
}

RealInterval encloseEllipticPiRealPiMultiple(
    const RealInterval& n,
    const Rational& piCoefficient,
    const RealInterval& m,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticPi precision must be at least one bit");
    const std::size_t bits = checkedAdd(
        precisionBits, 72, "elliptic Carlson working precision is too large");
    return reducedEllipticCarlsonWithReduction(
        EllipticSeriesKind::Pi, n,
        exactPiMultipleReduction(piCoefficient, bits), m, precisionBits);
}


RealInterval encloseExponentialIntegralEiReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Ei precision must be at least one bit");
    if (input.containsZero())
        throw std::domain_error("Ei interval contains the singular point zero");

    const Rational lower = input.lower().toRational();
    const Rational upper = input.upper().toRational();
    const RealInterval lo = pointExponentialIntegralEi(lower, precisionBits);
    const RealInterval hi = pointExponentialIntegralEi(upper, precisionBits);
    if (upper < rational(0))
        return RealInterval{hi.lower(), lo.upper()}; // x<0ではEi' = exp(x)/x < 0。
    return RealInterval{lo.lower(), hi.upper()};
}

RealInterval encloseSineIntegralSiReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 24, "Si interval precision is too large");
    const Rational lower = input.lower().toRational();
    const Rational upper = input.upper().toRational();
    RealInterval value = pointSineIntegralSi(lower, workBits);
    const Rational width = upper - lower;
    // |sin(x)/x|<=1（x=0では極限1）なので入力幅だけ膨らませれば十分。
    if (!width.isZero())
        value = add(value, symmetricError(width, workBits), workBits);
    return value.roundedOutward(precisionBits);
}

RealInterval encloseCosineIntegralCiPositive(
    const RealInterval& input,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (input.lower() <= zero)
        throw std::domain_error("Ci real certified backend requires a positive interval");
    const std::size_t workBits = checkedAdd(
        precisionBits, 24, "Ci interval precision is too large");
    const Rational lower = input.lower().toRational();
    const Rational upper = input.upper().toRational();
    RealInterval value = pointCosineIntegralCiPositive(lower, workBits);
    const Rational width = upper - lower;
    // |Ci'(x)|=|cos(x)/x|<=1/lower on a positive interval.
    if (!width.isZero())
        value = add(value,
            symmetricError(width / lower, workBits), workBits);
    return value.roundedOutward(precisionBits);
}

RealInterval encloseLogarithmicIntegralLiPositive(
    const RealInterval& input,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (input.lower() <= zero)
        throw std::domain_error("li real certified backend requires x > 0");
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "li working precision is too large");
    const RealInterval logarithm = encloseLogPositive(input.roundedOutward(workBits), workBits).interval;
    if (logarithm.containsZero())
        throw std::domain_error("li interval contains the singular point x=1");
    return encloseExponentialIntegralEiReal(logarithm, workBits).roundedOutward(precisionBits);
}

ComplexInterval encloseExponentialIntegralEiComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Ei precision must be at least one bit");
    return pointEiComplex(input, precisionBits);
}

ComplexInterval encloseSineIntegralSiComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Si precision must be at least one bit");
    return pointSiComplex(input, precisionBits);
}

ComplexInterval encloseCosineIntegralCiComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Ci precision must be at least one bit");
    return pointCiComplex(input, precisionBits);
}

ComplexInterval encloseLogarithmicIntegralLiComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("li precision must be at least one bit");
    const std::size_t workBits = checkedAdd(
        precisionBits, 32, "li complex working precision is too large");
    const ComplexInterval logarithm = enclosePrincipalComplexLog(
        input.roundedOutward(workBits), workBits).interval;
    // li(z)=Ei(Log(z)) はこのEiが最終演算である。Logの入力区間はガード桁で
    // 作るが，Eiの出力まで同じ余分な32bitを要求するとEi自身のガードと二重になる。
    return pointEiComplex(logarithm, precisionBits);
}

RealInterval enclosePolylogReal(
    std::uint64_t order,
    const Rational& z,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("polylog precision must be at least one bit");
    return pointPolylogPositiveOrder(order, z, precisionBits);
}

ComplexInterval enclosePolylogComplex(
    std::uint64_t order,
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("polylog precision must be at least one bit");
    return pointPolylogComplex(order, z, precisionBits);
}

RealInterval encloseBetaRational(
    const Rational& a,
    const Rational& b,
    std::size_t precisionBits) {
    return encloseBetaPositiveRational(a, b, precisionBits);
}

RealInterval encloseBetaLogRational(
    const Rational& a,
    const Rational& b,
    std::size_t precisionBits) {
    return encloseBetaLogPositiveRational(a, b, precisionBits);
}

RealInterval encloseBetaPositive(
    const RealInterval& a,
    const RealInterval& b,
    std::size_t precisionBits) {
    return encloseExp(encloseBetaLogPositive(a, b, precisionBits), precisionBits).interval;
}

RealInterval encloseZetaReal(
    const RealInterval& input,
    std::size_t precisionBits) {
    const Rational lower = input.lower().toRational();
    const Rational upper = input.upper().toRational();
    if (input.isPoint() && lower == rational(1))
        throw std::domain_error("zeta has a pole at s = 1");
    if (lower <= rational(1) && upper >= rational(1))
        throw PrecisionInsufficient{
            "zeta interval may contain the pole at s = 1",
            PrecisionInsufficientKind::InputInformation};

    // s>1では単調減少を使う旧real fast pathを維持する。その他は複素Euler-Maclaurin/
    // functional-equation backendへ送り，実軸上の結果の実部だけを返す。
    if (lower > rational(1)) {
        const RealInterval lowValue = pointZetaGreaterThanOne(upper, precisionBits);
        const RealInterval highValue = pointZetaGreaterThanOne(lower, precisionBits);
        return RealInterval{lowValue.lower(), highValue.upper()};
    }
    return pointZetaComplex(ComplexInterval::fromReal(input), precisionBits).real();
}

ComplexInterval encloseZetaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("zeta precision must be at least one bit");
    return pointZetaComplex(input, precisionBits);
}

ComplexInterval encloseDigammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("digamma precision must be at least one bit");
    return pointDigammaComplex(input, precisionBits);
}

ComplexInterval encloseTrigammaComplex(
    const ComplexInterval& input,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("trigamma precision must be at least one bit");
    return pointTrigammaComplex(input, precisionBits);
}

RealInterval encloseDigammaPositive(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (input.lower().toRational() <= rational(0))
        throw std::domain_error("digamma certified real backend currently requires x > 0");
    const RealInterval lower = pointDigammaPositive(input.lower().toRational(), precisionBits);
    const RealInterval upper = pointDigammaPositive(input.upper().toRational(), precisionBits);
    return RealInterval{lower.lower(), upper.upper()};
}

RealInterval encloseTrigammaPositive(
    const RealInterval& input,
    std::size_t precisionBits) {
    if (input.lower().toRational() <= rational(0))
        throw std::domain_error("trigamma certified real backend currently requires x > 0");
    const RealInterval lowerValue = pointTrigammaPositive(input.upper().toRational(), precisionBits);
    const RealInterval upperValue = pointTrigammaPositive(input.lower().toRational(), precisionBits);
    return RealInterval{lowerValue.lower(), upperValue.upper()};
}

RealInterval encloseIncompleteBetaRegularized(
    const RealInterval& a,
    const RealInterval& b,
    const RealInterval& x,
    std::size_t precisionBits) {
    const BigFloat zero;
    if (a.upper() <= zero || b.upper() <= zero)
        throw std::domain_error("ibeta certified backend requires positive a and b");
    if (a.lower() <= zero || b.lower() <= zero)
        throw PrecisionInsufficient{
            "ibeta parameter InformationEnclosure crosses the positive-real domain boundary",
            PrecisionInsufficientKind::InputInformation};

    const Rational xLower = x.lower().toRational();
    const Rational xUpper = x.upper().toRational();
    if (xUpper < rational(0) || xLower > rational(1))
        throw std::domain_error("ibeta certified backend requires x in [0,1]");
    if (xLower < rational(0) || xUpper > rational(1))
        throw PrecisionInsufficient{
            "ibeta x InformationEnclosure crosses the x in [0,1] domain boundary",
            PrecisionInsufficientKind::InputInformation};

    // a,b>0ではBeta分布のstochastic orderingから I_x(a,b) は
    // aに関して減少，bとxに関して増加する。したがってparameter interval全体は
    // [I_xL(aU,bL), I_xU(aL,bU)] で厳密に囲える。各endpointは既存exact-Rational
    // backendへ渡すため，有限precision入力からhidden center値を取り戻すこともない。
    const Rational aLower = a.lower().toRational();
    const Rational aUpper = a.upper().toRational();
    const Rational bLower = b.lower().toRational();
    const Rational bUpper = b.upper().toRational();
    const auto evaluatePoint = [&](const Rational& aa, const Rational& bb,
                                   const Rational& xx) {
        if (xx.isZero())
            return exactInterval(0, precisionBits);
        if (xx == rational(1))
            return exactInterval(1, precisionBits);
        const std::size_t normalizationBits = checkedAdd(
            precisionBits, xx > rational(1, 2) ? 80 : 40,
            "ibeta interval-parameter normalization precision is too large");
        const RealInterval betaNormalization = encloseBetaPositiveRational(
            aa, bb, normalizationBits);
        return pointIncompleteBetaRegularized(
            aa, bb, xx, betaNormalization, precisionBits);
    };

    const RealInterval lowValue = evaluatePoint(aUpper, bLower, xLower);
    const RealInterval highValue = evaluatePoint(aLower, bUpper, xUpper);
    return RealInterval{lowValue.lower(), highValue.upper()};
}

RealInterval encloseIncompleteBetaRegularized(
    const Rational& a,
    const Rational& b,
    const RealInterval& x,
    std::size_t precisionBits) {
    if (a <= rational(0) || b <= rational(0))
        throw std::domain_error("ibeta certified backend requires positive a and b");
    const Rational lower = x.lower().toRational();
    const Rational upper = x.upper().toRational();
    if (lower < rational(0) || upper > rational(1))
        throw std::domain_error("ibeta certified backend requires x in [0,1]");
    /*
    旧実装
    exact pointでも同一endpointを2回評価し，Beta(a,b)も2回再構築していた。

    const RealInterval lowValue = pointIncompleteBetaRegularized(a, b, lower, precisionBits);
    const RealInterval highValue = pointIncompleteBetaRegularized(a, b, upper, precisionBits);
    return RealInterval{lowValue.lower(), highValue.upper()};
    */

    if (lower.isZero() && upper.isZero())
        return exactInterval(0, precisionBits);
    if (lower == rational(1) && upper == rational(1))
        return exactInterval(1, precisionBits);

    if (lower == upper) {
        // x<=1/2のpointは再帰しないので+40 bitで十分。complementを使うpointだけ+80 bitを確保する。
        const std::size_t normalizationGuard = lower > rational(1, 2) ? 80 : 40;
        const std::size_t normalizationBits = checkedAdd(
            precisionBits, normalizationGuard, "ibeta normalization precision is too large");
        const RealInterval betaNormalization = encloseBetaPositiveRational(
            a, b, normalizationBits);
        return pointIncompleteBetaRegularized(
            a, b, lower, betaNormalization, precisionBits);
    }

    // interval endpointのどちらかがcomplementへ入っても共有できるよう+80 bitで一度だけnormalizationを作る。
    const std::size_t normalizationBits = checkedAdd(
        precisionBits, 80, "ibeta normalization precision is too large");
    const RealInterval betaNormalization = encloseBetaPositiveRational(
        a, b, normalizationBits);
    const RealInterval lowValue = pointIncompleteBetaRegularized(
        a, b, lower, betaNormalization, precisionBits);
    const RealInterval highValue = pointIncompleteBetaRegularized(
        a, b, upper, betaNormalization, precisionBits);
    return RealInterval{lowValue.lower(), highValue.upper()};
}

} // namespace mmcal::approximation
