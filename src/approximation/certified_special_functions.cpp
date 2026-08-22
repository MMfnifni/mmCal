// 特殊函数の保証付き評価
#include "certified_special_functions.hpp"

#include "certification_error.hpp"
#include "certified_constants.hpp"
#include "certified_complex_transcendental.hpp"
#include "certified_exponential.hpp"
#include "certified_logarithm.hpp"
#include "certified_sqrt.hpp"
#include "certified_trigonometry.hpp"
#include "interval_math.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
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
    if (absRational(z) > rational(160))
        throw CertifiedBackendUnsupported{
            "hypergeometric1F1 certified series currently requires |z| <= 160 for bounded work"};

    const auto absA = ceilAbsToUint64(a);
    const auto absB = ceilAbsToUint64(b);
    const auto absZ = ceilAbsToUint64(z);
    if (!absA || !absB || !absZ)
        throw CertifiedBackendUnsupported{"hypergeometric1F1 argument is too large for the series backend"};

    // j>=2|a|,2|b|,6|z| なら
    // |t_{j+1}/t_j| = |z||a+j|/(|b+j|(j+1)) <= 1/2。
    // 以後のtailは次項の2倍で厳密に上から押さえられる。
    const std::uint64_t ratioStart = std::max({
        2U * *absA + 2U,
        2U * *absB + 2U,
        6U * *absZ + 2U});
    constexpr std::uint64_t maximumTerms = 200000;
    if (ratioStart > maximumTerms)
        throw CertifiedBackendUnsupported{"hypergeometric1F1 requires too many series terms"};

    Rational term{BigInt{1}};
    Rational sum{BigInt{1}};
    const Rational tolerance = binaryThreshold(checkedAdd(
        precisionBits, 16, "hypergeometric1F1 precision is too large"));

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
    if (absZ > rational(9, 10))
        throw CertifiedBackendUnsupported{
            "hypergeometric2F1 certified series currently requires |z| <= 9/10 for bounded work"};

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
    if (absM > rational(9, 10))
        throw CertifiedBackendUnsupported{
            "elliptic certified series currently requires |m| <= 9/10 for bounded work"};
    if (kind == EllipticSeriesKind::Pi && absN > rational(9, 10))
        throw CertifiedBackendUnsupported{
            "ellipticPi certified series currently requires |n| <= 9/10 for bounded work"};

    const std::size_t bits = checkedAdd(
        precisionBits, 32, "elliptic working precision is too large");
    const auto sinBox = encloseSinRadian(phi, bits).interval;
    const auto cosBox = encloseCosRadian(phi, bits).interval;
    const RealInterval sinSquared = multiply(sinBox, sinBox, bits);
    RealInterval sineOdd = sinBox;
    RealInterval integral = exactInterval(phi, bits); // I_0(phi)=phi
    RealInterval sum = integral;

    Rational c{BigInt{1}};       // (1/2)_k/k!
    Rational e{BigInt{1}};       // coefficients of sqrt(1-x)
    Rational mPower{BigInt{1}};
    Rational q{BigInt{1}};       // Pi combined coefficient
    const Rational r = kind == EllipticSeriesKind::Pi
        ? (absM < absN ? absN : absM) : absM;
    const Rational tolerance = binaryThreshold(checkedAdd(
        precisionBits, 18, "elliptic precision is too large"));

    // 旧実装ではtail評価のたびにr^(k+1)を1から掛け直していたため、
    // 高精度ほど不要なO(k^2) Rational乗算が増えていた。級数本体と同様に
    // 冪を逐次更新し、保証境界は変えずにtail評価だけをO(k)へ落とす。
    Rational rPower = r;
    constexpr std::size_t maximumTerms = 200000;
    for (std::size_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        integral = evenSinePowerIntegral(
            k, cosBox, integral, sineOdd, sinSquared, bits);
        mPower *= m;
        c *= Rational{BigInt::fromUnsigned(2 * k - 1), BigInt::fromUnsigned(2 * k)};

        Rational coefficient;
        if (kind == EllipticSeriesKind::F) {
            coefficient = c * mPower;
        }
        else if (kind == EllipticSeriesKind::E) {
            e *= Rational{BigInt{static_cast<std::int64_t>(2 * k) - 3},
                BigInt::fromUnsigned(2 * k)};
            coefficient = e * Rational{BigInt::fromUnsigned(1)};
            // e already contains the sign and m^k is applied separately.
            coefficient *= mPower;
        }
        else {
            q = n * q + c * mPower;
            coefficient = q;
        }
        sum = add(sum, intervalTimesRational(integral, coefficient, bits), bits);

        Rational tail;
        rPower *= r;
        if (kind != EllipticSeriesKind::Pi) {
            tail = absRational(phi) * rPower / (rational(1) - r);
        }
        else {
            const Rational kp1{BigInt::fromUnsigned(k + 1)};
            const Rational kp2{BigInt::fromUnsigned(k + 2)};
            tail = absRational(phi) * rPower
                * (kp2 - kp1 * r)
                / ((rational(1) - r) * (rational(1) - r));
        }
        if (tail <= tolerance) {
            const RealInterval error = RealInterval::fromRationalBounds(-tail, tail, bits);
            return add(sum, error, bits).roundedOutward(precisionBits);
        }
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
        precisionBits, 40, "elliptic interval working precision is too large");
    const RealInterval nWork = n.roundedOutward(bits);
    const RealInterval phiWork = phi.roundedOutward(bits);
    const RealInterval mWork = m.roundedOutward(bits);

    const RealInterval absMInterval = absoluteInterval(mWork, bits);
    const RealInterval absNInterval = absoluteInterval(nWork, bits);
    const Rational absM = absMInterval.upper().toRational();
    const Rational absN = absNInterval.upper().toRational();
    const Rational limit = rational(9, 10);
    if (absM > limit) {
        if (absMInterval.lower().toRational() <= limit)
            throw PrecisionInsufficient{
                "Elliptic parameter interval straddles the |m| = 9/10 work boundary"};
        throw CertifiedBackendUnsupported{
            "elliptic certified series currently requires |m| <= 9/10 for bounded work"};
    }
    if (kind == EllipticSeriesKind::Pi && absN > limit) {
        if (absNInterval.lower().toRational() <= limit)
            throw PrecisionInsufficient{
                "Elliptic Pi parameter interval straddles the |n| = 9/10 work boundary"};
        throw CertifiedBackendUnsupported{
            "ellipticPi certified series currently requires |n| <= 9/10 for bounded work"};
    }

    const RealInterval sine = encloseSinRadianInterval(phiWork, bits).interval;
    const RealInterval cosine = encloseCosRadianInterval(phiWork, bits).interval;
    const RealInterval sineSquared = multiply(sine, sine, bits);
    RealInterval sineOdd = sine;
    RealInterval integral = phiWork;
    RealInterval sum = integral;

    Rational c{BigInt{1}};
    Rational e{BigInt{1}};
    RealInterval mPower = exactInterval(1, bits);
    RealInterval q = exactInterval(1, bits);
    const Rational r = kind == EllipticSeriesKind::Pi
        ? (absM < absN ? absN : absM) : absM;
    const Rational phiMagnitude = intervalAbsUpper(phiWork, bits);
    const Rational tolerance = binaryThreshold(checkedAdd(
        precisionBits, 18, "elliptic interval precision is too large"));

    Rational rPower = r;
    constexpr std::size_t maximumTerms = 200000;
    for (std::size_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        integral = evenSinePowerIntegral(
            k, cosine, integral, sineOdd, sineSquared, bits);
        mPower = multiply(mPower, mWork, bits);
        c *= Rational{BigInt::fromUnsigned(2 * k - 1), BigInt::fromUnsigned(2 * k)};

        RealInterval coefficient = exactInterval(0, bits);
        if (kind == EllipticSeriesKind::F) {
            coefficient = intervalTimesRational(mPower, c, bits);
        }
        else if (kind == EllipticSeriesKind::E) {
            e *= Rational{BigInt{static_cast<std::int64_t>(2 * k) - 3},
                BigInt::fromUnsigned(2 * k)};
            coefficient = intervalTimesRational(mPower, e, bits);
        }
        else {
            q = add(
                multiply(nWork, q, bits),
                intervalTimesRational(mPower, c, bits),
                bits);
            coefficient = q;
        }
        sum = add(sum, multiply(integral, coefficient, bits), bits);

        rPower *= r;
        Rational tail;
        if (kind != EllipticSeriesKind::Pi) {
            tail = r.isZero()
                ? Rational{}
                : phiMagnitude * rPower / (rational(1) - r);
        }
        else {
            const Rational kp1{BigInt::fromUnsigned(k + 1)};
            const Rational kp2{BigInt::fromUnsigned(k + 2)};
            tail = r.isZero()
                ? Rational{}
                : phiMagnitude * rPower
                    * (kp2 - kp1 * r)
                    / ((rational(1) - r) * (rational(1) - r));
        }
        if (tail <= tolerance)
            return add(sum, symmetricError(tail, bits), bits).roundedOutward(precisionBits);
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

[[nodiscard]] RealInterval pointExponentialIntegralEi(
    Rational x,
    std::size_t precisionBits) {
    if (x.isZero())
        throw std::domain_error("Ei is undefined at zero");
    const Rational absX = absRational(x);
    constexpr std::int64_t maximumCertifiedSeriesArgument = 96;
    if (absX > rational(maximumCertifiedSeriesArgument))
        throw CertifiedBackendUnsupported{"Ei certified series currently requires |x| <= 96"};

    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "Ei working precision is too large");
    RealInterval sum = add(
        encloseEulerGamma(workBits),
        encloseLogPositive(exactInterval(absX, workBits), workBits).interval,
        workBits);

    Rational factorialPower = x; // x^k/k!, k=1
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "Ei target precision is too large"));
    constexpr std::uint64_t maximumTerms = 200'000;
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const Rational term = factorialPower / Rational{BigInt::fromUnsigned(k)};
        sum = add(sum, exactInterval(term, workBits), workBits);

        if (absX < Rational{BigInt::fromUnsigned(k + 1)}) {
            const Rational ratio = absX * Rational{BigInt::fromUnsigned(k)}
                / Rational{numeric::pow(BigInt::fromUnsigned(k + 1), 2)};
            const Rational q = absX / Rational{BigInt::fromUnsigned(k + 1)};
            const Rational nextBound = absRational(term) * ratio;
            const Rational tail = nextBound / (rational(1) - q);
            if (tail <= target) {
                sum = add(sum, symmetricError(tail, workBits), workBits);
                return sum.roundedOutward(precisionBits);
            }
        }

        factorialPower *= x;
        factorialPower /= Rational{BigInt::fromUnsigned(k + 1)};
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
    constexpr std::int64_t maximumCertifiedSeriesArgument = 96;
    if (x > rational(maximumCertifiedSeriesArgument))
        throw CertifiedBackendUnsupported{"Si certified series currently requires |x| <= 96"};

    const auto magnitude = ceilAbsToUint64(x);
    if (!magnitude || *magnitude > std::numeric_limits<std::size_t>::max() / 2)
        throw CertifiedBackendUnsupported{"Si argument is too large for cancellation planning"};
    const std::size_t cancellationBits = static_cast<std::size_t>(*magnitude) * 2;
    const std::size_t workBits = checkedAdd(
        checkedAdd(precisionBits, 48, "Si working precision is too large"),
        cancellationBits, "Si cancellation precision is too large");
    const Rational x2 = x * x;
    Rational term = x;
    RealInterval sum = exactInterval(term, workBits);
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "Si target precision is too large"));

    constexpr std::uint64_t maximumTerms = 200'000;
    for (std::uint64_t k = 0; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const std::uint64_t a = 2 * k + 1;
        const std::uint64_t b = 2 * k + 2;
        const std::uint64_t c = 2 * k + 3;
        const Rational ratio = x2 * Rational{BigInt::fromUnsigned(a)}
            / Rational{BigInt::fromUnsigned(c) * BigInt::fromUnsigned(c) * BigInt::fromUnsigned(b)};
        const Rational nextBound = absRational(term) * ratio;
        const Rational q = x2
            / Rational{BigInt::fromUnsigned(b) * BigInt::fromUnsigned(c)};
        if (q < rational(1)) {
            const Rational tail = nextBound / (rational(1) - q);
            if (tail <= target) {
                sum = add(sum, symmetricError(tail, workBits), workBits);
                return sum.roundedOutward(precisionBits);
            }
        }
        term = -(term * ratio);
        sum = add(sum, exactInterval(term, workBits), workBits);
    }
    throw CertifiedBackendUnsupported{"Si series did not converge within the term limit"};
}

[[nodiscard]] RealInterval pointCosineIntegralCiPositive(
    const Rational& x,
    std::size_t precisionBits) {
    if (x <= rational(0))
        throw std::domain_error("Ci real certified backend requires x > 0");
    constexpr std::int64_t maximumCertifiedSeriesArgument = 96;
    if (x > rational(maximumCertifiedSeriesArgument))
        throw CertifiedBackendUnsupported{"Ci certified series currently requires 0 < x <= 96"};

    const auto magnitude = ceilAbsToUint64(x);
    if (!magnitude || *magnitude > std::numeric_limits<std::size_t>::max() / 2)
        throw CertifiedBackendUnsupported{"Ci argument is too large for cancellation planning"};
    const std::size_t cancellationBits = static_cast<std::size_t>(*magnitude) * 2;
    const std::size_t workBits = checkedAdd(
        checkedAdd(precisionBits, 56, "Ci working precision is too large"),
        cancellationBits, "Ci cancellation precision is too large");
    RealInterval sum = add(
        encloseEulerGamma(workBits),
        encloseLogPositive(exactInterval(x, workBits), workBits).interval,
        workBits);

    const Rational x2 = x * x;
    Rational term = -(x2 / rational(4)); // k=1: -x^2/(2*2!)
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "Ci target precision is too large"));
    constexpr std::uint64_t maximumTerms = 200'000;
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        sum = add(sum, exactInterval(term, workBits), workBits);
        const std::uint64_t a = 2 * k;
        const std::uint64_t b = 2 * k + 1;
        const std::uint64_t c = 2 * k + 2;
        const Rational ratio = x2 * Rational{BigInt::fromUnsigned(a)}
            / Rational{BigInt::fromUnsigned(c) * BigInt::fromUnsigned(c) * BigInt::fromUnsigned(b)};
        const Rational nextBound = absRational(term) * ratio;
        const Rational q = x2
            / Rational{BigInt::fromUnsigned(b) * BigInt::fromUnsigned(c)};
        if (q < rational(1)) {
            const Rational tail = nextBound / (rational(1) - q);
            if (tail <= target) {
                sum = add(sum, symmetricError(tail, workBits), workBits);
                return sum.roundedOutward(precisionBits);
            }
        }
        term = -(term * ratio);
    }
    throw CertifiedBackendUnsupported{"Ci series did not converge within the term limit"};
}

[[nodiscard]] RealInterval pointPolylogPositiveOrder(
    std::uint64_t order,
    Rational z,
    std::size_t precisionBits) {
    if (order == 0)
        throw std::invalid_argument("polylog certified series requires positive order");
    const Rational absZ = absRational(z);
    if (absZ > rational(49, 50))
        throw CertifiedBackendUnsupported{
            "polylog certified series currently requires |z| <= 49/50 for bounded work"};
    if (z.isZero())
        return exactInterval(0, precisionBits);

    const std::size_t workBits = checkedAdd(
        precisionBits, 40, "polylog working precision is too large");
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 18, "polylog target precision is too large"));
    RealInterval sum = exactInterval(0, workBits);
    Rational zPower{BigInt{1}};

    constexpr std::uint64_t maximumTerms = 1'000'000;
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        zPower *= z;
        const BigInt denominator = numeric::pow(BigInt::fromUnsigned(k), order);
        const Rational term = zPower / Rational{denominator};
        sum = add(sum, exactInterval(term, workBits), workBits);

        const BigInt nextDenominator = numeric::pow(BigInt::fromUnsigned(k + 1), order);
        const Rational nextBound = absRational(zPower * z) / Rational{nextDenominator};
        const Rational tail = nextBound / (rational(1) - absZ);
        if (tail <= target) {
            sum = add(sum, symmetricError(tail, workBits), workBits);
            return sum.roundedOutward(precisionBits);
        }
    }
    throw CertifiedBackendUnsupported{"polylog series did not converge within the term limit"};
}

[[nodiscard]] Rational rationalPowerInteger(Rational base, std::size_t exponent) {
    Rational result{BigInt{1}};
    while (exponent != 0) {
        consumeCertifiedWork();
        if ((exponent & 1U) != 0)
            result *= base;
        exponent >>= 1U;
        if (exponent != 0)
            base *= base;
    }
    return result;
}

[[nodiscard]] Rational risingRational(Rational value, std::size_t count) {
    Rational result{BigInt{1}};
    for (std::size_t i = 0; i < count; ++i) {
        consumeCertifiedWork();
        result *= value + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(i))};
    }
    return result;
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

    const RealInterval logarithm = encloseLogPositive(
        exactInterval(Rational{BigInt::fromUnsigned(base)}, precisionBits),
        precisionBits).interval;
    const RealInterval scaled = multiply(
        logarithm, exactInterval(exponent, precisionBits), precisionBits);
    return encloseExp(scaled, precisionBits).interval;
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
        for (std::size_t k = 2; k <= 64; ++k) {
            consumeCertifiedWork();
            const Rational rising = risingRational(s, 2 * k - 1);
            const BigInt factorial = numeric::factorial(static_cast<std::uint64_t>(2 * k));
            const Rational coefficient = absRational(bernoulliEven(k))
                * rising / Rational{factorial};
            const Rational exponent = -(s + Rational{BigInt::fromUnsigned(
                static_cast<std::uint64_t>(2 * k - 1))});
            const RealInterval power = positiveIntegerBasePower(n, exponent, workBits);
            const Rational bound = coefficient * power.upper().toRational();
            if (bound <= target) {
                chosenN = n;
                chosenK = k;
                chosenRising = rising;
                chosenPower = power;
                break;
            }
        }
    }
    if (chosenN == 0)
        throw CertifiedBackendUnsupported{"zeta Euler-Maclaurin budget is insufficient for the requested precision"};

    RealInterval result = exactInterval(0, workBits);
    for (std::uint64_t n = 1; n < chosenN; ++n) {
        consumeCertifiedWork();
        const RealInterval term = positiveIntegerBasePower(n, -s, workBits);
        result = add(result, term, workBits);
    }

    const RealInterval integralTail = divide(
        positiveIntegerBasePower(chosenN, rational(1) - s, workBits),
        exactInterval(s - rational(1), workBits), workBits);
    result = add(result, integralTail, workBits);
    result = add(result, multiply(
        exactInterval(rational(1, 2), workBits),
        positiveIntegerBasePower(chosenN, -s, workBits), workBits), workBits);

    for (std::size_t k = 1; k < chosenK; ++k) {
        consumeCertifiedWork();
        const Rational rising = risingRational(s, 2 * k - 1);
        const BigInt factorial = numeric::factorial(static_cast<std::uint64_t>(2 * k));
        const Rational coefficient = bernoulliEven(k) * rising / Rational{factorial};
        const Rational exponent = -(s + Rational{BigInt::fromUnsigned(
            static_cast<std::uint64_t>(2 * k - 1))});
        result = add(result, multiply(
            exactInterval(coefficient, workBits),
            positiveIntegerBasePower(chosenN, exponent, workBits), workBits), workBits);
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
    if (boundaryUpper > rational(512)) {
        if (complexAbsLower(z, boundaryBits) > rational(512))
            throw CertifiedBackendUnsupported{
                "complex Ei certified series currently requires |z| <= 512 for bounded work"};
        throw PrecisionInsufficient{
            "complex Ei magnitude interval straddles the |z| = 512 work boundary"};
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
    Rational majorant = q;

    constexpr std::uint64_t maximumTerms = 200000;
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const Rational ratio = q * Rational{BigInt::fromUnsigned(k)}
            / Rational{numeric::pow(BigInt::fromUnsigned(k + 1), 2)};
        const Rational nextMajorant = majorant * ratio;
        const Rational futureRatio = q / Rational{BigInt::fromUnsigned(k + 1)};
        if (futureRatio < rational(1)) {
            const Rational tail = nextMajorant / (rational(1) - futureRatio);
            if (tail <= target)
                return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
        }
        term = multiplyByRational(
            multiply(term, z, workBits),
            Rational{BigInt::fromUnsigned(k), numeric::pow(BigInt::fromUnsigned(k + 1), 2)},
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
    if (boundaryUpper > rational(128)) {
        if (complexAbsLower(z, boundaryBits) > rational(128))
            throw CertifiedBackendUnsupported{
                "complex Ci certified series currently requires |z| <= 128 for bounded work"};
        throw PrecisionInsufficient{
            "complex Ci magnitude interval straddles the |z| = 128 work boundary"};
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
    Rational majorant = qSquared / rational(4);

    constexpr std::uint64_t maximumTerms = 200000;
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        const BigInt a = BigInt::fromUnsigned(2 * k);
        const BigInt b = BigInt::fromUnsigned(2 * k + 1);
        const BigInt c = BigInt::fromUnsigned(2 * k + 2);
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
            multiply(term, zSquared, workBits),
            -Rational{a, b * c * c}, workBits);
        sum = add(sum, term, workBits);
        majorant = nextMajorant;
    }
    throw CertifiedBackendUnsupported{"complex Ci series did not converge within the term limit"};
}

[[nodiscard]] ComplexInterval pointFresnelComplex(
    const ComplexInterval& z,
    std::size_t precisionBits,
    bool sineVariant) {
    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "complex Fresnel working precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    const Rational qLower = complexAbsLower(z, workBits);
    // 複素Fresnelは現状Maclaurin級数のみで，大きい|z|では巨大中間項の
    // 相殺が性能の崖になる。実数経路には漸近backendがあるため，複素経路だけ
    // audited boundary |z|=8 でbounded-workを保証する。
    if (qLower >= rational(8))
        throw CertifiedBackendUnsupported{
            "complex Fresnel series currently requires |z| < 8 for bounded work"};
    if (q >= rational(8))
        throw PrecisionInsufficient{
            "complex Fresnel magnitude interval straddles the |z| = 8 work boundary"};
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

    ComplexInterval term = z.roundedOutward(workBits);
    Rational majorant = q;
    if (sineVariant) {
        term = multiply(term, zSquared, workBits);
        term = multiply(term, ComplexInterval::fromReal(pi), workBits);
        term = multiplyByRational(term, rational(1, 6), workBits);
        majorant = piUpper * q * qSquared / rational(6);
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

        ComplexInterval next = multiply(term, zFourth, workBits);
        next = multiply(next, ComplexInterval::fromReal(piSquaredQuarter), workBits);
        next = multiplyByRational(next, -Rational{numerator, denominator}, workBits);
        const Rational nextMajorant = majorant
            * baseMajorant * Rational{numerator, denominator};
        if (futureRatio < rational(1)) {
            const Rational tail = nextMajorant / (rational(1) - futureRatio);
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

    const std::size_t workBits = checkedAdd(
        precisionBits, 56, "complex hypergeometric1F1 working precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    if (q > rational(160)) {
        if (complexAbsLower(z, workBits) > rational(160))
            throw CertifiedBackendUnsupported{
                "complex hypergeometric1F1 certified series currently requires |z| <= 160 for bounded work"};
        throw PrecisionInsufficient{
            "complex hypergeometric1F1 magnitude interval straddles the |z| = 160 work boundary"};
    }
    const Rational absA = complexAbsUpper(a, workBits);
    const Rational absB = complexAbsUpper(b, workBits);
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 22, "complex hypergeometric1F1 target precision is too large"));

    ComplexInterval term = exactComplex(1, workBits);
    ComplexInterval sum = term;
    constexpr std::uint64_t maximumTerms = 250000;
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
        const Rational denominator = rationalPower(positiveLower, 2 * n - 1);
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

[[nodiscard]] ComplexInterval pointPolylogComplex(
    std::uint64_t order,
    const ComplexInterval& z,
    std::size_t precisionBits) {
    if (order == 0)
        throw std::invalid_argument("polylog certified complex series requires positive order");
    const std::size_t workBits = checkedAdd(
        precisionBits, 48, "complex polylog working precision is too large");
    const Rational q = complexAbsUpper(z, workBits);
    if (q > rational(49, 50)) {
        if (complexAbsLower(z, workBits) > rational(49, 50))
            throw CertifiedBackendUnsupported{
                "polylog certified complex series currently requires |z| <= 49/50 for bounded work"};
        throw PrecisionInsufficient{
            "complex polylog magnitude interval straddles the |z| = 49/50 work boundary"};
    }
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 20, "complex polylog target precision is too large"));
    ComplexInterval sum = exactComplex(0, workBits);
    ComplexInterval power = exactComplex(1, workBits);
    Rational qPower{BigInt{1}};
    constexpr std::uint64_t maximumTerms = 1'000'000;
    for (std::uint64_t k = 1; k < maximumTerms; ++k) {
        consumeCertifiedWork();
        power = multiply(power, z, workBits);
        qPower *= q;
        const BigInt denominator = numeric::pow(BigInt::fromUnsigned(k), order);
        sum = add(sum, multiplyByRational(
            power, Rational{BigInt{1}, denominator}, workBits), workBits);

        const BigInt nextDenominator = numeric::pow(BigInt::fromUnsigned(k + 1), order);
        const Rational nextBound = qPower * q / Rational{nextDenominator};
        const Rational tail = nextBound / (rational(1) - q);
        if (tail <= target)
            return inflateComplex(sum, tail, workBits).roundedOutward(precisionBits);
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

    if (q <= rational(9, 10)) {
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
    if (qLower <= rational(9, 10))
        throw PrecisionInsufficient{
            "complex hypergeometric2F1 magnitude interval straddles the |z| = 9/10 work boundary"};
    if (qLower < rational(1))
        throw CertifiedBackendUnsupported{
            "complex hypergeometric2F1 series currently requires |z| <= 9/10 for bounded work"};
    if (qLower <= rational(1))
        throw PrecisionInsufficient{
            "complex hypergeometric2F1 magnitude interval straddles the |z| = 1 continuation boundary"};

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
        const ComplexInterval gammaC = pointGammaComplex(c, workBits);
        const ComplexInterval gammaBMinusA = pointGammaComplex(subtract(b, a, workBits), workBits);
        const ComplexInterval gammaAMinusB = pointGammaComplex(subtract(a, b, workBits), workBits);
        const ComplexInterval gammaB = pointGammaComplex(b, workBits);
        const ComplexInterval gammaA = pointGammaComplex(a, workBits);
        const ComplexInterval gammaCMinusA = pointGammaComplex(subtract(c, a, workBits), workBits);
        const ComplexInterval gammaCMinusB = pointGammaComplex(subtract(c, b, workBits), workBits);

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

[[nodiscard]] ComplexInterval pointZetaComplex(
    const ComplexInterval& s,
    std::size_t precisionBits) {
    const Rational sigma = s.real().lower().toRational();
    if (sigma <= rational(1))
        throw CertifiedBackendUnsupported{"zeta certified complex backend currently requires Re(s) > 1"};
    const std::size_t workBits = checkedAdd(
        precisionBits, 72, "complex zeta working precision is too large");
    const Rational target = binaryThreshold(checkedAdd(
        precisionBits, 24, "complex zeta target precision is too large"));
    const Rational absS = complexAbsUpper(s, workBits);
    const Rational twoPiLower = rational(2) * enclosePi(workBits).interval.lower().toRational();

    std::uint64_t chosenN = 0;
    std::size_t chosenK = 0;
    Rational remainderBound;
    for (std::uint64_t N : {8ULL, 16ULL, 32ULL, 64ULL, 128ULL}) {
        for (std::size_t k = 2; k <= 40; ++k) {
            consumeCertifiedWork();
            Rational risingBound{BigInt{1}};
            for (std::size_t j = 0; j < 2 * k; ++j)
                risingBound *= absS + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(j))};
            const Rational fourierDenominator = rationalPower(twoPiLower, 2 * k);
            const Rational exponent = rational(1) - sigma
                - Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k))};
            const Rational nPower = positiveIntegerBasePower(N, exponent, workBits)
                .upper().toRational();
            // |B_2k({x})| <= 2(2k)! zeta(2k)/(2pi)^(2k), zeta(2k)<2.
            const Rational bound = rational(4) * risingBound * nPower
                / (fourierDenominator
                    * (sigma + Rational{BigInt::fromUnsigned(static_cast<std::uint64_t>(2 * k - 1))}));
            if (bound <= target) {
                chosenN = N;
                chosenK = k;
                remainderBound = bound;
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
    for (std::uint64_t n = 1; n < chosenN; ++n) {
        consumeCertifiedWork();
        result = add(result, positiveIntegerComplexPower(n, negativeS, workBits), workBits);
    }

    const ComplexInterval oneMinusS = subtract(exactComplex(1, workBits), s, workBits);
    const ComplexInterval integralNumerator = positiveIntegerComplexPower(
        chosenN, oneMinusS, workBits);
    const ComplexInterval sMinusOne = subtract(s, exactComplex(1, workBits), workBits);
    result = add(result, divide(integralNumerator, sMinusOne, workBits), workBits);
    result = add(result, multiplyByRational(
        positiveIntegerComplexPower(chosenN, negativeS, workBits), rational(1, 2), workBits), workBits);

    for (std::size_t k = 1; k < chosenK; ++k) {
        consumeCertifiedWork();
        ComplexInterval rising = exactComplex(1, workBits);
        for (std::size_t j = 0; j < 2 * k - 1; ++j)
            rising = multiply(rising,
                add(s, exactComplex(static_cast<std::int64_t>(j), workBits), workBits), workBits);
        const BigInt factorial = numeric::factorial(static_cast<std::uint64_t>(2 * k));
        const Rational coefficient = bernoulliEven(k) / Rational{factorial};
        const ComplexInterval exponent = negate(add(
            s, exactComplex(static_cast<std::int64_t>(2 * k - 1), workBits), workBits));
        ComplexInterval correction = multiply(
            rising, positiveIntegerComplexPower(chosenN, exponent, workBits), workBits);
        correction = multiplyByRational(correction, coefficient, workBits);
        result = add(result, correction, workBits);
    }
    return inflateComplex(result, remainderBound, workBits).roundedOutward(precisionBits);
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
    return pointEllipticSeries(EllipticSeriesKind::F, rational(0), phi, m, precisionBits);
}

RealInterval encloseEllipticEReal(
    const Rational& phi, const Rational& m, std::size_t precisionBits) {
    return pointEllipticSeries(EllipticSeriesKind::E, rational(0), phi, m, precisionBits);
}

RealInterval encloseEllipticPiReal(
    const Rational& n, const Rational& phi, const Rational& m, std::size_t precisionBits) {
    return pointEllipticSeries(EllipticSeriesKind::Pi, n, phi, m, precisionBits);
}

RealInterval encloseEllipticFReal(
    const RealInterval& phi, const RealInterval& m, std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticF precision must be at least one bit");
    return intervalEllipticSeries(
        EllipticSeriesKind::F, exactInterval(0, precisionBits), phi, m, precisionBits);
}

RealInterval encloseEllipticEReal(
    const RealInterval& phi, const RealInterval& m, std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticE precision must be at least one bit");
    return intervalEllipticSeries(
        EllipticSeriesKind::E, exactInterval(0, precisionBits), phi, m, precisionBits);
}

RealInterval encloseEllipticPiReal(
    const RealInterval& n, const RealInterval& phi, const RealInterval& m,
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("ellipticPi precision must be at least one bit");
    return intervalEllipticSeries(
        EllipticSeriesKind::Pi, n, phi, m, precisionBits);
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
    if (lower <= rational(1))
        throw CertifiedBackendUnsupported{
            "zeta certified real backend currently requires s > 1"};
    const RealInterval lowValue = pointZetaGreaterThanOne(upper, precisionBits);
    const RealInterval highValue = pointZetaGreaterThanOne(lower, precisionBits);
    return RealInterval{lowValue.lower(), highValue.upper()};
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
