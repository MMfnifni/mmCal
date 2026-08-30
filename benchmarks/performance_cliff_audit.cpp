#include "performance_cliff_audit.hpp"

#include "error/error_message.hpp"
#include "expression/expr.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "symbolic/algebraic_number.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::benchmarks {
namespace {

using Clock = std::chrono::steady_clock;
using mmcal::evaluation::EvaluationUsage;
using mmcal::numeric::BigInt;
using mmcal::numeric::Rational;

struct SweepPoint final {
    std::string axis;
    std::string expression;
};

struct Sweep final {
    std::string name;
    std::vector<SweepPoint> points;
    bool compareAdjacent = true;
};

enum class OutcomeKind {
    Numeric,
    Unevaluated,
    Error
};

struct Measurement final {
    OutcomeKind outcome = OutcomeKind::Error;
    std::string outcomeText;
    double firstMilliseconds = 0.0;
    double warmMilliseconds = 0.0;
    EvaluationUsage usage;
};

[[nodiscard]] bool isNumericResult(const mmcal::expression::Expr& value) {
    return value.isNumber()
        || value.isDecimalApproximation()
        || value.isComplexDecimalApproximation();
}

[[nodiscard]] std::string shortExpression(const mmcal::expression::Expr& value) {
    std::string text = mmcal::formatting::formatExpr(value);
    constexpr std::size_t maximumLength = 48;
    if (text.size() > maximumLength) {
        text.resize(maximumLength - 3);
        text += "...";
    }
    return text;
}

[[nodiscard]] Measurement measure(
    std::string_view expression,
    std::size_t iterations) {
    mmcal::kernel::KernelSession session;
    Measurement result;
    double warmTotal = 0.0;

    for (std::size_t iteration = 0; iteration < iterations; ++iteration) {
        session.resetForIndependentEvaluation();
        const auto start = Clock::now();
        try {
            const auto value = session.evaluate(expression);
            const auto end = Clock::now();
            const double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
            if (iteration == 0) {
                result.firstMilliseconds = elapsed;
                result.usage = session.lastEvaluationUsage();
                result.outcome = isNumericResult(value)
                    ? OutcomeKind::Numeric : OutcomeKind::Unevaluated;
                result.outcomeText = isNumericResult(value)
                    ? "numeric" : "held:" + shortExpression(value);
            }
            else {
                warmTotal += elapsed;
            }
        }
        catch (const mmcal::error::CalcError& error) {
            const auto end = Clock::now();
            const double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
            if (iteration == 0) {
                result.firstMilliseconds = elapsed;
                result.usage = session.lastEvaluationUsage();
                result.outcome = OutcomeKind::Error;
                result.outcomeText = std::string{mmcal::error::calcErrorTypeName(error.type())};
            }
            else {
                warmTotal += elapsed;
            }
        }
    }

    result.warmMilliseconds = iterations > 1
        ? warmTotal / static_cast<double>(iterations - 1)
        : result.firstMilliseconds;
    return result;
}

[[nodiscard]] double ratio(std::size_t current, std::size_t previous) {
    if (previous == 0)
        return current == 0 ? 1.0 : std::numeric_limits<double>::infinity();
    return static_cast<double>(current) / static_cast<double>(previous);
}

[[nodiscard]] std::string marker(
    const Measurement* previous,
    const Measurement& current,
    double cliffRatio) {
    if (!previous)
        return {};
    if (previous->outcome != current.outcome)
        return "BOUNDARY";
    if (current.outcome != OutcomeKind::Numeric)
        return {};

    const double timeRatio = current.warmMilliseconds
        / std::max(previous->warmMilliseconds, 0.001);
    const double refinementRatio = ratio(
        current.usage.certifiedRefinements,
        previous->usage.certifiedRefinements);

    if (timeRatio >= cliffRatio && current.warmMilliseconds >= 0.5)
        return "CLIFF(time x" + std::to_string(timeRatio).substr(0, 4) + ')';
    if (refinementRatio >= cliffRatio
        && current.usage.certifiedRefinements
            >= previous->usage.certifiedRefinements + 256)
        return "CLIFF(work x" + std::to_string(refinementRatio).substr(0, 4) + ')';
    return {};
}

void printSweep(
    const Sweep& sweep,
    std::size_t iterations,
    double cliffRatio) {
    std::cout << "\n[" << sweep.name << "]\n";
    std::cout << "  " << std::left << std::setw(12) << "axis"
              << std::right << std::setw(11) << "first-ms"
              << std::setw(11) << (iterations > 1 ? "warm-ms" : "time-ms")
              << std::setw(12) << "cert-ref"
              << std::setw(11) << "work"
              << std::setw(10) << "nodes"
              << "  outcome / marker\n";

    Measurement previous;
    bool hasPrevious = false;
    for (const auto& point : sweep.points) {
        const Measurement current = measure(point.expression, iterations);
        const std::string currentMarker = sweep.compareAdjacent
            ? marker(hasPrevious ? &previous : nullptr, current, cliffRatio)
            : std::string{};
        std::cout << "  " << std::left << std::setw(12) << point.axis
                  << std::right << std::fixed << std::setprecision(3)
                  << std::setw(11) << current.firstMilliseconds
                  << std::setw(11) << current.warmMilliseconds
                  << std::setw(12) << current.usage.certifiedRefinements
                  << std::setw(11) << current.usage.evaluationSteps
                  << std::setw(10) << current.usage.generatedNodes
                  << "  " << current.outcomeText;
        if (!currentMarker.empty())
            std::cout << "  << " << currentMarker;
        std::cout << '\n';
        previous = current;
        hasPrevious = true;
    }
}

[[nodiscard]] Sweep precisionSweep(
    std::string name,
    std::string expressionTemplate) {
    Sweep result{std::move(name), {}, true};
    for (const std::size_t digits : {10U, 20U, 40U, 80U}) {
        std::string expression = expressionTemplate;
        const std::string token{"{p}"};
        const auto position = expression.find(token);
        expression.replace(position, token.size(), std::to_string(digits));
        result.points.push_back({std::to_string(digits) + " digits", std::move(expression)});
    }
    return result;
}

[[nodiscard]] std::vector<Sweep> specialFunctionSweeps() {
    std::vector<Sweep> sweeps;

    // 高精度でbackend/guardの伸び方が変わる代表函数。
    sweeps.push_back(precisionSweep("gamma precision", "N[gamma[1/3],{p}]"));
    sweeps.push_back(precisionSweep("zeta precision", "N[zeta[3/2],{p}]"));
    sweeps.push_back(precisionSweep("complex zeta precision", "N[zeta[3/2+I],{p}]"));
    sweeps.push_back(precisionSweep("ibeta precision", "N[ibeta[1/3,2/3,1/4],{p}]"));
    sweeps.push_back(precisionSweep("2F1 precision", "N[hypergeometric2F1[1/2,1/3,5/4,4/5],{p}]"));
    sweeps.push_back(precisionSweep(
        "2F1 continuation precision",
        "N[hypergeometric2F1[3.4,5.6,4+I,4.6+2I],{p}]"));
    sweeps.push_back(precisionSweep("complex li precision", "N[li[-2+I],{p}]"));
    sweeps.push_back(precisionSweep("1F1 z=512 precision", "N[hypergeometric1F1[1/2,5/4,512],{p}]"));
    sweeps.push_back(precisionSweep("complex 1F1 z=256+I precision", "N[hypergeometric1F1[1/2,5/4,256+I],{p}]"));
    sweeps.push_back(precisionSweep("complex Fresnel 7+I precision", "N[fresnelc[7+I],{p}]"));
    sweeps.push_back(precisionSweep("complex Fresnel 32+I precision", "N[fresnelc[32+I],{p}]"));

    sweeps.push_back({"1F1 magnitude scaling", {
        {"z=64", "N[hypergeometric1F1[1/2,5/4,64],20]"},
        {"z=160", "N[hypergeometric1F1[1/2,5/4,160],20]"},
        {"z=161", "N[hypergeometric1F1[1/2,5/4,161],20]"},
        {"z=256", "N[hypergeometric1F1[1/2,5/4,256],20]"},
        {"z=512", "N[hypergeometric1F1[1/2,5/4,512],20]"},
        {"z=1000", "N[hypergeometric1F1[1/2,5/4,1000],20]"},
        {"z=-160", "N[hypergeometric1F1[1/2,5/4,-160],20]"},
        {"z=-256", "N[hypergeometric1F1[1/2,5/4,-256],20]"},
        {"z=-512", "N[hypergeometric1F1[1/2,5/4,-512],20]"},
        {"z=256+I", "N[hypergeometric1F1[1/2,5/4,256+I],20]"},
        {"z=-256+I", "N[hypergeometric1F1[1/2,5/4,-256+I],20]"},
    }});
    sweeps.push_back({"2F1 near unit circle", {
        {"z=1/2", "N[hypergeometric2F1[1/2,1/3,5/4,1/2],20]"},
        {"z=3/4", "N[hypergeometric2F1[1/2,1/3,5/4,3/4],20]"},
        {"z=4/5", "N[hypergeometric2F1[1/2,1/3,5/4,4/5],20]"},
        {"z=17/20", "N[hypergeometric2F1[1/2,1/3,5/4,17/20],20]"},
        {"z=9/10", "N[hypergeometric2F1[1/2,1/3,5/4,9/10],20]"},
        {"z=19/20", "N[hypergeometric2F1[1/2,1/3,5/4,19/20],20]"},
        {"z=49/50", "N[hypergeometric2F1[1/2,1/3,5/4,49/50],20]"},
        {"z=99/100", "N[hypergeometric2F1[1/2,1/3,5/4,99/100],20]"},
        {"z=1 Gauss", "N[hypergeometric2F1[1/2,1/3,5/4,1],20]"},
    }});
    sweeps.push_back({"polylog near unit circle", {
        {"z=1/2", "N[polylog[2,1/2],20]"},
        {"z=4/5", "N[polylog[2,4/5],20]"},
        {"z=9/10", "N[polylog[2,9/10],20]"},
        {"z=19/20", "N[polylog[2,19/20],20]"},
        {"z=49/50", "N[polylog[2,49/50],20]"},
        {"z=99/100", "N[polylog[2,99/100],20]"},
        {"z=199/200", "N[polylog[2,199/200],20]"},
        {"z=999/1000", "N[polylog[2,999/1000],20]"},
    }});
    sweeps.push_back(precisionSweep(
        "polylog z=99/100 precision", "N[polylog[2,99/100],{p}]"));
    sweeps.push_back(precisionSweep(
        "polylog z=999/1000 precision", "N[polylog[2,999/1000],{p}]"));
    sweeps.push_back({"polylog continuation", {
        {"1/2+I/4", "N[polylog[2,1/2+I/4],20]"},
        {"99/100+I/100", "N[polylog[2,99/100+I/100],20]"},
        {"I", "N[polylog[2,I],20]"},
        {"-2", "N[polylog[2,-2],20]"},
        {"2+I", "N[polylog[2,2+I],20]"},
    }});
    sweeps.push_back({"polylog higher-order near unit circle", {
        {"s=2", "N[polylog[2,999/1000],20]"},
        {"s=3", "N[polylog[3,999/1000],20]"},
        {"s=4", "N[polylog[4,999/1000],20]"},
        {"s=8", "N[polylog[8,999/1000],20]"},
    }});
    sweeps.push_back({"Ei real magnitude scaling", {
        {"x=32", "N[Ei[32],20]"},
        {"x=96", "N[Ei[96],20]"},
        {"x=97", "N[Ei[97],20]"},
        {"x=128", "N[Ei[128],20]"},
        {"x=256", "N[Ei[256],20]"},
        {"x=512", "N[Ei[512],20]"},
        {"x=1000", "N[Ei[1000],20]"},
    }});
    for (const std::string function : {"Si", "Ci"}) {
        sweeps.push_back({function + " real magnitude scaling", {
            {"x=32", "N[" + function + "[32],20]"},
            {"x=96", "N[" + function + "[96],20]"},
            {"x=97", "N[" + function + "[97],20]"},
            {"x=128", "N[" + function + "[128],20]"},
            {"x=256", "N[" + function + "[256],20]"},
            {"x=512", "N[" + function + "[512],20]"},
            {"x=1000", "N[" + function + "[1000],20]"},
            {"x=10000", "N[" + function + "[10000],20]"},
        }});
    }
    sweeps.push_back({"negative Ei magnitude scaling", {
        {"x=-16", "N[Ei[-16],20]"},
        {"x=-32", "N[Ei[-32],20]"},
        {"x=-64", "N[Ei[-64],20]"},
        {"x=-96", "N[Ei[-96],20]"},
        {"x=-256", "N[Ei[-256],20]"},
        {"x=-512", "N[Ei[-512],20]"},
    }});
    sweeps.push_back(precisionSweep("Ei x=256 precision", "N[Ei[256],{p}]"));
    sweeps.push_back(precisionSweep("Ei x=1000 precision", "N[Ei[1000],{p}]"));
    sweeps.push_back(precisionSweep("negative Ei x=-64 precision", "N[Ei[-64],{p}]"));
    sweeps.push_back(precisionSweep("Si x=512 precision", "N[Si[512],{p}]"));
    sweeps.push_back(precisionSweep("Ci x=512 precision", "N[Ci[512],{p}]"));
    sweeps.push_back({"complex Ei former magnitude boundary", {
        {"128+I", "N[Ei[128+I],20]"},
        {"512+I", "N[Ei[512+I],20]"},
        {"513+I", "N[Ei[513+I],20]"},
        {"1000+I normalized", "N[Ei[1000+I]/exp[1000+I],20]"},
        {"513I", "N[Ei[513I],20]"},
        {"-1000+I", "N[Ei[-1000+I],20]"},
    }});
    sweeps.push_back(precisionSweep(
        "complex Ei z=1000+I normalized precision",
        "N[Ei[1000+I]/exp[1000+I],{p}]"));
    sweeps.push_back({"complex Ci former magnitude boundary", {
        {"120+I", "N[Ci[120+I],20]"},
        {"128+I", "N[Ci[128+I],20]"},
        {"129+I", "N[Ci[129+I],20]"},
        {"140+I", "N[Ci[140+I],20]"},
        {"1000+I", "N[Ci[1000+I],20]"},
        {"513I", "N[Ci[513I],20]"},
        {"-1000+I", "N[Ci[-1000+I],20]"},
    }});
    sweeps.push_back(precisionSweep(
        "complex Ci z=1000+I precision", "N[Ci[1000+I],{p}]"));
    sweeps.push_back({"ellipticF parameter boundary", {
        {"m=1/3", "N[ellipticF[1/2,1/3],20]"},
        {"m=3/4", "N[ellipticF[1/2,3/4],20]"},
        {"m=4/5", "N[ellipticF[1/2,4/5],20]"},
        {"m=17/20", "N[ellipticF[1/2,17/20],20]"},
        {"m=9/10", "N[ellipticF[1/2,9/10],20]"},
        {"m=19/20", "N[ellipticF[1/2,19/20],20]"},
        {"m=49/50", "N[ellipticF[1/2,49/50],20]"},
        {"m=99/100", "N[ellipticF[1/2,99/100],20]"},
        {"m=999/1000", "N[ellipticF[1/2,999/1000],20]"},
    }});
    sweeps.push_back({"ellipticE parameter boundary", {
        {"m=1/3", "N[ellipticE[1/2,1/3],20]"},
        {"m=3/4", "N[ellipticE[1/2,3/4],20]"},
        {"m=4/5", "N[ellipticE[1/2,4/5],20]"},
        {"m=17/20", "N[ellipticE[1/2,17/20],20]"},
        {"m=9/10", "N[ellipticE[1/2,9/10],20]"},
        {"m=19/20", "N[ellipticE[1/2,19/20],20]"},
        {"m=49/50", "N[ellipticE[1/2,49/50],20]"},
        {"m=99/100", "N[ellipticE[1/2,99/100],20]"},
        {"m=999/1000", "N[ellipticE[1/2,999/1000],20]"},
    }});
    sweeps.push_back({"ellipticPi characteristic boundary", {
        {"n=1/5", "N[ellipticPi[1/5,1/2,1/3],20]"},
        {"n=3/4", "N[ellipticPi[3/4,1/2,1/3],20]"},
        {"n=4/5", "N[ellipticPi[4/5,1/2,1/3],20]"},
        {"n=17/20", "N[ellipticPi[17/20,1/2,1/3],20]"},
        {"n=9/10", "N[ellipticPi[9/10,1/2,1/3],20]"},
        {"n=19/20", "N[ellipticPi[19/20,1/2,1/3],20]"},
        {"n=49/50", "N[ellipticPi[49/50,1/2,1/3],20]"},
        {"n=99/100", "N[ellipticPi[99/100,1/2,1/3],20]"},
        {"n=999/1000", "N[ellipticPi[999/1000,1/2,1/3],20]"},
        {"m=n=99/100", "N[ellipticPi[99/100,1/2,99/100],20]"},
    }});
    sweeps.push_back({"ellipticF amplitude at m=9/10", {
        {"phi=1/2", "N[ellipticF[1/2,9/10],20]"},
        {"phi=1", "N[ellipticF[1,9/10],20]"},
        {"phi=4/3", "N[ellipticF[4/3,9/10],20]"},
        {"phi=3/2", "N[ellipticF[3/2,9/10],20]"},
        {"phi=2", "N[ellipticF[2,9/10],20]"},
    }});
    sweeps.push_back({"ellipticPi amplitude at n=9/10", {
        {"phi=1/2", "N[ellipticPi[9/10,1/2,1/3],20]"},
        {"phi=1", "N[ellipticPi[9/10,1,1/3],20]"},
        {"phi=4/3", "N[ellipticPi[9/10,4/3,1/3],20]"},
        {"phi=3/2", "N[ellipticPi[9/10,3/2,1/3],20]"},
        {"phi=2", "N[ellipticPi[9/10,2,1/3],20]"},
    }});
    sweeps.push_back(precisionSweep(
        "ellipticF m=9/10 precision", "N[ellipticF[1/2,9/10],{p}]"));
    sweeps.push_back(precisionSweep(
        "ellipticE m=9/10 precision", "N[ellipticE[1/2,9/10],{p}]"));
    sweeps.push_back(precisionSweep(
        "ellipticPi n=9/10 precision", "N[ellipticPi[9/10,1/2,1/3],{p}]"));
    sweeps.push_back(precisionSweep(
        "ellipticF near-one small-amplitude precision", "N[ellipticF[1/2,99/100],{p}]"));
    sweeps.push_back(precisionSweep(
        "ellipticE near-one small-amplitude precision", "N[ellipticE[1/2,99/100],{p}]"));
    sweeps.push_back(precisionSweep(
        "ellipticPi near-one small-amplitude precision", "N[ellipticPi[99/100,1/2,99/100],{p}]"));
    sweeps.push_back(precisionSweep(
        "ellipticF Carlson near-complete-amplitude precision", "N[ellipticF[3/2,99/100],{p}]"));
    sweeps.push_back(precisionSweep(
        "ellipticE Carlson near-complete-amplitude precision", "N[ellipticE[3/2,99/100],{p}]"));
    sweeps.push_back(precisionSweep(
        "ellipticPi Carlson near-pole precision", "N[ellipticPi[99/100,3/2,99/100],{p}]"));
    sweeps.push_back({"real zeta near pole", {
        {"s=3", "N[zeta[3],20]"},
        {"s=12/5", "N[zeta[12/5],20]"},
        {"s=9/5", "N[zeta[9/5],20]"},
        {"s=3/2", "N[zeta[3/2],20]"},
        {"s=6/5", "N[zeta[6/5],20]"},
        {"s=11/10", "N[zeta[11/10],20]"},
        {"s=101/100", "N[zeta[101/100],20]"},
    }});
    sweeps.push_back({"complex zeta near backend boundary", {
        {"3+I", "N[zeta[3+I],20]"},
        {"2+I", "N[zeta[2+I],20]"},
        {"3/2+I", "N[zeta[3/2+I],20]"},
        {"6/5+I", "N[zeta[6/5+I],20]"},
        {"11/10+I", "N[zeta[11/10+I],20]"},
    }});
    for (const std::string function : {"digamma", "trigamma"}) {
        sweeps.push_back({function + " near zero", {
            {"x=1", "N[" + function + "[1],20]"},
            {"x=1/10", "N[" + function + "[1/10],20]"},
            {"x=1/100", "N[" + function + "[1/100],20]"},
            {"x=1/1000", "N[" + function + "[1/1000],20]"},
        }});
    }
    sweeps.push_back({"ibeta x sweep", {
        {"x=1/4", "N[ibeta[1/3,2/3,1/4],20]"},
        {"x=1/2", "N[ibeta[1/3,2/3,1/2],20]"},
        {"x=3/4", "N[ibeta[1/3,2/3,3/4],20]"},
        {"x=9/10", "N[ibeta[1/3,2/3,9/10],20]"},
        {"x=99/100", "N[ibeta[1/3,2/3,99/100],20]"},
    }});
    sweeps.push_back({"complex Fresnel argument sweep", {
        {"1+I", "N[fresnelc[1+I],20]"},
        {"4+I", "N[fresnelc[4+I],20]"},
        {"7+I", "N[fresnelc[7+I],20]"},
        {"8+I", "N[fresnelc[8+I],20]"},
        {"20+I", "N[fresnelc[20+I],20]"},
        {"32+I", "N[fresnelc[32+I],20]"},
        {"1+32I", "N[fresnelc[1+32I],20]"},
        {"8+8I", "N[fresnelc[8+8I],20]"},
    }});
    sweeps.push_back({"Lambert W branch-point argument sweep", {
        {"W0 +i1e-2", "N[lambertw[-1/E+I/10^2],20]"},
        {"W0 +i1e-4", "N[lambertw[-1/E+I/10^4],20]"},
        {"W0 +i1e-8", "N[lambertw[-1/E+I/10^8],20]"},
        {"W0 -i1e-8", "N[lambertw[-1/E-I/10^8],20]"},
        {"W-1 +i1e-8", "N[lambertw[-1,-1/E+I/10^8],20]"},
        {"W1 -i1e-8", "N[lambertw[1,-1/E-I/10^8],20]"},
        {"W0 right 1e-12", "N[lambertw[-1/E+1/10^12],20]"},
        {"W-1 right 1e-12", "N[lambertw[-1,-1/E+1/10^12],20]"},
    }});
    sweeps.push_back({"Lambert W branch-point precision sweep", {
        {"20 digits", "N[lambertw[-1/E+I/10^12],20]"},
        {"40 digits", "N[lambertw[-1/E+I/10^12],40]"},
        {"80 digits", "N[lambertw[-1/E+I/10^12],80]"},
        {"100 digits", "N[lambertw[-1/E+I/10^12],100]"},
    }, false});
    sweeps.push_back({"Lambert W shallow branch-point precision", {
        {"20 digits", "N[lambertw[-1/E+I/10^2],20]"},
        {"40 digits", "N[lambertw[-1/E+I/10^2],40]"},
        {"80 digits", "N[lambertw[-1/E+I/10^2],80]"},
        {"100 digits", "N[lambertw[-1/E+I/10^2],100]"},
    }, false});
    sweeps.push_back({"Lambert W extreme branch-point offset", {
        {"W0 +i1e-160 20d", "N[lambertw[-1/E+I/10^160],20]"},
        {"W0 -i1e-160 20d", "N[lambertw[-1/E-I/10^160],20]"},
        {"W-1 +i1e-160 50d", "N[lambertw[-1,-1/E+I/10^160],50]"},
        {"W1 -i1e-160 50d", "N[lambertw[1,-1/E-I/10^160],50]"},
        {"W0 -i1e-160 100d", "N[lambertw[-1/E-I/10^160],100]"},
    }, false});
    sweeps.push_back({"real backend snapshot", {
        {"lgamma", "N[lgamma[1/3],20]"},
        {"erf", "N[erf[1],20]"},
        {"erfc", "N[erfc[1],20]"},
        {"lambertw", "N[lambertw[1],20]"},
        {"W-1", "N[lambertw[-1,-1/10],20]"},
        {"fresnelc", "N[fresnelc[1],20]"},
        {"fresnels", "N[fresnels[1],20]"},
        {"li+", "N[li[2],20]"},
        {"li-", "N[li[-2],20]"},
        {"beta", "N[beta[1/3,2/3],20]"},
        {"betaln", "N[betaln[1/3,2/3],20]"},
    }, false});
    sweeps.push_back({"complex backend snapshot", {
        {"gamma", "N[gamma[1+I],20]"},
        {"erf", "N[erf[1+I],20]"},
        {"Ei", "N[Ei[1+I],20]"},
        {"Si", "N[Si[1+I],20]"},
        {"Ci", "N[Ci[1+I],20]"},
        {"1F1", "N[hypergeometric1F1[1/2,5/4,1+I],20]"},
        {"2F1", "N[hypergeometric2F1[1/2,1/3,5/4,1/2+I/4],20]"},
        {"2F1 cont", "N[hypergeometric2F1[3.4,5.6,4+I,4.6+2I],20]"},
        {"polylog", "N[polylog[2,1/2+I/4],20]"},
        {"li", "N[li[-2+I],20]"},
    }, false});
    return sweeps;
}

struct AlgebraicRootCase final {
    std::string name;
    std::vector<Rational> polynomial;
};

[[nodiscard]] std::vector<Rational> binomialRootPolynomial(std::size_t degree) {
    std::vector<Rational> polynomial(degree + 1);
    polynomial.front() = Rational{BigInt{-2}};
    polynomial.back() = Rational{BigInt{1}};
    return polynomial;
}

[[nodiscard]] std::vector<Rational> scaleSeparatedRootPolynomial() {
    // (x^8-2^80)(x^8-2^-80): 根半径が2^-10と2^10の2群に分かれる。
    // 単一Cauchy半径seedとNewton-polygon multi-radius seedの差を監視する。
    const BigInt scale = BigInt{1} << 80;
    std::vector<Rational> polynomial(17);
    polynomial.front() = Rational{BigInt{1}};
    polynomial[8] = -(Rational{scale} + Rational{BigInt{1}, scale});
    polynomial.back() = Rational{BigInt{1}};
    return polynomial;
}

[[nodiscard]] std::vector<AlgebraicRootCase> algebraicRootCases() {
    std::vector<AlgebraicRootCase> cases;
    for (const std::size_t degree : {2U, 4U, 5U, 8U, 10U, 12U, 16U})
        cases.push_back({"x^" + std::to_string(degree) + "-2", binomialRootPolynomial(degree)});
    cases.push_back({"dense degree 5", {
        Rational{BigInt{1}}, Rational{BigInt{-1}}, Rational{BigInt{2}},
        Rational{BigInt{3}}, Rational{BigInt{-2}}, Rational{BigInt{1}}}});
    // (x^2+1)(x^2+1001/1000): 2組の純虚根が近いcertification stress case。
    cases.push_back({"clustered degree 4", {
        Rational{BigInt{1001}, BigInt{1000}}, Rational{},
        Rational{BigInt{2001}, BigInt{1000}}, Rational{}, Rational{BigInt{1}}}});
    cases.push_back({"two-radius degree 16", scaleSeparatedRootPolynomial()});
    return cases;
}

struct AlgebraicRootMeasurement final {
    bool succeeded = false;
    double createMilliseconds = 0.0;
    double refine80Milliseconds = 0.0;
    double refine320Milliseconds = 0.0;
};

struct AlgebraicSolveCase final {
    std::string name;
    std::string expression;
};

[[nodiscard]] const std::vector<AlgebraicSolveCase>& algebraicSolveCases() {
    static const std::vector<AlgebraicSolveCase> cases{
        {"degree 6 modular proof", "solve[x^6-3x^5-x^4+2x^3+2x^2-2x-1==0,x]"},
        {"degree 16 general", "solve[x^16+x+1==0,x]"},
        {"degree 32 sparse", "solve[x^32-x+1==0,x]"},
        {"degree 64 general", "solve[x^64+x+1==0,x]"},
        {"degree 65 cap crossing", "solve[x^65+x+1==0,x]"},
        {"reducible degree 8", "solve[(x^4-2)*(x^4+1)==0,x]"},
    };
    return cases;
}

[[nodiscard]] double measureAlgebraicSolve(
    const AlgebraicSolveCase& testCase,
    std::size_t iterations) {
    mmcal::kernel::KernelSession session;
    double milliseconds = 0.0;
    for (std::size_t iteration = 0; iteration < iterations; ++iteration) {
        session.resetForIndependentEvaluation();
        const auto start = Clock::now();
        const auto value = session.evaluate(testCase.expression);
        const auto end = Clock::now();
        if (mmcal::formatting::formatExpr(value).empty())
            std::abort();
        milliseconds += std::chrono::duration<double, std::milli>(end - start).count();
    }
    return milliseconds / static_cast<double>(iterations);
}

[[nodiscard]] AlgebraicRootMeasurement measureAlgebraicRoot(
    const AlgebraicRootCase& testCase,
    std::size_t iterations) {
    AlgebraicRootMeasurement result;
    for (std::size_t iteration = 0; iteration < iterations; ++iteration) {
        const auto createStart = Clock::now();
        const auto root = mmcal::symbolic::ComplexAlgebraicNumber::create(
            testCase.polynomial, 1);
        const auto createEnd = Clock::now();
        if (!root)
            return result;

        const auto refine80Start = Clock::now();
        (void)root->refined(80);
        const auto refine80End = Clock::now();
        const auto refine320Start = Clock::now();
        (void)root->refined(320);
        const auto refine320End = Clock::now();

        result.createMilliseconds += std::chrono::duration<double, std::milli>(
            createEnd - createStart).count();
        result.refine80Milliseconds += std::chrono::duration<double, std::milli>(
            refine80End - refine80Start).count();
        result.refine320Milliseconds += std::chrono::duration<double, std::milli>(
            refine320End - refine320Start).count();
    }
    const double divisor = static_cast<double>(iterations);
    result.succeeded = true;
    result.createMilliseconds /= divisor;
    result.refine80Milliseconds /= divisor;
    result.refine320Milliseconds /= divisor;
    return result;
}

} // namespace

void runSpecialFunctionPerformanceCliffAudit(
    const PerformanceCliffAuditOptions& options) {
    if (options.iterations == 0)
        throw std::invalid_argument("Performance cliff audit iterations must be positive");
    if (!(options.cliffRatio > 1.0))
        throw std::invalid_argument("Performance cliff ratio must be greater than one");

    std::cout << "special-function performance-cliff audit\n"
              << "iterations=" << options.iterations
              << " cliff-ratio=" << options.cliffRatio << "x\n"
              << "Numeric-to-held/error changes are BOUNDARY, not performance cliffs.\n";

    for (const auto& sweep : specialFunctionSweeps())
        printSweep(sweep, options.iterations, options.cliffRatio);
}

void runAlgebraicRootPerformanceCliffAudit(
    const PerformanceCliffAuditOptions& options) {
    if (options.iterations == 0)
        throw std::invalid_argument("Algebraic root audit iterations must be positive");

    std::cout << "complex algebraic Root performance-cliff audit\n"
              << "iterations=" << options.iterations << '\n'
              << "create = all-root isolation + ordering; refine = selected certified disk only.\n\n"
              << std::left << std::setw(22) << "case"
              << std::right << std::setw(12) << "create-ms"
              << std::setw(13) << "refine80-ms"
              << std::setw(14) << "refine320-ms"
              << "  outcome\n";

    for (const auto& testCase : algebraicRootCases()) {
        const auto measurement = measureAlgebraicRoot(testCase, options.iterations);
        std::cout << std::left << std::setw(22) << testCase.name;
        if (!measurement.succeeded) {
            std::cout << std::right << std::setw(39) << "" << "  unsupported\n";
            continue;
        }
        std::cout << std::right << std::fixed << std::setprecision(3)
                  << std::setw(12) << measurement.createMilliseconds
                  << std::setw(13) << measurement.refine80Milliseconds
                  << std::setw(14) << measurement.refine320Milliseconds
                  << "  numeric\n";
    }

    std::cout << "\npolynomial Solve canonicalization\n"
              << std::left << std::setw(26) << "case"
              << std::right << std::setw(12) << "solve-ms" << '\n';
    for (const auto& testCase : algebraicSolveCases())
        std::cout << std::left << std::setw(26) << testCase.name
                  << std::right << std::fixed << std::setprecision(3)
                  << std::setw(12) << measureAlgebraicSolve(testCase, options.iterations)
                  << '\n';
}

} // namespace mmcal::benchmarks
