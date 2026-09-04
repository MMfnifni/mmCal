// 全回帰テストのエントリポイント
#include "approximation/approximation_tests.hpp"
#include "approximation/certified_constants_tests.hpp"
#include "approximation/certified_complex_transcendental_tests.hpp"
#include "approximation/certified_sqrt_tests.hpp"
#include "approximation/certified_transcendental_tests.hpp"
#include "approximation/real_interval_tests.hpp"
#include "error/error_message_tests.hpp"
#include "builtins/linear_algebra_tests.hpp"
#include "builtins/numerical_calculus_tests.hpp"
#include "builtins/aggregate_array_elementary_tests.hpp"
#include "builtins/statistics_tests.hpp"
#include "builtins/combinatorics_signal_tests.hpp"
#include "builtins/special_function_extension_tests.hpp"
#include "builtins/random_function_tests.hpp"
#include "builtins/history_diagnostic_tests.hpp"
#include "builtins/session_approximation_tests.hpp"
#include "builtins/integration_tests.hpp"
#include "builtins/advanced_integration_tests.hpp"
#include "builtins/calculus_knowledge_tests.hpp"
#include "evaluation/builtin_registry_tests.hpp"
#include "evaluation/environment_tests.hpp"
#include "evaluation/user_function_registry_tests.hpp"
#include "evaluation/evaluator_tests.hpp"
#include "expression/expr_tests.hpp"
#include "numeric/big_int_tests.hpp"
#include "numeric/big_float_tests.hpp"
#include "numeric/decimal_approximation_tests.hpp"
#include "numeric/big_uint_tests.hpp"
#include "numeric/integer_algorithms_tests.hpp"
#include "numeric/number_tests.hpp"
#include "numeric/rational_tests.hpp"
#include "numeric/real_number_tests.hpp"
#include "kernel/kernel_session_tests.hpp"
#include "mathematics/math_registry_tests.hpp"
#include "mathematics/definedness_tests.hpp"
#include "mathematics/value_facts_tests.hpp"
#include "mathematics/knowledge_context_tests.hpp"
#include "solver/solution_set_tests.hpp"
#include "simplification/simplifier_tests.hpp"
#include "mathematics/exact_algebra_tests.hpp"
#include "mathematics/exact_trigonometry_tests.hpp"
#include "mathematics/exact_transcendental_tests.hpp"
#include "syntax/lexer_tests.hpp"
#include "syntax/lowerer_tests.hpp"
#include "syntax/parser_tests.hpp"
#include "symbols/symbol_registry_tests.hpp"
#include "symbols/symbol_table_tests.hpp"
#include "symbolic/differentiation_tests.hpp"
#include "symbolic/series_tests.hpp"
#include "test_framework.hpp"

#include <chrono>
#include <exception>
#include <iostream>
#include <string_view>

int main(int argc, char* argv[]) {
    try {
        bool timings = false;
        for (int i = 1; i < argc; ++i) {
            const std::string_view argument = argv[i];
            if (argument == "--timings") {
                timings = true;
                continue;
            }
            std::cerr << "Unknown argument: " << argument << '\n';
            return 2;
        }

        mmcal::tests::TestRunner tests;
        const auto run = [&](std::string_view name, auto&& function) {
            const auto started = std::chrono::steady_clock::now();
            function();
            if (timings) {
                const auto elapsed = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - started).count();
                std::cout << "[TIMING] " << name << ": " << elapsed << " ms\n";
            }
        };
        run("Approximation", [&] { mmcal::tests::runApproximationTests(tests); });
        run("RealInterval", [&] { mmcal::tests::runRealIntervalTests(tests); });
        run("CertifiedConstant", [&] { mmcal::tests::runCertifiedConstantTests(tests); });
        run("CertifiedSqrt", [&] { mmcal::tests::runCertifiedSqrtTests(tests); });
        run("CertifiedTranscendental", [&] { mmcal::tests::runCertifiedTranscendentalTests(tests); });
        run("CertifiedComplexTranscendental", [&] { mmcal::tests::runCertifiedComplexTranscendentalTests(tests); });
        run("BigUInt", [&] { mmcal::tests::runBigUIntTests(tests); });
        run("BigInt", [&] { mmcal::tests::runBigIntTests(tests); });
        run("BigFloat", [&] { mmcal::tests::runBigFloatTests(tests); });
        run("IntegerAlgorithm", [&] { mmcal::tests::runIntegerAlgorithmTests(tests); });
        run("DecimalApproximation", [&] { mmcal::tests::runDecimalApproximationTests(tests); });
        run("Rational", [&] { mmcal::tests::runRationalTests(tests); });
        run("RealNumber", [&] { mmcal::tests::runRealNumberTests(tests); });
        run("Number", [&] { mmcal::tests::runNumberTests(tests); });
        run("Expr", [&] { mmcal::tests::runExprTests(tests); });
        run("LinearAlgebra", [&] { mmcal::tests::runLinearAlgebraTests(tests); });
        run("NumericalCalculus", [&] { mmcal::tests::runNumericalCalculusTests(tests); });
        run("AggregateArrayElementary", [&] { mmcal::tests::runAggregateArrayElementaryTests(tests); });
        run("Statistics", [&] { mmcal::tests::runStatisticsTests(tests); });
        run("CombinatoricsSignal", [&] { mmcal::tests::runCombinatoricsSignalTests(tests); });
        run("SpecialFunctionExtension", [&] { mmcal::tests::runSpecialFunctionExtensionTests(tests); });
        run("RandomFunction", [&] { mmcal::tests::runRandomFunctionTests(tests); });
        run("HistoryDiagnostic", [&] { mmcal::tests::runHistoryDiagnosticTests(tests); });
        run("SessionApproximation", [&] { mmcal::tests::runSessionApproximationTests(tests); });
        run("Integration", [&] { mmcal::tests::runIntegrationTests(tests); });
        run("AdvancedIntegration", [&] { mmcal::tests::runAdvancedIntegrationTests(tests); });
        run("CalculusKnowledge", [&] { mmcal::tests::runCalculusKnowledgeTests(tests); });
        run("BuiltinRegistry", [&] { mmcal::tests::runBuiltinRegistryTests(tests); });
        run("SymbolTable", [&] { mmcal::tests::runSymbolTableTests(tests); });
        run("SymbolRegistry", [&] { mmcal::tests::runSymbolRegistryTests(tests); });
        run("MathRegistry", [&] { mmcal::tests::runMathRegistryTests(tests); });
        run("Definedness", [&] { mmcal::tests::runDefinednessTests(tests); });
        run("ValueFacts", [&] { mmcal::tests::runValueFactsTests(tests); });
        run("KnowledgeContext", [&] { mmcal::tests::runKnowledgeContextTests(tests); });
        run("SolutionSet", [&] { mmcal::tests::runSolutionSetTests(tests); });
        run("Simplifier", [&] { mmcal::tests::runSimplifierTests(tests); });
        run("ExactAlgebra", [&] { mmcal::tests::runExactAlgebraTests(tests); });
        run("ExactTrigonometry", [&] { mmcal::tests::runExactTrigonometryTests(tests); });
        run("ExactTranscendental", [&] { mmcal::tests::runExactTranscendentalTests(tests); });
        run("Differentiation", [&] { mmcal::tests::runDifferentiationTests(tests); });
        run("Series", [&] { mmcal::tests::runSeriesTests(tests); });
        run("Environment", [&] { mmcal::tests::runEnvironmentTests(tests); });
        run("UserFunctionRegistry", [&] { mmcal::tests::runUserFunctionRegistryTests(tests); });
        run("Evaluator", [&] { mmcal::tests::runEvaluatorTests(tests); });
        run("ErrorMessage", [&] { mmcal::tests::runErrorMessageTests(tests); });
        run("Lexer", [&] { mmcal::tests::runLexerTests(tests); });
        run("Parser", [&] { mmcal::tests::runParserTests(tests); });
        run("Lowerer", [&] { mmcal::tests::runLowererTests(tests); });
        run("KernelSession", [&] { mmcal::tests::runKernelSessionTests(tests); });
        return tests.result();
    }
    catch (const std::exception& error) {
        std::cerr << "Fatal error while running tests: " << error.what() << '\n';
        return 2;
    }
}
