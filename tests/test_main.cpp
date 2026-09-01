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

#include <exception>
#include <iostream>

int main() {
    try {
        mmcal::tests::TestRunner tests;
        mmcal::tests::runApproximationTests(tests);
        mmcal::tests::runRealIntervalTests(tests);
        mmcal::tests::runCertifiedConstantTests(tests);
        mmcal::tests::runCertifiedSqrtTests(tests);
        mmcal::tests::runCertifiedTranscendentalTests(tests);
        mmcal::tests::runCertifiedComplexTranscendentalTests(tests);
        mmcal::tests::runBigUIntTests(tests);
        mmcal::tests::runBigIntTests(tests);
        mmcal::tests::runBigFloatTests(tests);
        mmcal::tests::runIntegerAlgorithmTests(tests);
        mmcal::tests::runDecimalApproximationTests(tests);
        mmcal::tests::runRationalTests(tests);
        mmcal::tests::runRealNumberTests(tests);
        mmcal::tests::runNumberTests(tests);
        mmcal::tests::runExprTests(tests);
        mmcal::tests::runLinearAlgebraTests(tests);
        mmcal::tests::runNumericalCalculusTests(tests);
        mmcal::tests::runAggregateArrayElementaryTests(tests);
        mmcal::tests::runStatisticsTests(tests);
        mmcal::tests::runCombinatoricsSignalTests(tests);
        mmcal::tests::runSpecialFunctionExtensionTests(tests);
        mmcal::tests::runRandomFunctionTests(tests);
        mmcal::tests::runHistoryDiagnosticTests(tests);
        mmcal::tests::runSessionApproximationTests(tests);
        mmcal::tests::runIntegrationTests(tests);
        mmcal::tests::runAdvancedIntegrationTests(tests);
        mmcal::tests::runCalculusKnowledgeTests(tests);
        mmcal::tests::runBuiltinRegistryTests(tests);
        mmcal::tests::runSymbolTableTests(tests);
        mmcal::tests::runSymbolRegistryTests(tests);
        mmcal::tests::runMathRegistryTests(tests);
        mmcal::tests::runDefinednessTests(tests);
        mmcal::tests::runValueFactsTests(tests);
        mmcal::tests::runKnowledgeContextTests(tests);
        mmcal::tests::runSolutionSetTests(tests);
        mmcal::tests::runSimplifierTests(tests);
        mmcal::tests::runExactAlgebraTests(tests);
        mmcal::tests::runExactTrigonometryTests(tests);
        mmcal::tests::runExactTranscendentalTests(tests);
        mmcal::tests::runDifferentiationTests(tests);
        mmcal::tests::runSeriesTests(tests);
        mmcal::tests::runEnvironmentTests(tests);
        mmcal::tests::runUserFunctionRegistryTests(tests);
        mmcal::tests::runEvaluatorTests(tests);
        mmcal::tests::runErrorMessageTests(tests);
        mmcal::tests::runLexerTests(tests);
        mmcal::tests::runParserTests(tests);
        mmcal::tests::runLowererTests(tests);
        mmcal::tests::runKernelSessionTests(tests);
        return tests.result();
    }
    catch (const std::exception& error) {
        std::cerr << "Fatal error while running tests: " << error.what() << '\n';
        return 2;
    }
}
