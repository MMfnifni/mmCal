#include "core.hpp"
#include "evaluate.hpp"
#include "formatter.hpp"
#include "functions/functions.hpp"
#include "lexer_parser.hpp"
#include "repl_command.hpp"
#include "session_ops.hpp"

#include <iostream>
#include <string>

using namespace mm::cal;

/* ============================
  main
  ============================ */

int main(int argc, char *argv[]) {
 int base = 0;
 SystemConfig syscfg;
 std::vector<InputEntry> history;
 std::unordered_map<std::string, Value> globals;
 std::unordered_map<std::string, UserFunction> userFunctions;

 SessionState session{syscfg, history, base, globals, userFunctions};
 EvaluationContext ectx{session};

 initFunctions(syscfg);

 // --batchモード(tester用)
 if (argc > 1 && std::string(argv[1]) == "--batch") {
  std::cout << std::unitbuf;

  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  std::string line;

  while (std::getline(std::cin, line)) {
   if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

   try {
    EvalResult res = evalLine(line, syscfg, history, ectx);
    applySideEffects(ectx, res);

    if (res.explain != "") { std::cout << res.explain; }

    if (res.suppressDisplay) {
     std::cout << "evaluate success. display suppressed in silent mode" << '\n' << std::flush;
     continue;
    }

    if (res.hasDisplayOverride()) {
     std::cout << res.displayOverride << '\n' << std::flush;
    } else {
     std::cout << formatResult(res.value, syscfg) << '\n' << std::flush;
    }
   }

   catch (ControlRequest &e) {
    switch (e.kind) {
     case ControlRequest::Kind::Exit: return 0;

     case ControlRequest::Kind::Clear:
      resetSession(ectx, history);
      std::cout << "[cleared]" << '\n' << std::flush;
      break;
    }
   }

   catch (const CalcError &e) {
    std::cout << "Error: " << e.what() << '\n';
   }

   catch (const std::exception &e) {
    std::cout << "Error: " << e.what() << '\n';
   }

   catch (...) {
    std::cout << "Error: Unknown exception\n";
   }
  }

  return 0;
 }

 // ================================
 // CLI モード
 // ================================
 if (argc >= 2) {
  std::string line = argv[1];

  try {
   EvalResult res = evalLine(line, syscfg, history, ectx);
   applySideEffects(ectx, res);

   if (res.explain != "") { std::cout << res.explain; }

   if (res.suppressDisplay) {
    std::cout << "evaluate success. suppress display by silent" << '\n' << std::flush;
   } else if (res.hasDisplayOverride()) {
    std::cout << res.displayOverride << '\n' << std::flush;
   } else {
    std::cout << formatResult(res.value, syscfg) << '\n' << std::flush;
   }

  } catch (ControlRequest &e) {
   switch (e.kind) {
    case ControlRequest::Kind::Exit: return 0;

    case ControlRequest::Kind::Clear:
     resetSession(ectx, history);
     std::cout << "[cleared]\n";
     return 0;
   }

  } catch (const CalcError &e) {
   std::cout << "Error: " << e.what() << "\n" << line << "\n" << std::string(e.pos, ' ') << "^\n";
   return 1;

  } catch (const std::exception &e) {
   std::cout << "Error: " << e.what() << "\n";
   return 1;

  } catch (...) {
   std::cout << "Error: Unknown exception\n";
   return 1;
  }

  return 0;
 }

 // ================================
 // REPL モード
 // ================================
 std::cout << "================================\n"
              "  mm Calculator\n"
              "        mmKreutzef 2021-2026\n"
              "================================\n\n";

 std::string line;
 std::cout << std::unitbuf;
 while (true) {
  std::cout << "In [" << history.size() + 1 << "] := " << std::flush;
  if (!std::getline(std::cin, line)) break;

  if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

  try {
   EvalResult res = evalLine(line, syscfg, history, ectx);
   applySideEffects(ectx, res);
   history.push_back({line, res.value});

   if (res.explain != "") { std::cout << res.explain; }

   if (res.suppressDisplay) {
    std::cout << "Out[" << history.size() << "] := "
              << "evaluate success. suppress display by silent" << '\n'
              << std::flush;
    continue;
   }

   if (res.hasDisplayOverride()) {
    std::cout << "Out[" << history.size() << "] := " << res.displayOverride << "\n\n";
   } else {
    std::cout << "Out[" << history.size() << "] := " << formatResultMulti(res.value, syscfg) << "\n\n";
   }

  } catch (ControlRequest &e) {
   switch (e.kind) {
    case ControlRequest::Kind::Exit: std::cout << "bye...nara\n"; return 0;

    case ControlRequest::Kind::Clear:
     resetSession(ectx, history);
     std::cout << "variables and histories cleared\n";
     break;
   }

  } catch (const CalcError &e) { std::cout << "\nError: " << e.what() << "\n" << line << "\n" << std::string(e.pos, ' ') << "^\n\n"; } catch (const std::exception &e) {
   std::cout << "\nError: " << e.what() << "\n\n";

  } catch (...) { std::cout << "\nError: Unknown exception\n\n"; }
 }

 std::cout << "bye...nara\n";
 return 0;
}