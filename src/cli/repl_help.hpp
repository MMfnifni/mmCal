#pragma once

#include <iosfwd>
#include <string_view>

namespace mmcal::evaluation {
class BuiltinRegistry;
}

namespace mmcal::cli {

// Returns false only when the line is not a :help command. The caller handles
// it before KernelSession::evaluate(), so help never consumes an In[n] slot.
[[nodiscard]] bool handleReplHelpCommand(
    std::string_view line,
    const evaluation::BuiltinRegistry& registry,
    std::ostream& output);

} // namespace mmcal::cli
