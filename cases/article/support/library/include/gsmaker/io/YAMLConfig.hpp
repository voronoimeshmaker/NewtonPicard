#pragma once

#include <filesystem>

#include "gsmaker/core/Config.hpp"

namespace gss {

/// @brief YAML serialisation of configuration owned by GSMakerLib.
///
/// Case-specific YAML sections are ignored during load. They are handled
/// separately through ConfigTraits and ConfigFactory.
class YAMLConfig {
public:
    [[nodiscard]]
    static Config load(const std::filesystem::path& path);

    static void write(
        const Config& config,
        const std::filesystem::path& path);
};

} // namespace gss