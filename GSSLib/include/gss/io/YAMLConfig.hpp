#pragma once

#include <filesystem>

#include "gss/core/Config.hpp"

namespace gss {

/// @brief YAML configuration reader.
class YAMLConfig {
public:
    /// @brief Load and validate a full Config from disk.
    [[nodiscard]] static Config load(const std::filesystem::path& path);

    /// @brief Write a Config to YAML.
    static void write(const Config& config, const std::filesystem::path& path);
};

} // namespace gss
