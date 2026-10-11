#pragma once

#include <concepts>
#include <filesystem>

#include "gsmaker/core/Config.hpp"
#include "gsmaker/core/ConfigTraits.hpp"
#include "gsmaker/io/YAMLConfig.hpp"

namespace gss {

template<typename CaseTag>
concept ConfiguredCase =
    requires(const std::filesystem::path& path) {
        typename ConfigTraits<CaseTag>::SpecificConfig;

        {
            ConfigTraits<CaseTag>::load(path)
        } -> std::same_as<
            typename ConfigTraits<CaseTag>::SpecificConfig>;
    };


/// @brief Complete configuration seen by one program.
///
/// The common portion belongs to GSMakerLib.
/// The specific portion belongs to the scientific case/program.
template<ConfiguredCase CaseTag>
struct ProgramConfig {
    CommonConfig common{};
    typename ConfigTraits<CaseTag>::SpecificConfig specific{};
};


/// @brief Construct program configuration from its local configuration file.
///
/// The same local file may contain both the common GSMaker sections and
/// arbitrary case-specific sections. YAMLConfig reads only the common
/// portion; ConfigTraits<CaseTag> reads the case-specific portion.
template<ConfiguredCase CaseTag>
class ConfigFactory {
public:
    using SpecificConfig =
        typename ConfigTraits<CaseTag>::SpecificConfig;

    using Result = ProgramConfig<CaseTag>;

    [[nodiscard]]
    static Result load(const std::filesystem::path& path) {
        return Result{
            .common = YAMLConfig::load(path),
            .specific = ConfigTraits<CaseTag>::load(path)
        };
    }
};

} // namespace gss