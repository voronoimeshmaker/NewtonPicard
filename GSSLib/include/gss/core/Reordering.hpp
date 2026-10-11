#pragma once

#include <cstdint>

#include "gss/core/Types.hpp"

namespace gss {
namespace reorder {

/// @brief Number of quantization bits used by space-filling curves.
/// A side of 2^12 = 4096 cells per direction covers >16M cells unambiguously.
inline constexpr unsigned kSfcBits = 12;

/// @brief Quantize x in [xmin, xmax] to a Morton/Hilbert grid coordinate.
[[nodiscard]] inline std::uint32_t quantize(Real x, Real xmin, Real span) {
    const std::uint32_t side = std::uint32_t{1} << kSfcBits;
    if (!(span > Real{0})) {
        return 0;
    }
    Real t = (x - xmin) / span;
    t = t < Real{0} ? Real{0} : (t > Real{1} ? Real{1} : t);
    const Real scaled = t * static_cast<Real>(side);
    return static_cast<std::uint32_t>(scaled) & (side - 1);
}

/// @brief Morton (Z-order) key: bit-interleave of x and y.
[[nodiscard]] inline std::uint64_t morton_key(std::uint32_t x, std::uint32_t y) {
    std::uint64_t z = 0;
    for (unsigned b = 0; b < kSfcBits; ++b) {
        z |= static_cast<std::uint64_t>((x >> b) & 1u) << (2 * b);
        z |= static_cast<std::uint64_t>((y >> b) & 1u) << (2 * b + 1);
    }
    return z;
}

/// @brief Hilbert (XY -> d) key on a side x side grid, side a power of two.
/// Standard XY-to-Hilbert-distance algorithm (Wikipedia xy2d / Butz-Skilling).
[[nodiscard]] inline std::uint64_t hilbert_key(std::uint32_t x, std::uint32_t y, std::uint32_t side) {
    std::uint64_t d = 0;
    for (std::uint32_t s = side >> 1; s > 0; s >>= 1) {
        const std::uint32_t rx = (x & s) ? 1u : 0u;
        const std::uint32_t ry = (y & s) ? 1u : 0u;
        d += static_cast<std::uint64_t>(s) * static_cast<std::uint64_t>(s)
            * ((3 * rx) ^ ry);
        if (ry == 0) {
            if (rx == 1) {
                x = (s - 1) - x;
                y = (s - 1) - y;
            }
            x = x ^ y;
            y = x ^ y;
            x = x ^ y;
        }
    }
    return d;
}

} // namespace reorder
} // namespace gss