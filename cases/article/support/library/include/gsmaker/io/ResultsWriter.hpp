#pragma once

#include <filesystem>
#include <fstream>
#include <span>
#include <stdexcept>

#include "gss/core/MeshAdapter.hpp"

namespace gss {

/// @brief CSV/VTK result writer for cell-centered fields.
class ResultsWriter {
public:
    /// @brief Write a CSV file containing x, y, psi and source columns.
    static void write_csv(const std::filesystem::path& path, const MeshAdapter& mesh, std::span<const Real> psi) {
        std::ofstream out(path);
        if (!out) {
            throw std::runtime_error("failed to open results file: " + path.string());
        }
        out << "cell_id,x,y,psi,source,kappa\n";
        const auto& topo = mesh.topology();
        const auto& data = mesh.data();
        for (Index i = 0; i < data.n_cells; ++i) {
            const auto s = static_cast<std::size_t>(i);
            out << i << ',' << topo.coord_x[s] << ',' << topo.coord_y[s] << ','
                << psi[s] << ',' << data.source[s] << ',' << data.kappa[s] << '\n';
        }
    }
};

} // namespace gss
