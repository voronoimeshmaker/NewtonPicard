    #pragma once

    #include "gsmaker/core/Types.hpp"
    #include "gsmaker/verification/analytical/AnalyticalSolution.hpp"

    namespace gss {

        
    /**
     * @brief Rotated Fadeev analytical solution used for verification.
     *
     * The solution is
     *
     *   Psi(x,y) = log(alpha cos(x_tilde)
     *                  + sqrt(1 + alpha^2) cosh(y_tilde)),
     *
     * with
     *
     *   x_tilde = x cos(theta) - y sin(theta),
     *   y_tilde = x sin(theta) + y cos(theta).
     *
     * The associated Grad-Shafranov free function is
     *
     *   F(Psi) = -0.5 c exp(-2 Psi),
     *
     * so that dF/dPsi = c exp(-2 Psi), where c is the nonlinearity scale.
     *
     * The scale exists because the Fadeev family is spectrally narrow.  A
     * campaign of 860 draws over alpha and theta moved the Picard spectral
     * radius only between 0.078 and 0.202, and the Liouville scale covariance
     * Psi_L(x) = Psi(x/L) + log L leaves every eigenvalue of the fixed-point
     * operator unchanged, so the size of the domain cannot move it at all.
     * Scaling the free function does move it, and exactly:
     *
     *   B(c) = c B(1),   G(c) = c G(1),   rho(c) = c rho(1).
     *
     * Psi stays the exact solution for every c because the difference is
     * returned to the equation as a manufactured source,
     *
     *   Lap(Psi) = c exp(-2 Psi) + q,   q = (1 - c) exp(-2 Psi*),
     *
     * with Psi* the analytical field above.  At Psi = Psi* the two terms sum
     * to exp(-2 Psi*), which is the original equation, and q vanishes
     * identically at c = 1.  A campaign run at c = 1 is therefore the same
     * arithmetic as one run from before this parameter existed, digit for
     * digit, which is what makes it a usable regression.
     *
     * The analytical derivative is exposed only as an explicitly named
     * reference quantity.  Production source-term evaluation must continue to
     * obtain dF/dPsi numerically unless a verification case explicitly states
     * otherwise.
     *
     * Both F(Psi) and dF/dPsi grow without bound as Psi decreases, and
     * std::exp overflows to infinity for Psi below roughly -354.9.  The
     * evaluations below are therefore saturated so that they are total
     * functions and never return a non-finite value.  Saturation alone is not
     * a convergence safeguard: it keeps the arithmetic finite but does not
     * signal that the iterate is meaningless.  Callers that iterate on Psi
     * must test admissibility separately against psi_lower_bound().
     */
    class Fadeev {
    public:
        static constexpr double default_alpha = 0.225;
        static constexpr double default_theta_degrees = 5.7;
        static constexpr double default_nonlinearity_scale = 1.0;

        Fadeev();

        Fadeev(
            Real alpha,
            Real theta_degrees,
            Real nonlinearity_scale =
                static_cast<Real>(default_nonlinearity_scale));

        [[nodiscard]] Real alpha() const noexcept {
            return alpha_;
        }

        [[nodiscard]] Real theta_degrees() const noexcept {
            return theta_degrees_;
        }

        /**
         * @brief The multiplier c applied to F and to dF/dPsi.
         *
         * A study that builds a manufactured source must carry this value with
         * the field the source was built from.  Mixing a source made at one
         * scale with a field made at another is silently wrong arithmetic
         * rather than an error, so the study-side workspace cross-checks it.
         */
        [[nodiscard]] Real nonlinearity_scale() const noexcept {
            return nonlinearity_scale_;
        }

        /**
         * @brief Evaluate the analytical flux function Psi(x,y).
         */
        [[nodiscard]] Real psi(
            Real x,
            Real y) const noexcept;

        /**
         * @brief Evaluate F(Psi) from a supplied value of Psi.
         *
         * Saturated against exponential overflow; see psi_lower_bound().
         */
        [[nodiscard]] Real free_function(
            Real psi_value) const noexcept;

        /**
         * @brief Analytical dF/dPsi for verification only.
         *
         * This function exists to quantify the error of a numerically
         * reconstructed source term.  It is not the normal source-evaluation
         * path used by the scientific solver.
         *
         * Saturated against exponential overflow; see psi_lower_bound().
         */
        [[nodiscard]] Real dF_dpsi_reference(
            Real psi_value) const noexcept;

        /**
         * @brief Density q of the manufactured source at the exact field.
         *
         * Returns (1 - c) exp(-2 Psi*) at the supplied value of the analytical
         * solution, saturated by the same exponent clamp as free_function().
         * It is a density: a finite-volume caller multiplies it by the cell
         * volume using whatever sign its residual convention requires.  At
         * c = 1 it is identically zero.
         *
         * The argument must be the analytical Psi*, never the current iterate.
         * q is a fixed datum of the problem and not a function of the unknown;
         * a caller that passes the iterate is solving a different equation and
         * will not notice, because the result is still finite.
         */
        [[nodiscard]] Real manufactured_source_density(
            Real exact_psi_value) const noexcept;

        /**
         * @brief Greatest lower bound of Psi over the whole plane.
         *
         * The argument of the logarithm in psi() is
         *
         *   alpha cos(x_tilde) + sqrt(1 + alpha^2) cosh(y_tilde),
         *
         * which attains its minimum sqrt(1 + alpha^2) - |alpha| > 0 at
         * cos(x_tilde) = -1, y_tilde = 0.  Hence
         *
         *   Psi >= log(sqrt(1 + alpha^2) - |alpha|)
         *
         * everywhere on the admissible manifold, and correspondingly
         *
         *   |dF/dPsi| = |c| exp(-2 Psi) <= |c| (sqrt(1 + alpha^2) - |alpha|)^-2.
         *
         * For the default alpha = 0.225 and c = 1 this gives Psi >= log(0.8) =
         * -0.2231 and dF/dPsi <= 1.5625.  The bound on Psi is a property of the
         * analytical field alone and does not depend on c.
         *
         * A nonlinear iterate whose minimum falls appreciably below this bound
         * has left the set on which the problem is posed.  Detecting that
         * condition is several hundred iterations cheaper than waiting for
         * exponential overflow, and it is the test a globalised nonlinear
         * solver should use to reject a step.
         */
        [[nodiscard]] Real psi_lower_bound() const noexcept;

    private:
        Real alpha_{};
        Real theta_degrees_{};
        Real nonlinearity_scale_{};

        Real cos_theta_{};
        Real sin_theta_{};
        Real amplitude_{};
    };

    static_assert(AnalyticalSolution<Fadeev>);

    } // namespace gss