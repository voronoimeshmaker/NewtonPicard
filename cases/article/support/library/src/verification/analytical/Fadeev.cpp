    #include "gsmaker/verification/analytical/Fadeev.hpp"

    #include <algorithm>
    #include <cmath>
    #include <numbers>
    #include <stdexcept>

    namespace gss {

    namespace {

    /*
     * std::exp overflows to infinity for arguments above roughly 709.78 in
     * IEEE-754 double precision.  Clamping the exponent at 700 keeps
     * free_function() and dF_dpsi_reference() total: they return a large but
     * finite value instead of infinity, so a diverging nonlinear iterate can
     * no longer inject a non-finite field into the gradient reconstruction.
     *
     * This is a finiteness guard, not a physical bound.  The physically
     * admissible range of Psi is far narrower and is reported by
     * Fadeev::psi_lower_bound(); for the default alpha the overflow threshold
     * lies about 1600 times below it.
     */
    inline constexpr double max_exponent = 700.0;

    [[nodiscard]] bool is_finite(Real value) noexcept {

        return std::isfinite(
            static_cast<double>(value));
    }

    [[nodiscard]] Real degrees_to_radians(
        Real degrees) noexcept {

        const auto factor =
            static_cast<Real>(
                std::numbers::pi_v<double> / 180.0);

        return degrees * factor;
    }

    [[nodiscard]] double saturated_exponent(
        Real psi_value) noexcept {

        return std::min(
            -2.0 *
                static_cast<double>(
                    psi_value),
            max_exponent);
    }

    } // namespace


    Fadeev::Fadeev()
        : Fadeev(
              static_cast<Real>(default_alpha),
              static_cast<Real>(default_theta_degrees)) {}



              
    Fadeev::Fadeev(
        Real alpha,
        Real theta_degrees,
        Real nonlinearity_scale)
        : alpha_(alpha),
          theta_degrees_(theta_degrees),
          nonlinearity_scale_(nonlinearity_scale) {

        if (!is_finite(alpha_)) {
            throw std::invalid_argument(
                "Fadeev alpha must be finite");
        }

        if (!is_finite(theta_degrees_)) {
            throw std::invalid_argument(
                "Fadeev theta must be finite");
        }

        /*
         * Only finiteness is required of the scale.  Zero and negative values
         * are admitted deliberately and are not degenerate cases to guard
         * against:
         *
         *   c = 0 removes the nonlinearity entirely.  dF/dPsi vanishes, the
         *   manufactured source carries the whole right-hand side, B = 0 and
         *   rho = 0, so Picard converges in one step.  That is the exact lower
         *   endpoint of the ladder and a useful check that the manufactured
         *   term is assembled with the right sign.
         *
         *   c < 0 flips the sign of the nonlinearity.  It is admitted because
         *   rho(c) = |c| rho(1) reaches the same radii from the other side,
         *   which separates an effect of the spectral radius from an effect of
         *   the sign of the free function.
         *
         * Neither changes the admissible range of Psi, which psi_lower_bound()
         * derives from the analytical field alone.
         */
        if (!is_finite(nonlinearity_scale_)) {
            throw std::invalid_argument(
                "Fadeev nonlinearity scale must be finite");
        }

        const Real theta_radians =
            degrees_to_radians(theta_degrees_);

        cos_theta_ =
            static_cast<Real>(
                std::cos(
                    static_cast<double>(
                        theta_radians)));

        sin_theta_ =
            static_cast<Real>(
                std::sin(
                    static_cast<double>(
                        theta_radians)));

        amplitude_ =
            static_cast<Real>(
                std::sqrt(
                    1.0 +
                    static_cast<double>(
                        alpha_ * alpha_)));
    }


    Real Fadeev::psi(
        Real x,
        Real y) const noexcept {

        const Real x_tilde =
            x * cos_theta_ -
            y * sin_theta_;

        const Real y_tilde =
            x * sin_theta_ +
            y * cos_theta_;

        const Real argument =
            alpha_ *
                static_cast<Real>(
                    std::cos(
                        static_cast<double>(
                            x_tilde))) +
            amplitude_ *
                static_cast<Real>(
                    std::cosh(
                        static_cast<double>(
                            y_tilde)));

        return static_cast<Real>(
            std::log(
                static_cast<double>(
                    argument)));
    }


    Real Fadeev::free_function(
        Real psi_value) const noexcept {

        return static_cast<Real>(
            -0.5 *
            static_cast<double>(
                nonlinearity_scale_) *
            std::exp(
                saturated_exponent(
                    psi_value)));
    }


    Real Fadeev::dF_dpsi_reference(
        Real psi_value) const noexcept {

        return static_cast<Real>(
            static_cast<double>(
                nonlinearity_scale_) *
            std::exp(
                saturated_exponent(
                    psi_value)));
    }


    /*
     * q = (1 - c) exp(-2 Psi*).
     *
     * Adding this to the equation keeps Psi* as the exact solution for every
     * c, because at Psi = Psi* the scaled nonlinearity c exp(-2 Psi*) and q
     * sum to exp(-2 Psi*), which is what the unscaled problem had.  The same
     * exponent clamp as free_function() is used, so the two terms saturate
     * together rather than at different values of Psi.
     *
     * At c = 1 the factor (1 - c) is exactly zero in IEEE-754 arithmetic and
     * the product is exactly zero, not merely small.  That is what lets a run
     * at c = 1 reproduce a pre-existing campaign digit for digit.
     */
    Real Fadeev::manufactured_source_density(
        Real exact_psi_value) const noexcept {

        return static_cast<Real>(
            (1.0 -
             static_cast<double>(
                 nonlinearity_scale_)) *
            std::exp(
                saturated_exponent(
                    exact_psi_value)));
    }


    Real Fadeev::psi_lower_bound() const noexcept {

        /*
         * alpha cos(x_tilde) + sqrt(1 + alpha^2) cosh(y_tilde) attains its
         * minimum sqrt(1 + alpha^2) - |alpha| at cos(x_tilde) = -1 and
         * y_tilde = 0.  That minimum is strictly positive for every real
         * alpha, so the logarithm below is always well defined.
         */
        const double minimum_argument =
            static_cast<double>(amplitude_) -
            std::abs(
                static_cast<double>(alpha_));

        return static_cast<Real>(
            std::log(
                minimum_argument));
    }

    } // namespace gss