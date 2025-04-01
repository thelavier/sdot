#pragma once

#include "../Support/N.h"
#include "../Support/Assert.h"  // Make sure to include the assertions header.
#include <cmath>
#include <tuple>

namespace sdot {
    template <class TF> struct Point2;  // forward declaration
    template <class TF> struct Point3;

    // Manually define a trait that maps a point type to its dimension. 
    template<typename T> struct point_dimension; // leave undefined
    template<typename TF> struct point_dimension< Point2<TF> > { static constexpr std::size_t value = 2; };
    template<typename TF> struct point_dimension< Point3<TF> > { static constexpr std::size_t value = 3; };

namespace FunctionEnum {

template<class TS>
struct CompressibleFunc {
    // Evaluate the cost function at point p, given seed c and weight w.
    template<class PT, class TF>
    TF operator()( PT p, PT c, TF w ) const {
        constexpr std::size_t dim = sdot::point_dimension<PT>::value;
        if constexpr (dim == 2) {
            // 2D implementation:
            // Ensure the second coordinate of the seed is non-zero.
            ASSERT(c[1] != 0, "2D seed's 2nd coordinate must be non-zero");
            TF fc2 = static_cast<TF>(f_cor * f_cor);
            TF cexp1 = p[0] / fc2;
            TF cexp2 = (p[1] - (p[0] * p[0]) / static_cast<TF>(2) / fc2) / static_cast<TF>(g);
            TF term1 = fc2 / static_cast<TF>(2) * (cexp1 - c[0]) * (cexp1 - c[0]);
            TF term2 = static_cast<TF>(g) * cexp2;
            return (term1 + term2) / c[1] - static_cast<TF>(c_p * pi_0);
        } else if constexpr (dim == 3) {
            // 3D implementation:
            // Ensure the third coordinate of the seed is non-zero.
            ASSERT(c[2] != 0, "3D seed's 3rd coordinate must be non-zero");
            TF fc2 = static_cast<TF>(f_cor * f_cor);
            TF cexp1 = p[0] / fc2;
            TF cexp2 = p[1] / fc2;
            TF cexp3 = (p[2] - (p[0] * p[0] + p[1] * p[1]) / static_cast<TF>(2) / fc2) / static_cast<TF>(g);
            TF term1 = fc2 / static_cast<TF>(2) * (cexp1 - c[0]) * (cexp1 - c[0]);
            TF term2 = fc2 / static_cast<TF>(2) * (cexp2 - c[1]) * (cexp2 - c[1]);
            TF term3 = static_cast<TF>(g) * cexp3;
            return (term1 + term2 + term3) / c[2];
        } else {
            static_assert(dim == 2 || dim == 3, "Only 2D and 3D point types are supported.");
            return TF{}; // Fallback to suppress compiler warnings.
        }
    }

    const char *name() const {
        return "CompressibleFunc";
    }

    const CompressibleFunc &func_for_final_cp_integration() const {
        return *this;
    }

    N<0> need_ball_cut() const {
        return N<0>();
    }

    template<class TF>
    void span_for_viz( const TF &f, TS w ) const {
        // Optionally, define a span for visualization.
    }

    // Helper power functions:
    inline double fsd1(double x) const {
        double gammap = 1 / (1 - 1 / gamma);
        return (x > 0) ? std::pow(kappa * gamma, 1 - gammap) * std::pow(x, gammap - 1) : 0.0;
    }

    inline double fs(double x) const {
        double gammap = 1 / (1 - 1 / gamma);
        return (x > 0) ? (1 / gammap) * std::pow(kappa * gamma, 1 - gammap) * std::pow(x, gammap) : 0.0;
    }

    inline double fsa1(double x) const {
        double gammap = 1 / (1 - 1 / gamma);
        return (x > 0) ? (1 / (gammap * (gammap + 1))) * std::pow(kappa * gamma, 1 - gammap) * std::pow(x, gammap + 1) : 0.0;
    }

    inline double fsa2(double x) const {
        double gammap = 1 / (1 - 1 / gamma);
        return (x > 0) ? (1 / (gammap * (gammap + 1) * (gammap + 2))) * std::pow(kappa * gamma, 1 - gammap) * std::pow(x, gammap + 2) : 0.0;
    }

    // Inverse functions (2D version):
    inline sdot::Point2<TS> seed_inverse(const sdot::Point2<TS>& c) const {
        TS z2 = - TS(1) / (TS(2) * c[1]);
        TS z1 = - c[0] / c[1];
        return sdot::Point2<TS>(z1, z2);
    }

    inline double weight_inverse(double psi, const sdot::Point2<TS>& z) const {
        double term1 = (z[0] / (2 * z[1])) * (z[0] / (2 * z[1]));
        double term2 = (1 / (2 * z[1])) * (1 / (2 * z[1]));
        double term3 = (f_cor * f_cor / (2 * z[1])) * z[0] * z[0];
        return psi - (term1 + term2 - term3 + c_p * pi_0);
    }

    inline double ctd_c1_coeff(double s, const sdot::Point2<TS>& z, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, double w) const {

        // Compute differences between the p1 and p0 coordinates.
        double dp1 = p1[0] - p0[0]; // p^{j+1}_1 - p^j_1
        double dp2 = p1[1] - p0[1]; // p^{j+1}_2 - p^j_2
        
        double term1 = -2.0 * s * z[0] * (dp1 * dp1) + 4.0 * p0[0] * dp2;
        double term2 = dp1 * (
              - (f_cor * f_cor) * (z[0] * z[0])
              + 2.0 * z[1] * (c_p * pi_0 + w)
              - 2.0 * p0[0] * z[0]
              - 2.0 * p0[1] * (1.0 + s)
              + 2.0 * s * p1[1]
        );
        double term3 = gamma * (
              4.0 * s * z[0] * (dp1 * dp1)
              - 6.0 * p0[0] * dp2
              + dp1 * (
                    (f_cor * f_cor) * (z[0] * z[0])
                    - 2.0 * z[1] * (c_p * pi_0 + w)
                    + 4.0 * p0[0] * z[0]
                    + 2.0 * p0[1] * (1.0 + 2.0 * s)
                    - 4.0 * s * p1[1]
              )
        );

        // Sum all terms
        return term1 + term2 + term3;
    }

    inline double ctd_c2_coeff_1(double s, const sdot::Point2<TS>& z, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, double w) const {
        // Compute the differences between p1 and p0 coordinates.
        double dp1 = p1[0] - p0[0]; // p^{j+1}_1 - p^j_1
        double dp2 = p1[1] - p0[1]; // p^{j+1}_2 - p^j_2

        
        double term1 = 2.0 * s * (dp2 * dp2);
        double term2 = -4.0 * p0[1] * z[0] * dp1;
        double term3 = dp2 * (
              - (f_cor * f_cor) * (z[0] * z[0])
              + 2.0 * z[1] * (c_p * pi_0 + w)
              + 2.0 * p0[1]
              + 2.0 * z[0] * p0[0] * (1.0 + s)
              - 2.0 * s * z[0] * p1[0]
        );
        double term4 = gamma * (
              -4.0 * s * (dp2 * dp2)
              + 6.0 * p0[1] * z[0] * dp1
              + dp2 * (
                    (f_cor * f_cor) * (z[0] * z[0])
                    - 2.0 * z[1] * (c_p * pi_0 + w)
                    - 4.0 * p0[1]
                    - 2.0 * z[0] * p0[0] * (1.0 + 2.0 * s)
                    + 4.0 * s * z[0] * p1[0]
              )
        );

        // Sum all terms
        return term1 + term2 + term3 + term4;
    }

    inline double ctd_c2_coeff_2(double s, const sdot::Point2<TS>& z, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, double w) const {
        // Compute common subexpressions:
        double A = (s - 1.0) * p0[0] - s * p1[0];
        double B = p1[1] - p0[1] + z[0] * (p0[0] - p1[0]);
        double cfuncp0 = this->operator()(p0, z, w);  
        double cfuncp1 = this->operator()(p1, z, w);  
        double X = p1[0] * (w - cfuncp0) - p0[0] * (w - cfuncp1);

         // Compute the terms of the formula:
        double term1 = 4.0 * A * A * A;
        double term2 = (4.0 * z[1] * gamma * A * A * X) / ((3.0 * gamma - 2.0) * B);
        double term3 = (4.0 * z[1] * gamma * (gamma - 1.0) * A * X * X) / ((6.0 * gamma * gamma - 7.0 * gamma + 2.0) * (B * B));
        double term4 = (2.0 * z[1] * (gamma - 1.0) * (gamma - 1.0) * X * X * X) / ((6.0 * gamma * gamma - 7.0 * gamma + 2.0) * (B * B * B));

        return term1 + term2 + term3 + term4;
        
    }

    // Parameters for the cost function and the fstar functions:
    TS kappa, gamma, g, f_cor, pi_0, c_p;
};

} // namespace FunctionEnum
} // namespace sdot