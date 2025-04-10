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

    // F Star Functions
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

    // Inverse Transform Functions
    template<class PT>
    inline auto seed_inverse(PT c) const {
        constexpr std::size_t dim = sdot::point_dimension<PT>::value;
        if constexpr (dim == 2) {
            TS z1 = - c[0] / c[1];
            TS z2 = - TS(1) / (TS(2) * c[1]);
            return sdot::Point2<TS>(z1, z2);
        } else if constexpr (dim == 3) {
            TS z1 = - c[0] / c[2];
            TS z2 = - c[1] / c[2];
            TS z3 = - TS(1) / (TS(2) * c[2]);
            return sdot::Point3<TS>(z1, z2, z3);
        } else {
            static_assert(dim == 2 || dim == 3, "Only 2D and 3D point types are supported.");
        }
    }

    template<class PT>
    inline double weight_inverse(double psi, PT z) const {
        constexpr std::size_t dim = sdot::point_dimension<PT>::value;
        if constexpr (dim == 2) {
            double term1 = (z[0] / (2 * z[1])) * (z[0] / (2 * z[1]));
            double term2 = (1 / (2 * z[1])) * (1 / (2 * z[1]));
            double term3 = (f_cor * f_cor / (2 * z[1])) * z[0] * z[0];
            return psi - (term1 + term2 - term3 + c_p * pi_0);
        } else if constexpr (dim == 3) {
            double term1 = (z[0] / (2 * z[2])) * (z[0] / (2 * z[2]));
            double term2 = (z[1] / (2 * z[2])) * (z[1] / (2 * z[2]));
            double term3 = (1 / (2 * z[2])) * (1 / (2 * z[2]));
            double term4 = (f_cor * f_cor / (2 * z[2])) * (z[0] * z[0] + z[1] * z[1]);
            return psi - (term1 + term2 + term3 - term4);
        } else {
            static_assert(dim == 2 || dim == 3, "Only 2D and 3D point types are supported.");
        }
    }

    // Centroid Coefficients
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
        // Precomputed constants
        double f2    = f_cor * f_cor;      // f^2

        // Non-gamma terms:
        double term1 = -2.0 * p0[1] * p1[1]
            - 2.0 * (p0[1] * p0[1]) * (-1.0 + s)
            + 4.0 * p0[1] * p1[1] * s
            - 2.0 * (p1[1] * p1[1]) * s
            + 4.0 * p0[1] * p1[0] * z[0]
            - 2.0 * p0[0] * p1[1] * z[0]
            + 2.0 * p0[0] * p0[1] * (-1.0 + s) * z[0]
            - 2.0 * p0[1] * p1[0] * s * z[0]
            - 2.0 * (p0[0] - p1[0]) * p1[1] * s * z[0]
            - f2 * p0[1] * z[0] * z[0]
            + f2 * p1[1] * z[0] * z[0]
            + 2.0 * p0[1] * (c_p * pi_0 + w) * z[1]
            - 2.0 * p1[1] * (c_p * pi_0 + w) * z[1];

        // Gamma-multiplied terms:
        double term2 = p0[1] * p1[1] * (4.0 - 8.0 * s)
            + 4.0 * (p0[1] * p0[1]) * (-1.0 + s)
            + 4.0 * (p1[1] * p1[1]) * s
            + 2.0 * p0[0] * p1[1] * z[0]
            - 4.0 * p0[0] * p0[1] * (-1.0 + s) * z[0]
            + 4.0 * (p0[0] - p1[0]) * p1[1] * s * z[0]
            + p0[1] * p1[0] * (-6.0 + 4.0 * s) * z[0]
            + f2 * p0[1] * z[0] * z[0]
            - f2 * p1[1] * z[0] * z[0]
            - 2.0 * p0[1] * (c_p * pi_0 + w) * z[1]
            + 2.0 * p1[1] * (c_p * pi_0 + w) * z[1];

        // Sum both parts and return the result
        return term1 + gamma * term2;
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
        double term3 = (8.0 * z[1] * z[1] * gamma * (gamma - 1.0) * A * X * X) / ((6.0 * gamma * gamma - 7.0 * gamma + 2.0) * (B * B));
        double term4 = (8.0 * z[1] * z[1] * z[1] * (gamma - 1.0) * (gamma - 1.0) * X * X * X) / ((6.0 * gamma * gamma - 7.0 * gamma + 2.0) * (B * B * B));

        return term1 + term2 - term3 + term4;
        
    }

    //AppellF1 Hypergeometric Function
    template <typename T>
    static T appell_F1(T a, T b1, T b2, T c, T x, T y) {
        const T tol = static_cast<T>(1e-12);
        const int max_iter = 1000;
        T sum = 0;
        for (int k = 0; k < max_iter; ++k) {
            T term_k = 0;
            for (int m = 0; m <= k; ++m) {
                int n = k - m;
                T term = ( std::tgamma(a + m + n) / std::tgamma(a) )
                         / ( std::tgamma(c + m + n) / std::tgamma(c) );
                term *= ( std::tgamma(b1 + m) / std::tgamma(b1) ) / std::tgamma(m + 1);
                term *= ( std::tgamma(b2 + n) / std::tgamma(b2) ) / std::tgamma(n + 1);
                term *= std::pow(x, m) * std::pow(y, n);
                term_k += term;
            }
            sum += term_k;
            if (k > 0 && std::abs(term_k) < tol * std::abs(sum)) {
                break;
            }
        }
        return sum;
    }

    // Hessian Boundary Integrand
    inline TS hess_bdry_integrand(TS t, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk, TS w) const
    {
        // Compute the point on the edge for parameter t.
        sdot::Point2<TS> pt = t * p1 + (1 - t) * p0;

        // Evaluate the cost function at pt.
        TS cfunc_pt = (*this)(pt, zi, w);

        // Evaluate the norm of the gradient of the cost function at pt. 
        TS grad_norm = std::sqrt(hess_bdry_A1(p0, p1, zi, zk) * t * t + hess_bdry_A2(p0, p1, zi, zk) * t + hess_bdry_A3(p0, p1, zi, zk));

        return fsd1(w - cfunc_pt) / grad_norm;
    }

    // Hessian Boundary Coefficients
    inline double hess_bdry_A1 (const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk) const {
        double dp = p1[0] - p0[0];
        double dz = zi[1] - zk[1];
        return (dp * dp * dz * dz) / ((zi[1] * zi[1]) * (zk[1] * zk[1]));

    }

    inline double hess_bdry_A2 (const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk) const {
        double dp = p1[0] - p0[0];
        double dz = zi[1] - zk[1];
        double term = p0[0] * dz + (f_cor * f_cor) * (zi[0] * zk[1] - zi[1] * zk[0]);
        return (2.0 * dp * dz * term) / ((zi[1] * zi[1]) * (zk[1] * zk[1]));
        
    }

    inline double hess_bdry_A3 (const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk) const {
        double dz = zi[1] - zk[1];
        double term = p0[0] * dz + (f_cor * f_cor) * (zi[0] * zk[1] - zi[1] * zk[0]);
        return ((g * g * dz * dz) + (term * term)) / ((zi[1] * zi[1]) * (zk[1] * zk[1]));
        
    }

    inline double hess_bdry_K_plus(double s, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk, double w) const {
        double A1 = hess_bdry_A1(p0, p1, zi, zk);
        double A2 = hess_bdry_A2(p0, p1, zi, zk);
        double A3 = hess_bdry_A3(p0, p1, zi, zk);
        double cfuncp0 = this->operator()(p0, zi, w);
        double D = std::sqrt(A2*A2 - 4 * A1 * A3);
        double commonB = (p1[1] - p0[1]) + zi[0]*(p0[0] - p1[0]);
        double num = (D + 2*A1*s + A2) * commonB;
        double den = (D + A2) * commonB + 2*A1 * zi[1] * (w - cfuncp0);
        return std::sqrt(num / den);
    }

    inline double hess_bdry_K_minus(double s, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk, double w) const {
        double A1 = hess_bdry_A1(p0, p1, zi, zk);
        double A2 = hess_bdry_A2(p0, p1, zi, zk);
        double A3 = hess_bdry_A3(p0, p1, zi, zk);
        double cfuncp0 = this->operator()(p0, zi, w);
        double D = std::sqrt(A2*A2 - 4 * A1 * A3);
        double commonB = (p1[1] - p0[1]) + zi[0]*(p0[0] - p1[0]);
        double num = (D - 2*A1*s - A2) * commonB;
        double den = (D - A2) * commonB - 2*A1 * zi[1] * (w - cfuncp0);
        return std::sqrt(num / den);
    }

    inline double hess_bdry_C_plus(double s, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk, double w) const {
        double A1 = hess_bdry_A1(p0, p1, zi, zk);
        double A2 = hess_bdry_A2(p0, p1, zi, zk);
        double A3 = hess_bdry_A3(p0, p1, zi, zk);
        double D = std::sqrt(A2 * A2 - 4 * A1 * A3);
        double commonB = (p1[1] - p0[1]) + zi[0] * (p0[0] - p1[0]);
    
        // Compute the interpolated point p(s)
        sdot::Point2<TS> p_s = p0;
        p_s[0] = p0[0] + s * (p1[0] - p0[0]);
        p_s[1] = p0[1] + s * (p1[1] - p0[1]);
    
        // Evaluate the cost function at s and at s = 0.
        double cfuncps = this->operator()(p_s, zi, w); // c(Φ(γ(s)), zᵢ)
        double cfuncp0 = this->operator()(p0, zi, w);    // c(Φ(γ(0)), zᵢ)
    
        double num = -2.0 * A1 * zi[1] * (w - cfuncps);
        double den = ((D - A2) * commonB) - 2.0 * A1 * zi[1] * (w - cfuncp0);
    
        return num / den;
    }
    
    
    inline double hess_bdry_C_minus(double s, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk, double w) const {
        double A1 = hess_bdry_A1(p0, p1, zi, zk);
        double A2 = hess_bdry_A2(p0, p1, zi, zk);
        double A3 = hess_bdry_A3(p0, p1, zi, zk);
        double D = std::sqrt(A2 * A2 - 4 * A1 * A3);
        double commonB = (p1[1] - p0[1]) + zi[0] * (p0[0] - p1[0]);
    
        // Compute the interpolated point p(s)
        sdot::Point2<TS> p_s = p0;
        p_s[0] = p0[0] + s * (p1[0] - p0[0]);
        p_s[1] = p0[1] + s * (p1[1] - p0[1]);
    
        // Evaluate the cost function at s and at s = 0.
        double cfuncps = this->operator()(p_s, zi, w); // c(Φ(γ(s)), zᵢ)
        double cfuncp0 = this->operator()(p0, zi, w);    // c(Φ(γ(0)), zᵢ)
    
        double num = -2.0 * A1 * zi[1] * (w - cfuncps);
        double den = ((-D - A2) * commonB) - 2.0 * A1 * zi[1] * (w - cfuncp0);
    
        return num / den;
    }    

    // Parameters for the cost function and the fstar functions:
    TS kappa, gamma, g, f_cor, pi_0, c_p;
};

} // namespace FunctionEnum
} // namespace sdot