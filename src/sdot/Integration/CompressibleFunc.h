#pragma once

#include "../Support/N.h"
#include "../Support/Assert.h"  // Make sure to include the assertions header.
#include <cmath>
#include <tuple>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

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

    inline bool IntegralType() const {
        return Int; // True for line integrals, False for vertex evaluations
    }

    inline int IntegralResolution() const {
        return Int_res;
    }

    N<0> need_ball_cut() const {
        return N<0>();
    }

    template<class TF>
    void span_for_viz( const TF &f, TS w ) const {
        // Optionally, define a span for visualization.
    }

    // Construct and immediately build the GL table
    CompressibleFunc(TS kappa_, TS gamma_, TS g_, TS f_cor_, TS pi_0_, TS c_p_, TS Int_, int Int_res_) : kappa(kappa_), gamma(gamma_), g(g_), f_cor(f_cor_), pi_0(pi_0_), c_p(c_p_), Int(Int_), Int_res(Int_res_) {
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

    inline double fsd1_sf(double x) const {
        double gammap = 1 / (1 - 1 / gamma);
        return (x > 0) ? std::pow(x, gammap - 1) : 0.0;
    }

    inline double fs_sf(double x) const {
        double gammap = 1 / (1 - 1 / gamma);
        return (x > 0) ? std::pow(x, gammap) : 0.0;
    }

    // Inverse Transform Functions
    template<class PT>
    inline auto seed_inverse(PT y) const {
        constexpr std::size_t dim = sdot::point_dimension<PT>::value;
        if constexpr (dim == 2) {
            TS z1 = - y[0] / y[1];
            TS z2 = - TS(1) / (TS(2) * y[1]);
            return sdot::Point2<TS>(z1, z2);
        } else if constexpr (dim == 3) {
            TS z1 = - y[0] / y[2];
            TS z2 = - y[1] / y[2];
            TS z3 = - TS(1) / (TS(2) * y[2]);
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

    inline double hess_bdry_coeff (const sdot::Point2<TS>& p, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk) const {
        double dz = zi[1] - zk[1];
        double x1 = p[0];
        double cross = -zi[1] * zk[0] + zi[0] * zk[1];
        double term = x1 * dz + f_cor * f_cor * cross;
        double denom = std::abs(zi[1] * zk[1]); 
        return std::sqrt((g * g * dz * dz + term * term)) / denom;
    }   

    // Integrands
    // Volume Integrand
    inline TS volume_integrand(TS t, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, TS w) const
    {
        // Compute the point on the edge for parameter t.
        sdot::Point2<TS> pt = t * p1 + (1 - t) * p0;

        // Evaluate the cost function at pt.
        TS cfunc_pt = (*this)(pt, zi, w);

        return fs_sf(w - cfunc_pt);
    }
    // Internal Energy Inegrand
    inline TS ie_integrand(TS t, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, TS w) const
    {
        // Compute the point on the edge for parameter t.
        sdot::Point2<TS> pt = t * p1 + (1 - t) * p0;

        // Evaluate the cost function at pt.
        TS cfunc_pt = (*this)(pt, zi, w);

        return (w - cfunc_pt) * std::pow(fsd1_sf(w - cfunc_pt), gamma);
    }
    // Centroid Component 1
    inline TS ctd_0_integrand(TS t, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, TS w) const
    {
        // Compute the point on the edge for parameter t.
        sdot::Point2<TS> pt = t * p1 + (1 - t) * p0;

        // Evaluate the cost function at pt.
        TS cfunc_pt = (*this)(pt, zi, w);

        return pt[0] * fs_sf(w - cfunc_pt);
    }
    // Centroid Component 2 Part 1
    inline TS ctd_1_1_integrand(TS t, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, TS w) const
    {
        // Compute the point on the edge for parameter t.
        sdot::Point2<TS> pt = t * p1 + (1 - t) * p0;

        // Evaluate the cost function at pt.
        TS cfunc_pt = (*this)(pt, zi, w);

        return pt[1] * fs_sf(w - cfunc_pt);
    }
    // Centroid Component 2 Part 2
    inline TS ctd_1_2_integrand(TS t, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, TS w) const
    {
        // Compute the point on the edge for parameter t.
        sdot::Point2<TS> pt = t * p1 + (1 - t) * p0;

        // Evaluate the cost function at pt.
        TS cfunc_pt = (*this)(pt, zi, w);

        return pt[0] * pt[0] * fs_sf(w - cfunc_pt);
    }
    // Hessian Volume Integrand
    inline TS hess_volume_integrand(TS t, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, TS w) const
    {
        // Compute the point on the edge for parameter t.
        sdot::Point2<TS> pt = t * p1 + (1 - t) * p0;

        // Evaluate the cost function at pt.
        TS cfunc_pt = (*this)(pt, zi, w);

        return fsd1_sf(w - cfunc_pt);
    }
    // Hessian Boundary Integrand
    inline TS hess_bdry_integrand(TS t, const sdot::Point2<TS>& p0, const sdot::Point2<TS>& p1, const sdot::Point2<TS>& zi, const sdot::Point2<TS>& zk, TS w) const
    {
        // Compute the point on the edge for parameter t.
        sdot::Point2<TS> pt = t * p1 + (1 - t) * p0;

        // Evaluate the cost function at pt.
        TS cfunc_pt = (*this)(pt, zi, w);
        double t_value = w - cfunc_pt;
        double f_star_prime_value = fsd1_sf(t_value);

        // Evaluate the norm of the gradient of the cost function at pt. 
        TS grad_norm = hess_bdry_coeff(pt, zi, zk);

        // build dp = p1 - p0
        TS dxdt = p1[0] - p0[0];
        TS dydt = p1[1] - p0[1];

        // Jacobian D Phi at pt:  [[1/f^2, 0], [–pt.x/f^2, 1/g]]
        TS inv_f2 = TS(1) / (f_cor * f_cor);
        TS inv_g  = TS(1) / g;

        // push‐forward tangent v = DPhi·dp
        TS v1 = inv_f2 * dxdt;
        
        // --- Deconstruct the v2 calculation for analysis ---
        TS v2_term1 = -pt[0] * inv_f2 * dxdt;
        TS v2_term2 = inv_g * dydt;
        TS v2 = v2_term1 + v2_term2;

        // arclength element
        TS ds = std::hypot(v1, v2);

        double result = (f_star_prime_value / grad_norm) * ds;

        // std::ostringstream debugStream;
        // debugStream << " w: " << w << "\n";
        // debugStream << " c(pt,z^i): " << cfunc_pt << "\n";
        // debugStream << " f_star_prime_value: " << f_star_prime_value << "\n";
        // debugStream << " grad_norm: " << grad_norm << "\n";
        // debugStream << " ds: " << ds << "\n";
        // debugStream << " result: " << result << "\n";
        // debugStream << "------------------------------------" << std::endl;

        // std::cout << debugStream.str() << std::flush;

        return result;
    }

    // storage for the table
    std::vector<TS> quad_nodes, quad_weights;

    // Parameters for the cost function class:
    TS kappa, gamma, g, f_cor, pi_0, c_p, Int;
    int Int_res;
};

} // namespace FunctionEnum
} // namespace sdot