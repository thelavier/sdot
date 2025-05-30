#pragma once

#include "../Support/N.h"
#include <cmath>

namespace sdot {
namespace FunctionEnum {

/**
*/
template<class TS>
struct ExpWmR2db {
    template<class PT,class TF>
    auto operator()( PT p, PT c, TF w ) const {
        using std::exp; return exp( ( w - norm_2_p2( p - c ) ) / eps );
    }

    const char *name() const {
        return "ExpWmR2db";
    }

    const auto &func_for_final_cp_integration() const {
        return *this;
    }

    N<0> need_ball_cut() const {
        return {};
    }

    template<class TF>
    void span_for_viz( const TF& f, TS w ) const {
        using std::sqrt;
        using std::exp;
        using std::pow;

        TS dr = 0.1 * sqrt( eps );
        f( 5 * sqrt( eps ), 50 * sqrt( eps ), exp( ( w - pow( 27.5 * sqrt( eps ), 2 ) ) / eps ) );
        for( TS r = 5 * sqrt( eps ); r > 0; r -= dr )
            f( r, r + dr, exp( ( w - pow( r + 0.5 * dr, 2 ) ) / eps ) );
    }

    TS eps;
};

}
}
