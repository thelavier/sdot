#pragma once

#include "../Integration/SpaceFunctions/Constant.h"
#include "../Integration/FunctionEnum.h"
#include <vector>

namespace sdot {

/**
   We assume that grid has already been initialized by diracs
*/
template<class TF,class Grid,class Bounds,class Pt,class Func>
void get_internal_energy( TF *res, Grid &grid, Bounds &bounds, const Pt *positions, const TF *weights, std::size_t nb_diracs, const Func &radial_func ) {
    grid.for_each_laguerre_cell( [&]( auto &lc, auto num_dirac, int ) {
        TF ie = 0;
        bounds.for_each_intersection( lc, [&]( auto &cp, const auto &space_func ) {
            ie +=  cp.internal_energy( space_func, radial_func.func_for_final_cp_integration(), weights[ num_dirac ] );
        } );
        res[ num_dirac ] = ie;
    }, bounds.englobing_convex_polyhedron(), positions, weights, nb_diracs, false, radial_func.need_ball_cut() );
}

template<class TF,class Grid,class Bounds,class Pt>
void get_internal_energy( TF *res, Grid &grid, Bounds &bounds, const Pt *positions, const TF *weights, std::size_t nb_diracs ) {
    get_internal_energy( res, grid, bounds, positions, weights, nb_diracs, FunctionEnum::Unit() );
}

} // namespace sdot
