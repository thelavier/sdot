#ifndef SDOT_ConvexPolyhedron2_H
#define SDOT_ConvexPolyhedron2_H

//// nsmake avoid_inc boost/

#include "../Integration/SpaceFunctions/Polynomial.h"
#include "../Integration/FunctionEnum.h"
#include "../Display/VtkOutput.h"
#include <boost/dynamic_bitset.hpp>
#include <functional>
#include <algorithm>
#include <bitset>

#define _USE_MATH_DEFINES
#include <math.h>

namespace sdot {

/**
  Pc must contain
    - dim (2, 3, ...)
    - TI (std::size_t, ...) => index type
    - TF (double, ...) => floating point type
  CI = Cut id type

  Beware: ball_cuts must be done AFTER the plane_cuts.
*/
template<class Pc>
class ConvexPolyhedron2 {
public:
    static constexpr bool          keep_min_max_coords       = false;
    static constexpr bool          allow_ball_cut            = Pc::allow_ball_cut;
    using                          VecBool                   = boost::dynamic_bitset<std::uint64_t>;
    using                          Node                      = std::size_t;
    using                          TI                        = typename Pc::TI; ///< index type
    using                          TF                        = typename Pc::TF; ///< scalar type
    using                          Pt                        = Point2<TF>;  ///< 3D point
    using                          CI                        = typename Pc::CI; ///< cut info

    //
    static constexpr bool          store_the_normals         = true;

    // types for simd
    // using                       AF                        = std::array<TF,64>;
    // using                       AC                        = std::array<CI,64>;
    // using                       AB                        = std::bitset<64>;

    // types for the ctor
    struct                         EnglobingSimplex          { Pt p; TF r; };
    struct                         BoundaryItem              { std::array<Pt,2> points; TF measure, a0, a1; CI id; Pt pos_integral, momentum[ Pc::dim ]; };
    struct                         Simplex                   { Pt pts[ 3 ]; };
    struct                         Box                       { Pt p0, p1; };

    /// we start from a triangle that includes the circle defined by englobing_center and englobing_radius (but this sphere is not used, it's just here to construct a triangle)
    /**/                           ConvexPolyhedron2         ( const EnglobingSimplex &es, CI cut_id = {} );
    /**/                           ConvexPolyhedron2         ( const Simplex &simplex, CI cut_id = {} );
    /**/                           ConvexPolyhedron2         ( const Box &box, CI cut_id = {} );
    /**/                           ConvexPolyhedron2         ( const ConvexPolyhedron2 &that );
    /**/                           ConvexPolyhedron2         ();

    // traversal
    template<class S,class R,class Grid> void for_each_boundary_measure ( const S &sf, const R &rf, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( TF boundary_measure, CI id)> &f, TF weight = 0, Pt home_position = {}  ) const;

    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::CompressibleFunc<TF> &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::ExpWmR2db<TF>        &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Arfd                 &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::WmR2                 &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Unit                 &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::R2                   &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
       
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::CompressibleFunc<TF> &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::ExpWmR2db<TF>        &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Arfd                 &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::WmR2                 &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Unit                 &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;
    template<class Grid> void      for_each_boundary_item    ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::R2                   &f, const Grid &grid, const std::size_t nb_diracs, const Pt* positions, const std::function<void( const BoundaryItem &boundary_item )> &cb, TF weight = 0, Pt home_position = {} ) const;

    void                           for_each_approx_seg       ( const std::function<void( Pt )> &f, TF max_ratio_area_error = 1e-1 ) const; ///<
    void                           for_each_simplex          ( const std::function<void( CI num_0, CI num_1 )> &f ) const;
    void                           for_each_bound            ( const std::function<void( Pt p0, Pt p1, CI id )> &f ) const;
    void                           for_each_node             ( const std::function<void( Pt v )> &f ) const;

    // information
    void                           display_html_canvas       ( std::ostream &os, TF weight, bool ext = false ) const;
    void                           write_to_stream           ( std::ostream &os ) const;
    Pt                             min_position              () const;
    Pt                             max_position              () const;
    void                           display_asy               ( std::ostream &os, const std::string &draw_info = "", const std::string &fill_info = "", bool fill = false, bool avoid_bounds = false, bool want_line = true ) const; ///< ouput asymptote format
    int                            is_a_ball                 () const { return _nb_points ? 0 : ( sphere_radius > 0 ? 1 : -1 ); } // 0 is not a ball. -1 if void, 1 if not void
    std::size_t                    nb_points                 () const { return _nb_points; }
    template<class V> void         display                   ( V &vo, const typename V::CV &cell_data = {}, bool filled = true, TF max_ratio_area_error = 1e-1, bool display_tangents = false ) const;
    template<class F> bool         all_pos                   ( const F &f ) const;
    Pt                             normal                    ( std::size_t n ) const { return { normals[ 0 ][ n ], normals[ 1 ][ n ] }; }
    Pt                             point                     ( std::size_t n ) const { return { points[ 0 ][ n ], points[ 1 ][ n ] }; }
    bool                           empty                     () const { return _nb_points == 0 && sphere_radius <= 0; }

    // modifications
    void                           intersect_with            ( const ConvexPolyhedron2 &cp );
    void                           set_cut_ids               ( CI cut_id ); ///< replace all the cut_ids
    template<int no> bool          plane_cut                 ( Pt origin, Pt normal, CI cut_id, N<no> normal_is_normalized ); ///< return true if effective cut
    bool                           plane_cut                 ( Pt origin, Pt normal, CI cut_id = {} ); ///< return true if effective cut
    void                           ball_cut                  ( Pt center, TF radius, CI cut_id = {} ); ///< beware: only one sphere cut is authorized, and it must be done after all the plane cuts.
    void                           clear                     ( Pt englobing_center, TF englobing_radius, CI englobing_cut_id = {} );

    // tests
    bool                           is_a_cutting_plane        ( Pt origin, Pt normal ) const;
    TF                             distance                  ( const Pt &pos, bool count_domain_boundaries = false ) const; // signed distance to the boundaries (<0 means inside)
    bool                           contains                  ( const Pt &pos ) const;

    // computations
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::CompressibleFunc<TF> &f, TF weight = 0 ) const { TODO; }
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::ExpWmR2db<TF>        &f, TF weight = 0 ) const { TODO; }
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Arfd                 &f, TF weight = 0 ) const { TODO; }
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::WmR2                 &f, TF weight = 0 ) const { TODO; }
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Unit                 &f, TF weight = 0 ) const;
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::R2                   &f, TF weight = 0 ) const { TODO; }

    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::CompressibleFunc<TF> &f, TF weight = 0 ) const;
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::ExpWmR2db<TF>        &f, TF weight = 0 ) const;
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Arfd                 &f, TF weight = 0 ) const;
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::WmR2                 &f, TF weight = 0 ) const;
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Unit                 &f, TF weight = 0 ) const;
    void                           add_centroid_contrib      ( Pt &ctd, TF &vol, const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::R2                   &f, TF weight = 0 ) const;

    void                           add_centroid_contrib      ( Pt &ctd, TF &vol ) const;

    //    TF                       boundary_measure          ( FunctionEnum::ExpWmR2db<TF> ) const;
    //    TF                       boundary_measure          ( FunctionEnum::Unit          ) const;
    //    TF                       boundary_measure          () const;

    Pt                             random_point              () const;

    template<class S,class R> Pt   centroid                  ( const S &sf, const R &rf, TF w = 0 ) const;
    Pt                             centroid                  () const;
    TF                             measure                   () const;

    TF                             internal_energy           ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::CompressibleFunc<TF> &f, TF weight = 0 ) const;
    TF                             internal_energy           ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::ExpWmR2db<TF>        &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::WmR2                 &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Unit                 &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Arfd                 &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Arf                  &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::R2                   &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::R4                   &f, TF weight = 0 ) const { TODO; return 0; }

    TF                             internal_energy           ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::CompressibleFunc<TF> &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::ExpWmR2db<TF>        &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::WmR2                 &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Unit                 &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Arfd                 &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Arf                  &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::R2                   &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             internal_energy           ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::R4                   &f, TF weight = 0 ) const { TODO; return 0; }

    TF                             integration               ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::CompressibleFunc<TF> &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             integration               ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::ExpWmR2db<TF>        &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             integration               ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::WmR2                 &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             integration               ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Unit                 &f, TF weight = 0 ) const;
    TF                             integration               ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Arfd                 &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             integration               ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::Arf                  &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             integration               ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::R2                   &f, TF weight = 0 ) const { TODO; return 0; }
    TF                             integration               ( const SpaceFunctions::Polynomial<TF,6> &sf, const FunctionEnum::R4                   &f, TF weight = 0 ) const { TODO; return 0; }

    TF                             integration               ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::CompressibleFunc<TF> &f, TF weight = 0 ) const;
    TF                             integration               ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::ExpWmR2db<TF>        &f, TF weight = 0 ) const;
    TF                             integration               ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::WmR2                 &f, TF weight = 0 ) const;
    TF                             integration               ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Unit                 &f, TF weight = 0 ) const;
    TF                             integration               ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Arfd                 &f, TF weight = 0 ) const;
    TF                             integration               ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::Arf                  &f, TF weight = 0 ) const;
    TF                             integration               ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::R2                   &f, TF weight = 0 ) const;
    TF                             integration               ( const SpaceFunctions::Constant<TF>     &sf, const FunctionEnum::R4                   &f, TF weight = 0 ) const;

    template<class Sf> TF          integration_der_wrt_weight( const Sf &sf, const FunctionEnum::CompressibleFunc<TF> &f, TF weight ) const;
    template<class Sf> TF          integration_der_wrt_weight( const Sf &sf, const FunctionEnum::ExpWmR2db<TF>        &f, TF weight ) const;
    template<class Sf> TF          integration_der_wrt_weight( const Sf &sf, const FunctionEnum::Arfd                 &f, TF weight ) const;
    template<class Sf> TF          integration_der_wrt_weight( const Sf &sf, const FunctionEnum::WmR2                 &f, TF weight ) const;
    template<class Sf,class FU> TF integration_der_wrt_weight( const Sf &sf, const FU                                 &f, TF weight ) const;

    // approximate computations
    TF                             boundary_measure_ap       ( TF max_ratio_area_error = 1e-4 ) const; ///<
    template<class Fu> TF          integration_ap            ( const Fu &func, TF weight = 0, std::size_t n = 1e6 ) const;
    template<class Fu> Pt          centroid_ap               ( const Fu &func, TF weight = 0, std::size_t n = 1e6 ) const; ///<
    TF                             measure_ap                ( TF max_ratio_area_error = 1e-4 ) const; ///<
    static TF                      pi                        () { if ( std::is_same<TF,double>::value ) return M_PI; using std::atan; return 4 * atan( TF( 1 ) ); }

    std::vector<TF>                normals[ 2 ];
    std::vector<TF>                points[ 2 ];
    std::vector<TF>                distances;
    std::vector<CI>                cut_ids;
    VecBool                        outside;
    VecBool                        arcs;

    Pt                             sphere_center;
    TF                             sphere_radius;
    CI                             sphere_cut_id;

    Pt                             min_coord;
    Pt                             max_coord;


private:
    enum                           CutType                   { LINE = 0, ARC = 1 };
    struct                         Cut                       { int cut_type; CI cut_id; Pt normal; Pt point; };
    template<class Coeffs> TF      _r_polynomials_integration( const Coeffs &coeffs, TF scaling = 1 ) const;
    template<class Coef> void      _r_centroid_integration   ( TF &r_x, TF &r_y, const Coef &coeffs, TF scale = 1 ) const;
    void                           _centroid_arc             ( Pt &ctd, TF &mea, Pt p0, Pt p1, TF coeff ) const;
    TF                             _arc_length               ( Pt p0, Pt p1 ) const;
    TF                             _arc_area                 ( Pt p0, Pt p1 ) const;

    void                           set_nb_points             ( std::size_t nb_points );

    std::size_t                    _nb_points;
    std::vector<Cut>               _tmp_cuts;
};

} // namespace sdot

#include "ConvexPolyhedron2.tcc"

#endif // SDOT_ConvexPolyhedron2_H
