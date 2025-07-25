#pragma once

#include "../../Geometry/ConvexPolyhedron2.h"
#include "../../Geometry/ConvexPolyhedron3.h"
#include <functional>

namespace sdot {

/**
  Pc is expected to contain
  - static int dim
  - TF => floating point type
  - TI => index type

*/
template<class Pc>
class SpZGrid {
public:
    // data from Pc
    static constexpr bool   allow_translations    = Pc::allow_translations;
    static constexpr int    dim                   = Pc::dim;
    using                   TF                    = typename Pc::TF;
    using                   TI                    = typename Pc::TI;

    // parameters
    static constexpr int    degree_w_approx       = 0;
    static constexpr bool   allow_mpi             = true;

    // static definitions
    static constexpr int    nb_coeffs_w_approx    = 1 + dim * ( degree_w_approx >= 1 ) + dim * ( dim + 1 ) / 2 * ( degree_w_approx >= 2 );
    using                   CP2                   = ConvexPolyhedron2<Pc>;
    using                   CP3                   = ConvexPolyhedron3<Pc>;
    using                   CP                    = typename std::conditional<dim==3,CP3,CP2>::type;
    using                   Pt                    = typename CP::Pt;

    // methods
    /* ctor */              SpZGrid               ( TI max_diracs_per_cell = 11 );

    void                    update                ( const Pt *positions, const TF *weights, TI nb_diracs, bool positions_have_changed = true, bool weights_have_changed = true, bool clip_at_sqrt_weight = false );

    int                     for_each_laguerre_cell( const std::function<void( CP &lc, TI num )> &f, const CP &starting_lc, const Pt *positions, const TF *weights, TI nb_diracs, bool stop_if_void_lc = false ); ///< starting_lc can be a polygonal bound
    int                     for_each_laguerre_cell( const std::function<void( CP &lc, TI num, int num_thread )> &f, const CP &starting_lc, const Pt *positions, const TF *weights, TI nb_diracs, bool stop_if_void_lc = false, bool ball_cut = false ); ///< version with num_thread
    template<int bc> int    for_each_laguerre_cell( const std::function<void( CP &lc, TI num, int num_thread )> &f, const CP &starting_lc, const Pt *positions, const TF *weights, TI nb_diracs, bool stop_if_void_lc, N<bc> ball_cut ); ///< version with num_thread

    template<class V> void  display               ( V &vtk_output, TF z = 0 ) const; ///< for debug purpose

    Pt                      inv_sym               ( Pt pt, int num_sym ) const { if ( allow_translations && num_sym >= 0 ) return pt - translations[ num_sym ]; return pt; };
    Pt                      sym                   ( Pt pt, int num_sym ) const { if ( allow_translations && num_sym >= 0 ) return pt + translations[ num_sym ]; return pt; }

    // values used by update
    TI                      max_diracs_per_cell;
    int                     depth_initial_send;
    std::vector<Pt>         translations;         ///< symetries

    struct CustomReplica {
        Pt position;       // The pre-calculated y_rep position (Eigen vector)
        TF weight;         // The pre-calculated psi_rep weight
        TI original_index; // The global index (0, 1, ...) of the seed this is a replica of
        TI replica_id;     // A unique ID for this replica to be used by the integral routine
    };
    std::vector<std::vector<CustomReplica>> custom_replicas_for_seed;

private:
    using                   CoeffsWApprox         = std::array<TF,nb_coeffs_w_approx>;
    using                   TFIsStd               = N<std::is_same<TF,double>::value>;

    struct                  PWI {
        TI                  num_dirac;
        Pt                  position;
        TF                  weight;
    };

    struct                  Box {
        float               dist_2                ( Pt p, TF w ) const;

        TI                  last_child_index;
        TI                  sibling_index;

        CoeffsWApprox       coeffs_w_approx;
        TI                  beg_indices;          ///< only for local boxes
        TI                  end_indices;          ///< only for local boxes
        Box*                last_child;
        Box*                sibling;
        std::vector<PWI>    ext_pwi;              ///< only for external boxes
        Pt                  min_pt;
        Pt                  max_pt;
        TI                  depth;
        int                 rank;
    };
    struct                  Neighbor {
        int                 mpi_rank;
        Box*                root;
    };

    Box*                    deserialize_rec       ( const std::vector<char> &dst, int ext_rank, N<0> );
    Box*                    deserialize_rec       ( const std::vector<char> &dst, int ext_rank, N<1> );
    std::vector<char>       serialize_rec         ( const Pt *positions, const TF *weights, std::vector<Box *> front, TI max_depth, N<0> );
    std::vector<char>       serialize_rec         ( const Pt *positions, const TF *weights, std::vector<Box *> front, TI max_depth, N<1> );
    void                    initial_send          ( const Pt *positions, const TF *weights );
    void                    update_box            ( const Pt *positions, const TF *weights, Box *box, TI beg_indices, TI end_indices, TI depth );
    static TI               nb_diracs             ( Box *box );
    std::string             ext_info              () const;

    template                <class Node>
    bool                    can_be_evicted        ( CP &lc, Pt c0, TF w0, Box *box, int num_sym, std::vector<Node *> &front );

    template<class TA> TF   w_approx              ( const TA &c, Pt x ) const;

    std::vector<TI>         dirac_indices;
    std::vector<Neighbor>   neighbors;
    std::deque<Box>         boxes;
    Box*                    root;
};

} // namespace sdot

#include "SpZGrid.tcc" // IWYU pragma: export

