#include "../src/sdot/Support/Stream.h"
#include "../src/sdot/Geometry/ConvexPolyhedron2.h"
#include "../src/sdot/Display/VtkOutput.h"
#include "./catch_main.h"
#include <iostream>

//// nsmake cpp_flag -march=native

using namespace sdot;

TEST_CASE( "Arf", "Arf" ) {
    struct Pc { enum { dim = 2, allow_ball_cut = 1 }; using TI = std::size_t; using TF = double; using CI = std::string; };
    using  LC = ConvexPolyhedron2<Pc>;
    using  TF = LC::TF;
    using  std::sqrt;
    using  std::pow;

    // box { 0, 0 }, { 2, 1 }
    LC icp( LC::Box{ { 0, 0 }, { 2, 1 } } );
    icp.sphere_center = { 0, 0 };

    // full disc radius 2
    LC sph( LC::Box{ { 0, 0 }, { 8, 8 } } );
    sph.sphere_center = { 4, 4 };

    //
    FunctionEnum::Arfd arf;
    arf.values = []( double r ) -> double {
        return r < 1 ? 1 - r * r : 0;
    };
    arf.inp_scaling = []( TF w ) -> TF {
        return TF( 1 ) / sqrt( w );
    };
    arf.out_scaling = []( TF w ) -> TF {
        return w;
    };
    arf.stops = { 1.0 };

    // integrations
    // Integrate[ Integrate[ 100 - ( x * x + y * y ), { x, 0, 2 } ], { y, 0, 1 } ] => 590.0 / 3.0
    // Integrate[ Integrate[ 10 - ( x * x + y * y ), { x, 0, 2 } ], { y, 0, 1 } ] => 50.0 / 3.0
    // Integrate[ Integrate[ ( 1 - ( x * x +  y * y ) ) * UnitStep[ 1^2 - x^2 - y^2 ], { x, 0, 1 } ], { y, 0, 1 } ] => 50.0 / 3.0
    // Integrate[ Integrate[ ( 1/2 - ( x * x + y * y ) ) * UnitStep[ 1/2 - x^2 - y^2 ], { x, 0, 1 } ], { y, 0, 1 } ] => M_PI / 32
    CHECK_THAT( icp.integration( arf, 100 ), WithinAbs( 590.0 /  3, 1e-5 ) );
    CHECK_THAT( icp.integration( arf,  10 ), WithinAbs(  50.0 /  3, 1e-5 ) );
    CHECK_THAT( icp.integration( arf,   1 ), WithinAbs(  M_PI /  8, 1e-5 ) );
    CHECK_THAT( icp.integration( arf, 0.5 ), WithinAbs(  M_PI / 32, 1e-5 ) );

    // centroids
    // Integrate[ Integrate[ x * ( 100 - ( x * x + y * y ) ), { x, 0, 2 } ], { y, 0, 1 } ] / Integrate[ Integrate[ 100 - ( x * x + y * y ), { x, 0, 2 } ], { y, 0, 1 } ]
    // Integrate[ Integrate[ x * ( 10 - ( x * x + y * y ) ), { x, 0, 2 } ], { y, 0, 1 } ] / Integrate[ Integrate[ 10 - ( x * x + y * y ), { x, 0, 2 } ], { y, 0, 1 } ]
    P( icp.centroid( arf, 100 ), 293.0 / 295 );
    P( icp.centroid( arf,  10 ),  23.0 /  25 );
}
