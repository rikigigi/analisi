//#define BOOST_TEST_MODULE lammps2020_msd_tests
#include <boost/test/included/unit_test.hpp>
#include "test_fixtures.h"
#include "neighbour.h"


struct TrajTest{
    TrajTest(bool pbc=false,bool l2020=false):t{pbc,l2020}{};
    TrajSetup t;
    DataRegression<double> data;
};

struct NeighTest{
    NeighTest(bool l2020=true):t{true,l2020},n{&t.traj,{{99,1.0,1.0},{67,1.0,1.0}}} {};
        TrajSetup t;
        DataRegression<double> data;
        Neighbours<Trajectory,double> n;
};


BOOST_FIXTURE_TEST_SUITE(traj_min_image,TrajTest)
BOOST_AUTO_TEST_CASE(min_image){

    size_t n_atoms=t.traj.get_natoms();
    double * distances=new double[4*n_atoms*n_atoms];

    for (size_t i=0;i<n_atoms;++i){
        for (size_t j=0;j<n_atoms;++j) {
            distances[i*4*n_atoms+j*4+3]=t.traj.d2_minImage(i,j,0,0,distances+i*4*n_atoms+j*4);
            for (size_t idim=0;idim<3;idim++) {
                BOOST_TEST (fabs(distances[i*4*n_atoms+j*4+idim])<=t.traj.scatola(0)[3+idim]);
            }
        }
    }

    BOOST_TEST(data.test_regression("min_image",distances,4*n_atoms*n_atoms));

    delete [] distances;

}
BOOST_AUTO_TEST_SUITE_END()




BOOST_AUTO_TEST_CASE(pbc){
    {
        auto tt=TrajTest(true,false);
        size_t n_atoms=tt.t.traj.get_natoms();
        BOOST_TEST(tt.data.test_regression("pbc_1",tt.t.traj.posizioni(0,0),3*n_atoms));
    }
    {
        auto tt=TrajTest(true,true);
        size_t n_atoms=tt.t.traj.get_natoms();
        BOOST_TEST(tt.data.test_regression("pbc_2",tt.t.traj.posizioni(0,0),3*n_atoms));
    }


}

BOOST_FIXTURE_TEST_SUITE(traj_neigh,NeighTest)
BOOST_AUTO_TEST_CASE(neigh_test) {
    n.update_neigh(0,true);
    n.get_sann(0,0);
}

BOOST_AUTO_TEST_SUITE_END()
