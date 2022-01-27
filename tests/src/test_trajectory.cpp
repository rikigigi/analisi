#define BOOST_TEST_MODULE lammps2020_msd_tests
#include <boost/test/included/unit_test.hpp>
#include "test_fixtures.h"


struct TrajTest{
    TrajTest(bool pbc=false,bool l2020=false):t{pbc,l2020}{};
    TrajSetup t;
    DataRegression<double> data;
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
