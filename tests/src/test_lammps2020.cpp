#define BOOST_TEST_MODULE lammps2020_msd_tests
#include <boost/test/included/unit_test.hpp>
#include "msd.h"
#include "test_fixtures.h"


struct TrajTest{
    TrajTest():t{false,false}{};
    TrajSetup t;
    DataRegression<double> data;
};


BOOST_FIXTURE_TEST_SUITE(traj,TrajTest)
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

template<int NTH>
struct MsdFixture {
    MsdFixture() : traj{false,true},
        msd{&traj.traj, 1,0,NTH,true}
    {}
    TrajSetup traj;
    MSD<Traiettoria> msd;
    DataRegression<double> data;
    MSD<Traiettoria> & calc(int primo) {
        msd.reset(75-primo);
        msd.calcola(primo);
        return msd;
    }
    size_t size(){auto s= msd.get_shape(); return s[0]*s[1]*s[2];}
};


#include "gofrt.h"

template<int NTH, int TMAX>
struct GofrFixture {
    GofrFixture() : traj{false,true},
        gofr{&traj.traj,
             0.9, //MIN r
             2.0, //MAX r
             20, //number of bins
             TMAX, // t max
             NTH, //number of threads
             70 // skip
}
    {}
    TrajSetup traj;
    Gofrt<double, Traiettoria> gofr;
    DataRegression<double> data;
    Gofrt<double, Traiettoria> & calc(int primo) {
        gofr.reset(75-primo);
        gofr.calcola(primo);
        return gofr;
    }
    size_t size(){auto s= gofr.get_shape(); return s[0]*s[1]*s[2];}
};

BOOST_AUTO_TEST_CASE(test_cell_permutator){
    double d[6]={1,2,3,4,5,6};
    double e[6]={1,2,3,4,5,6};
    Traiettoria::lammps_to_internal(d);
    Traiettoria::internal_to_lammps(d);
    for (int i=0;i<6;++i){
        BOOST_TEST(d[i]==e[i]);
    }
}

#define TEST_MULTIT(N)\
using MsdFixture_ ## N = MsdFixture<N>;\
BOOST_FIXTURE_TEST_SUITE(test_msd##N, MsdFixture_ ## N)\
BOOST_AUTO_TEST_CASE(test_msd_##N){\
    calc(0);\
    BOOST_TEST(data.test_regression("msd_"#N"a",msd.accesso_lista(),size()));\
    calc(13);\
    BOOST_TEST(data.test_regression("msd_"#N"b",msd.accesso_lista(),size()));\
}\
BOOST_AUTO_TEST_SUITE_END()

TEST_MULTIT(1)
TEST_MULTIT(3)
TEST_MULTIT(4)


#define TEST_MULTITGOFR(N,T)\
using GofrFixture_ ## N ## _ ## T = GofrFixture<N,T>;\
BOOST_FIXTURE_TEST_SUITE(test_gofr##N ## _ ## T, GofrFixture_ ## N ## _ ## T)\
BOOST_AUTO_TEST_CASE(test_gofr_##N ## _ ## T){\
    calc(0);\
    BOOST_TEST(data.test_regression("gofr_"#N"a"#T,gofr.accesso_lista(),size()));\
    calc(13);\
    BOOST_TEST(data.test_regression("gofr_"#N"b"#T,gofr.accesso_lista(),size()));\
}\
BOOST_AUTO_TEST_SUITE_END()

TEST_MULTITGOFR(1,1)
TEST_MULTITGOFR(3,1)
TEST_MULTITGOFR(3,4)
TEST_MULTITGOFR(1,4)
