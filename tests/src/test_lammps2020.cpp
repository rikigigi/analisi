#define BOOST_TEST_MODULE lammps2020_msd_tests
#include <boost/test/included/unit_test.hpp>
#include "msd.h"
#include "test_fixtures.h"

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

