//#define BOOST_TEST_MODULE gk_tests
#include <boost/test/included/unit_test.hpp>

#include "test_fixtures.h"
#include "greenkuboNcomponentionicfluid.h"
#include "readlog.h"

template<unsigned int nth>
struct GkFixture{
    GkFixture() : gk{&traj.traj,".test",1,{"c_flux[1]","c_vcm[1][1]"},false,0,nth} {
        data.path=data.path+"gk/";
        gk.reset(5000);
    }
    LogSetup traj;
    DataRegression<double> data;
    GreenKuboNComponentIonicFluid<ReadLog<>, double,double> gk;
};


BOOST_FIXTURE_TEST_SUITE(gk_traj_1, GkFixture<1>)
BOOST_AUTO_TEST_CASE(test_gk_2comp_nth1)
{
    gk.calculate(0);
    BOOST_TEST(data.test_regression("gk_2comp",gk.access_vdata(),gk.lunghezza()));
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_FIXTURE_TEST_SUITE(gk_traj_2, GkFixture<2>)
BOOST_AUTO_TEST_CASE(test_gk_2comp_nth2)
{
    gk.calculate(0);
    BOOST_TEST(data.test_regression("gk_2comp",gk.access_vdata(),gk.lunghezza()));
}

BOOST_AUTO_TEST_SUITE_END()
