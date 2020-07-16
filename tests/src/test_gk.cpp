#define BOOST_TEST_MODULE gk_tests
#include <boost/test/included/unit_test.hpp>

#include "test_fixtures.h"
#include "greenkuboNcomponentionicfluid.h"

struct GkFixture{
    GkFixture() : gk{&traj.traj,".test",1,{"c_flux[1]","c_vcm[1][1]"}} {
        data.path=data.path+"gk/";
        gk.reset(5000);
    }
    LogSetup traj;
    DataRegression data;
    GreenKuboNComponentIonicFluid<double,double> gk;
};


BOOST_FIXTURE_TEST_SUITE(gk_traj, GkFixture)
BOOST_AUTO_TEST_CASE(test_gk_2comp)
{
    gk.calcola(0);
    BOOST_TEST(data.test_regression("gk_2comp",gk.accesso_lista(),gk.lunghezza()));
}

BOOST_AUTO_TEST_SUITE_END()
