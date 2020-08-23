#define BOOST_TEST_MODULE ph_tests
#include <boost/test/included/unit_test.hpp>
#include "atomicdensity.h"

#include "test_fixtures.h"

template <int NTH,typename H>
struct PhFixture {
    PhFixture():traj{true},ret{&traj.traj,{20,20,20},NTH,1}{
        data.path+="ad/";
    }
    TrajSetup traj;
    using AD = AtomicDensity<Traiettoria,H>;
    AD ret;
    DataRegression<H> data;

    AD & calc(int primo){
        ret.reset(150-primo);
        ret.calcola(primo);
        return ret;
    }
    size_t size(){auto s= ret.get_shape(); return s[0]*s[1]*s[2];}
};

using PhFixture_1 = PhFixture<1,long>;

BOOST_FIXTURE_TEST_SUITE(test_atomic_density, PhFixture_1 )
BOOST_AUTO_TEST_CASE(test_position_histogram)
{
    calc(0);
    BOOST_TEST(data.test_regression("test_calcola_0",ret.accesso_lista(),size()));
    calc(12);
    BOOST_TEST(data.test_regression("test_calcola_1",ret.accesso_lista(),size()));
}
BOOST_AUTO_TEST_SUITE_END()
