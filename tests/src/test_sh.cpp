//#define BOOST_TEST_MODULE sh_tests
#include <boost/test/included/unit_test.hpp>
#include "sphericalcorrelations.h"
#include "steinhardt.h"
#include <vector>
#include <algorithm>

#include "test_fixtures.h"

template <int l,int NTH>
struct ShFixture{
    ShFixture() : nbin{4}, natoms{traj.traj.get_natoms()}, ntypes{traj.traj.get_ntypes()}, sh{&(traj.traj), {{0.5, 3.0},{0.5, 3.0},{0.5, 3.0},{0.5, 3.0}}, nbin, 17, NTH, 13, 10,false} {
        sh.reset(100);
        data.path=data.path+"sh/";
    }
    ~ShFixture(){}
    TrajSetup traj;
    const unsigned int nbin;
    const size_t natoms;
    const size_t ntypes;
    SphericalCorrelations<l,double,Traiettoria> sh;
    DataRegression<double> data;
    double  workspace[(l+1)*(l+1)], cheby[(l+1)*2];
    void calc(int timestep, double * res){
        sh.calc(timestep, res, workspace,cheby);
    }
    size_t new_res_array_size(){
	return (l+1)*(l+1)*natoms*nbin*ntypes;
    }
    double* new_res_array() {
	return new double[new_res_array_size()];
    }
};

template<int l,int NTH, int skip=17>
struct SteinhardtFixture{
    SteinhardtFixture(): nbin{4}, natoms{traj.traj.get_natoms()}, ntypes{traj.traj.get_ntypes()},
        sh{&(traj.traj),{{0.5, 3.0},{0.5, 3.0},{0.5, 3.0},{0.5, 3.0}},nbin,100,{4,6},NTH,skip,true,{}} {
        sh.reset(traj.traj.get_nloaded_timesteps());
        data.path=data.path+"sh_stein/";
    }
    TrajSetup traj;
    const unsigned int nbin;
    const size_t natoms;
    const size_t ntypes;
    Steinhardt<l,double,Traiettoria> sh;
    DataRegression<double> data;
};
template<int l,int NTH>
struct SteinhardtFixtureNeigh{
    SteinhardtFixtureNeigh(): nbin{4}, natoms{traj.traj.get_natoms()}, ntypes{traj.traj.get_ntypes()},
        sh{&(traj.traj),{{0.5, 3.0},{0.5, 3.0},{0.5, 3.0},{0.5, 3.0}},nbin,100,{4,6},NTH,13,false,{{19,4.0,4.0},{17,4.0,4.0}}} {
        sh.reset(traj.traj.get_nloaded_timesteps());
        data.path=data.path+"sh_stein/";
    }
    TrajSetup traj;
    const unsigned int nbin;
    const size_t natoms;
    const size_t ntypes;
    Steinhardt<l,double,Traiettoria> sh;
    DataRegression<double> data;
};

typedef ShFixture<10,1> ShFix10_1 ;
typedef ShFixture<10,2> ShFix10_2 ;
typedef ShFixture<10,3> ShFix10_3 ;

typedef SteinhardtFixture<6,1,17> StFix6_1;
typedef SteinhardtFixture<6,2,17> StFix6_2;
typedef SteinhardtFixture<6,3,17> StFix6_3;
typedef SteinhardtFixture<6,1,53> StFix6_sk_1;
typedef SteinhardtFixture<6,2,53> StFix6_sk_2;
typedef SteinhardtFixture<6,3,53> StFix6_sk_3;

typedef SteinhardtFixtureNeigh<6,1> StFix6_tr_1;
typedef SteinhardtFixtureNeigh<6,2> StFix6_tr_2;
typedef SteinhardtFixtureNeigh<6,3> StFix6_tr_3;

#define TESTS(T,suff)\
BOOST_FIXTURE_TEST_SUITE(sh ## T, T )\
BOOST_AUTO_TEST_CASE(test_calcola)\
{\
    sh.calcola(0);\
    BOOST_TEST(data.test_regression("test_calcola" # suff ,sh.accesso_lista(),sh.lunghezza()));\
}\
BOOST_AUTO_TEST_SUITE_END()

TESTS(ShFix10_1,)
TESTS(ShFix10_2,)
TESTS(ShFix10_3,)

TESTS(StFix6_1,)
TESTS(StFix6_2,)
TESTS(StFix6_3,)
TESTS(StFix6_tr_1,_tr_neigh)
TESTS(StFix6_tr_2,_tr_neigh)
TESTS(StFix6_tr_3,_tr_neigh)
TESTS(StFix6_sk_1,_skip)
TESTS(StFix6_sk_2,_skip)
TESTS(StFix6_sk_3,_skip)

BOOST_FIXTURE_TEST_SUITE(sh_snapshots, ShFix10_1)
BOOST_AUTO_TEST_CASE(test_single_snapshot)
{
    BOOST_TEST(new_res_array_size()==sh.get_snap_size());
    double * res = new_res_array();
    calc(0,res);
    BOOST_TEST(data.test_regression("test_single_snapshot",res,new_res_array_size()));
    calc(0,res);
    BOOST_TEST(data.test_regression("test_single_snapshot",res,new_res_array_size()));
    delete [] res;
}
BOOST_AUTO_TEST_CASE(test_corr_consistency){
    double *sh1= new_res_array();
    double *sh2= new_res_array();
    double *sh3= new_res_array();
    double *res1= new double[sh.get_final_snap_size()];
    double *res2= new double[sh.get_final_snap_size()];
    int * cont=new int[traj.traj.get_ntypes()];
    data.fuzz(sh1,sh.get_snap_size());
    data.fuzz(sh2,sh.get_snap_size());
    data.fuzz(sh3,sh.get_snap_size());
    data.fuzz(res1,sh.get_final_snap_size());
    data.fuzz(res2,sh.get_final_snap_size());
    calc(2,sh1);
    data.fuzz(sh3,sh.get_snap_size());
    data.fuzz(res1,sh.get_final_snap_size());
    data.fuzz(res2,sh.get_final_snap_size());
    calc(17,sh2);
    data.fuzz(sh3,sh.get_snap_size());
    data.fuzz(res1,sh.get_final_snap_size());
    data.fuzz(res2,sh.get_final_snap_size());
    sh.corr_sh_calc(sh1,sh2,res1,sh3,sh.get_snap_size(),sh.get_final_snap_size(),cont);
    data.fuzz(sh3,sh.get_snap_size());
    data.fuzz(res2,sh.get_final_snap_size());
    sh.corr_sh_calc(sh1,sh2,res2,sh3,sh.get_snap_size(),sh.get_final_snap_size(),cont);
    BOOST_TEST(data.is_same(res1,res2,sh.get_final_snap_size()));
    data.fuzz(sh3,sh.get_snap_size());
    data.fuzz(res2,sh.get_final_snap_size());
    sh.corr_sh_calc(sh1,sh2,res2,sh1,sh.get_snap_size(),sh.get_final_snap_size(),cont);
    BOOST_TEST(data.is_same(res1,res2,sh.get_final_snap_size()));

    delete [] sh1;
    delete [] sh2;
    delete [] sh3;
    delete [] res1;
    delete [] res2;
    delete [] cont;

}
BOOST_AUTO_TEST_SUITE_END()


