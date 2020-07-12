#define BOOST_TEST_MODULE sh_tests
#include <boost/test/included/unit_test.hpp>
#include "config.h"
#include "traiettoria.h"
#include "sphericalcorrelations.h"
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

struct TestPath{
	TestPath(): path(PROJ_DIR "/tests/"){}
	std::string path;
};

struct TrajSetup{
    TrajSetup(): traj{path.path+"lammps.bin"}{
	traj.imposta_dimensione_finestra_accesso(150);
	traj.imposta_inizio_accesso(0);
    }
    TestPath path;
    Traiettoria traj;

};

struct DataRegression{
    DataRegression(): path{base_path.path + "cpp_regression_data/"} {}
    TestPath base_path;
    std::string path;
    bool test_regression(std::string name, double * data, size_t size) {
	if (size==0) {
	   BOOST_TEST_MESSAGE("zero size data not supported");
           return false;
	}
	std::ifstream from_fs(path+name, std::ios::binary);
        std::vector<unsigned char> buffer(std::istreambuf_iterator<char>(from_fs), {});
	if (buffer.size()==0) {
	    BOOST_TEST_MESSAGE("no data in file. Writing it");
        BOOST_TEST_MESSAGE(path+name);
	    from_fs.close();
	    std::ofstream output(path+name, std::ios::binary );
	    std::copy(
                (char*)data,(char*)data+size*sizeof(double),
        std::ostreambuf_iterator<char>(output)
			    );
	    return false;
	}
    if (buffer.size() < size*sizeof(double)) {
        BOOST_TEST_MESSAGE("wrong size of data");
        BOOST_TEST_MESSAGE(size*sizeof(double));
        BOOST_TEST_MESSAGE(buffer.size());
	    return false;
	}
	double * data_fs= (double*) buffer.data();
	for (size_t i=0;i<size;++i){
	    if (data_fs[i]!=data[i]) return false;
	}
	return true;
    }
};

template <int l>
struct ShFixture{
    ShFixture() : nbin{4}, natoms{traj.traj.get_natoms()}, ntypes{traj.traj.get_ntypes()}, sh{&(traj.traj), 0.5, 3.0, nbin, 17, 2, 13, false} {
	    sh.reset(100);
            data.path=data.path+"sh/";
    }
    ~ShFixture(){}
    TrajSetup traj;
    unsigned int nbin;
    int natoms,ntypes;
    SphericalCorrelations<l,double,Traiettoria> sh;
    DataRegression data;
    double  workspace[(l+1)*(l+1)], cheby[(l+1)*2];
    void calc(int timestep, double * res){
	double l_[3]={traj.traj.scatola(timestep)[1]-traj.traj.scatola(timestep)[0],
                     traj.traj.scatola(timestep)[3]-traj.traj.scatola(timestep)[2],
                     traj.traj.scatola(timestep)[5]-traj.traj.scatola(timestep)[4]};
	    sh.calc(timestep, res, workspace,cheby,l_);
    }
    size_t new_res_array_size(){
	return (l+1)*(l+1)*natoms*nbin*ntypes;
    }
    double* new_res_array() {
	return new double[new_res_array_size()];
    }
};


BOOST_FIXTURE_TEST_SUITE(sh, ShFixture<10>)
    
    BOOST_AUTO_TEST_CASE(test_single_snapshot)
    {
	double * res = new_res_array();
        calc(0,res);
	BOOST_TEST(data.test_regression("test_single_snapshot",res,new_res_array_size()));

    }

BOOST_AUTO_TEST_SUITE_END()


