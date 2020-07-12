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
    DataRegression():path{base_path.path + "cpp_regression_data/"}, max_double_relative_error{1e-10} {}
    TestPath base_path;
    std::string path;
    double max_double_relative_error;
    bool is_same(double a, double b) {
        double max=fabs(a)>fabs(b) ? a : b;
        double min=fabs(a)>fabs(b) ? b : a;
        if (a==0 and b==0) return true;
        if ((max-min)/max > max_double_relative_error) return false;
        return true;
    }
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
        if (! is_same(data_fs[i],data[i])) {
            BOOST_TEST_MESSAGE("data differs: " <<data_fs[i] << " " << data[i]<<" " <<data_fs[i]-data[i] );
            for (size_t j=0;j<size;j++){
                BOOST_TEST_MESSAGE(data_fs[j] << " " << data[j]<<" " <<data_fs[j]-data[j] );
            }
            return false;
        }
	}
	return true;
    }
};

template <int l,int NTH>
struct ShFixture{
    ShFixture() : nbin{4}, natoms{traj.traj.get_natoms()}, ntypes{traj.traj.get_ntypes()}, sh{&(traj.traj), 0.5, 3.0, nbin, 17, NTH, 13, false} {
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

typedef ShFixture<10,1> ShFix10_1 ;
typedef ShFixture<10,2> ShFix10_2 ;
typedef ShFixture<10,3> ShFix10_3 ;

#define TESTS(T)\
BOOST_FIXTURE_TEST_SUITE(sh ## T, T )\
BOOST_AUTO_TEST_CASE(test_single_snapshot)\
{\
    double * res = new_res_array();\
    calc(0,res);\
    BOOST_TEST(data.test_regression("test_single_snapshot",res,new_res_array_size()));\
    delete [] res;\
}\
BOOST_AUTO_TEST_CASE(test_calcola)\
{\
    sh.calcola(0);\
    BOOST_TEST(data.test_regression("test_calcola",sh.accesso_lista(),sh.lunghezza()));\
}\
BOOST_AUTO_TEST_SUITE_END()

TESTS(ShFix10_1)
TESTS(ShFix10_2)
TESTS(ShFix10_3)


#include "calc_buffer.h"

struct FakeCalc{
    FakeCalc(): n_eval{0},n_check{0} {}
    void calc(int k, int* res) {
        *res=k*42;
        n_eval++;
    }
    bool check(int k, int* res) {
        n_check++;
        return *res == k*42;
    }
    int n_eval,n_check;
};



BOOST_AUTO_TEST_CASE(test_buffer){
    CalcBuffer<int> test(29,1);
    FakeCalc calculator;
    for (int i=35;i<97;++i){
        for (int j=i-1;j<i+33;++j){
            BOOST_TEST(calculator.check(i,test.buffer_calc(calculator,i)));
            BOOST_TEST(calculator.check(j,test.buffer_calc(calculator,j)));
        }
    }
    BOOST_TEST_MESSAGE("times used: "<<calculator.n_check<<"; time calculated: "<< calculator.n_eval);
}

