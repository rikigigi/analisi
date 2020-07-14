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
    const std::string path;
};

struct TrajSetup{
    TrajSetup(): traj{path.path+"data/lammps.bin"}{
	traj.imposta_dimensione_finestra_accesso(150);
	traj.imposta_inizio_accesso(0);
    }
    const TestPath path;
    Traiettoria traj;

};

struct DataRegression{
    DataRegression():path{base_path.path + "cpp_regression_data/"}, max_double_relative_error{1e-10} {}
    const TestPath base_path;
    std::string path;
    double max_double_relative_error;
    void fuzz(double *a,size_t size) {
        static double s=0.76;
        for (size_t i=0;i<size;++i) {
            s=3.78*s*(1.0-s);
            a[i]=s;
        }
    }
    bool is_same(double a, double b) {
        double max=fabs(a)>fabs(b) ? a : b;
        double min=fabs(a)>fabs(b) ? b : a;
        if (a==0 and b==0) return true;
        if ((max-min)/max > max_double_relative_error) return false;
        return true;
    }
    bool is_same(double *a, double *b, size_t size) {
        for (size_t i=0;i<size;++i){
            if (! is_same(a[i],b[i])) {
                BOOST_TEST_MESSAGE("data differs: " <<a[i] << " " << b<<" " <<a[i]-b[i] );
                for (size_t j=0;j<size;j++){
                    BOOST_TEST_MESSAGE(a[j] << " " << b[j]<<" " <<a[j]-b[j] );
                }
                return false;
            }
        }
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
    return is_same(data_fs,data,size);
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
BOOST_AUTO_TEST_CASE(test_calcola)\
{\
    sh.calcola(0);\
    BOOST_TEST(data.test_regression("test_calcola",sh.accesso_lista(),sh.lunghezza()));\
}\
BOOST_AUTO_TEST_SUITE_END()

TESTS(ShFix10_1)
TESTS(ShFix10_2)
TESTS(ShFix10_3)

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


//calc_buffer test -- eheh, this took 48 hours to work...

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
            int * a=test.buffer_calc(calculator,i),
                * b=test.buffer_calc(calculator,j);
            BOOST_TEST(calculator.check(i,a));
            BOOST_TEST(calculator.check(j,b));
        }
    }
    BOOST_TEST_MESSAGE("times used: "<<calculator.n_check<<"; time calculated: "<< calculator.n_eval);
}

//double loop splitter test

#include "twoloopsplit.h"
template <class T>
bool test_double_loop(T nworkers_, const T size1_, const T skip1_,const T block1_, const T size2_, const T skip2_, const T block2_) {

    TwoLoopSplit test(nworkers_,size1_,skip1_,block1_,size2_,skip2_,block2_);
    std::vector<std::pair<T,T> >elements,missing;
    for (T iw=0;iw<nworkers_;++iw) {
        bool end=false;
        while (!end){
            T idx1,idx2;
            test.get_next_idx_pair(iw,idx1,idx2,end);
            elements.push_back({idx1,idx2});
        }
    }
    T counter=0;
    for (T idx1=0;idx1<size1_;idx1+=skip1_)
        for (T idx2=idx1;idx2<size2_+idx1;idx2+=skip2_) {
            counter++;
            std::pair<T,T> e{idx1,idx2};
            auto res=std::find(elements.begin(),elements.end(),e);
            if (res == elements.end()) {
                missing.push_back(e);
            } else {
                elements.erase(res);
            }
        }
    if (0 != elements.size() || missing.size()>0) {
        return false;
    } else {
        return true;
    }
}

BOOST_AUTO_TEST_CASE(twoloopsplit) {
    BOOST_TEST(test_double_loop<size_t>(2,10,1,5,10,1,5));
    BOOST_TEST(test_double_loop<size_t>(2,10,2,5,10,1,5));
    BOOST_TEST(test_double_loop<size_t>(2,10,1,5,10,2,5));
    BOOST_TEST(test_double_loop<size_t>(2,10,2,5,10,2,5));
    BOOST_TEST(test_double_loop<size_t>(3,17,3,5,29,2,7));
    BOOST_TEST(test_double_loop<size_t>(7,200,3,5,70,2,9));
    BOOST_TEST(test_double_loop<size_t>(7,200,3,5,70,1,9));
}

