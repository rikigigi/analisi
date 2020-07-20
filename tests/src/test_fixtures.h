#ifndef TEST_FIXTURES_H
#define TEST_FIXTURES_H

#include <string>
#include <fstream>
#include <boost/test/included/unit_test.hpp>
#include "config.h"
#include "traiettoria.h"
#include "readlog.h"

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

struct LogSetup{
    LogSetup(): traj{path.path+"data/gk_integral.dat"} {}
    const TestPath path;
    ReadLog<double> traj;
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
        if (std::isnan(a)) {
            BOOST_TEST_WARN("Found NaN while testing 'a' array!");
        }
        if (std::isnan(b)) {
            BOOST_TEST_WARN("Found NaN while testing 'b' array!");
        }
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
                double ave_diff=0.0,ave_violated_diff=0.0;
                double max_diff=fabs(a[i]-b[i]);
                double min_diff=fabs(a[i]-b[i]);
                unsigned n_violated=0,n=0;
                for (size_t j=0;j<size;j++){
                    double t=fabs(a[j]-b[j]);
                    if (! is_same(a[j],b[j])) {
                        ave_violated_diff+=(t-ave_violated_diff)/(++n_violated);
                    }
                    //BOOST_TEST_MESSAGE(a[j] << " " << b[j]<<" " <<a[j]-b[j] );
                    if (std::isnan(t)){
                        BOOST_TEST_WARN("found NaN ");
                    } else {
                        ave_diff+=(t - ave_diff)/(++n);
                    }
                    if (t<min_diff) min_diff=t;
                    if (t>max_diff) max_diff=t;
                }
                BOOST_TEST_ERROR("\naverage difference: "<<ave_diff <<"\navereage difference when !is_same(a,b): "<<ave_violated_diff<<"\nnumber of !is_same(a,b): "<<n_violated<<"\nmin,max diff :"<<min_diff<<" "<<max_diff<<"\n");
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

#endif
