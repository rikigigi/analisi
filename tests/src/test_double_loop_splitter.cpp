//#define BOOST_TEST_MODULE twoloopsplit_tests
#include <boost/test/included/unit_test.hpp>

//double loop splitter test
#include <vector>
#include <algorithm>

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

