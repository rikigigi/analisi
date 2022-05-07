//#define BOOST_TEST_MODULE calc_buffer_tests
#include <boost/test/included/unit_test.hpp>



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
BOOST_AUTO_TEST_CASE(test_buffer_discard){
    CalcBuffer<int> test(29,1);
    FakeCalc calculator;
    for (int i=35;i<97;++i){
        for (int j=i-1;j<i+33;++j){
            if (j==i+3)
                test.discard(42);
            int * a=test.buffer_calc(calculator,i),
                * b=test.buffer_calc(calculator,j);
            BOOST_TEST(calculator.check(i,a));
            BOOST_TEST(calculator.check(j,b));
        }
    }
    BOOST_TEST_MESSAGE("times used: "<<calculator.n_check<<"; time calculated: "<< calculator.n_eval);
}
BOOST_AUTO_TEST_CASE(test_buffer_discard_all){
    CalcBuffer<int> test(29,1);
    FakeCalc calculator;
    for (int i=35;i<97;++i){
        for (int j=i-1;j<i+33;++j){
            if (j==42)
                test.discard();
            int * a=test.buffer_calc(calculator,i),
                * b=test.buffer_calc(calculator,j);
            BOOST_TEST(calculator.check(i,a));
            BOOST_TEST(calculator.check(j,b));
        }
    }
    BOOST_TEST_MESSAGE("times used: "<<calculator.n_check<<"; time calculated: "<< calculator.n_eval);
}

BOOST_AUTO_TEST_CASE(test_buffer_hit_miss) {

    int t2_stop=20,t2_start=0;
    int t1_start=0,t1_stop=10;
    int overlap=t2_start< t1_stop ?t1_stop- t2_start : ( t1_start < t2_stop ?t2_stop -t1_start : 0 );
    for (int istart=0;istart<10;istart++){
        int start=istart*12345;
        CalcBuffer<int> test(t2_stop-t2_start+1,1);
        FakeCalc calculator;
        for (int t1=t1_start;t1<t1_stop;++t1){
            for (int t2=t2_start;t2<t2_stop;++t2){

                int * a=test.buffer_calc(calculator,start+t1),
                    * b=test.buffer_calc(calculator,start+t2);
                BOOST_TEST(calculator.check(start+t1,a));
                BOOST_TEST(calculator.check(start+t2,b));
                if (t2==t2_stop-1){
                    test.discard(start+t1);
                }
            }
        }
        BOOST_TEST(test.get_miss()==calculator.n_eval);
        BOOST_TEST(test.get_hit()==calculator.n_check-calculator.n_eval);
        BOOST_TEST(test.get_miss()==t1_stop-t1_start + t2_stop-t2_start - overlap);
    }
}


