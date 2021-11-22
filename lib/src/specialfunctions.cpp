#include "specialfunctions.h"

using namespace SpecialFunctions;



void test() {

    //init data
    double theta,phi;
    std::cin >> theta >> phi;
    double cost=cos(theta),sinp=sin(phi),cosp=cos(phi);
    constexpr int lmax=10;
    double cheby[2*(lmax+1)];
    double result[(lmax+1)*(lmax+1)];

    {
        std::cout << "+++++++++++++++++++++++"<<std::endl;
        SphericalHarmonics<lmax,double,false> test(cost,sinp,cosp,cheby,nullptr);
        test.calc();
        test.get_val().print(std::cout);
    }


    {
        std::cout << "+++++++++++++++++++++++"<<std::endl;
        SphericalHarmonics<lmax,double,true> test(cost,sinp,cosp,cheby,result);
        test.calc();
        test.get_val().print(std::cout);
    }
}

