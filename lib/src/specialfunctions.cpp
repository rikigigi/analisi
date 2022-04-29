#include "specialfunctions.h"

using namespace SpecialFunctions;




void print_realSpericalHarmonics_coeff(int l,std::ostream & out) {
    for (int j=0;j<=l;++j) {
        for (int i=-j;i<=j;++i){
            out << "(l,m) = ("<<j<<","<<i <<"): "<<realSpericalHarmonics_coeff<double>(j,i,j-abs(i)) <<std::endl;
        }
    }
}

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

