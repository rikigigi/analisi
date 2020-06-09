#include "associatedlegendrepoly.h"


double mixall(int l, int m, double val, double * cheby){
    if (m>=0){
        return val*realSpericalHarmonics_coeff<double>(l,m,l-m+1)*cheby[2*m];
    } else {
        return val*realSpericalHarmonics_coeff<double>(l,m,l-abs(m)+1)*cheby[2*abs(m)+1];
    }
}

int main() {

    //init data
    double theta,phi;
    std::cin >> theta >> phi;
    double cost=cos(theta),sinp=sin(phi),cosp=cos(phi);
    constexpr int lmax=10;
    double cheby[2*(lmax+1)];
    MultiVal<lmax,double> res;

    //legendre polynomials
    AssociatedLegendrePoly<lmax,0,0,double>::calc(cost,res);
    res.copy_mplus_mminus();
    res.print(std::cout);
    //chebyshev recursive angle multiplication
    chebyshevRecursiveAngleMultiplication<lmax>(cosp,sinp,cheby);
    print_chebyshevRecursiveAngleMultiplication<lmax>(cheby,std::cout);
    //mix all toghether with the correct coefficients to get spherical harmonics (yee)
    print_realSpericalHarmonics_coeff(lmax,std::cout);
    res.apply_f(mixall,cheby);


    res.print(std::cout);


}
