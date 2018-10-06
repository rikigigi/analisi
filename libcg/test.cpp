#include "cg.h"
#include "config.h"

#ifdef HAVEeigen3EigenDense
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif


template <unsigned int N,unsigned int DIM,int flags >
class MultiPair : public Function <Eigen::Matrix<double,N*DIM,1>, double, Eigen::Matrix<double,N*DIM,N*DIM> > {
public:
    virtual double operator() (const Eigen::Matrix<double,N*DIM,1> & x ) {
        double res=0.0;
        for (unsigned int i=0;i<N;i++) {
            for (unsigned int j=i+1;j<N;j++) {
                ///TODO: aggiungere pbc
                 Eigen::Matrix<double,DIM,1> dx=(x.block<DIM>(i*DIM)-x.block<DIM>(j*DIM));
                 double r2=dx.transpose()*dx;
                 res+=pair(r2);
            }
        }
        return res;
    }
    virtual Eigen::Matrix<double,N*DIM,1> deriv(const Eigen::Matrix<double,N*DIM,1> & x) {
        Eigen::Matrix<double,N*DIM,1> res(0.0);
        for (unsigned int i=0;i<N;i++) {
            for (unsigned int j=i+1;j<N;j++) {
                ///TODO: aggiungere pbc
                 Eigen::Matrix<double,DIM,1> dx=(x.block<DIM>(i*DIM)-x.block<DIM>(j*DIM));
                 double r2=dx.transpose()*dx;
                 double dp_dr=pair_deriv_r2(r2);
                 for (unsigned int i0=0;i0<DIM;i0++) {
                    res(i*DIM+i0)+=dp_dr*dx(i0)*2;
                 }
            }
        }
        return res;
    }

    virtual Eigen::Matrix<double,N*DIM,N*DIM> hessian(const Eigen::Matrix<double,N*DIM,1> & x) {

    }

protected:
    constexpr bool want_pair_r()       {return (flags & 1) == 1;}
    constexpr bool want_pair_deriv_r() {return (flags & 4) == 4;}

    virtual double pair       (const double & r2)=0;
    virtual double pair_deriv_r2 (const double & r2)=0;
    virtual double pair_deriv2_r2(const double & r2);
};

template <unsigned int N,unsigned int DIM>
class LJPair : public MultiPair<N,DIM,0> {
protected:
    virtual double pair       (const double & r2) {
        double r6=r2*r2*r2;
        return 1.0/(r6*r6)-1.0/r6;
    }
    virtual double pair_deriv_r2 (const double & r2) {
        double r6=r2*r2*r2;
        return -6.0/(r6*r6*r2)+3.0/(r6*r2);
    }

    virtual double pair_deriv2_r2(const double & r2){
        double r4=r2*r2;
        double r6=r4*r2;
        double r12=r6*r6;
        return 7.0*6.0/(r12*r4)-4.0*3.0/(r6*r4);
    }

};

int main() {
    ///... testare il cj
}
