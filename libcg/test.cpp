#include "cg.h"
#include "config.h"
#include <iostream>

#ifdef HAVEeigen3EigenDense
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif


template <int N,int DIM,int flags >
class MultiPair : public Function <Eigen::Matrix<double,N*DIM,1>, double, Eigen::Matrix<double,N*DIM,N*DIM> > {
public:
    virtual double operator() (const Eigen::Matrix<double,N*DIM,1> & x ) final {
        double res=0.0;
        for (unsigned int i=0;i<N;i++) {
            for (unsigned int j=i+1;j<N;j++) {
                double r2;
                Eigen::Matrix<double,DIM,1> dx;
                if (want_pbc()){
                    dx=pbc(x.template segment<DIM>(i*DIM,DIM),x.template segment<DIM>(j*DIM,DIM),r2);
                }else{

                    dx=(x.template segment<DIM>(i*DIM,DIM)-x.template segment<DIM>(j*DIM,DIM));
                    r2=dx.squaredNorm();
                }
                res+=pair(r2);
            }
        }
        return res;
    }
    virtual Eigen::Matrix<double,N*DIM,1> deriv(const Eigen::Matrix<double,N*DIM,1> & x) final {
        Eigen::Matrix<double,N*DIM,1> res=Eigen::Matrix<double,N*DIM,1>::Zero();
        for (unsigned int i=0;i<N;i++) {
            for (unsigned int j=0;j<N;j++) {
                if (i==j)
                    continue;
                double r2;
                Eigen::Matrix<double,DIM,1> dx;
                if (want_pbc()){
                    dx=pbc(x.template segment<DIM>(i*DIM,DIM),x.template segment<DIM>(j*DIM,DIM),r2);
                }
                else {
                    dx=(x.template segment<DIM>(i*DIM,DIM)-x.template segment<DIM>(j*DIM,DIM));
                    r2=dx.squaredNorm();
                }
                double dp_dr=pair_deriv_r2(r2);
                for (unsigned int i0=0;i0<DIM;i0++) {
                    res(i*DIM+i0)+=dp_dr*dx(i0)*2;
                }
            }
        }
        return res;
    }

    virtual Eigen::Matrix<double,N*DIM,N*DIM> hessian(const Eigen::Matrix<double,N*DIM,1> & x) final {
        Eigen::Matrix<double,N*DIM,N*DIM> res;
        if (has_deriv2()) {
            for (unsigned int i1=0;i1<N;i1++) {
                for (unsigned int i2=0;i2<N;i2++) {
                    for (unsigned int d0=0;d0<DIM;d0++) {
                        for (unsigned int d1=0;d1<DIM;d1++) {

                        }
                    }
                }
            }
        } else {
            throw std::runtime_error("Not implemented!");
        }

    }
    virtual void init_pbc(const Eigen::Matrix<double,DIM,DIM> &t){
        T=t;
        Tinv=t.inverse();
    }

private:
    constexpr bool want_pbc()       {return (flags & 1) == 1;}
    constexpr bool has_deriv2()     {return (flags & 2) == 2;}
    constexpr bool newton_forces()  {return (flags & 4) == 4;}
    Eigen::Matrix<double,DIM,DIM> T,Tinv;
    template <typename D1,typename D2>
    inline
    Eigen::Matrix<typename D1::Scalar,DIM,1>  pbc(const Eigen::MatrixBase<D1> &x1, const Eigen::MatrixBase<D2> &x0,
                       typename D1::Scalar & r2min) {
        Eigen::Matrix<typename D1::Scalar,DIM,1> u1,u0,dx;
        u1=Tinv*x1;
        u0=Tinv*x0;
        u1=u1-Eigen::floor(u1.array()).matrix();
        u0=u0-Eigen::floor(u0.array()).matrix();
        dx=T*(u1-u0);
        r2min=dx.squaredNorm();
        for (unsigned int i=0;i<DIM;i++) {
            for (int i0=-1;i0<2;i0+=2){
                Eigen::Matrix<typename D1::Scalar,DIM,1> dxn=dx+double(i0)*T.col(i);
                typename D1::Scalar r2=dxn.squaredNorm();
                if (r2<r2min){
                    r2min=r2;
                    dx=dxn;
                }
            }
        }
        return dx;
    }
    virtual double pair       (const double & r2)=0;
    virtual double pair_deriv_r2 (const double & r2)=0;
    virtual double pair_deriv2_r2(const double & r2)=0;
};

template <unsigned int N,unsigned int DIM>
class LJPair : public MultiPair<N,DIM,1 | 2 > {
protected:
    virtual double pair       (const double & r2) override {
        double r6=r2*r2*r2;
        return 1.0/(r6*r6)-1.0/r6;
    }
    virtual double pair_deriv_r2 (const double & r2) override {
        double r6=r2*r2*r2;
        return -6.0/(r6*r6*r2)+3.0/(r6*r2);
    }

    virtual double pair_deriv2_r2(const double & r2) override{
        double r4=r2*r2;
        double r6=r4*r2;
        double r12=r6*r6;
        return 7.0*6.0/(r12*r4)-4.0*3.0/(r6*r4);
    }

};

int main() {
    ///... testare il cj

    LJPair<25,3> test;
    Eigen::Matrix3d cel;
    cel << 6.0 , 0.0 , 0.0
        , 0.0 , 6.0 , 0.0
        , 0.0 , 0.0 , 6.0;
    test.init_pbc( cel);
    Eigen::Matrix<double,25*3,1> x=Eigen::Matrix<double,25*3,1>::Random()*6.0;

    ParabolaLineMinimization <Eigen::Matrix<double,25*3,1>,double,LJPair<25,3> > lineMin(0.02,0.1,4,3);
    Cg<Eigen::Matrix<double,25*3,1>,double,LJPair<25,3> > testcg(test,x,test(x),lineMin,8);
    for (unsigned int i=0;i<2000;i++) {
        std::cout << i << " " << testcg.get_fx()<< "\n";
        testcg.iteration();
    }
    std::cout <<  "Final: " << testcg.get_fx()<< "\n";

}
