#include "cg.h"
#include "config.h"
#include <iostream>

#ifdef HAVEeigen3EigenDense
#include <eigen3/Eigen/Dense>
#else
#include <Eigen/Dense>
#endif

template <int D>
class QuadraticForm : public Function <Eigen::Matrix<double,D,1>, double> {
  public:
    virtual double operator () (const Eigen::Matrix<double,D,1> & x) final {
        return (0.5*x.transpose()*A*x-b.transpose()*x)(0,0);
    }

    virtual Eigen::Matrix<double,D,1> deriv (const Eigen::Matrix<double,D,1> & x) final {
        return (A*x - b).eval();
    }

    void set_A_b(const Eigen::Matrix<double,D,D> & A_, const Eigen::Matrix<double,D,1> &b_) {
        A=A_;
        b=b_;
        x0=A.colPivHouseholderQr().solve(b);
        std::cout<<"sol:" << x0<<"\n";
    }
    double get_solution_distance(const Eigen::Matrix<double,D,1> & x) {
        return (x-x0).norm();
    }
protected:
    Eigen::Matrix<double,D,D> A;
    Eigen::Matrix<double,D,1> b,x0;
};


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
                for (unsigned int i2=i1+1;i2<N;i2++) { //hessian is symmetric
                    // note: here I calculate also the diagonal term, that has a different form with one more sum
                    // If one makes the calculation,
                    // the diagonal term is simply the sum of the other terms in the row. Because the matrix is symmetric, at the end of the day,
                    // for every out of diagonal term (i1,i2) that I calculate, I add the same quantity to the diagonal terms (i1,i1) and (i2,i2)
                    // the factor in front of the 3x3 matrix is the same because of the symmetry of the potential

                    //get pbc distance
                    Eigen::Matrix<double,DIM,1> dx;
                    double r2;
                    if (want_pbc()){
                        dx=pbc(x.template segment<DIM>(i*DIM,DIM),x.template segment<DIM>(j*DIM,DIM),r2);
                    }
                    else {
                        dx=(x.template segment<DIM>(i*DIM,DIM)-x.template segment<DIM>(j*DIM,DIM));
                        r2=dx.squaredNorm();
                    }

                    double f2=4*pair_deriv2_r2(r2); //second derivative with respect to r^2_ij
                    double f1=2*pair_deriv_r2(r2); //first derivative with respect to r^2_ij
                    double sub[DIM*DIM]={0.0};
                    for (unsigned int d0=0;d0<DIM;d0++) {
                        // diagonal only term
                        sub[d0*DIM+d0]+=f1;
                        for (unsigned int d1=d0;d1<DIM;d1++) {
                            sub[d0*DIM+d1]+=f2*dx(d0)*dx(d1);
                        }
                    }

                    //add the matrix contribution to diagonal and off diagonal term (matrix is symmetric)
                    for (unsigned int d0=0;d0<DIM;d0++) {
                        //off diagonal
                        res(i1*DIM+d0,i2*DIM+d0)=sub[d0*DIM+d0];
                        //off diagonal symmetric
                        res(i2*DIM+d0,i1*DIM+d0)=sub[d0*DIM+d0];
                        //diagonal 1
                        res(i1*DIM+d0,i1*DIM+d0)+=sub[d0*DIM+d0];
                        //diagonal 2
                        res(i2*DIM+d0,i2*DIM+d0)+=sub[d0*DIM+d0];
                        for (unsigned int d1=d0+1;d1<DIM;d1++) {
                            // sub[d0*DIM+d1]=sub[d1*DIM+d0];

                            //off diagonal
                            res(i1*DIM+d0,i2*DIM+d1)=sub[d1*DIM+d0];
                            res(i1*DIM+d1,i2*DIM+d0)=sub[d1*DIM+d0];

                            //off diagonal symmetric
                            res(i2*DIM+d0,i1*DIM+d1)=sub[d1*DIM+d0];
                            res(i2*DIM+d1,i1*DIM+d0)=sub[d1*DIM+d0];

                            //diagonal 1
                            res(i1*DIM+d0,i1*DIM+d1)=sub[d1*DIM+d0];
                            res(i1*DIM+d1,i1*DIM+d0)=sub[d1*DIM+d0];
                            //diagonal 2
                            res(i2*DIM+d0,i2*DIM+d1)=sub[d1*DIM+d0];
                            res(i2*DIM+d1,i2*DIM+d0)=sub[d1*DIM+d0];
                        }
                    }
                }
            }
            return res;
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
    /// test con una forma quadratica:

    Eigen::Matrix4d m=Eigen::Matrix4d::Random();
    m = Eigen::Matrix4d::Random();
    Eigen::Matrix<double,4,1> xq=2*Eigen::Matrix<double,4,1>::Random(), b=Eigen::Matrix<double,4,1>::Random();

    QuadraticForm<4> testq;
    testq.set_A_b(0.5*(m.transpose()+m)+Eigen::Matrix4d::Identity()*4,b);
    ParabolaLineMinimization <Eigen::Matrix<double,4,1>,double,QuadraticForm<4> > lineMinq(5,4,20,3);
    Cg<Eigen::Matrix<double,4,1>,double,QuadraticForm<4> > testcgq(testq,xq,testq(xq),lineMinq,8);
    std::cout << 0 << " " << testq.get_solution_distance(testcgq.get_x())<< "\n";
    for (unsigned int i=0;i<4;i++){
        testcgq.iteration();
        std::cout << i+1 << " " << testq.get_solution_distance(testcgq.get_x())<< "\n";
    }

    std::cout << "======================\n";

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
        if (!testcg.iteration())
            break;
    }
    std::cout <<  "Final: " << testcg.get_fx()<< "\n";

}
