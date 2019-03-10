#include "cg.h"
#include "config.h"
#include "../rnd.h"
#include "../lib/json.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <functional>

#ifdef HAVEeigen3EigenDense
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#else
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#endif

template<typename Derived>
std::istream & operator >>
(std::istream & s,
 Eigen::MatrixBase<Derived> & m)
{
    for (int i = 0; i < m.rows(); ++i)
    for (int j = 0; j < m.cols(); j++)
      s >> m(i,j);

  return s;
}


#define NATOMS -1
#define DIMs 3
using json = nlohmann::json;

template <int D>
class QuadraticForm : public Function <
        Eigen::Matrix<double,D,1>,
        double,
        Eigen::Matrix<double,D,D>,
        Eigen::Ref<Eigen::Matrix<double,D,1> >
      > {
  public:
    virtual double operator () (const Eigen::Matrix<double,D,1> & x) final {
        return (0.5*x.transpose()*A*x-b.transpose()*x)(0,0);
    }

    virtual void deriv (const Eigen::Matrix<double,D,1> & x, Eigen::Ref<Eigen::Matrix<double,D,1> > res) final {
        res= (A*x - b).eval();
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


constexpr int eigen_matrix_dim(const int N,const int DIM) {
    return N>0? (N*DIM):(Eigen::Dynamic);
}

template <int N,int DIM,int flags >
class MultiPair : public Function <
        Eigen::Ref< const Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> >,
        double,
        Eigen::Ref <Eigen::Matrix<double,eigen_matrix_dim(N,DIM),eigen_matrix_dim(N,DIM)> >,
        Eigen::Ref <Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> >,
        Eigen::Ref<Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> >,
        Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1>
> {
public:
    MultiPair(double cut) : cutoff2(cut*cut),cutoff(cut){

    }
    virtual double operator() (const Eigen::Ref< const Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > & x ) final {
        double res=0.0;
        for (unsigned int i=0;i<x.rows()/DIM;i++) {
            for (unsigned int j=i+1;j<x.rows()/DIM;j++) {
                double r2;
                Eigen::Matrix<double,DIM,1> dx;
                if (want_pbc()){
                    dx=pbc(x.template segment<DIM>(i*DIM,DIM),x.template segment<DIM>(j*DIM,DIM),r2);
                }else{

                    dx=(x.template segment<DIM>(i*DIM,DIM)-x.template segment<DIM>(j*DIM,DIM));
                    r2=dx.squaredNorm();
                }
                if (r2>cutoff2) continue;
                res+=pair(r2);
            }
        }
        return res;
    }
    virtual void deriv(const Eigen::Ref< const Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > & x,Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > res) final {
        res.setZero();
        zero_virial();
        zero_energy();
        for (unsigned int i=0;i<x.rows()/DIM;i++) {
            for (unsigned int j=i+1;j<x.rows()/DIM;j++) {
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
                if (r2>cutoff2) continue;
                double dp_dr=2*pair_deriv_r2(r2);
                virial_pair(dp_dr,dx.data());
                energy_pair(pair(r2));
                for (unsigned int i0=0;i0<DIM;i0++) {
                    res(i*DIM+i0)+=dp_dr*dx(i0);
                    res(j*DIM+i0)-=dp_dr*dx(i0);
                }
            }
        }
    }

    virtual void hessian_deriv(const Eigen::Ref< const Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > & x,
                               Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > res1,
                               Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),eigen_matrix_dim(N,DIM)> > res
                               ) final {
        res.setZero();
        res1.setZero();
        zero_virial();
        zero_energy();
//        Eigen::Matrix<double,N*DIM,N*DIM> res1;
        if (has_deriv2()) {
            for (unsigned int i1=0;i1<x.rows()/DIM;i1++) {
                for (unsigned int i2=i1+1;i2<x.rows()/DIM;i2++) { //hessian is symmetric
                    // note: here I calculate also the diagonal term, that has a different form with one more sum
                    // If one makes the calculation,
                    // the diagonal term is simply the sum of the other terms in the row. Because the matrix is symmetric, at the end of the day,
                    // for every out of diagonal term (i1,i2) that I calculate, I add the same quantity to the diagonal terms (i1,i1) and (i2,i2)
                    // the factor in front of the 3x3 matrix is the same because of the symmetry of the potential

                    //get pbc distance
                    Eigen::Matrix<double,DIM,1> dx;
                    double r2;
                    if (want_pbc()){
                        dx=pbc(x.template segment<DIM>(i1*DIM,DIM),x.template segment<DIM>(i2*DIM,DIM),r2);
                    }
                    else {
                        dx=(x.template segment<DIM>(i1*DIM,DIM)-x.template segment<DIM>(i2*DIM,DIM));
                        r2=dx.squaredNorm();
                    }

                    if (r2>cutoff2) continue;
                    double f2=4*pair_deriv2_r2(r2); //second derivative with respect to r^2_ij
                    double f1=2*pair_deriv_r2(r2); //first derivative with respect to r^2_ij
                    double sub[DIM*DIM]={0.0};

                    virial_pair(f1,dx.data());
                    energy_pair(pair(r2));

                    for (unsigned int d0=0;d0<DIM;d0++) {
                        // calculate also force (it is free, at this point)
                        res1(i1*DIM+d0)+=f1*dx(d0);
                        res1(i2*DIM+d0)-=f1*dx(d0);

                        // diagonal only term
                        sub[d0*DIM+d0]+=f1;
                        for (unsigned int d1=d0;d1<DIM;d1++) {
                            sub[d0*DIM+d1]+=f2*dx(d0)*dx(d1);
                        }
                    }
                    /*
                    for (unsigned int d0=0;d0<DIM;d0++) {
                        for (unsigned int d1=0;d1<DIM;d1++)
                            std::cout << sub[d0*DIM+d1]<<" ";
                        std::cout << "\n";
                    }
                    std::cout << "\n\n";
*/
                    //add the matrix contribution to diagonal and off diagonal term (matrix is symmetric)
                    for (unsigned int d0=0;d0<DIM;d0++) {
                        //off diagonal
                        res(i1*DIM+d0,i2*DIM+d0)=-sub[d0*DIM+d0];
                        //off diagonal symmetric
                        res(i2*DIM+d0,i1*DIM+d0)=-sub[d0*DIM+d0];
                        //diagonal 1
                        res(i1*DIM+d0,i1*DIM+d0)+=sub[d0*DIM+d0];
                        //diagonal 2
                        res(i2*DIM+d0,i2*DIM+d0)+=sub[d0*DIM+d0];
                        for (unsigned int d1=0;d1<d0;d1++) {
                            // sub[d0*DIM+d1]=sub[d1*DIM+d0];

                            //off diagonal
                            res(i1*DIM+d0,i2*DIM+d1)=-sub[d1*DIM+d0];
                            res(i1*DIM+d1,i2*DIM+d0)=-sub[d1*DIM+d0];

                            //off diagonal symmetric
                            res(i2*DIM+d0,i1*DIM+d1)=-sub[d1*DIM+d0];
                            res(i2*DIM+d1,i1*DIM+d0)=-sub[d1*DIM+d0];

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
            res1=res1;
        } else {
            throw std::runtime_error("Not implemented!");
        }

    }
    virtual void init_pbc(const Eigen::Matrix<double,DIM,DIM> &t){
        T=t;
        Tinv=t.inverse();
    }

    bool check_forces(const Eigen::Ref< const Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > & x, const double & dx_over_x, const double & max_error=0.001 ) {
        Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> force,x_;
        force.resize(x.rows(),1);
        x_.resize(x.rows(),1);
        deriv(x,force);
        bool res=true;
        for (unsigned int i=0;i<x.rows();i++) {
            double dvdr=0.0;
            x_=x;
            x_(i)=x(i)+dx_over_x;
            double vp=operator ()(x_);
            x_(i)=x(i)-dx_over_x;
            double vm=operator ()(x_);
            dvdr=(vp-vm)/(2*dx_over_x);
            std::cerr << "Index "<< i<< ": \tnumerical = " << dvdr << " calculated = "<< force(i)<<"\n";
            if (fabs((dvdr-force(i)))/force(i)>max_error) {
                res=false;
            }
        }
        return res;
    }
    bool check_hessian_forces(const Eigen::Ref< const Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > & x, const double & dx_over_x, const double & max_error=0.001 ) {
        bool res=false;
        Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> force,x_;
        force.resize(x.rows(),1);
        x_.resize(x.rows(),1);
        Eigen::Matrix<double,eigen_matrix_dim(N,DIM),eigen_matrix_dim(N,DIM)> H;
        H.resize(x.rows(),x.rows());
        hessian_deriv(x,force,H);
        for (unsigned int i=0;i<x.rows();i++) {
            {
                double dvdr=0.0;
                x_=x;
                x_(i)=x(i)+dx_over_x;
                double vp=operator ()(x_);
                x_(i)=x(i)-dx_over_x;
                double vm=operator ()(x_);
                dvdr=(vp-vm)/(2*dx_over_x);
                std::cerr << "Force "<< i<< ": \t numerical = " << dvdr << " calculated = "<< force(i)<<"\n";
                if (fabs((dvdr-force(i)))/force(i)>max_error) {
                    res=false;
                }
            }
            for (unsigned int j=0;j<x.rows();j++) {
                double dp,dm,d1p,d1m;
                x_=x;
                x_(j)=x_(j)+dx_over_x;
                x_(i)=x_(i)+dx_over_x;
                dp=operator ()(x_);
                x_=x;
                x_(j)=x_(j)+dx_over_x;
                x_(i)=x_(i)-dx_over_x;
                dm=operator ()(x_);
                d1p=(dp-dm)/(2*dx_over_x);

                x_=x;
                x_(j)=x_(j)-dx_over_x;
                x_(i)=x_(i)+dx_over_x;
                dp=operator ()(x_);
                x_=x;
                x_(j)=x_(j)-dx_over_x;
                x_(i)=x_(i)-dx_over_x;
                dm=operator ()(x_);
                d1m=(dp-dm)/(2*dx_over_x);

                double d2=(d1p-d1m)/(2*dx_over_x);

                std::cerr << i << " "<<j << " " << d2 << " " << H(i,j) << "\n";
                if (fabs((d2-H(i,j)))/H(i,j)>max_error) {
                    res=false;
                }
            }

        }

        std::cout << "\n\n\n"<<H<<"\n";
        return res;
    }

    void pbc_wrap(Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > x ,double L=0.0) {
        if (L==0.0) L=T(0,0);
        x=x-L*Eigen::floor(x.array()/L).matrix();
    }

    Eigen::Matrix<double,DIM,DIM> get_virial() {
        return virial;
    }

    double get_energy() {
        return energy;
    }

private:
    constexpr bool want_pbc()     const  {return (flags & 1) == 1;}
    constexpr bool has_deriv2()     const {return (flags & 2) == 2;}
    constexpr bool newton_forces() const  {return (flags & 4) == 4;}
    constexpr bool orthorombic_box() const {return (flags & 8) == 8;}
    Eigen::Matrix<double,DIM,DIM> T,Tinv,virial;
    double energy,cutoff2,cutoff;
    template <typename D1,typename D2>
    inline
    Eigen::Matrix<typename D1::Scalar,DIM,1>  pbc(const Eigen::MatrixBase<D1> &x1, const Eigen::MatrixBase<D2> &x0,
                       typename D1::Scalar & r2min) {
        Eigen::Matrix<typename D1::Scalar,DIM,1> u1,u0,dx;
        u1=Tinv*(x1-x0);
        u0=u1-Eigen::round(u1.array()).matrix();
        dx=T*u0;
        r2min=dx.squaredNorm();
        if (! orthorombic_box()){
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
        }
        return dx;
    }


    void zero_virial() {
        virial.setZero();
    }

    void zero_energy() {
        energy=0;
    }

    void energy_pair(double e_ij){
        energy+=e_ij;
    }

    void virial_pair(const double f_ij,double * const  x_ij) {
        for (unsigned int x=0;x<DIM;x++) {
            virial(x,x)+=f_ij*x_ij[x]*x_ij[x];
            for (unsigned int y=x+1;y<DIM;y++) {
                virial(x,y)+=f_ij*x_ij[x]*x_ij[y];
                virial(y,x)+=f_ij*x_ij[x]*x_ij[y];
            }
        }
    }

    virtual double pair       (const double & r2)=0;
    virtual double pair_deriv_r2 (const double & r2)=0;
    virtual double pair_deriv2_r2(const double & r2)=0;
};

template <int N,unsigned int DIM>
class LJPair : public MultiPair<N,DIM,1 | 2 | 8 > {
public:
    LJPair(double cut) : MultiPair<N,DIM,1 | 2 | 8 >(cut) {}
protected:
    virtual double pair       (const double & r2) override {
        double r6=r2*r2*r2;
        return 4.0/(r6*r6)-4.0/r6;
    }
    virtual double pair_deriv_r2 (const double & r2) override {
        double r6=r2*r2*r2;
        return 4.0*(-6.0/(r6*r6*r2)+3.0/(r6*r2));
    }

    virtual double pair_deriv2_r2(const double & r2) override{
        double r4=r2*r2;
        double r6=r4*r2;
        double r12=r6*r6;
        return 4.0*(7.0*6.0/(r12*r4)-4.0*3.0/(r6*r4));
    }

};

template <int N,unsigned int DIM,unsigned int FLAGS>
class Integrator {
public:
    Integrator (MultiPair<N,DIM,FLAGS> * p,unsigned int natoms) : p(p) {
        if (N<=0) {
            pos_m.resize(DIM*natoms,1);
            deriv_m.resize(DIM*natoms,1);
            hessian_m.resize(DIM*natoms,DIM*natoms);
        }
    }

    virtual void step(Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > pos ) {
        std::cerr << "Error: not implemented\n";
        abort();
    }
    virtual void step(Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > pos,Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > vel  ){
        std::cerr << "Error: not implemented\n";
        abort();
    }

    void dump(std::ostream & out, const Eigen::Ref< const Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > & x) {

        out << x.rows()/DIM <<"\n\n";
        for (unsigned int i=0;i<x.rows()/DIM;i++) {
            out  <<"1 ";
             for (unsigned int j=0;j<DIM;j++)
                out << x(i*DIM+j) << " ";
             out << "\n";
        }
       // out << "\n;

    }


protected:
    MultiPair<N,DIM,FLAGS> *p;
    double energy;
    Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> pos_m;
    Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> deriv_m;
    Eigen::Matrix<double,eigen_matrix_dim(N,DIM),eigen_matrix_dim(N,DIM)> hessian_m;
};
template <int N, unsigned int DIM,unsigned int FLAGS>
class IntegratorAcceleratedLangevin : public Integrator<N,DIM,FLAGS> {
public:
    static double regularizer(const double & e,double * const & p) {
        /* eq. (45) of Mazzola, Sorella, "Accelerated molecular dynamics for ab-initio Electronic Simulations""
         * e    --> lambda
         * p[0] --> epsilon
         * p[1] --> tau
         * p[2] --> delta
        */
        return 1.0/(1.0+exp((e-p[0])/p[1]))/p[2]+e*(1.0-1.0/(1.0+exp((e-p[0])/p[1])));
    }
    IntegratorAcceleratedLangevin (MultiPair<N,DIM,FLAGS> *p,
                                   Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > pos,
                                   double delta,
                                   double T,
                                   unsigned int natoms,
                                   bool accelerated=true
            ) : Integrator<N,DIM,FLAGS>(p,natoms), delta(delta), T(T),accelerated(accelerated),d(1.0) {
        pos_m=pos;
        first=true;
        c=sqrt(2*T*delta);


        eigenvalue.resize(DIM*natoms,1);
        z.resize(DIM*natoms,1);
        pos_p.resize(DIM*natoms,1);
        H.resize(DIM*natoms,DIM*natoms);
        H_inv.resize(DIM*natoms,DIM*natoms);
        eigenvectors.resize(DIM*natoms,DIM*natoms);

    }

    void set_T(double T_) {
        T=T_;
        c=sqrt(2*T*delta);
    }

    void set_d(double d_) {
        d=d_;
    }

    virtual void step(Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > pos) final {

        double regularizer_parameters[]={1.0,2/d,d};

        for (unsigned int i=0;i<z.rows();i++) {
            z(i)=normal_gauss();
        }

        if (accelerated){
            p->hessian_deriv(pos,deriv_m,H);
            //regularize the matrix: diagonalize it
            eigensolver.compute(H);
            eigenvalue=eigensolver.eigenvalues();
            eigenvectors=eigensolver.eigenvectors();
            //regularize the eigenvalues with regularizer function
            eigenvalue=eigenvalue.unaryExpr(std::bind(&regularizer,std::placeholders::_1,regularizer_parameters));
            H=eigenvectors.transpose()*eigenvalue.asDiagonal()*eigenvectors;
            H_inv=(eigenvectors.transpose()*( Eigen::pow(eigenvalue.array(),-1) ).matrix().asDiagonal()*eigenvectors);
            //trasforma le variabili secondo la matrice di covarianze
            z=( (Eigen::pow(eigenvalue.array(),-0.5) ).matrix().asDiagonal()*eigenvectors*z).eval();
        } else{
            p->deriv(pos,deriv_m);
        }
//ok fino a qui
        /*
        auto eigenvalues=H.eigenvalues();
        std::cout << eigenvalues<<"\n";
        //if (eigenvalues.col(0)[0]<0){
        //    H=H+eigenvalues(0)*1.1*Eigen::Matrix<double,N*DIM,N*DIM>::Identity();
        //}
        */
        //mette in z un vettore di variabili normali distribuite come N(0,1)


        if (accelerated ){
            if (!first)
                pos_p=pos-c*z-delta*H_inv*deriv_m-
                    H_inv*(hessian_m-H)*(pos_m-pos)/2.0;
            else
                pos_p=pos-c*z-H_inv*deriv_m;
            hessian_m=H;
        }
        else
            pos_p=pos-c*z -delta*deriv_m;
/*
        this->dump(std::cerr,deriv_m);
 //       this->dump(std::cerr,deriv_m);
        std::cerr <<"z:"<<z<< "\n\nH: " << H << "\n\nH_inv"<<H_inv<<"\n\n\n";
*/
        pos_m=pos;
        pos=pos_p;

        first=false;
    }

    virtual void step(Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > pos,Eigen::Ref < Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> > vel  ) {
        std::cerr << "Warning: called method for second order dynamics, but this is first order\n";
        step(pos);
    }

    void init_global_rnd(unsigned int seed,unsigned int thermalization_steps=10000) {
        set_SHR3_jsr((seed+1)*123);
        for (unsigned int i=0;i<thermalization_steps;i++) {
            rnd_shr3();
        }
        init_cmwc4096();
        for (unsigned int i=0;i<thermalization_steps;i++) {
            cmwc4096();
        }
    }

    Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1> get_eigenvalue(){return eigenvalue;}

private:
    double delta,T,c,d;
    bool accelerated,first;

    Eigen::SelfAdjointEigenSolver< Eigen::Matrix<double,eigen_matrix_dim(N,DIM),eigen_matrix_dim(N,DIM),1> > eigensolver;

    Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1>                         eigenvalue;
    Eigen::Matrix<double,eigen_matrix_dim(N,DIM),1>                         z,pos_p;
    Eigen::Matrix<double,eigen_matrix_dim(N,DIM),eigen_matrix_dim(N,DIM) > H,H_inv,eigenvectors;

    using Integrator<N,DIM,FLAGS>::p;
    using Integrator<N,DIM,FLAGS>::pos_m;
    using Integrator<N,DIM,FLAGS>::deriv_m;
    using Integrator<N,DIM,FLAGS>::hessian_m;


};

int main(int argc,char *argv[]) {

    bool accelerated=false, test_hessian=false,plot_eigenvalue=false,eoutput=false,random_state=true;
    unsigned int natoms=0;
    std::string eigen_out,energy_out,file_in,file_out,json_in="input.json",prefix="";
    double d=1.0,cutoff=10.0;
    if (argc != 1) {
        json_in.assign(argv[1]);
    }
    std::cerr << "Reading input from \""<<json_in<<"\"\n";

    std::ifstream inputjs(json_in);
    json js;
    inputjs >> js;
    inputjs.close();

    /// test con una forma quadratica:
    if (false){
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
    }

    if (NATOMS<=0) {
        if (js.count("natoms")==0) {
            std::cerr << "Error: you must specify \"natoms\": [#atoms]\n";
            return -1;
        } else {
            natoms=js["natoms"];
        }
    } else {
        natoms=NATOMS;
    }

    if (js.count("cutoff")==0) {
        std::cerr << "Error: you must specify \"cutoff\": [potential cutoff]\n";
        return -1;
    } else {
        cutoff=js["cutoff"];
    }

    if (js.count("prefix")!=0) {
        prefix=js["prefix"];
        std::cerr << "Using the prefix for all filenames: \""<<prefix<<"\"\n";
    }


    if (js.count("initial_state")==0) {
        std::cerr << "Error: you must specify \"initial_state\": [true|false]\n";
        return -1;
    } else {
        random_state=!js["initial_state"];
        if (!random_state) {
            if (js.count("file_in")==0){
                std::cerr <<"Error: you must specify the input file for the initial coordinates with \"file_in\": \"[inputfile]\".\n";
                return -1;
            } else {
                file_in=js["file_in"];
            }
        }
    }

    if (js.count("output_file")==0) {
        std::cerr << "Error: please specify \"output_file\" in input file.\n";
        return -1;
    } else {
        file_out=js["output_file"];
    }

    if (js.count("cell_size")==0 && js.count("rho")==0 || js.count("cell_size")>0 && js.count("rho")>0) {
        std::cerr << "Error: you must specify only one of \"cell_size\": [cellsize, for example 5.0]  or \"rho\": [reduced density, for example 0.5]\n";
        return -1;
    }




    if (js.count("test_hessian")!=0) {
        test_hessian=js["test_hessian"];
        std::cerr << "test_hessian: " << (test_hessian?"true":"false") << "\n";
    }



    if (js.count("dynamics")==0) {
        std::cerr << "Error: you must specify the \"dynamics\": {} section\n";
        return -1;
    } else {
        if (js["dynamics"].count("nsteps")==0) {
            std::cerr << "Error: you must specify \"dynamics\":{\"nsteps\": [#number of dynamic steps]}\n";
            return -1;
        }
        if (js["dynamics"].count("output")==0) {
            std::cerr << "Error: you must specify \"dynamics\":{\"output\": [output file for the dynamics]}\n";
            return -1;
        }
        if (js["dynamics"].count("write_energy")>0) {
            eoutput=js["dynamics"]["write_energy"];
            if (eoutput){
                if (js["dynamics"].count("energy_out")==0){
                    std::cerr << "Error: you must specify \"dynamics\":{\"energy_out\": [output file for the thermodynamic quantities]}\n";
                    return -1;
                }else {
                    energy_out=js["dynamics"]["energy_out"];
                }
            }
        }
        if (js["dynamics"].count("T")==0) {
            std::cerr << "Error: you must specify \"dynamics\":{\"T\": [temperature of the system]}\n";
            return -1;
        }
        if (js["dynamics"].count("dt")==0) {
            std::cerr << "Error: you must specify \"dynamics\":{\"dt\": [timestep of the dynamic]}\n";
            return -1;
        }
        if (js["dynamics"].count("accelerated")>0) {
            std::cerr << "Accelerated dynamics: ";
            accelerated=js["dynamics"]["accelerated"];
            std::cerr << (accelerated?"true":"false")<<"\n";
            if (js["dynamics"].count("delta")>0) {
                d=js["dynamics"]["delta"];
            }
            if (js["dynamics"].count("plot_eigenvalue")>0) {
                plot_eigenvalue=js["dynamics"]["plot_eigenvalue"];
                if (plot_eigenvalue){
                    if (js["dynamics"].count("eigen_out")==0){
                        std::cerr << "Error: you must specify \"dynamics\":{\"eigen_out\": [output file for the eigenvalues of the hessian]}\n";
                        return -1;
                    }
                    eigen_out=js["dynamics"]["eigen_out"];
                }
            }
        }
    }
    if (js.count("minimization")==0) {
        std::cerr << "Error: you must specify \"minimization\": {}\n";
        return -1;
        if (js["minimization"].count("nsteps")==0) {
            std::cerr << "Error: you must specify \"minimization\":{\"nsteps\": [#number of cg steps]}\n";
            return -1;
        }
    }

    unsigned int nsteps_d=js["dynamics"]["nsteps"];
    unsigned int nsteps_cg=js["minimization"]["nsteps"];
    unsigned int dump_mod=1;
    double cell_size;//
    if (js.count("rho")==0) {
        cell_size=js["cell_size"];
    } else {
        double rho=js["rho"];
        cell_size=std::pow(double(natoms)/rho,1.0/3.0);
    }
    double temperature=js["dynamics"]["T"];
    double dt=js["dynamics"]["dt"];
    double Tfinal;
    std::string outname=js["dynamics"]["output"];
    if (js["dynamics"].count("print")==1){
        dump_mod=js["dynamics"]["print"];
    }

    if (js["dynamics"].count("Tfinal")==0) {
        Tfinal=temperature;
    } else {
        Tfinal=js["dynamics"]["Tfinal"];
    }



    LJPair<NATOMS,DIMs> test(cutoff);
    Eigen::Matrix<double,DIMs,DIMs> cel;

    cel = cell_size * Eigen::Matrix<double,DIMs,DIMs>::Identity();
    double volume=std::pow(cell_size,DIMs);
    test.init_pbc( cel);

    Eigen::Matrix<double,eigen_matrix_dim(NATOMS,DIMs),1>x;
    if (NATOMS<=0) {
        x.resize(natoms*DIMs,1);
    }
    if (random_state){
        std::cerr << "Initializing atoms in random positions.\n";
        x=Eigen::Matrix<double,eigen_matrix_dim(NATOMS,DIMs),1>::Random(natoms*DIMs,1)*cell_size;
    } else {
        std::cerr << "Reading initial positions from file \""<<file_in<<"\".\n";
        std::ifstream infile(file_in);
        infile >> x;
        if (!infile.good()) {
            std::cerr << "An error occured while reading file (wrong number of atoms? Wrong file?).\n";
            return -1;
        }
        infile.close();
    }
    if (nsteps_cg>0){

        ParabolaLineMinimization <Eigen::Matrix<double,eigen_matrix_dim(NATOMS,DIMs),1>,double,LJPair<NATOMS,DIMs> > lineMin(0.02,0.1,4,3);
        Cg<Eigen::Matrix<double,eigen_matrix_dim(NATOMS,DIMs),1>,double,LJPair<NATOMS,DIMs> > testcg(test,x,test(x),lineMin,8);
        std::cerr << "Performing cg minimization of the initial state.\n";
        for (unsigned int i=0;i<nsteps_cg;i++) {
            if (i%100==0)
                std::cout << i << " " << testcg.get_fx()<< "\n";
            if (!testcg.iteration())
                break;
        }
        std::cout <<  "Final: " << testcg.get_fx()<< "\n";
        x=testcg.get_x();
    }

    if (test_hessian){
        std::cerr << "Testing hessian and forces algorithm..";
        std::cout << "Forces calculation test: " << test.check_forces(x,0.0001,0.001)<<"\n";
        std::cerr << ".";
        std::cout << "Hessian calculation test: " << test.check_hessian_forces(x,0.0001,0.001)<<"\n";
        std::cerr << " Done!\n";
    }

    if (nsteps_d>0){
        std::ofstream output_cmd(prefix+".vmd_cmd");
        output_cmd <<"#mol delete 0\n\
             #mol addrep 0\n\
             #display resetview\n\
                 mol new {"<< prefix+outname<<"} type {xyz} first 0 last -1 step 1 waitfor -1\n"
                 << "set cell [pbc set { "<< cell_size <<" " <<cell_size<<" " <<cell_size<< " } -all]\npbc wrap -all\npbc box\nmol modstyle 0 0 VDW 0.100000 12.000000\n";
        output_cmd.close();
        std::cout << "\n\n========================\n| Begin of the dynamic |\n========================\n\n";
        std::ofstream output(prefix+outname,std::ios_base::app);
        std::ofstream output_energy(prefix+energy_out,std::ios_base::app);
        std::ofstream output_eigen(prefix+eigen_out,std::ios_base::app);
        IntegratorAcceleratedLangevin<NATOMS,DIMs,11> firstOrderAcceleratedLangevin(&test,x,dt,temperature,natoms,accelerated);
        if(accelerated)
            firstOrderAcceleratedLangevin.set_d(d);
        //forse ok in c++17
        //    IntegratorAcceleratedLangevin<> firstOrderAcceleratedLangevin(&test,x,0.001,0.5);
        firstOrderAcceleratedLangevin.init_global_rnd(67857);
        //    test.pbc_wrap(x);


        std::cerr <<"natoms = "<< natoms << ";  DIMs = " << DIMs<<";  volume = " << volume<<"\n";

        for (unsigned int istep=0;istep<nsteps_d;istep++) {
            double T=temperature+ (Tfinal-temperature)*double(istep)/double(nsteps_d-1);
            firstOrderAcceleratedLangevin.set_T(T);
            firstOrderAcceleratedLangevin.step(x);
            if (istep%dump_mod==0){
                firstOrderAcceleratedLangevin.dump(output,x);
                if (plot_eigenvalue){
                    output_eigen << firstOrderAcceleratedLangevin.get_eigenvalue().transpose()<<"\n";
                }
            }
            if (eoutput &&  istep%dump_mod==0){
                std::cout <<istep<<" "<<T<<" " << test.get_energy() << " "<< ( T*natoms- test.get_virial().trace()/DIMs)/volume<< " "<<test.get_virial().trace()/DIMs/volume<< "\n";
                output_energy <<istep<<" "<<T<<" " << test.get_energy() << " "<< ( T*natoms- test.get_virial().trace()/DIMs)/volume<< " "<<test.get_virial().trace()/DIMs/volume<<"\n";
            }
        }
        std::cout << "Finished!\nVMD:\n\n vmd -e .vmd_cmd\n\nClean:\n\n rm "<<(eoutput?"\""+prefix+energy_out+"\" ":" ") <<(plot_eigenvalue?"\""+prefix+eigen_out+"\" ":" ")<< "\"" << prefix+outname << "\""<< prefix+".vmd_cmd \n\n" ;




        std::cerr << "End of run. Writing output in \""<<prefix+file_out<<"\"...\n";
    }

    std::ofstream output_final(prefix+file_out);
    output_final << x<<"\n";
    if (output_final.good()) {
        std::cerr<< "Ok.\n";
    } else {
        std::cerr << "Fail!!\n";
    }
}
