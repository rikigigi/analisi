#ifndef CG_H
#define CG_H

#include "lineminimizations.h"




template <class X,class K,class F>
class Cg
{
public:
    Cg(F &f_,LineMinimization<X,K,F> &lm_) : f(f_),lm(lm_) {}
    Cg(F &f_,X x0_,K  f_x0, LineMinimization<X,K,F> &lm_,unsigned int min_iteration_without_restarts):
        f(f_),lm(lm_) {init(x0_,f_x0,min_iteration_without_restarts);}
    void init(X x0_,K  f_x0, unsigned int min_iteration_without_restarts ) {
        min_iteration_without_restart=min_iteration_without_restarts;
        x=x0_;
        f_x=f_x0;

        //first iteration
        n_iterations=0;
        restart_iteration();
    }

    void restart_iteration() {
        std::cout << "Restarting CG\n";
        d=-f.deriv(x);
        r=d;
        lm(f,x,d,f_x);
        restart=false;
        last_restart_iteration=n_iterations;
        n_iterations++;
    }

    void iteration() {
        if (restart) {
            restart_iteration();
        } else {
            rp=r;
            r=-f.deriv(x);
            beta=(r.transpose()*(r-rp)/(rp.squaredNorm()))(0);
            if (beta<0) beta=0;
            d=r+beta*d;
            auto res = lm(f,x,d,f_x);
            if (res!=  LineMinimization<X, K, F>::Result::MinimumFound &&
                    n_iterations - last_restart_iteration > min_iteration_without_restart) {
                restart=true;
            }
            n_iterations++;
        }
    }
    X & get_x() {return x;}
    K & get_fx() {return f_x;}
    unsigned int get_iterations() {return n_iterations;}

private:
    F &f;
    LineMinimization<X,K,F> & lm;
    X x,d,r,rp;
    K beta,f_x;
    unsigned int n_iterations,last_restart_iteration,min_iteration_without_restart;
    bool restart;
};

#endif // CG_H
