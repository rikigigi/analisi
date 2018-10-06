#ifndef CG_H
#define CG_H

#include "lineminimizations.h"




template <class X,class K>
class Cg
{
public:
    Cg() {}
    void init(Function<X,K> &f_,X x0_,K & f_x0, LineMinimization<X,K> &lm_) {
        f=f_;
        lm=lm_;
        x=x0_;
        f_x=f_x0;

        //first iteration
        n_iterations=0;
        restart_iteration();
    }

    void restart_iteration() {
        d=-f.deriv(x);
        r=d;
        lm(f,x,d,f_x);
        n_iterations++;
    }

    void iteration() {
        rp=r;
        r=-f.deriv(x);
        beta=r.transpose()*(r-rp)/(rp.transpose()*rp);
        if (beta<0) beta=0;
        d=r+beta*d;
        lm(f,x,d,f_x);
        n_iterations++;
    }
    X & get_x() {return x;}
    unsigned int get_iterations() {return n_iterations;}

private:
    Function<X,K> &f;
    LineMinimization<X,K> & lm;
    X x,d,r,rp;
    K beta,f_x;
    unsigned int n_iterations;
};

#endif // CG_H
