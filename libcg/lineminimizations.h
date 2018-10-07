#ifndef LINEMINIMIZATIONS_H
#define LINEMINIMIZATIONS_H

#include "function.h"
#include <cmath>
#include <iostream>


template <class X, class K,class F, class XX = void> class LineMinimization {
  public:
    enum Result {MissedMinimum,MinimumFound,MaximumIterationsExceeded};
    /**
      * Minimizza la funzione f nella direzione di d:  Min_a f(x + a*d)
      * Deve aggiornare x al valore del minimo trovato
      * Viene fornito il valore di f in x f_x, che poi viene aggiornato al nuovo valore.
    **/
    virtual Result operator () (F f,X & x, const X & d, K & f_x)=0;
};



/**
  * deve essere definito sqrt(K). Trova un minimo approssimato supponendo che la funzione sia localmente una parabola
**/
template <class X, class K, class F> class ParabolaLineMinimization : public LineMinimization<X,K,F> {
public:
    typedef typename LineMinimization<X, K, F>::Result Result;
    ParabolaLineMinimization(const K dx, const K xmax,
                             unsigned int max_search=3,unsigned int raise_rate=7,
                             unsigned int stop_raise_epoch=42,
                             unsigned int outer_iter=1,
                             double raise_coeff    =1.1,
                             double decr_coeff     =0.75,
                             double decr_coeff_cost=0.998
            ) : outer_iter(outer_iter), dx(dx), xmax(xmax),max_search(max_search),
    stop_raise_epoch(stop_raise_epoch),raise_rate(raise_rate),tot_iter(0),
      raise_coeff(raise_coeff),decr_coeff(decr_coeff),decr_coeff_cost(decr_coeff_cost)
    {
        if (outer_iter==0) outer_iter=1;
        fail_count=0;
    }
    virtual Result operator () (F f,X & x0, const X & d, K & f_a) override {
        Result res=Result::MinimumFound;
        if (calls_without_fail%raise_rate==0 && calls_without_fail>0 &&tot_iter<stop_raise_epoch ){
            dx=dx*raise_coeff ;
            std::cout << "Setting new dx_0="<<dx<<"\n";
        }else {
            dx=dx*decr_coeff_cost;
        }
        calls_without_fail++;
        tot_iter++;

        X du=d/d.norm();
        K dx_o=dx;

        unsigned int retries=0;

        for (int i1=0;i1<outer_iter;i1++){
            K f_b,f_c,b,c,a_min; //a=0.0
            X x1,x2;
            b=dx;
            x1=x0+du*b;
            f_b=f(x1);
            for (unsigned int i=0;i<max_search && f_b>=f_a;i++){
                b=b*(i+1);
                x1=x0+du*b;
                f_b=f(x1);
            }
            if (f_b>=f_a){
                res= Result::MaximumIterationsExceeded;
                dx/=2.0;
                retries++;
                if (retries<max_search) {
                    i1--;
                } else {
                    x0=x1;
                    f_a=f_b;
                    std::cout << "Failed to find a minimum\ndx="<<dx<< " b=" << b << " i1=" <<i1<<"\n";
                    calls_without_fail=0;
                    fail_count++;
                    if (fail_count%max_search==0 && fail_count>0) {
                        dx_o=dx_o*decr_coeff;
                        std::cout << "Setting new dx_0="<<dx_o<<"\n";
                    }
                }
                continue;
            }
            c=b+dx;
            x2=x0+du*c;
            f_c=f(x2);
            for (unsigned int i=0;i<max_search && f_c<=f_b;i++){
                c=b+dx*(i+1);
                x2=x0+du*c;
                f_c=f(x2);
            }
            if (f_c<=f_b){
                res= Result::MaximumIterationsExceeded;
            }

            if (res==Result::MaximumIterationsExceeded){
                x0=x2;
                f_a=f_c;
            } else {
                a_min=    b     +  0.5  *     (   (f_a-f_b)*(c-b)*(c-b)  -  (f_c-f_b)*b*b   )
       //                                   ---------------------------------------------------
        /                                     (   (f_a-f_b)*(c-b)        +  (f_c-f_b)*b     )        ;
                x0=x0+du*a_min;
                f_a=f(x0);
            }
        }

        dx=dx_o;

        return res;

    }
private:
    K dx,xmax;
    unsigned int max_search,outer_iter,fail_count,calls_without_fail,
    raise_rate,stop_raise_epoch,tot_iter;
    double raise_coeff, decr_coeff, decr_coeff_cost;
};





#endif // LINEMINIMIZATIONS_H
