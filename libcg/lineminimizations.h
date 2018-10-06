#ifndef LINEMINIMIZATIONS_H
#define LINEMINIMIZATIONS_H

#include "function.h"


template <class X, class K, class XX = void> class LineMinimization {
  public:
    enum Result {MissedMinimum,MinimumFound,MaximumIterationsExceeded};
    /**
      * Minimizza la funzione f nella direzione di d:  Min_a f(x + a*d)
      * Deve aggiornare x al valore del minimo trovato
      * Viene fornito il valore di f in x f_x, che poi viene aggiornato al nuovo valore.
    **/
    virtual Result operator () (Function<X,K,XX> & f,X & x, const X & d, K & f_x)=0;
};



/**
  * deve essere definito sqrt(K). Trova un minimo approssimato supponendo che la funzione sia localmente una parabola
**/
template <class X, class K> class ParabolaLineMinimization : public LineMinimization<X,K> {
public:
    typedef typename LineMinimization<X,K>::Result Result;
    ParabolaLineMinimization(K & dx, K& xmax,unsigned int max_search=3,unsigned int outer_iter=1) : outer_iter(outer_iter), dx(dx), xmax(xmax),max_search(max_search) {
        if (outer_iter==0) outer_iter=1;
    }
    virtual Result operator () (Function<X,K> & f,X & x0, const X & d, K & f_a) {
        Result res=Result::MinimumFound;

        K du=d/sqrt(d.transpose()*d);

        for (unsigned int i1=0;i1<outer_iter;i1++){
            K f_b,f_c,x1,x2,b,c,a_min; //a=0.0
            b=dx;
            x1=x0+du*b;
            f_b=f(x1);
            for (unsigned int i=0;i<max_search && f_b>=f_a;i++){
                b=b*(i+1);
                x1=x0+du*b;
                f_b=f(x1);
            }
            if (f_b>=f_a)
                res= Result::MaximumIterationsExcedeed;
            c=b+dx;
            x2=x0+du*c;
            f_c=f(x2);
            for (unsigned int i=0;i<max_search && f_c<=f_b;i++){
                c=b+dx*(i+1);
                x2=x0+du*c;
                f_c=f(x2);
            }
            if (f_c<=f_b)
                res= Result::MaximumIterationsExcedeed;

            if (res==Result::MaximumIterationsExcedeed){
                x0=x0+du*xmax;
                f_a=f(x0);
            } else {
                a_min=    b     +  0.5  *     (   (f_a-f_b)*(c-b)*(c-b)  -  (f_c-f_b)*b*b   )
       //                                   ---------------------------------------------------
        /                                     (   (f_a-f_b)*(c-b)        +  (f_c-f_b)*b     )        ;
                x0=x0+du*a_min;
                f_a=f(x0);
                res=Result::MinimumFound;
            }
        }



        return res;

    }
private:
    K dx,xmax;
    unsigned int max_search,outer_iter;
};





#endif // LINEMINIMIZATIONS_H
