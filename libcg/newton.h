#ifndef NEWTON_H
#define NEWTON_H

template <class X,class K,class F>
class Newton
{
public:
    void init ( X x0_,K  f_x0 );
    bool iteration();
    X & get_x() {return x;}
    K & get_fx() {return f_x;}
private:
   F &f;
   X x;
   K f_x;
   unsigned int n_iterations;
};

#endif // NEWTON_H
