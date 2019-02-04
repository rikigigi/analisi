#ifndef FUNCTION_H
#define FUNCTION_H
/**
  * Classe da derivare per definire una funzione minimizzabile con il CG
  * X deve avere definito il prodotto e .transpose()
  * K deve avere funzionare con < 0 e = 0, e deve essere restituito da X::transpose() * X
**/
template <class X,
          class K,
          class XX = X,
          class X6=X,
          class X7=const X & ,
          class X8=X6,
          class X2=X,
          class X3=X,
          class X4=X,
          class X5=X,
          class K2=K
          > class Function {
  public:
    /**
      * Restituisce il valore della funzione in x (scalare)
    **/
    virtual K  operator() (const X & x)=0;
    /**
      * Restituisce il valore della derivata rispetto all'argomento in x (vettore)
    **/
    virtual void deriv(const X2 & x, X6 y)=0; //// nuova interfaccia
    /**
      * Restituisce il valore della derivata direzionale nella direzione d ( il comportamento di default è grad(f)*d )
      * calcolata in x
    **/
    /*
    virtual K2 dderiv(const X3 & x, const X4 & d) {
        X8 y;
        deriv(x,y);
        return d.transpose()*y;
    }*/
    /**
      * Restituisce la matrice hessiana calcolata nel punto x
    **/
    /*
    virtual void hessian(const X5 & x, XX y) {
        X8 x2;
        hessian_deriv(x,x2,y);
    }
    */
    /**
      * Calcola la matrice hessiana e la derivata prima assiem
      *  (la derivata prima a questo punto è gratis)
    **/
    virtual void hessian_deriv(const X5 & x, X7 x2, XX y)  {}  //// nuova interfaccia
};

#endif // FUNCTION_H
