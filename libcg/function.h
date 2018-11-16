#ifndef FUNCTION_H
#define FUNCTION_H
/**
  * Classe da derivare per definire una funzione minimizzabile con il CG
  * X deve avere definito il prodotto e .transpose()
  * K deve avere funzionare con < 0 e = 0, e deve essere restituito da X::transpose() * X
**/
template <class X, class K, class XX = void,
          class X6=X,class X7=const X & ,class X2=X,
          class X3=X,class X4=X,class X5=X,class K2=K
          > class Function {
  public:
    /**
      * Restituisce il valore della funzione in x (scalare)
    **/
    virtual K  operator() (const X & x)=0;
    /**
      * Restituisce il valore della derivata rispetto all'argomento in x (vettore)
    **/
    virtual X6 deriv(const X2 & x)=0;
    /**
      * Restituisce il valore della derivata direzionale nella direzione d ( il comportamento di default è grad(f)*d )
      * calcolata in x
    **/
    virtual K2 deriv(const X3 & x, const X4 & d) {
        return d.transpose()*deriv(x);
    }
    /**
      * Restituisce la matrice hessiana calcolata nel punto x
    **/
    virtual XX hessian(const X5 & x) {
        X6 x2;
        return hessian_deriv(x,x2);
    }
    /**
      * Calcola la matrice hessiana e la derivata prima assiem
      *  (la derivata prima a questo punto è gratis)
    **/
    virtual XX hessian_deriv(const X5 & x, X7 x2) {}
};

#endif // FUNCTION_H
