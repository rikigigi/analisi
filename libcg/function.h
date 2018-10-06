#ifndef FUNCTION_H
#define FUNCTION_H
/**
  * Classe da derivare per definire una funzione minimizzabile con il CG
  * X deve avere definito il prodotto e .transpose()
  * K deve avere funzionare con < 0 e = 0, e deve essere restituito da X::transpose() * X
**/
template <class X, class K, class XX = void> class Function {
  public:
    /**
      * Restituisce il valore della funzione in x (scalare)
    **/
    virtual K  operator() (const X & x)=0;
    /**
      * Restituisce il valore della derivata rispetto all'argomento in x (vettore)
    **/
    virtual X deriv(const X & x)=0;
    /**
      * Restituisce il valore della derivata direzionale nella direzione d ( il comportamento di default è grad(f)*d )
      * calcolata in x
    **/
    virtual K deriv(const X & x, const X & d) {
        return d.transpose()*deriv(x);
    }
    /**
      * Restituisce la matrice hessiana calcolata nel punto x
    **/
    virtual XX hessian(const X & x) {}
};

#endif // FUNCTION_H
