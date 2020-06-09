/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef OPERAZIONISULISTA_H
#define OPERAZIONISULISTA_H

#include "mediablocchi.h"
//#include <memory>

/**
  * Definisce le operazioni, eseguite membro a membro, delle liste.
  * Definisce le operazioni, eseguite membro a membro, con uno scalare.
**/

template <class T,class TFLOAT=double> class OperazioniSuLista
{
public:
    OperazioniSuLista<T,TFLOAT> &operator =(const OperazioniSuLista<T,TFLOAT> &);
    OperazioniSuLista<T,TFLOAT> &operator =(OperazioniSuLista<T,TFLOAT> &&) = default;
    T & operator+= (const T&);
    T & operator-= (const T&);
    T & operator*= (const T&);
    T & operator/= (const T&);

    /**TODO
     * QUESTI NON FUNZIONANO
     * sistemare l'operatore di copia ed risolvere i problemi con la memoria
    **/
    const T operator+ (const T&) const ;
    const T operator- (const T&) const ;
    const T operator* (const T&) const ;
    const T operator/ (const T&) const ;

    T & operator+= (const TFLOAT&);
    T & operator-= (const TFLOAT&);
    T & operator*= (const TFLOAT&);
    T & operator/= (const TFLOAT&);

    /**TODO
     * QUESTI NON SONO MAI STATI TESTATI
     * E PROBABILMENTE NON FUNZIONANO
    **/

    const T operator+ (const TFLOAT &) const ;
    const T operator- (const TFLOAT&) const ;
    const T operator* (const TFLOAT&) const ;
    const T operator/ (const TFLOAT&) const ;
    unsigned int lunghezza() const;
    TFLOAT elemento(unsigned int i)const;
    TFLOAT * accesso_lista(){return lista;}
    void azzera();
protected:
    OperazioniSuLista<T,TFLOAT> (const OperazioniSuLista<T,TFLOAT> & other);
    TFLOAT * lista;
    OperazioniSuLista();
    ~OperazioniSuLista();
    unsigned int lunghezza_lista;
};

#endif // OPERAZIONISULISTA_H
