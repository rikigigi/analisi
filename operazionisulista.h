/**
  *
  * (c) Riccardo Bertossa, 2017
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy to receive a copy
  *   of the good modified code, with comments, at
  *    riccardo dot bertossa at gmail dot com
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

template <class T> class OperazioniSuLista
{
public:
    OperazioniSuLista<T> &operator =(const OperazioniSuLista<T> &);
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

    T & operator+= (const double&);
    T & operator-= (const double&);
    T & operator*= (const double&);
    T & operator/= (const double&);

    /**TODO
     * QUESTI NON SONO MAI STATI TESTATI
     * E PROBABILMENTE NON FUNZIONANO
    **/

    const T operator+ (const double &) const ;
    const T operator- (const double&) const ;
    const T operator* (const double&) const ;
    const T operator/ (const double&) const ;
    unsigned int lunghezza() const;
    double elemento(unsigned int i)const;
    double * accesso_lista(){return lista;}
    void azzera();
protected:
    OperazioniSuLista<T> (const OperazioniSuLista<T> & other);
    double * lista;
    OperazioniSuLista();
    ~OperazioniSuLista();
    unsigned int lunghezza_lista;
};

#endif // OPERAZIONISULISTA_H
