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
#include "config.h"
//#include <memory>

/**
  * Definisce le operazioni, eseguite membro a membro, delle liste.
  * Definisce le operazioni, eseguite membro a membro, con uno scalare.
**/

template <class T,class TFLOAT=double> class OperazioniSuLista
{
public:
    OperazioniSuLista<T,TFLOAT> &operator =(const OperazioniSuLista<T, TFLOAT> &destra) {
#ifdef DEBUG2
        std::cerr << "chiamato OperazioniSuLista<T, TFLOAT>::operator = " __FILE__ ":"<<__LINE__<<"\n";
#endif
        if (lunghezza_lista!=destra.lunghezza_lista) { //rialloca la memoria
#ifdef DEBUG2
            std::cerr << "delete [] "<<lista<<"\n";
#endif
            delete  [] lista;
            lunghezza_lista=destra.lunghezza_lista;
            lista = new TFLOAT[lunghezza_lista];
#ifdef DEBUG2
            std::cerr << "new TFLOAT [] "<<lista<<"\n";
#endif
        }
        if (lunghezza_lista==destra.lunghezza_lista){
            for (unsigned int i=0;i<lunghezza_lista;i++) {
                lista[i]=destra.lista[i];
            }
        } else {
            abort();
        }

        return *this;
    }
    OperazioniSuLista<T,TFLOAT> &operator =(OperazioniSuLista<T,TFLOAT> &&) = default;
    T & operator+= (const T&destra){
        if (destra.lunghezza()== lunghezza_lista) {
            for (unsigned int i=0;i<lunghezza_lista;i++) {
                lista[i]+=destra.elemento(i);
            }
        } else {
            std::cerr << "Errore: lunghezza delle liste diverse! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
            abort();
        }
        return static_cast<T&>(*this);
    }
    T & operator-= (const T & destra) {
        if (destra.lunghezza()== lunghezza_lista) {
            for (unsigned int i=0;i<lunghezza_lista;i++) {
                lista[i]-=destra.elemento(i);
            }
        } else {
            std::cerr << "Errore: lunghezza delle liste diverse! ("<< __FILE__ << ":"<< __LINE__ <<")\n";
            abort();
        }
        return static_cast<T&>(*this);
    }
    T & operator*= (const T & destra) {
        if (destra.lunghezza()== lunghezza_lista) {
#ifdef DEBUG2
            std::cerr << "chiamato OperazioniSuLista<T>::operator *= " __FILE__ ":"<<__LINE__<<"\n";
            std::cerr << "lista = "<<lista<<", destra.lista = " << destra.lista << "\n";
#endif
            for (unsigned int i=0;i<lunghezza_lista;i++) {
                lista[i]*=destra.elemento(i);
            }
        } else {
            std::cerr << "Errore: lunghezza delle liste diverse! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
            abort();
        }
        return static_cast<T&>(*this);
    }
    T & operator/= (const T & destra) {
        if (destra.lunghezza()== lunghezza_lista) {
            for (unsigned int i=0;i<lunghezza_lista;i++) {
                lista[i]/=destra.elemento(i);
            }
        } else {
            std::cerr << "Errore: lunghezza delle liste diverse! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
            abort();
        }
        return static_cast<T&>(*this);
    }

    T & operator+= (const TFLOAT & destra) {
        for (unsigned int i=0;i<lunghezza_lista;i++) {
            lista[i]+=destra;
        }
        return static_cast<T&>(*this);
    }
    T & operator-= (const TFLOAT & destra) {
        for (unsigned int i=0;i<lunghezza_lista;i++) {
            lista[i]-=destra;
        }
        return static_cast<T&>(*this);
    }
    T & operator*= (const TFLOAT & destra) {
        for (unsigned int i=0;i<lunghezza_lista;i++) {
            lista[i]*=destra;
        }
        return static_cast<T&>(*this);
    }
    T & operator/= (const TFLOAT & destra) {
        for (unsigned int i=0;i<lunghezza_lista;i++) {
            lista[i]/=destra;
        }
        return static_cast<T&>(*this);
    }


//    /**TODO
//     * QUESTI NON FUNZIONANO
//     * sistemare l'operatore di copia ed risolvere i problemi con la memoria
//    **/
//    const T operator+ (const T&) const ;
//    const T operator- (const T&) const ;
//    const T operator* (const T&) const ;
//    const T operator/ (const T&) const ;

//    /**TODO
//     * QUESTI NON SONO MAI STATI TESTATI
//     * E PROBABILMENTE NON FUNZIONANO
//    **/

//    const T operator+ (const TFLOAT &) const ;
//    const T operator- (const TFLOAT&) const ;
//    const T operator* (const TFLOAT&) const ;
//    const T operator/ (const TFLOAT&) const ;
    unsigned int lunghezza() const{
        return lunghezza_lista;
    }

    TFLOAT elemento(unsigned int i)const{
        if (i<lunghezza_lista) {
            return lista[i];
        } else {
            std::cerr << "Errore: fuori dal range! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
            abort();
        }
    }
    TFLOAT * accesso_lista(){return lista;}
    void azzera(){
        for (unsigned int i=0;i<lunghezza_lista;i++) {
            lista[i]=0;
        }
    }
    inline void azzera(int start, int stop){
        for (unsigned int i=start;i<stop;i++) {
            lista[i]=0;
        }
    }
protected:
    OperazioniSuLista<T,TFLOAT> (const OperazioniSuLista<T, TFLOAT> & copiare) {
        lista=0;
        lunghezza_lista=0;
        operator=(copiare);
    #ifdef DEBUG
        std::cerr << "OperazioniSuLista<T, TFLOAT>::OperazioniSuLista(const OperazioniSuLista<T, TFLOAT> & copiare) "<<lista<<"\n";
    #endif
    }

    TFLOAT * lista;
    OperazioniSuLista(){
        lista=0;
        lunghezza_lista=0;
    }
    ~OperazioniSuLista(){
#ifdef DEBUG2
    std::cerr << "chiamato OperazioniSuLista<T, TFLOAT>::~OperazioniSuLista " __FILE__ ":"<<__LINE__<<"\n";
    std::cerr << "delete [] "<<lista<<"\n";
#endif
    delete [] lista;

}
    unsigned int lunghezza_lista;
};

#endif // OPERAZIONISULISTA_H
