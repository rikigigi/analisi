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



#include "operazionisulista.h"
#include <cstdlib>
#include <iostream>
#include "spettrovibrazionale.h"
#include "posizioniequilibrio.h"
#include "modivibrazionali.h"
#include "istogrammavelocita.h"
#include "config.h"
#include "greenkubo2componentionicfluid.h"

template <class T>  OperazioniSuLista<T>::OperazioniSuLista()
{
    lista=0;
    lunghezza_lista=0;
}

template <class T> double OperazioniSuLista<T>::elemento(unsigned int i) const{
    if (i<lunghezza_lista) {
        return lista[i];
    } else {
        std::cerr << "Errore: fuori dal range! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
        abort();
    }
}
template <class T> unsigned int  OperazioniSuLista<T>::lunghezza() const{
    return lunghezza_lista;
}

template <class T> T &  OperazioniSuLista<T>::operator += (const T & destra) {
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

template <class T> T &  OperazioniSuLista<T>::operator -= (const T & destra) {
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

template <class T> T &  OperazioniSuLista<T>::operator *= (const T & destra) {
    if (destra.lunghezza()== lunghezza_lista) {
#ifdef DEBUG
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

template <class T> T &  OperazioniSuLista<T>::operator /= (const T & destra) {
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

template <class T> const T  OperazioniSuLista<T>::operator + (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res+=destra;
    return res;
}
template <class T> const T  OperazioniSuLista<T>::operator - (const T & destra) const {
    T res =  static_cast<const T&>(*this);
    res-=destra;
    return res;
}
template <class T> const T  OperazioniSuLista<T>::operator * (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res*=destra;
    return res;
}
template <class T> const T  OperazioniSuLista<T>::operator / (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res/=destra;
    return res;
}

template <class T> T &  OperazioniSuLista<T>::operator += (const double & destra) {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]+=destra;
    }
    return static_cast<T&>(*this);
}
template <class T> T &  OperazioniSuLista<T>::operator -= (const double & destra) {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]-=destra;
    }
    return static_cast<T&>(*this);
}
template <class T> T &  OperazioniSuLista<T>::operator *= (const double & destra) {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]*=destra;
    }
    return static_cast<T&>(*this);
}
template <class T> T &  OperazioniSuLista<T>::operator /= (const double & destra) {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]/=destra;
    }
    return static_cast<T&>(*this);
}

template <class T> const T  OperazioniSuLista<T>::operator + (const double & destra) const {
    T res = static_cast<const T&>(*this);
    res+=destra;
    return res;
}
template <class T> const T  OperazioniSuLista<T>::operator - (const double & destra) const {
    T res = static_cast<const T&>(*this);
    res-=destra;
    return res;
}
template <class T> const T  OperazioniSuLista<T>::operator * (const double & destra) const {
    T res = static_cast<const T&>(*this);
    res*=destra;
    return res;
}
template <class T> const T  OperazioniSuLista<T>::operator / (const double & destra) const {
    T res = static_cast<const T&>(*this);
    res/=destra;
    return res;
}


template <class T> OperazioniSuLista<T> & OperazioniSuLista<T>::operator =(const OperazioniSuLista<T> &destra) {
#ifdef DEBUG
    std::cerr << "chiamato OperazioniSuLista<T>::operator = " __FILE__ ":"<<__LINE__<<"\n";
#endif
    if (lunghezza_lista!=destra.lunghezza_lista) { //rialloca la memoria
#ifdef DEBUG
    std::cerr << "delete [] "<<lista<<"\n";
#endif
        delete  [] lista;
        lunghezza_lista=destra.lunghezza_lista;
        lista = new double[lunghezza_lista];
#ifdef DEBUG
    std::cerr << "new double [] "<<lista<<"\n";
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

template <class T> OperazioniSuLista<T>::OperazioniSuLista(const OperazioniSuLista<T> & copiare) {
    lista=0;
    lunghezza_lista=0;
    operator=(copiare);
#ifdef DEBUG
    std::cerr << "OperazioniSuLista<T>::OperazioniSuLista(const OperazioniSuLista<T> & copiare) "<<lista<<"\n";
#endif
}

template <class T> OperazioniSuLista<T>::~OperazioniSuLista() {
#ifdef DEBUG
    std::cerr << "chiamato OperazioniSuLista<T>::~OperazioniSuLista " __FILE__ ":"<<__LINE__<<"\n";
    std::cerr << "delete [] "<<lista<<"\n";
#endif
    delete [] lista;

}

template <class T> void OperazioniSuLista<T>::azzera() {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]=0;
    }
}

template class OperazioniSuLista<SpettroVibrazionale>;
template class OperazioniSuLista<PosizioniEquilibrio>;
template class OperazioniSuLista<ModiVibrazionali>;
template class OperazioniSuLista<IstogrammaVelocita>;
template class OperazioniSuLista<GreenKubo2ComponentIonicFluid>;
