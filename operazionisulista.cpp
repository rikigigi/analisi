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
#include "greenkuboNcomponentionicfluid.h"
#include "gofrt.h"
#include "msd.h"

template < class T, class TFLOAT >  OperazioniSuLista < T, TFLOAT >::OperazioniSuLista()
{
    lista=0;
    lunghezza_lista=0;
}

template < class T, class TFLOAT > TFLOAT OperazioniSuLista < T, TFLOAT >::elemento(unsigned int i) const{
    if (i<lunghezza_lista) {
        return lista[i];
    } else {
        std::cerr << "Errore: fuori dal range! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
        abort();
    }
}
template < class T, class TFLOAT > unsigned int OperazioniSuLista < T, TFLOAT >::lunghezza() const{
    return lunghezza_lista;
}

template < class T, class TFLOAT > T &   OperazioniSuLista < T, TFLOAT >::operator += (const T & destra) {
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

template < class T, class TFLOAT > T &   OperazioniSuLista < T, TFLOAT >::operator -= (const T & destra) {
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

template < class T, class TFLOAT > T &   OperazioniSuLista < T, TFLOAT >::operator *= (const T & destra) {
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

template < class T, class TFLOAT > T &   OperazioniSuLista < T, TFLOAT >::operator /= (const T & destra) {
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

template < class T, class TFLOAT > const T  OperazioniSuLista < T, TFLOAT >::operator + (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res+=destra;
    return res;
}
template < class T, class TFLOAT > const T  OperazioniSuLista < T, TFLOAT >::operator - (const T & destra) const {
    T res =  static_cast<const T&>(*this);
    res-=destra;
    return res;
}
template < class T, class TFLOAT > const T  OperazioniSuLista < T, TFLOAT >::operator * (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res*=destra;
    return res;
}
template < class T, class TFLOAT > const T  OperazioniSuLista < T, TFLOAT >::operator / (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res/=destra;
    return res;
}

template < class T, class TFLOAT > T &   OperazioniSuLista < T, TFLOAT >::operator += (const TFLOAT & destra) {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]+=destra;
    }
    return static_cast<T&>(*this);
}
template < class T, class TFLOAT > T &   OperazioniSuLista < T, TFLOAT >::operator -= (const TFLOAT & destra) {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]-=destra;
    }
    return static_cast<T&>(*this);
}
template < class T, class TFLOAT > T &   OperazioniSuLista < T, TFLOAT >::operator *= (const TFLOAT & destra) {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]*=destra;
    }
    return static_cast<T&>(*this);
}
template < class T, class TFLOAT > T &   OperazioniSuLista < T, TFLOAT >::operator /= (const TFLOAT & destra) {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]/=destra;
    }
    return static_cast<T&>(*this);
}

template < class T, class TFLOAT > const T  OperazioniSuLista < T, TFLOAT >::operator + (const TFLOAT & destra) const {
    T res = static_cast<const T&>(*this);
    res+=destra;
    return res;
}
template < class T, class TFLOAT > const T  OperazioniSuLista < T, TFLOAT >::operator - (const TFLOAT & destra) const {
    T res = static_cast<const T&>(*this);
    res-=destra;
    return res;
}
template < class T, class TFLOAT > const T  OperazioniSuLista < T, TFLOAT >::operator * (const TFLOAT & destra) const {
    T res = static_cast<const T&>(*this);
    res*=destra;
    return res;
}
template < class T, class TFLOAT > const T  OperazioniSuLista < T, TFLOAT >::operator / (const TFLOAT & destra) const {
    T res = static_cast<const T&>(*this);
    res/=destra;
    return res;
}


template < class T, class TFLOAT >  OperazioniSuLista < T, TFLOAT > & OperazioniSuLista<T, TFLOAT>::operator =(const OperazioniSuLista<T, TFLOAT> &destra) {
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

template < class T, class TFLOAT >  OperazioniSuLista < T, TFLOAT >::OperazioniSuLista(const OperazioniSuLista<T, TFLOAT> & copiare) {
    lista=0;
    lunghezza_lista=0;
    operator=(copiare);
#ifdef DEBUG
    std::cerr << "OperazioniSuLista<T, TFLOAT>::OperazioniSuLista(const OperazioniSuLista<T, TFLOAT> & copiare) "<<lista<<"\n";
#endif
}

template < class T, class TFLOAT >  OperazioniSuLista < T, TFLOAT >::~OperazioniSuLista() {
#ifdef DEBUG2
    std::cerr << "chiamato OperazioniSuLista<T, TFLOAT>::~OperazioniSuLista " __FILE__ ":"<<__LINE__<<"\n";
    std::cerr << "delete [] "<<lista<<"\n";
#endif
    delete [] lista;

}

template < class T, class TFLOAT > void OperazioniSuLista < T, TFLOAT >::azzera() {
    for (unsigned int i=0;i<lunghezza_lista;i++) {
        lista[i]=0;
    }
}

template class OperazioniSuLista<SpettroVibrazionale>;
template class OperazioniSuLista<PosizioniEquilibrio>;
template class OperazioniSuLista<ModiVibrazionali>;
template class OperazioniSuLista<IstogrammaVelocita>;
template class OperazioniSuLista<GreenKubo2ComponentIonicFluid>;
template class OperazioniSuLista<GreenKuboNComponentIonicFluid<> >;
template class OperazioniSuLista<GreenKuboNComponentIonicFluid<long double,double>,long double >;
template class OperazioniSuLista<GreenKuboNComponentIonicFluid<long double,long double>,long double >;
template class OperazioniSuLista<MSD>;
template class OperazioniSuLista<Gofrt<double>,double>;
template class OperazioniSuLista<Gofrt<long double>,long double>;
