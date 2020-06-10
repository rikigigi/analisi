/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/

/*

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
#include "correlatorespaziale.h"
#include "heatc.h"
#include "calcolamultithread.h"
#include "centerdiff.h"
#include "centerofmassdiff.h"


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
template class OperazioniSuLista<HeatC, double>;
template class OperazioniSuLista<CorrelatoreSpaziale,double>;
template class OperazioniSuLista<CenterDiff,double>;
template class OperazioniSuLista<CenterOfMassDiff,double>;
//template class OperazioniSuLista<CalcolaMultiThread,double>;
*/
