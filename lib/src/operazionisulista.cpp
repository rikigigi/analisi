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

template < class T, class TFLOAT > const T  VectorOp < T, TFLOAT >::operator + (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res+=destra;
    return res;
}
template < class T, class TFLOAT > const T  VectorOp < T, TFLOAT >::operator - (const T & destra) const {
    T res =  static_cast<const T&>(*this);
    res-=destra;
    return res;
}
template < class T, class TFLOAT > const T  VectorOp < T, TFLOAT >::operator * (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res*=destra;
    return res;
}
template < class T, class TFLOAT > const T  VectorOp < T, TFLOAT >::operator / (const T & destra) const {
    T res = static_cast<const T&>(*this);
    res/=destra;
    return res;
}


template < class T, class TFLOAT > const T  VectorOp < T, TFLOAT >::operator + (const TFLOAT & destra) const {
    T res = static_cast<const T&>(*this);
    res+=destra;
    return res;
}
template < class T, class TFLOAT > const T  VectorOp < T, TFLOAT >::operator - (const TFLOAT & destra) const {
    T res = static_cast<const T&>(*this);
    res-=destra;
    return res;
}
template < class T, class TFLOAT > const T  VectorOp < T, TFLOAT >::operator * (const TFLOAT & destra) const {
    T res = static_cast<const T&>(*this);
    res*=destra;
    return res;
}
template < class T, class TFLOAT > const T  VectorOp < T, TFLOAT >::operator / (const TFLOAT & destra) const {
    T res = static_cast<const T&>(*this);
    res/=destra;
    return res;
}




template class VectorOp<SpettroVibrazionale>;
template class VectorOp<PosizioniEquilibrio>;
template class VectorOp<ModiVibrazionali>;
template class VectorOp<IstogrammaVelocita>;
template class VectorOp<GreenKubo2ComponentIonicFluid>;
template class VectorOp<GreenKuboNComponentIonicFluid<> >;
template class VectorOp<GreenKuboNComponentIonicFluid<long double,double>,long double >;
template class VectorOp<GreenKuboNComponentIonicFluid<long double,long double>,long double >;
template class VectorOp<MSD>;
template class VectorOp<Gofrt<double>,double>;
template class VectorOp<Gofrt<long double>,long double>;
template class VectorOp<CorrelatoreSpaziale,double>;
*/
