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
#include "compiler_types.h"
//#include <memory>

/**
  * Definisce le operazioni, eseguite membro a membro, delle liste.
  * Definisce le operazioni, eseguite membro a membro, con uno scalare.
**/

template <class T,class TFLOAT=double> class VectorOp
{
public:
    VectorOp<T,TFLOAT> &operator =(const VectorOp<T, TFLOAT> &destra) {
#ifdef DEBUG2
        std::cerr << "chiamato VectorOp<T, TFLOAT>::operator = " __FILE__ ":"<<__LINE__<<"\n";
#endif
        if (data_length!=destra.data_length) { //rialloca la memoria
#ifdef DEBUG2
            std::cerr << "delete [] "<<vdata<<"\n";
#endif
            delete  [] vdata;
            data_length=destra.data_length;
            vdata = new TFLOAT[data_length];
#ifdef DEBUG2
            std::cerr << "new TFLOAT [] "<<vdata<<"\n";
#endif
        }
        if (data_length==destra.data_length){
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]=destra.vdata[i];
            }
        } else {
            abort();
        }

        return *this;
    }
    VectorOp<T,TFLOAT> &operator =(VectorOp<T,TFLOAT> &&) = default;
    T & operator+= (const T&destra){
        if (destra.lunghezza()== data_length) {
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]+=destra.elemento(i);
            }
        } else {
            std::cerr << "Errore: lunghezza delle liste diverse! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
            abort();
        }
        return static_cast<T&>(*this);
    }
    T & operator-= (const T & destra) {
        if (destra.lunghezza()== data_length) {
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]-=destra.elemento(i);
            }
        } else {
            std::cerr << "Errore: lunghezza delle liste diverse! ("<< __FILE__ << ":"<< __LINE__ <<")\n";
            abort();
        }
        return static_cast<T&>(*this);
    }
    T & operator*= (const T & destra) {
        if (destra.lunghezza()== data_length) {
#ifdef DEBUG2
            std::cerr << "chiamato VectorOp<T>::operator *= " __FILE__ ":"<<__LINE__<<"\n";
            std::cerr << "vdata = "<<vdata<<", destra.vdata = " << destra.vdata << "\n";
#endif
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]*=destra.elemento(i);
            }
        } else {
            std::cerr << "Errore: lunghezza delle liste diverse! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
            abort();
        }
        return static_cast<T&>(*this);
    }
    T & operator/= (const T & destra) {
        if (destra.lunghezza()== data_length) {
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]/=destra.elemento(i);
            }
        } else {
            std::cerr << "Errore: lunghezza delle liste diverse! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
            abort();
        }
        return static_cast<T&>(*this);
    }

    T & operator+= (const TFLOAT & destra) {
        for (unsigned int i=0;i<data_length;i++) {
            vdata[i]+=destra;
        }
        return static_cast<T&>(*this);
    }
    T & operator-= (const TFLOAT & destra) {
        for (unsigned int i=0;i<data_length;i++) {
            vdata[i]-=destra;
        }
        return static_cast<T&>(*this);
    }
    T & operator*= (const TFLOAT & destra) {
        for (unsigned int i=0;i<data_length;i++) {
            vdata[i]*=destra;
        }
        return static_cast<T&>(*this);
    }
    T & operator/= (const TFLOAT & destra) {
        for (unsigned int i=0;i<data_length;i++) {
            vdata[i]/=destra;
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
        return data_length;
    }

    TFLOAT elemento(unsigned int i)const{
        if (i<data_length) {
            return vdata[i];
        } else {
            std::cerr << "Errore: fuori dal range! ("<< __FILE__<< ":"<< __LINE__ <<")\n";
            abort();
        }
    }
    TFLOAT * access_vdata(){return vdata;}
    void azzera(){
        for (unsigned int i=0;i<data_length;i++) {
            vdata[i]=0;
        }
    }
    inline void azzera(int start, int stop){
        for (unsigned int i=start;i<stop;i++) {
            vdata[i]=0;
        }
    }
protected:
    VectorOp<T,TFLOAT> (const VectorOp<T, TFLOAT> & copiare) {
        vdata=0;
        data_length=0;
        operator=(copiare);
    #ifdef DEBUG
        std::cerr << "VectorOp<T, TFLOAT>::VectorOp(const VectorOp<T, TFLOAT> & copiare) "<<vdata<<"\n";
    #endif
    }

    TFLOAT * vdata;
    VectorOp(){
        vdata=0;
        data_length=0;
    }
    ~VectorOp(){
#ifdef DEBUG2
    std::cerr << "chiamato VectorOp<T, TFLOAT>::~VectorOp " __FILE__ ":"<<__LINE__<<"\n";
    std::cerr << "delete [] "<<vdata<<"\n";
#endif
    delete [] vdata;

}
    unsigned int data_length;
};

#endif // OPERAZIONISULISTA_H
