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

#include "config.h"
#include "compiler_types.h"
#include <stdexcept>

/**
 * Simple class to copy and operate on vectors on the heap
 * Implements vector-vector elementwise + - * /, and operations with scalars
**/

template <class T,class TFLOAT=double> class VectorOp
{
public:
    VectorOp<T,TFLOAT> &operator =(const VectorOp<T, TFLOAT> &destra) {
        if (data_length!=destra.data_length) { //rialloca la memoria
            delete  [] vdata;
            data_length=destra.data_length;
            vdata = new TFLOAT[data_length];
        }
        if (data_length==destra.data_length){
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]=destra.vdata[i];
            }
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
            throw std::runtime_error("Trying to operate on VectorOp of different sizes!");
        }
        return static_cast<T&>(*this);
    }
    T & operator-= (const T & destra) {
        if (destra.lunghezza()== data_length) {
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]-=destra.elemento(i);
            }
        } else {
            throw std::runtime_error("Trying to operate on VectorOp of different sizes!");
        }
        return static_cast<T&>(*this);
    }
    T & operator*= (const T & destra) {
        if (destra.lunghezza()== data_length) {
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]*=destra.elemento(i);
            }
        } else {
            throw std::runtime_error("Trying to operate on VectorOp of different sizes!");
        }
        return static_cast<T&>(*this);
    }
    T & operator/= (const T & destra) {
        if (destra.lunghezza()== data_length) {
            for (unsigned int i=0;i<data_length;i++) {
                vdata[i]/=destra.elemento(i);
            }
        } else {
            throw std::runtime_error("Trying to operate on VectorOp of different sizes!");
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


    unsigned int lunghezza() const{
        return data_length;
    }

    TFLOAT elemento(unsigned int i)const{
        if (i<data_length) {
            return vdata[i];
        } else {
	    throw std::runtime_error("Out of range index");
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
    }

    TFLOAT * vdata;
    VectorOp(){
        vdata=0;
        data_length=0;
    }
    ~VectorOp(){
    delete [] vdata;

}
    unsigned int data_length;
};

#endif // OPERAZIONISULISTA_H
