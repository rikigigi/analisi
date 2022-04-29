#ifndef ASSOCIATEDLEGENDREPOLY_H
#define ASSOCIATEDLEGENDREPOLY_H
#include <iostream>
#include <cmath>

#define PI 3.14159265358979323846264338327950288419716939937510

namespace SpecialFunctions {


template <int m,class T>
void inline chebyshevRecursiveAngleMultiplication(const T &x,const T &y, T * res) noexcept{
    res[0]=1.0; //cos(0*teta)
    res[1]=0.0; //sin(0*teta)
    res[2]=x; //cos(1*teta)
    res[3]=y; //sin(1*teta)
    for (int i=2;i<=m;++i){ //recursion for sin,cos (i*teta)
        int j=0;
        res[i*2+j]=2*x*res[(i-1)*2+j]-res[(i-2)*2+j];
        j=1;
        res[i*2+j]=2*x*res[(i-1)*2+j]-res[(i-2)*2+j];
    }
}
template <int m,class T>
void print_chebyshevRecursiveAngleMultiplication( T * res,std::ostream & out) {
    for (int i=0;i<=m;++i){
        out << "[cos("<<i<<"*teta)"<< ",sin("<<i<<"*teta)] = ["<<res[i*2]<<","<<res[i*2+1]<<"]"<<std::endl;
    }
}


template <class T>
constexpr T realSpericalHarmonics_coeff(long int l, long int m, long int cur_v, T value=T(1)) {
    if (cur_v==0) {
        return realSpericalHarmonics_coeff(l,m,1,value);
    }
    if (m==0) {
        return sqrt(T(2*l+1)/(4*PI));
    } else if (m<0) {
        return realSpericalHarmonics_coeff(l,-m,cur_v,value);
    }
    if (cur_v<=l+m) {
        return realSpericalHarmonics_coeff<T>(l,m,cur_v+1,value*cur_v);
    } else {
        return sqrt(T(2*l+1)/(4*PI)*2 / value)*pow(-1,m);
    }
}


namespace MultiValStatic {

template <int l, class T>
struct MultiVal {
    T val[l+1]; // only positive m + zero
    T val_minus[l]; //negative values of m
    MultiVal<l-1,T> val_;

    template <typename ...Args>
    void inline apply_f ( T f(int,int,T,Args...) , Args...args) noexcept {
        for (int m=0;m<l+1;++m) {
            val[m]=f(l,m,val[m],args...);
        }
        for (int m=0;m<l;++m) {
            val_minus[m]=f(l,-1-m,val_minus[m],args...);
        }
        val_.apply_f(f,args...);
    }

    void inline running_average_m_zero_vector(T* average, const T & counter) const {
        average[0]+=(val[0]-average[0])/counter;
        val_.running_average_m_zero_vector(average+1,counter);
    }

    void inline sum_in_m_zero() noexcept {
        for (int m=0;m<l;++m) { //negative m
            val[0]+=val_minus[m];
        }
        for (int m=0;m<l;++m) { //positives m
            val[0]+=val[m+1];
        }
        val_.sum_in_m_zero();
    }
    void inline copy_mplus_mminus() noexcept{
        for (int m=1;m<l+1;++m) {
            val_minus[m-1]=val[m];
        }
        val_.copy_mplus_mminus();
    }

    void print(std::ostream & out) {
        val_.print(out);
        for (int m=0;m<l;m++){
            out << "(l,m) = ("<<l<<","<<-(l-m) <<"): "<< val_minus[l-m-1]<<std::endl;
        }
        for (int m=0;m<l+1;m++){
            out << "(l,m) = ("<<l<<","<<m <<"): "<< val[m]<<std::endl;
        }
    }

    template<int j>
    constexpr MultiVal<j,T> & valm() {if constexpr (j==l) return *this; else return val_.template valm<j>();}

    T get_l_m0(int j) {
        if (l==j) {
            return val[0];
        } else {
            return val_.get_l_m0(j);
        }
    }

};

template <class T>
struct MultiVal<0,T> {
    T val[1];

    template <typename ...Args>
    inline void apply_f ( T f(int,int,T,Args...) , Args...args)  noexcept{
        val[0]=f(0,0,val[0],args...);
    }
    void inline running_average_m_zero_vector(T* average, const T & counter) const {
        average[0]+=(val[0]-average[0])/counter;
    }
    void inline sum_in_m_zero() noexcept {}
    void copy_mplus_mminus(){}
    void print(std::ostream & out) {
        out << "(l,m) = ("<<0<<","<<0 <<"): "<< val[0]<<std::endl;
    }
    template<int j>
    constexpr MultiVal<j,T> & valm()  noexcept{
        static_assert (j==0, "You specified a wrong value of j!" );
        return *this;
    }
    T get_l_m0(int j) {
        if (j==0) {
            return val[0];
        } else {
            throw std::runtime_error("Wrong value of l!");
        }
    }
};

}



namespace MultiValDynamic {

template <int l, class T>
struct MultiVal {
    T *val; // only positive m + zero
    T *val_minus; //negative values of m
    MultiVal<l-1,T> val_;

    void inline init(T * array) noexcept{
        val_minus=array;
        val=array+l;
        val_.init(val+l+1);
    }

    template <typename ...Args>
    void inline apply_f ( T f(int,int,T,Args...) , Args...args) noexcept {
        for (int m=0;m<l+1;++m) {
            val[m]=f(l,m,val[m],args...);
        }
        for (int m=0;m<l;++m) {
            val_minus[m]=f(l,-1-m,val_minus[m],args...);
        }
        val_.apply_f(f,args...);
    }

    void inline running_average_m_zero_vector(T* average, const T & counter) const {
        average[0]+=(val[0]-average[0])/counter;
        val_.running_average_m_zero_vector(average+1,counter);
    }

    void inline copy_mplus_mminus() noexcept{
        for (int m=1;m<l+1;++m) {
            val_minus[m-1]=val[m];
        }
        val_.copy_mplus_mminus();
    }

    void inline sum_in_m_zero() noexcept {
        for (int m=0;m<l;++m) { //negative m
            val[0]+=val_minus[m];
        }
        for (int m=0;m<l;++m) { //positives m
            val[0]+=val[m+1];
        }
        val_.sum_in_m_zero();
    }

    void print(std::ostream & out) {
        val_.print(out);
        for (int m=0;m<l;m++){
            out << "(l,m) = ("<<l<<","<<-(l-m) <<"): "<< val_minus[l-m-1]<<std::endl;
        }
        for (int m=0;m<l+1;m++){
            out << "(l,m) = ("<<l<<","<<m <<"): "<< val[m]<<std::endl;
        }
    }

    template<int j>
    constexpr MultiVal<j,T> & valm() {if constexpr (j==l) return *this; else return val_.template valm<j>();}

    T get_l_m0(int j) {
        if (l==j) {
            return val[0];
        } else {
            return val_.get_l_m0(j);
        }
    }


};

template <class T>
struct MultiVal<0,T> {
    T *val;

    void inline init(T * array) noexcept{
        val=array;
    }
    template <typename ...Args>
    inline void apply_f ( T f(int,int,T,Args...) , Args...args)  noexcept{
        val[0]=f(0,0,val[0],args...);
    }
    void inline running_average_m_zero_vector(T* average, const T & counter) const {
        average[0]+=(val[0]-average[0])/counter;
    }
    void inline sum_in_m_zero() noexcept {}
    void copy_mplus_mminus()noexcept{}
    void print(std::ostream & out) {
        out << "(l,m) = ("<<0<<","<<0 <<"): "<< val[0]<<std::endl;
    }
    template<int j>
    constexpr MultiVal<j,T> & valm()  noexcept{
        static_assert (j==0, "You specified a wrong value of j!" );
        return *this;
    }

    T get_l_m0(int j) {
        if (j==0) {
            return val[0];
        } else {
            throw std::runtime_error("Wrong value of l!");
        }
    }
};

}

//selector for memory model
template <int l, class T, bool dynamic>
struct MultiVal;
template <int l, class T>
struct MultiVal<l,T,false> : public MultiValStatic::MultiVal<l,T> {};
template <int l, class T>
struct MultiVal<l,T,true> : public MultiValDynamic::MultiVal<l,T> {};


template <int lmax,int l, int m, class T, bool dynamic>
struct AssociatedLegendrePoly
{
    static inline void calc1(const T &x, MultiVal<lmax,T,dynamic> & res) noexcept {
        /* l+1,l step:
         * Implements the following:
         * P^{l}_{l+1}=(2l+1) x P^l_l(x)
        */
        res.template valm<l>().val[m]=
                res.template valm<l-1>().val[m] *x*(2*m+1);
        if constexpr (l<lmax){
            AssociatedLegendrePoly<lmax,l+1,m,T,dynamic>::calc2(x,res); //start l+1,m recursion
        }
    }
    static inline void calc2(const  T &x, MultiVal<lmax,T,dynamic> & res) noexcept{
        /*
         * l+1,m recursion (that can go forever)
         * P^m_{l} = (x(2l - 1)P^m_{l-1} - (l + m - 1)P^m_{l-2}) / (l - m)
        */
        res.template valm<l>().val[m]=(res.template valm<l-1>().val[m] *x*(2*l-1) - res.template valm<l-2>().val[m] *(l+m-1))/
                (l-m)                                    ;
        if constexpr(l<lmax) {
            AssociatedLegendrePoly<lmax,l+1,m,T,dynamic>::calc2(x,res); //continue l+1,m recursion
        }

    }
};

template <int lmax, int l,class T, bool dynamic>
struct AssociatedLegendrePoly<lmax,l,l,T,dynamic> {
    static inline void calc (const T&x,  MultiVal<lmax,T,dynamic> & res) noexcept {
        /* l,l recursion:
         * Implements the following:
         * P^{l}_{l}=-(2l-1)\sqrt{1-x^2}P^{l-1}_{l-1}(x)
        */
        res.template valm<l>().val[l]= -(2*l-1)* sqrt(1-x*x) * res.template valm<l-1>().val[l-1];

        if constexpr (lmax>l){
            AssociatedLegendrePoly<lmax,l+1,l+1,T,dynamic>::calc(x,res);//continue l,l recursion
            AssociatedLegendrePoly<lmax,l+1,l,T,dynamic>::calc1(x,res);//start a new l+1,l recursion
        }
    }
};

template <int lmax, class T, bool dynamic>
struct AssociatedLegendrePoly<lmax,0,0,T,dynamic> {
    static inline void calc (const T&x,  MultiVal<lmax,T,dynamic> & res) noexcept {
        //begin with l=0 and go up till lmax
        //store the result while recursing in res
        res.template valm<0>().val[0]=1.0;
        if constexpr (lmax>0){
            AssociatedLegendrePoly<lmax,1,1,T,dynamic>::calc(x,res);//start l,l recursion
            AssociatedLegendrePoly<lmax,1,0,T,dynamic>::calc1(x,res);//start l+1,l recursion
        }
    }
};


template<int L, class T,bool plus=true,int l=L, bool dynamic>
inline T* get_multival (int j, MultiVal<L,T,dynamic> & val) noexcept {
    if (j==l) {
        if constexpr (plus)
                return val.template valm<l>().val;
        else{
            if constexpr (l>0){
                return val.template valm<l>().val_minus;
            } else {
                return nullptr;
            }
        }
    } else {
        if constexpr (l>0){
            return get_multival<L,T,plus,l-1>(j,val);
        } else {
            return nullptr;
        }
    }
}

template <int l,class T, bool dynamic,bool cartesian=false>
class SphericalHarmonics{
public:
    SphericalHarmonics(T x, T y, T z, T* cheby, T* result) noexcept : cheby{cheby} {
        if constexpr (dynamic) {
            val.init(result);
        }
        if constexpr(cartesian) {
            T rxy=sqrt(x*x+y*y);
            T r=sqrt(x*x+y*y+z*z);
            cost=z/r;
            sinp=y/rxy;
            cosp=x/rxy;
        } else {
            cost=x; sinp=y; cosp=z;
        }
    }
    inline
    T* get_l_mplus(int j) const noexcept{
        return get_multival<l,T,true>(j,val);
    }
    inline
    T* get_l_mminus(int j) const noexcept{
        return get_multival<l,T,false>(j,val);
    }
    inline
    MultiVal<l,T,dynamic> & get_val() noexcept {return val;}

    void inline calc() noexcept{
        AssociatedLegendrePoly<l,0,0,T,dynamic>::calc(cost,val);
        val.copy_mplus_mminus();
        chebyshevRecursiveAngleMultiplication<l>(cosp,sinp,cheby);
        val.apply_f(&mixall,cheby);
    }
private:
    T static inline mixall(int ll, int m, T val, T * cheby) noexcept {
        if (m>=0){
            return val*realSpericalHarmonics_coeff<T>(ll,m,ll-m+1)*cheby[2*m];
        } else {
            return val*realSpericalHarmonics_coeff<T>(ll,m,ll-abs(m)+1)*cheby[2*abs(m)+1];
        }
    }
    T * cheby;
    T cost,sinp,cosp;
    MultiVal<l,T,dynamic> val;
};

}




#endif // ASSOCIATEDLEGENDREPOLY_H
