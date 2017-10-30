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



#ifndef MEDIABLOCCHI_H
#define MEDIABLOCCHI_H

#include "traiettoria.h"
#include "assert.h"
#include <type_traits>
#include <iostream>
#include "config.h"
#include "cronometro.h"
#include "readlog.h"



template <class TR> class TraiettoriaF {
public:
    static void imposta_dimensione_finestra_accesso(unsigned int s_,TR* traiettoria ){}
    static void imposta_inizio_accesso(unsigned int s_,TR* traiettoria  ){}
    static unsigned int get_ntimesteps(TR* traiettoria ){}
};

template <> class TraiettoriaF<Traiettoria> {
public:
static void imposta_dimensione_finestra_accesso(unsigned int s_,Traiettoria* traiettoria ){
    traiettoria->imposta_dimensione_finestra_accesso(s_);
}
static void imposta_inizio_accesso(unsigned int s_,Traiettoria* traiettoria  ){
    traiettoria->imposta_inizio_accesso(s_);
}
static unsigned int get_ntimesteps(Traiettoria* traiettoria ){
    return traiettoria->get_ntimesteps();
}

};


template <> class TraiettoriaF<ReadLog> {
public:
static void imposta_dimensione_finestra_accesso(unsigned int s_,ReadLog* traiettoria ){

}
static void imposta_inizio_accesso(unsigned int s_,ReadLog* traiettoria  ){

}
static unsigned int get_ntimesteps(ReadLog* traiettoria ){
    return traiettoria->n_timestep();
}

};


template <class TR,class T, typename ... Args > class MediaBlocchiG
{
public:
    MediaBlocchiG(TR * t,
                 const unsigned int & numero_blocchi,
                 Args ... arg
                 ) {
//        static_assert(std::is_base_of<Calcolo,T>::value,"T deve essere derivato da Calcolo!");
        traiettoria=t;
        n_b=numero_blocchi;
        ok=false;

        Tmedio = new T(traiettoria, arg...);
        Tvar = new T(traiettoria,arg...);
        delta = new T(traiettoria,arg...);
        tmp=new T(traiettoria,arg...);
        int timestepsPerBlocco=(TraiettoriaF<TR>::get_ntimesteps(traiettoria)-Tmedio->numeroTimestepsOltreFineBlocco(n_b))/n_b;
        if(timestepsPerBlocco>0){
            s=timestepsPerBlocco;
            ok=true;
        } else {
            std::cerr<< "Impossibile dividere la traiettoria in "<<n_b<<"blocchi!\n";
        }

        calcolo = new T (traiettoria,arg...);
    }

    ~MediaBlocchiG() {
#ifdef DEBUG
        std::cerr << "~MediaBlocchiG(): Tmedio="<<Tmedio<<", Tvar="<<Tvar<< ", delta="<<delta<<", tmp="<<tmp<<".\n";
#endif
        delete Tmedio;
        delete Tvar;
        delete delta;
        delete tmp;
        delete calcolo;
    }

    void calcola() {
        if (!ok) return;
        Tmedio->reset(s);
        Tmedio->azzera();
        Tvar->reset(s);
        Tvar->azzera();
        delta->reset(s);
        tmp->reset(s);
        TraiettoriaF<TR>::imposta_dimensione_finestra_accesso(s+Tmedio->numeroTimestepsOltreFineBlocco(n_b),traiettoria);
        for (unsigned int iblock=0;iblock<n_b;iblock++) {
#ifdef DEBUG
            std::cerr << "calcolo->calcola(iblock*s);\n";
#endif
            cronometro cron;
            cron.start();
            calcolo->reset(s);
            TraiettoriaF<TR>::imposta_inizio_accesso(iblock*s,traiettoria);
            calcolo->calcola(iblock*s);
#ifdef DEBUG2
            std::cerr << "*delta = *calcolo - *Tmedio;\n";
#endif
            //algoritmo media e varianza online
            *delta = *calcolo;
            *delta -= *Tmedio;
            //*delta = *calcolo - *Tmedio
#ifdef DEBUG2
            std::cerr << "*delta = *calcolo - *Tmedio;\n";
#endif
            *tmp = *delta;
            *tmp /= (iblock+1);
            *Tmedio += *tmp;
            //*Tmedio += (*delta)/(iblock+1); // calcolo / double , calcolo+=calcolo
#ifdef DEBUG2
            std::cerr << "*Tvar += (*delta)*(*calcolo-*Tmedio);\n";
#endif
            *tmp = *calcolo;
            *tmp -= *Tmedio;
            *tmp *= *delta;
            *Tvar += *tmp;
            //*Tvar += (*delta)*(*calcolo-*Tmedio); // calcolo * calcolo
            cron.stop();
#ifdef DEBUG
            std::cerr << "Tempo cpu per il calcolo del blocco "<<iblock+1<<" su "<<n_b <<": "<< cron.time()<<"s.\n";
#endif
        }
        *Tvar/=((n_b-1)*n_b);


    }

    T * media() {return Tmedio;}
    T * varianza() {return Tvar;}
    T * puntatoreCalcolo() {return calcolo;}

private:
    unsigned int n_b,s; //numero di blocchi e dimensione
    T *Tmedio,*Tvar,*calcolo,*delta,*tmp;
    TR * traiettoria;
    bool ok;

};


template<class T, typename ... Args> using MediaBlocchi = MediaBlocchiG<Traiettoria,T,Args...>;



#endif // MEDIABLOCCHI_H
