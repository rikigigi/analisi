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
#include <vector>
#include <functional>
#include "calcoliblocchi.h"
#ifdef USE_MPI
#include "mp.h"
#endif



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


template <> class TraiettoriaF<ReadLog<> > {
public:
static void imposta_dimensione_finestra_accesso(unsigned int s_,ReadLog<>* traiettoria ){

}
static void imposta_inizio_accesso(unsigned int s_,ReadLog<>* traiettoria  ){

}
static unsigned int get_ntimesteps(ReadLog<>* traiettoria ){
    return traiettoria->n_timestep();
}

};

template <> class TraiettoriaF<ReadLog<long double> > {
public:
static void imposta_dimensione_finestra_accesso(unsigned int s_,ReadLog<long double>* traiettoria ){

}
static void imposta_inizio_accesso(unsigned int s_,ReadLog<long double>* traiettoria  ){

}
static unsigned int get_ntimesteps(ReadLog<long double>* traiettoria ){
    return traiettoria->n_timestep();
}

};


template <class TR,class T, typename ... Args > class MediaBlocchiG
{
public:
    MediaBlocchiG(TR * t,
                 const unsigned int & numero_blocchi
                 ) {
//        static_assert(std::is_base_of<Calcolo,T>::value,"T deve essere derivato da Calcolo!");
        traiettoria=t;
        n_b=numero_blocchi;
        ok=false;
    Tmedio=0;
    Tvar=0;
    delta=0;
    tmp=0;
    calcolo=0;

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
    template<class Calcolo> void calcola_custom (Calcolo * calc,
                                            Args ... arg) {


        calcolo = new T (traiettoria,arg...);
        int timestepsPerBlocco=(TraiettoriaF<TR>::get_ntimesteps(traiettoria)-calcolo->numeroTimestepsOltreFineBlocco(n_b))/n_b;
        if(timestepsPerBlocco>0){
            s=timestepsPerBlocco;
            ok=true;
        } else {
            std::cerr<< "Impossibile dividere la traiettoria in "<<n_b<<"blocchi!\n";
            abort();
        }
        calcolo->reset(s);
        calc->calcola_begin(s,calcolo);

        TraiettoriaF<TR>::imposta_dimensione_finestra_accesso(s+calcolo->numeroTimestepsOltreFineBlocco(n_b),traiettoria);
        cronometro cron;
#ifndef USE_MPI
        cron.set_expected(1.0/double(n_b));
        cron.start();
        for (unsigned int iblock=0;iblock<n_b;iblock++){

#ifdef DEBUG
            std::cerr << "calcolo->calcola(iblock*s);\n";
#endif

            calcolo->reset(s);
            TraiettoriaF<TR>::imposta_inizio_accesso(iblock*s,traiettoria);
            calcolo->calcola(iblock*s);

            calc->calcola(calcolo);
            cron.stop();
#ifdef DEBUG
            std::cerr << "Tempo per il calcolo del blocco "<<iblock+1<<" su "<<n_b <<": "<< cron.time_last()<<
                         "s. Tempo trascorso e rimanente: "<<cron.time()<<"s "<<cron.expected()<<"s.\n";
#endif
        }
#else // MPI
        int mpime=Mp::mpi().me(),mpisize=Mp::mpi().size();
        std::cerr << "Inizio del calcolo con " << mpisize << " processi MPI...\n";

        if (n_b%mpisize!=0){
            cron.set_expected(1.0/double(n_b/mpisize+1));
            std::cerr << "Attenzione: sto usando un numero di blocchi che non è multiplo del numero di processi MPI. Dei processi saranno inattivi per un certo periodo di tempo.\n";
        } else
            cron.set_expected(1.0/double(n_b/mpisize));
        cron.start();
        for (unsigned int iblock=0;iblock<n_b;iblock+=mpisize) {
            unsigned int ib=iblock+mpime;
            if (ib>=n_b){
                break;
            }
            calcolo->reset(s);
            TraiettoriaF<TR>::imposta_inizio_accesso(ib*s,traiettoria);
            calcolo->calcola(ib*s);
            cronometro cronmpi;
            unsigned int mpirecv=1;
            if (mpime==0){
                //ricevi i risultati degli altri e calcola ciò che è da calcolare.
                calc->calcola(calcolo);
                if (iblock+mpisize<n_b){
                    mpirecv=mpisize;
                }else{
                    mpirecv=n_b-iblock;
                    std::cerr << mpisize-mpirecv << " processi MPI inattivi per questa iterazione.\n";
                }
                for (unsigned int i=1;i<mpirecv;i++) {
                    cronmpi.start();
                    Mp::mpi().recv_root(calcolo,i);
                    cronmpi.stop();
                    calc->calcola(calcolo);
                }
            } else {
                Mp::mpi().send_to_root(calcolo);
            }
            cron.stop();
            std::cerr << "Tempo trascorso per i blocchi "<<iblock<<"-"<<iblock+mpirecv-1<<": "<<cron.time()<<"s ;  Tempo rimanente: "<<cron.expected()<<"s.\nTempo mpi recv: "<<cronmpi.time()<<"s\n";
            std::cerr.flush();
        }

#endif // MPI
        calc->calcola_end(n_b);
    }



    void calcola(Args ... arg) {
        Tmedio = new T(traiettoria, arg...);
        Tvar = new T(traiettoria,arg...);
        delta = new T(traiettoria,arg...);
        tmp=new T(traiettoria,arg...);

        MediaVar<T> media_var(Tmedio,Tvar,delta,tmp);

        calcola_custom<MediaVar<T> >(&media_var,arg...);

    }


    T * media() {if (ok) return Tmedio; else abort();}
    T * varianza() {if (ok) return Tvar; else abort();}
    T * puntatoreCalcolo() {if (ok) return calcolo; else abort();}

private:
    unsigned int n_b,s; //numero di blocchi e dimensione
    T *Tmedio,*Tvar,*calcolo,*delta,*tmp;
    TR * traiettoria;
    std::tuple<Args...> arg;
    bool ok;

};


template<class T, typename ... Args> using MediaBlocchi = MediaBlocchiG<Traiettoria,T,Args...>;



#endif // MEDIABLOCCHI_H
