/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include <cmath>
#include<thread>
#include <vector>
#include <fstream>
#include <sstream>
#include "msd.h"
#include "config.h"
#include "traiettoria.h"
#include "floating_exceptions.h"

#ifdef USE_MPI
#include "mp.h"
#endif
template <class T,bool FPE>
MSD<T,FPE>::MSD(T *t, unsigned int skip, unsigned int tmax, unsigned int nthreads, bool calcola_msd_centro_di_massa, bool calcola_msd_nel_sistema_del_centro_di_massa, bool debug) :
    traiettoria(t), lmax(tmax), cm_msd(calcola_msd_centro_di_massa),debug(debug),cm_self(calcola_msd_nel_sistema_del_centro_di_massa),ntypes{0} ,
    CMT{nthreads,skip,0}
{
    if (calcola_msd_centro_di_massa)
        f_cm=2;
    else
        f_cm=1;
}
template <class T,bool FPE>
unsigned int MSD<T,FPE>::numeroTimestepsOltreFineBlocco(unsigned int n_b){
    return (traiettoria->get_ntimesteps()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->get_ntimesteps()/(n_b+1)+1 : lmax;
}
template <class T,bool FPE>
void MSD<T,FPE>::reset(const unsigned int numeroTimestepsPerBlocco) {

    leff=(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    ntypes=traiettoria->get_ntypes();
    lunghezza_lista=leff*ntypes*f_cm;

    ntimesteps=numeroTimestepsPerBlocco;
    delete [] lista;
    lista=new double [lunghezza_lista];
}

template <class T, bool FPE>
void MSD<T,FPE>::calc_init(size_t primo) {

    if (leff+ntimesteps+primo > traiettoria->get_ntimesteps()){
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the average or the lenght of the time lag");
    } else if (leff+ntimesteps > traiettoria->get_nloaded_timesteps()) {
        std::stringstream ss;
        ss <<"there are not enough loaded timesteps inside the trajectory object. I need at least " << leff+ntimesteps << " timesteps to do the requested calculation";
        throw std::runtime_error(ss.str());
    }
}

template <class T,bool FPE>
void MSD<T,FPE>::calc_single_th(size_t begin, size_t ultimo, size_t primo, unsigned int ith) {
    uint64_t *cont=new uint64_t [ntypes*f_cm];

    for (size_t t=begin;t<ultimo;t++) {
        for (size_t i=0;i<ntypes*f_cm;i++){
            lista[ntypes*t*f_cm+i]=0.0;
            cont[i]=0;
        }
        for (size_t imedia=0;imedia<ntimesteps;imedia+=skip){

            auto fpem = FloatingPointExceptionManager<void,FPE>([=](int fpe)->void{
                std::cerr<<fpe << " Thread "<<ith<< ", timesteps "<<primo+imedia<< " and " <<primo+imedia+t<< " raised a floating point exception or a NaN was found"<<std::endl;
            });
            if (cm_self){
                for (size_t iatom=0;iatom<traiettoria->get_natoms();iatom++) {
                    size_t itype=traiettoria->get_type(iatom);
                    double delta=(pow(
                                   traiettoria->template posizioni<false>(primo+imedia,iatom)[0]-traiettoria->template posizioni<false>(primo+imedia+t,iatom)[0]
                                  -(traiettoria->template posizioni_cm<false>(primo+imedia,itype)[0]-traiettoria->template posizioni_cm<false>(primo+imedia+t,itype)[0])
                            ,2)+
                            pow(
                                traiettoria->template posizioni<false>(primo+imedia,iatom)[1]-traiettoria->template posizioni<false>(primo+imedia+t,iatom)[1]
                                -(traiettoria->template posizioni_cm<false>(primo+imedia,itype)[1]-traiettoria->template posizioni_cm<false>(primo+imedia+t,itype)[1])
                            ,2)+
                            pow(
                                traiettoria->template posizioni<false>(primo+imedia,iatom)[2]-traiettoria->template posizioni<false>(primo+imedia+t,iatom)[2]
                                -(traiettoria->template posizioni_cm<false>(primo+imedia,itype)[2]-traiettoria->template posizioni_cm<false>(primo+imedia+t,itype)[2])
                            ,2))
                            -lista[ntypes*t*f_cm+traiettoria->get_type(iatom)];
                    lista[ntypes*t*f_cm+traiettoria->get_type(iatom)]+=delta/(++cont[traiettoria->get_type(iatom)]);

                }
            }else{
                for (size_t iatom=0;iatom<traiettoria->get_natoms();iatom++) {
                    double delta=(pow(traiettoria->template posizioni<false>(primo+imedia,iatom)[0]-traiettoria->template posizioni<false>(primo+imedia+t,iatom)[0],2)+
                            pow(traiettoria->template posizioni<false>(primo+imedia,iatom)[1]-traiettoria->template posizioni<false>(primo+imedia+t,iatom)[1],2)+
                            pow(traiettoria->template posizioni<false>(primo+imedia,iatom)[2]-traiettoria->template posizioni<false>(primo+imedia+t,iatom)[2],2))
                            -lista[ntypes*t*f_cm+traiettoria->get_type(iatom)];
                    lista[ntypes*t*f_cm+traiettoria->get_type(iatom)]+=delta/(++cont[traiettoria->get_type(iatom)]);
                    if constexpr (FPE) fpem.check_nan(lista[ntypes*t*f_cm+traiettoria->get_type(iatom)]);

                }
            }
            if constexpr (FPE) {
                for (size_t itype=0;itype<ntypes;++itype){
                    fpem.check_nan(lista[ntypes*t*f_cm+itype]);
                }
            }
            if (cm_msd) {
                for (size_t itype=0; itype < ntypes; itype++) {
                double delta=(pow(traiettoria->template posizioni_cm<false>(primo+imedia,itype)[0]-traiettoria->template posizioni_cm<false>(primo+imedia+t,itype)[0],2)+
                        pow(traiettoria->template posizioni_cm<false>(primo+imedia,itype)[1]-traiettoria->template posizioni_cm<false>(primo+imedia+t,itype)[1],2)+
                        pow(traiettoria->template posizioni_cm<false>(primo+imedia,itype)[2]-traiettoria->template posizioni_cm<false>(primo+imedia+t,itype)[2],2))
                        -lista[ntypes*t*f_cm+ntypes+itype];
                    lista[ntypes*t*f_cm+ntypes+itype]+=delta/(++cont[ntypes+itype]);
                    if constexpr (FPE) fpem.check_nan(lista[ntypes*t*f_cm+ntypes+itype]);
               }
            }
        }
    }
    delete [] cont;
}


template <class T,bool FPE>
void MSD<T,FPE>::calc_end() {

    if (debug) {
#ifndef USE_MPI
        std::ofstream out("msd.dump",std::ios::app);
#else
        std::ofstream out(Mp::mpi().outname("msd.dump"));
#endif
        for (size_t ts=0;ts<leff;ts++) {
            out << ts;
            for (unsigned int itype=0;itype<ntypes*f_cm;itype++){
                out <<" "<<lista[ntypes*ts*f_cm+itype];
            }
            out << "\n";
        }
        out << "\n\n";
    }
}

template <class T,bool FPE>
MSD<T,FPE> & MSD<T,FPE>::operator=(const MSD<T,FPE> &destra) {
    OperazioniSuLista<MSD<T,FPE> >::operator =( destra);
    return *this;
}

#ifdef BUILD_MMAP
template class MSD<Traiettoria,true>;
template class MSD<Traiettoria,false>;
#endif
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class MSD<Traiettoria_numpy,true>;
template class MSD<Traiettoria_numpy,false>;
#endif
