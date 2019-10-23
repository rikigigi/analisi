/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/


#include "heatc.h"
#include <cmath>
#include<thread>

HeatC::HeatC(Traiettoria *t, double sigma, unsigned int nthreads, unsigned int skip) :
    sigma{sigma}, nthreads{nthreads},skip{skip},t{t},ntimesteps{0}
{
    if (nthreads==0)
        nthreads=1;
    if (skip==0)
        skip=1;
}

void HeatC::reset(const unsigned int numeroTimestepsPerBlocco) {
    if (ntimesteps<= numeroTimestepsPerBlocco && t->get_ntypes()<=tipi_atomi){
        ntimesteps=numeroTimestepsPerBlocco;
        return;
    }
    delete [] lista;
    tipi_atomi=t->get_ntypes();
    ntimesteps=numeroTimestepsPerBlocco;
    lunghezza_lista=ntimesteps*3*3*tipi_atomi/skip;
    lista = new double[lunghezza_lista];
}

double calc0(const double & pre, const double & A,const double & B) {
    return pre*(B*B-A*A);
}

template <int C>
void calc(double * box,double * pos, const double * val, const double & sigma, double * res) {
    const double sigma3=sigma*sigma*sigma*2;
    for (unsigned int ival=0;ival<3;++ival){
        for (unsigned int icoord=0;icoord<3;++icoord) {
            if (pos[icoord]+sigma/2.0>box[icoord*2+1]) {
                res[3*ival+icoord]+=C*calc0(val[ival]/sigma3,pos[icoord]-sigma/2.0, box[icoord*2+1]);
                res[3*ival+icoord]+=C*calc0(val[ival]/sigma3,box[icoord*2+0],       pos[icoord]+sigma/2.0-box[icoord*2+1]+box[icoord*2+0]);
            } else if (pos[icoord]-sigma/2.0<box[icoord*2+0]){
                res[3*ival+icoord]+=C*calc0(val[ival]/sigma3,box[icoord*2+0],                                       pos[icoord]+sigma/2.0);
                res[3*ival+icoord]+=C*calc0(val[ival]/sigma3,pos[icoord]-sigma/2.0+box[icoord*2+1]-box[icoord*2+0], box[icoord*2+1]);
            } else { // no overlapping beetween atomic energy and box boundaries
                res[3*ival+icoord]+=C*calc0(val[ival]/sigma3,pos[icoord]-sigma/2.0,pos[icoord]+sigma/2.0);
            }
        }
    }
}

void HeatC::calcola(unsigned int primo) {

    unsigned int npassith=ntimesteps/skip/nthreads;
    std::vector<std::thread> threads;
    for (unsigned int ith=0;ith<nthreads;++ith){
        threads.push_back(std::thread([&,ith](){
            unsigned int start=primo+npassith*ith*skip;
            unsigned int stop=start+npassith*skip;
            if (ith==nthreads-1)
                stop=primo+ntimesteps;
            int ilista=(start-primo)/skip;
            for (unsigned int itimestep=start;itimestep<stop;itimestep+=skip){

                //more specific code here
                for (unsigned int i=0;i<3*3*tipi_atomi;++i) {
                    lista[ilista*tipi_atomi*3*3+i]=0.0;
                }

                //loop over atoms
                double *cellm = t->scatola(itimestep);
                double *cellp = t->scatola(itimestep+1);
                for (unsigned int iatom=0;iatom<t->get_natoms();++iatom){
                    double *rm = t->posizioni(itimestep,iatom);
                    double *rp = t->posizioni(itimestep+1,iatom);
                    double *vm = t->velocita(itimestep,iatom);
                    double *vp = t->velocita(itimestep+1,iatom);
                    calc<-1>(cellm,rm,vm,sigma,&lista[ilista*tipi_atomi*3*3+t->get_type(iatom)*3*3]);
                    calc<1>(cellp,rp,vp,sigma,&lista[ilista*tipi_atomi*3*3+t->get_type(iatom)*3*3]);
                }

                ilista++;
            }
        }));
    }

    for (auto & t : threads){
        t.join();
    }
    threads.clear();
}

std::vector<ssize_t> HeatC::get_shape() const{
    return {ntimesteps/skip,tipi_atomi,3,3};
}
std::vector<ssize_t> HeatC::get_stride() const{
    return {static_cast<long>(sizeof(double)*tipi_atomi*3*3),
                static_cast<long>(sizeof(double)*3*3),
                sizeof(double)*3,
                sizeof(double)};
}

HeatC::~HeatC(){

}

HeatC & HeatC::operator =(const HeatC & destra) {
    OperazioniSuLista<HeatC>::operator =(destra);

}
