#include "msd.h"
#include <cmath>

MSD::MSD(Traiettoria *t, unsigned int skip, unsigned int tmax) :
    traiettoria(t), lmax(tmax),skip(skip)
{
}

unsigned int MSD::numeroTimestepsOltreFineBlocco(unsigned int n_b){
    return (traiettoria->get_ntimesteps()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->get_ntimesteps()/(n_b+1)+1 : lmax;
}

void MSD::reset(const unsigned int numeroTimestepsPerBlocco) {

    leff=(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    lunghezza_lista=leff*traiettoria->get_ntypes();
    ntimesteps=numeroTimestepsPerBlocco;
    delete [] lista;
    lista=new double [lunghezza_lista];
}


void MSD::calcola(unsigned int primo) {

    unsigned int *cont=new unsigned int [traiettoria->get_ntypes()];

    for (unsigned int t=0;t<leff;t++) {
        for (unsigned int i=0;i<traiettoria->get_ntypes();i++){
            lista[traiettoria->get_ntypes()*t+i]=0.0;
            cont[i]=0;
        }
        for (unsigned int imedia=0;imedia<ntimesteps;imedia+=skip){
            for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
                double delta=(pow(traiettoria->posizioni(primo+imedia,iatom)[0]-traiettoria->posizioni(primo+imedia+t,iatom)[0],2)+
                             pow(traiettoria->posizioni(primo+imedia,iatom)[1]-traiettoria->posizioni(primo+imedia+t,iatom)[1],2)+
                             pow(traiettoria->posizioni(primo+imedia,iatom)[2]-traiettoria->posizioni(primo+imedia+t,iatom)[2],2))/3.0
                        -lista[traiettoria->get_ntypes()*t+traiettoria->get_type(iatom)];
                lista[traiettoria->get_ntypes()*t+traiettoria->get_type(iatom)]+=delta/(++cont[traiettoria->get_type(iatom)]);

            }
        }
    }
    delete [] cont;
}

MSD & MSD::operator=(const MSD & destra) {
    OperazioniSuLista<MSD>::operator =( destra);
    return *this;
}
