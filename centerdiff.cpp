#include "centerdiff.h"
#include <array>
#include <cmath>

CenterDiff::CenterDiff(Traiettoria *t, unsigned int nthreads, unsigned int skip, unsigned int nit) : CalcolaMultiThread (nthreads,skip), t{t}, nit{nit},lista_alloc{0}
{
    if (nit==0) nit=1;
}

void CenterDiff::reset(const unsigned int numeroTimestepsPerBlocco) {
    lunghezza_lista=numeroTimestepsPerBlocco*3*3*nit/skip;
    ntimesteps=numeroTimestepsPerBlocco;
    if (lunghezza_lista> lista_alloc){
        delete [] lista;
        lista = new double[lunghezza_lista];
        lista_alloc=lunghezza_lista;
    }
}



void CenterDiff::calc_single_th(const unsigned int &start, const unsigned int &stop, const unsigned int &primo, const unsigned int &ith) {
    unsigned int ilista=(start-primo)/skip;
    std::array<double,3*3> dx,cn;
    std::array<double,3> sums;
    dx.fill(0.0);
    for (unsigned int itimestep=start;itimestep<stop;++itimestep) {
        double *box = t->scatola(itimestep);
        for (unsigned int i=0;i<nit;++i) {
            sums.fill(0.0);
            cn.fill(0.0);
            for (unsigned int iatom=0;iatom<t->get_natoms();++iatom) {
                double * coord = t->posizioni(itimestep,iatom);
                double * vel = t->velocita(itimestep,iatom);
                for (unsigned int idim1=0;idim1<3;idim1++){
                    double l=box[idim1*2+1]-box[idim1*2+0];
                    sums[idim1]+=vel[idim1];
                    for (unsigned int idim2=0;idim2<3;idim2++) {
                        double new_coord=coord[idim1]-dx[3*idim2+idim1];
                        double wrapped_coord=new_coord-std::round(new_coord/l)*l;
                        cn[idim2*3+idim1]+=wrapped_coord*vel[idim2];
                    }
                }
            }
            for (unsigned int idim1=0;idim1<3;++idim1){
                for (unsigned int idim2=0;idim2<3;++idim2){
                    dx[idim2*3+idim1]=dx[idim2*3+idim1]+cn[idim2*3+idim1]/sums[idim2];
                    lista[ilista*3*3*nit+i*3*3+idim2*3+idim1]=dx[idim2*3+idim1];
                }
            }
        }
        ++ilista;
    }
}

CenterDiff::~CenterDiff() {

}
