#include "centerofmassdiff.h"

CenterOfMassDiff::CenterOfMassDiff(Traiettoria *t, unsigned int nthreads, unsigned int skip) :
    CalcolaMultiThread (nthreads,skip), t{t},lista_alloc{0},starting_center{},ntype{0},zero{1.0e-10},maxiter{100},error{false}
{

}

void CenterOfMassDiff::reset(const unsigned int numeroTimestepsPerBlocco) {
    ntype=t->get_ntypes();
    ntimesteps=numeroTimestepsPerBlocco;
    lunghezza_lista=ntimesteps*3*ntype/skip;
    starting_center.resize(3*ntype);
    if (lunghezza_lista>lista_alloc) {
        delete [] lista;
        lista = new double [lunghezza_lista];
        lista_alloc=lunghezza_lista;
    }
}

void CenterOfMassDiff::calc_single_th(const unsigned int &start, const unsigned int &stop, const unsigned int &primo, const unsigned int &ith) noexcept {
    unsigned int ilista=(start-primo)/skip;
    std::valarray<double> dx,cn;
    std::valarray<unsigned int> sums;
    dx=starting_center;
    cn.resize((starting_center.size()));
    sums.resize(ntype);
    for (unsigned int itimestep=start;itimestep<stop;++itimestep) {
        double *box = t->scatola(itimestep);
        double cns=0.0;
        int iter=0;
        do {
            sums=0.0;
            cn=0.0;
            for (unsigned int iatom=0;iatom<t->get_natoms();++iatom) {
                double * coord = t->posizioni(itimestep,iatom);
                unsigned int itype=t->get_type(iatom);
                sums[itype]++;
                for (unsigned int idim1=0;idim1<3;idim1++){
                    double l=box[idim1*2+1]-box[idim1*2+0];
                    double new_coord=coord[idim1]-dx[3*itype+idim1]; //shift of the origin, according to the old center
                    double wrapped_coord=new_coord-std::round(new_coord/l)*l; // center of the cell is in zero
                    cn[itype*3+idim1]+=wrapped_coord;
                }
            }
            cns=0.0;
            for (unsigned int itype=0;itype<ntype;++itype){
                for (unsigned int idim1=0;idim1<3;++idim1){
                    dx[itype*3+idim1]=dx[itype*3+idim1]+cn[itype*3+idim1]/sums[itype];
                    cns+=cn[itype*3+idim1]/sums[itype];
                }
            }

        } while (cns>zero && ++iter<maxiter);
        if (iter>=maxiter) {
            error=true;
        }
        for (unsigned int itype=0;itype<ntype;++itype){
            for (unsigned int idim1=0;idim1<3;++idim1){
                lista[ilista*3*ntype+3*itype+idim1]=dx[itype*3+idim1];
            }
        }
        ++ilista;
    }
}

CenterOfMassDiff::~CenterOfMassDiff(){

}
