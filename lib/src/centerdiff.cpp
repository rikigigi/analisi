#include "centerdiff.h"
#include <array>
#include <cmath>

CenterDiff::CenterDiff(Traiettoria *t, unsigned int nthreads, unsigned int skip, unsigned int nit, bool sum_first_two_and_ignore_vz,bool sum_1) :
    CalcolaMultiThread (nthreads,skip), t{t}, nit{nit},lista_alloc{0},starting_center{},sum_first_two_and_ignore_vz{sum_first_two_and_ignore_vz},sum_1{sum_1}
{
    if (nit==0) nit=1;
}

void CenterDiff::reset(const unsigned int numeroTimestepsPerBlocco) {
    data_length=numeroTimestepsPerBlocco*3*3*nit/skip;
    ntimesteps=numeroTimestepsPerBlocco;
    if (data_length> lista_alloc){
        delete [] vdata;
        vdata = new double[data_length];
        lista_alloc=data_length;
    }
}



void CenterDiff::calc_single_th(const unsigned int &start, const unsigned int &stop, const unsigned int &primo, const unsigned int &ith) noexcept {
    unsigned int ilista=(start-primo)/skip;
    std::array<double,3*3> dx,cn;
    std::array<double,3> sums;
    dx=starting_center;
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
                    if (sum_first_two_and_ignore_vz && idim1==2){
                        if (sum_1){
                            sums[idim1]+=1.0;
                        } else {
                            sums[idim1]+=vel[0]+vel[1];
                        }
                    } else {
                        sums[idim1]+=vel[idim1];
                    }
                    if (!sum_first_two_and_ignore_vz){
                        for (unsigned int idim2=0;idim2<3;idim2++) {
                            double new_coord=coord[idim1]-dx[3*idim2+idim1]; //shift of the origin
                            double wrapped_coord=new_coord-std::round(new_coord/l)*l; // center of the cell is in zero
                            cn[idim2*3+idim1]+=wrapped_coord*vel[idim2];
                        }
                    } else {
                        for (unsigned int idim2=0;idim2<2;idim2++) {
                            double new_coord=coord[idim1]-dx[3*idim2+idim1]; //shift of the origin
                            double wrapped_coord=new_coord-std::round(new_coord/l)*l; // center of the cell is in zero
                            cn[idim2*3+idim1]+=wrapped_coord*vel[idim2];
                        }
                        unsigned int idim2=2;
                        double new_coord=coord[idim1]-dx[3*idim2+idim1]; //shift of the origin
                        double wrapped_coord=new_coord-std::round(new_coord/l)*l; // center of the cell is in zero
                        if (sum_1){
                            cn[idim2*3+idim1]+=wrapped_coord;
                        } else {
                            cn[idim2*3+idim1]+=wrapped_coord*(vel[0]+vel[1]);
                        }

                    }
                }
            }
            for (unsigned int idim1=0;idim1<3;++idim1){
                for (unsigned int idim2=0;idim2<3;++idim2){
                    dx[idim2*3+idim1]=dx[idim2*3+idim1]+cn[idim2*3+idim1]/sums[idim2];
                    vdata[ilista*3*3*nit+i*3*3+idim2*3+idim1]=dx[idim2*3+idim1];
                }
            }
        }
        ++ilista;
    }
}

CenterDiff::~CenterDiff() {

}
