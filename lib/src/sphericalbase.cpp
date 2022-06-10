#include "sphericalbase.h"
#include "specialfunctions.h"
#include <sstream>

template <int lmax, class TFLOAT, class T>
SphericalBase<lmax,TFLOAT,T>::SphericalBase(T *t, const size_t nbin, const Rminmax_t rminmax):
    t{*t}, natoms{t->get_natoms()},ntypes{static_cast<size_t>(t->get_ntypes())},nbin{nbin},
    detail{rminmax,(TFLOAT) nbin}
{
    if (ntypes*ntypes != rminmax.size()) {
        std::stringstream ss;
        ss << "you must provide a radial range for each pair of atomic types, in total ntypes*ntypes pair of numbers. You provided "<<
              rminmax.size() << " elements while ntypes is " << ntypes <<" .";
        throw std::runtime_error(ss.str());
    }
}


template <int lmax, class TFLOAT, class T>
void SphericalBase<lmax,TFLOAT,T>::calc(int timestep,
                                                TFLOAT *result,
                                                TFLOAT *workspace,//workspace array
                                                TFLOAT * cheby,//workspace array,
                                                int * counter, // counter of #of atoms in the bin, for every atom
                                                Neighbours_T * nnl
                                                ) const {
    //zero result
    for (int i=0;i<(lmax+1)*(lmax+1)*natoms*nbin*ntypes;++i) {
        result[i]=0;
    }
    if (counter){
        for (int i=0;i<natoms*ntypes*nbin;++i){
            counter[i]=0;
        }
    }

    if (nnl == nullptr){
        for (unsigned int iatom=0;iatom<natoms;iatom++) {
            //other atom loop
            const unsigned int itype=t.get_type(iatom);
            for (unsigned int jatom=0;jatom<natoms;jatom++) {
                const unsigned int jtype=t.get_type(jatom);
                const unsigned int pairidx=ntypes*itype+jtype;
                if (iatom==jatom)
                    continue;

                //minimum image distance
                double x[3];
                double d=sqrt(t.d2_minImage(iatom,jatom,timestep,timestep,x));
                //bin index
                int idx=(int)floorf((d-detail.rminmax[pairidx].first)/detail.dr[pairidx]);

                if (idx<nbin && idx >= 0){
                    //calculate sin and cos
                    //calculate spherical harmonics and add to the correct average
                    SpecialFunctions::SphericalHarmonics<lmax,TFLOAT,true,true> sh(x[0],x[1],x[2],cheby,workspace);
                    sh.calc();
                    //now in workspace you have all the spherical harmonics components of the density of the current jatom around iatom
                    //add to the sh density of the current iatom of the current bin of the type jtype
                    for (int ll=0;ll<(lmax+1)*(lmax+1);++ll) {
                        result[index_wrk(iatom,jtype,idx)+ll]+=workspace[ll];
                    }
                    if (counter){
                        counter[index_wrk_counter(iatom,jtype,idx)]++;
                    }
                }
            }
        }
    } else {
        nnl->update_neigh(timestep,true); //sorted list of neighbours, to use sann algorithm
        for (unsigned int iatom=0;iatom<natoms;iatom++) {
            const unsigned int itype=t.get_type(iatom);
            for (unsigned int jtype=0;jtype<ntypes;jtype++) {
                //other atom loop
                auto position_iterator = nnl->get_sann_r(iatom,jtype);
                const unsigned int pairidx=ntypes*itype+jtype;
                for (const auto & x : position_iterator){

                    //minimum image distance is already present in the neighbour list:
                    // x is: {r, r_x, r_y, r_z}

                    //calculate sin and cos
                    //calculate spherical harmonics and add to the correct average
                    SpecialFunctions::SphericalHarmonics<lmax,TFLOAT,true,true> sh(x[1],x[2],x[3],cheby,workspace);
                    sh.calc();
                    //now in workspace you have all the spherical harmonics components of the density of the current jatom around iatom
                    //add to the sh density of the current iatom of the current bin of the type jtype
                    for (int ll=0;ll<(lmax+1)*(lmax+1);++ll) {
                        result[index_wrk(iatom,jtype,0)+ll]+=workspace[ll];
                    }
                    if (counter){
                        counter[index_wrk_counter(iatom,jtype,0)]++;
                    }

                }
            }
        }
    }
}
#include "config.h"
#ifdef BUILD_MMAP
#include "traiettoria.h"
template class SphericalBase<6,double,Traiettoria>;
template class SphericalBase<8,double,Traiettoria>;
template class SphericalBase<10,double,Traiettoria>;
#endif
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class SphericalBase<6,double,Traiettoria_numpy>;
template class SphericalBase<8,double,Traiettoria_numpy>;
template class SphericalBase<10,double,Traiettoria_numpy>;
#endif
