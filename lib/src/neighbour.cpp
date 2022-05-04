#include "neighbour.h"
#include <algorithm>
#include <cstring>
#include <cmath>
#include "config.h"

template <class T, class TType>
void Neighbours<T,TType>::update_neigh(const size_t timestep, bool sort) {
    for (size_t i=0;i<natoms*(max_neigh+1);i+=max_neigh+1){
       list[i]=0;
    }
    for (size_t iatom=0;iatom<natoms;++iatom) {
       size_t itype = t.get_type(iatom);
       for (size_t jatom=iatom+1; jatom<natoms; ++jatom) {
          size_t jtype = t.get_type(jatom);
          TType x[3];
          TType min_img_dist2=t.d2_minImage(iatom,jatom,timestep,timestep,x); // xi - xj
          if (min_img_dist2 <= cutoff2(jtype)) {
              // update neigh list of iatom and jatom
              size_t idx=info.list_offset[jtype] + iatom*(nneigh(jtype)+1); //i atom in the center, j of type jtype around
              size_t jdx=info.list_offset[itype] + jatom*(nneigh(itype)+1); // the j atom in the center, i of type itype around
              size_t ni=list[idx], nj=list[jdx];
              if (ni < nneigh(jtype) && nj < nneigh(itype)) {
                  TType min_img_dist = sqrt(min_img_dist2);
                  list[idx+ni] = jatom;
                  list[jdx+nj] = iatom;
                  //save distance
                  size_t pidx=info.rpos_offset[jtype] + (iatom*nneigh(jtype)+list[idx])*4;
                  size_t pjdx=info.rpos_offset[itype] + (jatom*nneigh(itype)+list[jdx])*4;
                  rpos[pidx]=min_img_dist;
                  rpos[pjdx]=min_img_dist;
                  for (int i=0;i<3;++i){
                      rpos[pidx+1+i]=x[i];
                      rpos[pjdx+1+i]=-x[i];
                  }
                  list[idx]++;
                  list[jdx]++;
              } else {
                  throw std::runtime_error("Too many neighbours in shell!");
              }
          }
       }
    }
    if (sort){
        //sort neighs by distance
        for (size_t jtype=0;jtype<ntypes;++jtype){
            for (size_t iatom=0;iatom<natoms;++iatom){
                auto iter = get_neigh(iatom,jtype);
                size_t nneigh_i=iter.size();
                for (size_t i=0;i<nneigh_i;++i){
                    tmp_idxs[i]=i;
                }
                const size_t offset=info.rpos_offset[jtype]+iatom*4;
                std::sort(tmp_idxs,tmp_idxs+nneigh_i,
                          [&](const size_t &a, const size_t &b){
                    return rpos[offset+a] < rpos[offset+b];
                }
                );
                //apply the permutation
                for (size_t i=0;i<nneigh_i;++i){
                    for (size_t j=0;j<4;++j){
                        tmp_rpos[i*4+j]=rpos[offset+tmp_idxs[i]];
                    }
                }
                std::memcpy(rpos+offset,tmp_rpos,nneigh_i*4*sizeof(TType));
                for (size_t i=0;i<nneigh_i;++i){
                    tmp_idxs[i]=*(iter.begin()+tmp_idxs[i]);
                }
                std::memcpy(iter.begin_w(),tmp_idxs,nneigh_i*sizeof(size_t));
            }
        }
    }
    sorted=true;
}

template <class T, class TType>
typename Neighbours<T,TType>::template NeighIterator<TType> Neighbours<T,TType>::get_sann(const size_t iatom, const size_t jtype) const  {
    size_t nneigh_i=get_neigh(iatom,jtype).size();
    if (nneigh_i <3) return NeighIterator<TType>{rpos,0};
    size_t sann_n=3;
    TType s = 0.0;
    for (size_t i=0;i<3;++i){
        s += rpos[info.rpos_offset[jtype]+(iatom*nneigh(jtype)+i)*4];
    }
    while (sann_n < nneigh_i) {
       if (s <= rpos[info.rpos_offset[jtype]+(iatom*nneigh(jtype)+sann_n)*4]*(sann_n-2)) break;
    }
    return NeighIterator<TType>{rpos+info.rpos_offset[jtype]+iatom*nneigh(jtype)*4,sann_n*4};
}


#ifdef BUILD_MMAP
#include "traiettoria.h"
template class Neighbours<Traiettoria,double>;
#endif
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class Neighbours<Traiettoria_numpy,double>;
#endif
