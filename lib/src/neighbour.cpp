#include "neighbour.h"
#include <algorithm>
#include <cstring>

template <class T, class TType>
void Neighbours<T,TType>::update_neigh(const size_t timestep, bool sort) {
    TType cutoff2=cutoff*cutoff;
    for (size_t i=0;i<natoms*(max_neigh+1);i+=max_neigh+1){
       list[i]=0;
    }
    for (size_t iatom=0;iatom<natoms;++iatom) {
       for (size_t jatom=iatom+1; jatom<natoms; ++jatom) {
          TType x[3];
          TType min_img_dist2=t.d2_minImage(iatom,jatom,timestep,timestep,x);
          if (min_img_dist2 <= cutoff2) {
              // update neigh list of iatom and jatom
              size_t idx=iatom*(max_neigh+1);
              size_t jdx=jatom*(max_neigh+1);
              size_t ni=list[idx], nj=list[jdx];
              if (ni < max_neigh && nj < max_neigh) {
                  TType min_img_dist = sqrt(min_img_dist2);
                  list[idx+ni] = jatom;
                  list[jdx+nj] = iatom;
                  //save distance
                  size_t pidx=(iatom*max_neigh+list[idx])*4;
                  size_t pjdx=(jatom*max_neigh+list[jdx])*4;
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
        for (size_t iatom=0;iatom<natoms;++iatom){
           auto iter = get_neigh(iatom);
           size_t nneigh=iter.size();
           for (size_t i=0;i<nneigh;++i){
               tmp_idxs[i]=i;
           }
           const size_t offset=iatom*4;
           std::sort(tmp_idxs,tmp_idxs+nneigh,
                       [&](const size_t &a, const size_t &b){
                          return rpos[offset+a] < rpos[offset+b];
                       }
                     );
           //apply the permutation
           for (size_t i=0;i<nneigh;++i){
               for (size_t j=0;j<4;++j){
                   tmp_rpos[i*4+j]=rpos[offset+tmp_idxs[i]];
               }
           }
           std::memcpy(rpos+offset,tmp_rpos,nneigh*4*sizeof(TType));
           for (size_t i=0;i<nneigh;++i){
               tmp_idxs[i]=*(iter.begin()+tmp_idxs[i]);
           }
           std::memcpy(iter.begin(),tmp_idxs,nneigh*sizeof(size_t));
        }
        sorted=true;
    }
}

template <class T, class TType>
typename Neighbours<T,TType>::template NeighIterator<TType> Neighbours<T,TType>::get_sann(const size_t iatom) {
    size_t nneigh=get_neigh(iatom).size();
    if (nneigh <3) return NeighIterator<TType>{rpos,0};
    size_t sann_n=3;
    TType s = 0.0;
    for (size_t i=0;i<3;++i){
        s += rpos[(iatom*max_neigh+i)*4];
    }
    while (sann_n < nneigh) {
       if (s <= rpos[(iatom*max_neigh+sann_n)*4]*(sann_n-2)) break;
    }
    return NeighIterator<TType>{rpos+iatom*4,sann_n};
}

