#include "neighbour.h"
#include <algorithm>
#include <cstring>
#include <cmath>
#include "config.h"
#include <sstream>

template <class T, class TType>
void Neighbours<T,TType>::update_neigh(const size_t timestep, bool sort) {
    for (size_t i=0;i<natoms*(info.nneigh+ntypes);++i){
       list[i]=0;
    }
    for (size_t iatom=0;iatom<natoms;++iatom) {
       size_t itype = t.get_type(iatom);
       for (size_t jatom=iatom+1; jatom<natoms; ++jatom) {
          size_t jtype = t.get_type(jatom);
          TType x[3];
          TType min_img_dist2=t.d2_minImage(iatom,jatom,timestep,timestep,x); // xi - xj
          if (min_img_dist2 <= cutoff2(jtype) || min_img_dist2 <= cutoff2(itype)) {
              // update neigh list of iatom and jatom
              auto iiter = get_neigh(iatom,jtype);
              auto jiter = get_neigh(jatom,itype);
              const size_t ni=iiter.size(), nj=jiter.size();
              if (ni < nneigh(jtype) && nj < nneigh(itype)) {
                  TType min_img_dist = sqrt(min_img_dist2);
                  iiter.begin_w()[ni] = jatom;
                  jiter.begin_w()[nj] = iatom;
                  //save distance
                  auto ipos = get_neigh_r(iatom,jtype);
                  auto jpos = get_neigh_r(jatom,itype);
                  ipos.begin_w()[ni][0]=min_img_dist;
                  jpos.begin_w()[nj][0]=min_img_dist;
                  for (int i=0;i<3;++i){
                      ipos.begin_w()[ni][1+i]=x[i];
                      jpos.begin_w()[nj][1+i]=-x[i];
                  }
                  if (min_img_dist2 <= cutoff2(jtype))
                    iiter.begin_w()[-1]++;
                  if( min_img_dist2 <= cutoff2(itype))
                    jiter.begin_w()[-1]++;
              } else {
                  std::stringstream ss;
                  ss <<"Too many neighbours in shell! (iatom="<< iatom << " jatom=" << jatom << " ni="<<ni << " nj="<<nj <<" itype="<<itype<<" jtype="<<jtype<<")";
                  throw std::runtime_error(ss.str());
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
                if (nneigh_i==0) continue;
                for (size_t i=0;i<nneigh_i;++i){
                    tmp_idxs[i]=i;
                }
                auto riter=get_neigh_r(iatom,jtype);
                std::sort(tmp_idxs,tmp_idxs+nneigh_i,
                          [&](const size_t &a, const size_t &b){
                    return riter.begin()[a][0] < riter.begin()[b][0];
                }
                );
                //apply the permutation
                for (size_t i=0;i<nneigh_i;++i){
                    TType * _t = (TType*) (riter.begin()[tmp_idxs[i]]);
                    for (size_t j=0;j<4;++j){
                        tmp_rpos[i*4+j]=_t[j];
                    }
                }
                std::memcpy(riter.begin_w(),tmp_rpos,nneigh_i*4*sizeof(TType));
                for (size_t i=0;i<nneigh_i;++i){
                    tmp_idxs[i]=iter.begin()[tmp_idxs[i]];
                }
                std::memcpy(iter.begin_w(),tmp_idxs,nneigh_i*sizeof(size_t));
            }
        }
    }
    sorted=true;
}


template <class T, class TType>
size_t Neighbours<T,TType>::get_sann_n(const size_t iatom, const size_t jtype) const  {
    auto nneigh_i=get_neigh_r(iatom,jtype);
    if (nneigh_i.size() <3) return 0;
    size_t sann_n=3;
    TType s = 0.0;
    for (size_t i=0;i<sann_n;++i){
        s += nneigh_i.begin()[i][0];
    }
    while (sann_n < nneigh_i.size()) {
       if (s <= nneigh_i.begin()[sann_n][0]*(sann_n-2)) break;
       s += nneigh_i.begin()[sann_n][0];
       sann_n++;
    }
    return sann_n;
}

#ifdef BUILD_MMAP
#include "traiettoria.h"
template class Neighbours<Traiettoria,double>;
#endif
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class Neighbours<Traiettoria_numpy,double>;
#endif
