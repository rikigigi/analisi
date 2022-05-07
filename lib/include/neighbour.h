#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <iterator>
#include <vector>
#include <tuple>

#include "macros.h"

template <class T, class TType = double>
class Neighbours {
public:
    using ListSpec = std::vector<std::tuple<size_t, TType, TType> >;
    using TType4 = TType[4];
    template<class TN>
    class NeighIterator {
    public:
        NeighIterator(TN *idxs, const size_t len) : detail{idxs,len} {}
        const TN * begin() {return detail.idxs;}
        TN * begin_w() {return detail.idxs;}
        const TN * end() {return detail.idxs+detail.len;}
        size_t size() const {return detail.len;}
    private:
        struct {
            TN * idxs;
            const size_t len;
        } detail;
   
    };

    Neighbours(T *t, const ListSpec nneigh_cut2_skin2  ) :
        t{*t}, nneigh_cut_skin{nneigh_cut2_skin2} ,sorted{false},
        natoms{Neighbours::t.get_natoms()},
        ntypes{Neighbours::t.get_ntypes()},
        info{nneigh_cut_skin,natoms,ntypes}
    {
        if (ntypes != nneigh_cut_skin.size()) {
            throw std::runtime_error("In neighbours list you must specify parameters for each atomic type" AT);
        }
        list = new size_t[natoms*(info.nneigh+ntypes)];
        rpos = new TType[natoms*4*info.nneigh];
        tmp_idxs = new size_t[info.max_nneigh];
        tmp_rpos = new TType[info.max_nneigh*4];
    }
    ~Neighbours(){
       delete [] list;
       delete [] rpos;
       delete [] tmp_idxs;
       delete [] tmp_rpos;
    }

    void update_neigh(const size_t timestep,bool sort);
    NeighIterator<size_t> get_neigh(const size_t iatom, const size_t jtype) const{
        return NeighIterator<size_t>{list+ info.list_offset[jtype]+ iatom*(nneigh(jtype)+1)+1, list[info.list_offset[jtype]+ iatom*(nneigh(jtype)+1)]};
    }
    NeighIterator<TType4> get_neigh_r(const size_t iatom, const size_t jtype) const{
        return NeighIterator<TType4>{ (TType4 *) (rpos+ info.rpos_offset[jtype]+ iatom*nneigh(jtype)*4), list[info.list_offset[jtype]+ iatom*(nneigh(jtype)+1)]};
    }
    size_t get_sann_n(const size_t iatom, const size_t jtype) const;

    NeighIterator<size_t> get_sann(const size_t iatom, const size_t jtype) const{
        return NeighIterator<size_t>{list+ info.list_offset[jtype]+ iatom*(nneigh(jtype)+1)+1,get_sann_n(iatom,jtype)};
    }
    NeighIterator<TType4> get_sann_r(const size_t iatom, const size_t jtype) const{
    return NeighIterator<TType4>{(TType4*)(rpos+info.rpos_offset[jtype]+iatom*nneigh(jtype)*4),get_sann_n(iatom,jtype)};
}

private:


    T & t;
    const ListSpec nneigh_cut_skin;

    TType cutoff,skin, *rpos,*tmp_rpos;
    bool sorted;
    const size_t natoms,ntypes;
    size_t max_neigh, *tmp_idxs, *list; // first #of neigh, then index list for each atom
    //constants: strides, sizes, etc calculated once and left there forever
    const struct _tot_nneigh {
        _tot_nneigh(const ListSpec ncss, const size_t natoms, const size_t ntypes): nneigh{0},max_nneigh{0} {
            //total sizes
            size_t offset=0;
            size_t offset_r=0;
            for (const auto & ncs : ncss) {
                list_offset.push_back(offset);
                rpos_offset.push_back(offset_r);
                nneigh += std::get<0>(ncs);
                if (std::get<0>(ncs) > max_nneigh) max_nneigh=std::get<0>(ncs);
                offset += (std::get<0>(ncs)+1)*natoms;
                offset_r += std::get<0>(ncs)*natoms*4;
            }
        }
        size_t nneigh; // total number of neighbours per atom, for all types
        size_t max_nneigh; //maximum of the number of neighbours per atom over for all types
        std::vector<size_t> list_offset; // where the list starts for type i
        std::vector<size_t> rpos_offset; // where the rpos starts for type i
    } info;
    TType cutoff2(size_t itype) const {return std::get<1>(nneigh_cut_skin[itype]);}
    size_t nneigh(size_t itype) const {return std::get<0>(nneigh_cut_skin[itype]);}
};



#endif // MEIGHBOUR_H
