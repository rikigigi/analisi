#ifndef NEIGHBOUR_H
#define NEIGHBOUR_H

#include <iterator>

template <class T, class TType = double>
class Neighbours {
public:
    template<class TN>
    class NeighIterator {
    public:
        NeighIterator(const TN *idxs, const size_t len) : detail{idxs,len} {}
        const TN * begin() {return detail.idxs;}
        const TN * end() {return detail.idxs+detail.len;}
        size_t size() const {return detail.len;}
    private:
        struct {
            const TN * idxs;
            const size_t len;
        } detail;
   
    };

    Neighbours(T *t, const TType cutoff, size_t max_neigh, const TType skin  ) :
        t{*t}, cutoff{cutoff}, skin{skin}, sorted{false} {
        natoms = Neighbours::t.get_natoms();
        list = new size_t[natoms*(max_neigh+1)];
        rpos = new TType[natoms*4*max_neigh];
        tmp_idxs = new size_t[max_neigh];
        tmp_rpos = new TType[max_neigh*4];
    }
    ~Neighbours(){
       delete [] list;
       delete [] rpos;
       delete [] tmp_idxs;
       delete [] tmp_rpos;
    }

    void update_neigh(const size_t timestep,bool sort);
    NeighIterator<size_t> get_neigh(const size_t iatom){
        return NeighIterator{list+iatom*(max_neigh+1)+1, list[iatom*(max_neigh+1)]};
    }
    NeighIterator<TType> get_sann(const size_t iatom);

private:
    T & t;
    TType cutoff,skin, *rpos,*tmp_rpos;
    bool sorted;
    size_t max_neigh,natoms, *tmp_idxs, *list; // first #of neigh, then index list for each atom

};



#endif // MEIGHBOUR_H
