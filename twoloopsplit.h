#ifndef TWOLOOPSPLIT_H
#define TWOLOOPSPLIT_H

#include <vector>
#include <stdexcept>
#include <cmath>
template <class T>
struct WorkerState{
    T idx1_begin;
    T idx1;
    T idx1_end;
    T idx2_begin;
    T idx2;
    T idx2_end;
};

template <class T>
class TwoLoopSplit {
public:
    /**
      * WARNING: nworkers can be modified by this constructor if the block size is too big
    **/
    TwoLoopSplit(T & nworkers_, const T size1_, const T skip1_,const T block1_, const T size2_, const T skip2_, const T block2_ ) : size1{size1_}, skip1{skip1_}, block1{block1_}, size2{size2_}, skip2{skip2_}, block2{block2_},nworkers{nworkers_} {
        if (block1<skip1 || block2<skip2) {
            throw std::runtime_error("cannot have a block smaller than skip");
        }
        T tot_points=0,points_per_worker=0;
        for (T idx2=0; idx2<size2+size1; idx2+=block2) {
            T min_begin=size1,max_end=0,tot_elements=0;
            for (T i=idx2;i<idx2+block2;++i){
                T begin1,end1,nelements;
                row_begin_end(i,begin1,end1,nelements);
                if (nelements == 0) continue;
                tot_elements+=nelements;
                if (begin1<min_begin){
                    min_begin=begin1;
                }
                if (end1>max_end){
                    max_end=end1;
                }
            }
            row_begin_end1.push_back({min_begin,max_end});
            n_elements.push_back(tot_elements);
            tot_points+=tot_elements;
        }
        points_per_worker=tot_points/nworkers;
        //greedy division of the work
        work.resize(nworkers);
        state.resize(nworkers);
        T batch_size=0,cur_worker_idx=0;
        for (T i = 0; i<n_elements.size();++i){
            //see if by adding the next number of elements the counter is nearer and greater to the limit than the current one
            T cur_size=n_elements[i];
            if (cur_size==0) continue;
            if (batch_size>=points_per_worker && cur_size+batch_size>points_per_worker && cur_worker_idx!=nworkers-1) {
                work[cur_worker_idx].second=batch_size;
                batch_size=0;
                cur_worker_idx++;
            }
            work[cur_worker_idx].first.push_back(i);
            batch_size+=cur_size;
        }
        work[cur_worker_idx].second=batch_size;
        if (cur_worker_idx<nworkers_-1) {
            nworkers=cur_worker_idx+1;
            nworkers_=nworkers;
        }
        //init worker state
        for (T i=0;i<nworkers;++i) {
            init_worker_state(i);
        }
    }


    /**
      * get next pair of points
      * if false the block is finished and the next result will be in the next block
      * if end becomes true, there is nothing more to do here, and this worker is done
    **/
    bool get_next_idx_pair(const T iworker, T & idx1, T &idx2,bool &end){
        auto & state_ = state[iworker];
        end=false;
        idx1=state_.idx1;
        idx2=state_.idx2;
        state_.idx2+=skip2;
        if (state_.idx2>=state_.idx2_end || state_.idx2 >= size2 + state_.idx1) {
            state_.idx1+=skip1;
            if (state_.idx1>state_.idx1_end) {
                end=!init_worker_state(iworker);
                return false;
            }
            if (state_.idx1>state_.idx2_begin) {
                state_.idx2=round_high2(state_.idx1,state_.idx1);
            } else {
                state_.idx2=round_high2(state_.idx1,state_.idx2_begin);
            }
        }
        return idx1/block1 == state_.idx1/block1;
    }
private:
    bool init_worker_state(const T iworker){
        if (work[iworker].first.size()==0) return false;
        WorkerState<T> state_;
        //get and pop blocks row idx
        auto idx=work[iworker].first.back();
        work[iworker].first.pop_back();
        state_.idx1_begin=row_begin_end1[idx].first;
        state_.idx1_end=row_begin_end1[idx].second;
        state_.idx1=state_.idx1_begin;
        state_.idx2_begin=idx*block2;
        state_.idx2_end=(idx+1)*block2;
        state_.idx2=round_high2(state_.idx1,state_.idx2_begin);
        state[iworker]=state_;
        return true;
    }

    T round_high2(const T idx1, const T idx2_begin){
        return idx2_begin + (skip2-(idx2_begin-idx1)%skip2)%skip2;
    }

    T round_low1(const T idx, const T skip) {
        return idx-idx%skip;
    }
    T round_high1(const T idx, const T skip) {
        return idx + (skip-idx%skip)%skip;
    }

    void row_begin_end(const T idx2, T & begin1, T & end1,T & nelements) {
        begin1=size1;
        end1=0;
        nelements=0;
        for (T i=0;i<size1;++i){
            if (element_in_range(i,idx2)){
                if (begin1>i) begin1=i;
                end1=i;
                ++nelements;
            }
        }
    }

    bool element_in_range(const T idx1, const T idx2) {
        if (idx1%skip1 != 0 || (idx2-idx1)%skip2 != 0) return false;
        return idx2 < size2+idx1 && idx2>=idx1;
    }

    const T size1,skip1,block1, size2,skip2,block2;
    T nworkers;
    std::vector<std::pair<T,T>> row_begin_end1; // both begin and end are included in the range!
    std::vector<T> n_elements;
    std::vector<std::pair<std::vector<T>,T > > work;
    std::vector<WorkerState<T>> state;
};

#endif // TWOLOOPSPLIT_H
