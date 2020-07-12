#ifndef CALC_BUFFER_H
#define CALC_BUFFER_H

#include <unordered_map>
#include <vector>
#include <stack>

template <class Data, class K = int, class Container=std::unordered_map<K,size_t> >
class CalcBuffer {
public:
    CalcBuffer(const size_t & max_size_, const size_t & chunk_size_):data{nullptr} {
        set_max_size(max_size_, chunk_size_);
    }
    ~CalcBuffer(){
        delete [] data;
    }
    template< class Calculator, class ...Args >
    Data * buffer_calc(Calculator & c, const K & key, Args ... args){
        auto element=buffer.find(key);
        if (element!=buffer.end()) {
            return data+element->second;
        } else {
            if (buffer.size()>=max_size){
                //delete something
                if (!remove_last_discarded() && buffer.size()>0)
                    remove_item(buffer.begin()->first);
            }
            //get space for result
            size_t idx=add_item(key);
            c.calc(key,data+idx, args...);
            return data+idx;
        }
    }
    void discard(const K & key){
        discard_first.push_back(key);
    }
    void discard(const std::vector<K> & keys){
        discard_first.insert(discard_first.end(),keys.begin(),keys.end());
    }
    void set_max_size(const size_t &max_s, const size_t & chunk_s){
        max_size=max_s;
        chunk_size=chunk_s;
        delete [] data;
        data=new Data[max_size*chunk_size];
        free_idxs.resize(max_size);
        for (size_t i=0;i<max_size;++i){
            free_idxs[i]=i*chunk_size;
        }
    }
    bool remove_last_discarded(){
        if (discard_first.size()>0){
            remove_item(discard_first.back());
            discard_first.pop_back();
            return true;
        }
        return false;
    }
    size_t add_item(const K & key) {
        auto idx=free_idxs.back();
        free_idxs.pop_back();
        buffer.insert({key,idx});
        return idx;
    }
    void remove_item(const K & key) {
        free_idxs.push_back(buffer[key]);
        buffer.erase(key);
    }
private:
    size_t max_size,chunk_size;
    Container buffer;
    Data* data;
    std::vector<size_t> free_idxs;
    std::vector<K> discard_first;
};

#endif // CALC_BUFFER_H
