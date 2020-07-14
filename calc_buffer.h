#ifndef CALC_BUFFER_H
#define CALC_BUFFER_H

#include <unordered_map>
#include <vector>
#include <deque>
#include <algorithm>

template <class Data, class K = int, class Container=std::unordered_map<K,size_t> >
class CalcBuffer {
public:
    CalcBuffer(const size_t & max_size_, const size_t & chunk_size_, const size_t & n_valid_last=2):data{nullptr},n_valid_last{n_valid_last} {
        if (n_valid_last>=max_size_) {
            throw std::runtime_error("cannot ask for more elements than the size of the buffer");
        }
        set_max_size(max_size_, chunk_size_);
    }
    ~CalcBuffer(){
        delete [] data;
    }
    template< class Calculator, class ...Args >
    Data * buffer_calc(Calculator & c, const K & key, Args ... args){
        requested_keys.push_back(key);
        if (requested_keys.size()>n_valid_last){
            requested_keys.pop_front();
        }
        auto element=buffer.find(key);
        if (element!=buffer.end()) {
            return data+element->second;
        } else {
            if (buffer.size()>=max_size){
                //delete something
                if (!remove_last_discarded() && buffer.size()>0){
                    auto rkey=buffer.begin();
                    while ( std::find(requested_keys.begin(),requested_keys.end(),rkey->first ) != requested_keys.end()){
                        ++rkey;
                        if (rkey==buffer.end()){
                            throw std::runtime_error("error: cannot find a element to delete in the buffer");
                        }
                    }
                    remove_item(rkey->first);
                }
            }
            //get space for result
            size_t idx=add_item(key);
            c.calc(key,data+idx, args...);
            return data+idx;
        }
    }
    void discard(){
        discard_first.clear();
        for (auto & b : buffer) {
            discard_first.push_back(b.first);
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
            while (discard_first.size()>0 && (buffer.find(discard_first.back())==buffer.end() || std::find(requested_keys.begin(),requested_keys.end(),discard_first.back() ) != requested_keys.end())){
                discard_first.pop_back();
            }

            if (discard_first.size()==0) return false;
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
    size_t max_size,chunk_size,n_valid_last;
    Container buffer;
    Data* data;
    std::vector<size_t> free_idxs;
    std::vector<K> discard_first;
    std::deque<K> requested_keys;
};

#endif // CALC_BUFFER_H
