#ifndef READLOG_H
#define READLOG_H
#include <string>
#include "traiettoria.h"
#include <vector>

template <class TFLOAT=double>
class ReadLog
{
public:
    ReadLog(std::string filename,Traiettoria * t=0,unsigned int skip=1);
    ~ReadLog();
    TFLOAT * line(unsigned int index);
    unsigned int n_timestep(){return data.size()/data_size;}
    unsigned int n_data(){return headers.size()-1;}
    unsigned int timestep(unsigned int index);
    unsigned int get_natoms(){std::cerr << "!!Attenzione: per pigrizia il codice assume in alcune parti che ci sono 1728 atomi! Per favore modificami. (" << __FILE__ <<":" << __LINE__<<")\n";
        return 1728;}
    /**
     *  introduco una sintassi strana per calcolare al volo la corrente di carica dal file binario della traiettoria:
     * "#traj:JZ N q1 ... qN" --> calcola la corrente dalla traiettoria utilizzando le N cariche per gli atomi q1 ... qN
     * NOT IMPLEMENTED
    **/
    std::pair<unsigned int,bool> get_index_of(std::string header);
    bool need_binary(std::vector<std::string> headers);
    void set_traj(Traiettoria * t);
private:
    Traiettoria * traiettoria;
    std::vector<std::string> headers;
    std::vector<TFLOAT > data;
    std::vector<unsigned int > timesteps;
    unsigned int skip,size,step_index;
    size_t data_size;
    bool if_only_numbers(std::string str);
    unsigned int natoms;
    std::vector<double> qs(std::string header);

};

#endif // READLOG_H
