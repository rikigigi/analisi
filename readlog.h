#ifndef READLOG_H
#define READLOG_H
#include <string>
#include "traiettoria.h"
#include <vector>

template <class TFLOAT=double>
class ReadLog
{
public:
    ReadLog(std::string filename, Traiettoria * t=0, unsigned int skip=1,unsigned int nthreads=0, unsigned int nbatch=200,
            std::vector<std::string> req_headers=std::vector<std::string>());
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
    int need_binary(std::vector<std::string> headers);
    void calc_currents(Traiettoria * t,unsigned int blocks);
private:
    Traiettoria * traiettoria;
    std::vector<std::string> headers;
    std::vector<TFLOAT > data;
    std::vector<unsigned int > timesteps;
    unsigned int skip,size,step_index,nbatch;
    size_t data_size,data_size_from_binary;
    bool if_only_numbers(std::string str);
    unsigned int natoms,nthreads;

    ///qui vengono memorizzati le stringhe delle correnti da calcolare e le cariche da utilizzare per calcolarle
    std::vector< std::pair<std::string,std::vector<TFLOAT> > > q_current_type;
    /// analizza la stringa che definisce la corrente da calcolare
    std::pair<std::string,std::vector<TFLOAT> >qs(std::string header);
    /// ritorna l'indice della corrente calcolata (da 0 a #correnti calcolate dalla traiettoria binaria)
    unsigned int get_calc_j_index(std::string header);

};

#endif // READLOG_H
