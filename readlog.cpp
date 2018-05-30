#include "readlog.h"
#include <fstream>
#include <sstream>
#include <limits>

template <class TFLOAT> ReadLog<TFLOAT>::ReadLog(std::string filename, Traiettoria *t, unsigned int skip):
    traiettoria(t),skip(skip)
{
    unsigned int lcont=0;
    std::string tmp,header;
    std::ifstream log(filename);
    step_index=std::numeric_limits< unsigned int>::max();

    //trova la prima linea fatta di soli numeri. Quella prima

    while (log.good()) {
        header=tmp;
        std::getline(log,tmp);
        lcont++;
        if (if_only_numbers(tmp)&&tmp.length()>1)
            break;
    }

    //trova i nomi delle colonne e la colonna del timestep
    std::string delim = " ";

     size_t start = 0;
     size_t end = header.find(delim);
     while (end != std::string::npos)
     {
         if (end!=start&& start<header.length()) headers.push_back(header.substr(start, end - start));
         start = end + delim.length();
         end = header.find(delim, start);
         if (headers.back()=="Step") step_index=headers.size()-1;
     }
     if (end!=start && start<header.length()) headers.push_back(header.substr(start, end - start));
     if (headers.back()=="Step") step_index=headers.size()-1;
#ifdef DEBUG
     std::cerr << "Indice degli header:\n";
     for(int i=0;i<headers.size();i++){
         std::cerr << i << " " << headers.at(i)<<"\n";
     }
#endif

    unsigned int cont=0;
    data_size= (step_index==std::numeric_limits< unsigned int>::max())?headers.size():(headers.size()-1);
    TFLOAT * reads=new TFLOAT[data_size];
    while (log.good()) {

        if(if_only_numbers(tmp)&&tmp.length()>1){

            unsigned int TS=cont;
            std::stringstream sstr(tmp);

            if (cont++ % skip==0)
            for(unsigned int i=0;i<headers.size();i++){
                if (i<step_index){
                    sstr >> reads[i];
                } else if (i>step_index) {
                    sstr >> reads[i-1];
                } else {
                    sstr >> TS;
                }
            }
            if (sstr.fail()){
                std::cerr << "Ignoro la linea "<<lcont<<": '"<<tmp<<"'\n";
                sstr.clear();
                cont--;
            } else {
                for (unsigned int i=0;i<data_size;i++)
                    data.push_back(reads[i]);
                timesteps.push_back(TS);
            }
        }

        std::getline(log,tmp);
        lcont++;
    }
    delete [] reads;
#ifdef DEBUG
    std::cerr << "Numero di dati letti dal file '"<<filename<<"': "<<data.size()<<"\n";
#endif

}

//questo analizza la stringa speciale "#traj:JZ N q1 ... qN" e ritorna le cariche
template <class TFLOAT> std::vector<double> ReadLog<TFLOAT>::qs(std::string header) {

}

template <class TFLOAT> bool ReadLog<TFLOAT>::need_binary(std::vector<std::string> headers) {
    for (auto it=headers.begin();it!=headers.end();++it) {
        //usa la funzione qs per verificare se serve la traiettoria NOT IMPLEMENTED
    }

    return false;

}

template <class TFLOAT> void ReadLog<TFLOAT>::set_traj(Traiettoria * t){
    traiettoria=t;
}

template <class TFLOAT> std::pair<unsigned int, bool> ReadLog<TFLOAT>::get_index_of(std::string header) {
    //qui devo controllare la sintassi. Se "#traj:JZ N q1 ... qN", allora calcolo le correnti
    //partendo dalla traiettoria binaria NOT IMPLEMENTED

    unsigned int idx=0;
    for (unsigned int i=0;i<headers.size();i++){
        if (headers.at(i)!="Step"){
            if (headers.at(i)==header){
                return std::pair<unsigned int,bool>(idx,true);
            } else {
                idx++;
            }
        }
    }
    return std::pair<unsigned int ,bool>(idx,false);
}

template <class TFLOAT> unsigned int ReadLog<TFLOAT>::timestep(unsigned int index){
    return timesteps[index];
}

template <class TFLOAT> bool ReadLog<TFLOAT>::if_only_numbers(std::string str){
    return str.find_first_not_of("Ee0123456789 \t.-+")==std::string::npos;
}

template <class TFLOAT> ReadLog<TFLOAT>::~ReadLog(){


}

template <class TFLOAT> TFLOAT * ReadLog<TFLOAT>::line(unsigned int index){
    return &data[index*data_size];
}

template class ReadLog<double>;
template class ReadLog<long double>;
