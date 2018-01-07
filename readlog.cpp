#include "readlog.h"
#include <fstream>
#include <sstream>

ReadLog::ReadLog(std::string filename, Traiettoria *t, unsigned int skip):
    traiettoria(t),skip(skip),step_index(0)
{
    unsigned int lcont=0;
    std::string tmp,header;
    std::ifstream log(filename);

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
    while (log.good()) {

        if(if_only_numbers(tmp)&&tmp.length()>1){

            double * reads=new double[headers.size()];
            unsigned int TS=0;
            std::stringstream sstr(tmp);

            if (cont++ % skip==0)
            for(int i=0;i<headers.size();i++){
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
                delete [] reads;
            } else {
                data.push_back(std::pair<unsigned int ,double*>(TS,reads));
            }
        }

        std::getline(log,tmp);
        lcont++;
    }
#ifdef DEBUG
    std::cerr << "Numero di dati letti dal file '"<<filename<<"': "<<data.size()<<"\n";
#endif

}


std::pair<unsigned int, bool> ReadLog::get_index_of(std::string header) {
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

unsigned int ReadLog::timestep(unsigned int index){
    return data.at(index).first;
}

bool ReadLog::if_only_numbers(std::string str){
    return str.find_first_not_of("e0123456789    .-+")==std::string::npos;
}

ReadLog::~ReadLog(){

    for (unsigned int i=0;i<data.size();i++){
        delete [] data.at(i).second;
    }

}

double * ReadLog::line(unsigned int index){
    return data.at(index).second;
}
