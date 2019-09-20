/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "readlog.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <cerrno>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <thread>
#include "cronometro.h"

#include "chargefluxts.h"

template <class TFLOAT> ReadLog<TFLOAT>::ReadLog(std::string filename, Traiettoria *t, unsigned int skip, unsigned int nthreads, unsigned int nbatch, std::vector<std::string> req_headers):
    traiettoria(t),skip(skip),nthreads(nthreads),nbatch(nbatch)
{
    cronometro cron;
    cron.start();

    if (nthreads<0)
        nthreads=1;
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

    if (header.size()==0) {
        std::cerr << "Errore: non riesco a trovare l'intestazione delle colonne!\n";
        abort();
    }

    //trova i nomi delle colonne e la colonna del timestep
    std::string delim = " \t";

     size_t start = 0;
     size_t end = header.find_first_of(delim);
     while (end != std::string::npos)
     {
         if (end!=start&& start<header.length()) {
             headers.push_back(header.substr(start, end - start));
             if (headers.back().find_first_of(delim)!=std::string::npos) headers.pop_back();
         }
         start = end + 1;
         end = header.find_first_of(delim, start);
         if (headers.size()>0 && headers.back()=="Step") step_index=headers.size()-1;
     }
     if (end!=start && start<header.length()) headers.push_back(header.substr(start, end - start));
     if (headers.size()<=0) {
         std::cerr << "Errore: non ho trovato nessun header nel file '"<<filename<< "'\n";
         abort();
     }
     if (headers.back()=="Step") step_index=headers.size()-1;
#ifdef DEBUG
     std::cerr << "Indice degli header:\n";
     for(int i=0;i<headers.size();i++){
         std::cerr << i << " " << headers.at(i)<<"\n";
     }
#endif

     //determina il numero di colonne da allocare, e comprende le eventuali cariche
     unsigned int n_columns_from_binary=0;
     data_size_from_binary=0;
     if (req_headers.size()>0){
         for (auto it = req_headers.begin();it!= req_headers.end();++it) {
             if (*it == "Step"){
                 std::cerr << "Errore: non puoi chiedere di fare un analisi con al posto della corrente i timestep (colonna Step)\n";
                 abort();
             }
             auto Qs=qs(*it);
             if (Qs.first !=""){
                 n_columns_from_binary++; //aggiungi i valori delle cariche letti
                 q_current_type.push_back(Qs);
             }
         }
         data_size_from_binary=n_columns_from_binary*3;
     }

    unsigned int cont=0;
    data_size= (step_index==std::numeric_limits< unsigned int>::max())?(headers.size()+data_size_from_binary):(headers.size()-1+data_size_from_binary);



    TFLOAT * reads=new TFLOAT[data_size*nbatch*nthreads];
    unsigned int *TS=new unsigned int [nbatch*nthreads];
    bool * readok=new bool[nbatch*nthreads];
    for (unsigned int i=0;i<data_size*nbatch*nthreads;++i)
        reads[i]=0.0;
    for (unsigned int ib=0;ib<nbatch*nthreads;ib++){
        readok[ib]=false;
    }

    std::string * buffer=new std::string [nbatch*nthreads];
    buffer[0]=tmp;
    unsigned int ib=0;
    for (ib=1;ib<nbatch*nthreads && log.good();ib++){
        std::getline(log,buffer[ib]);
        lcont++;
    }


    while (ib>0) {


        std::vector<std::thread> threads;

        for (unsigned int ith=0;ith<nthreads;ith++){
            threads.push_back(std::thread([&,ith](){

                for (unsigned int il=0;il<nbatch;il++){
                    if (ith*nbatch+il>=ib){
                        continue;
                    }
                    std::string &tmp=buffer[ith*nbatch+il];

                    if(if_only_numbers(tmp)&&tmp.length()>1){
                        TS[nbatch*ith+il]=0;
                        std::stringstream sstr(tmp);

                            for(unsigned int i=0;i<headers.size();i++){
                                if (i<step_index){
                                    sstr >> reads[data_size*nbatch*ith+il*data_size+i];
                                } else if (i>step_index) {
                                    sstr >> reads[data_size*nbatch*ith+il*data_size+i-1];
                                } else {
                                    sstr >> TS[nbatch*ith+il];
                                }
                            }
                        if (sstr.fail()){
                            readok[nbatch*ith+il]=false;
                        } else {
                            readok[nbatch*ith+il]=true;
                        }
                    }
                }

            }));
        }


        for (unsigned int ith=0;ith<nthreads;ith++)
            threads[ith].join();
        threads.clear();

        //copia i dati nell'array grande

        for (unsigned int i=0;i<nthreads*nbatch;i++) {
            if (readok[i] ) {
                readok[i]=false;
                if (cont++ % skip == 0){
                    data.insert(data.end(),&reads[i*data_size],&reads[(i+1)*data_size]);
                    if (TS[i]>0)
                        timesteps.push_back(TS[i]);
                    else
                        timesteps.push_back(cont);
                }
            }
        }


        //carica un altro gruppo di linee
        ib=0;
        for (ib=0;ib<nbatch*nthreads && log.good();ib++){
            std::getline(log,buffer[ib]);
            lcont++;
        }

    }
    delete [] reads;
    delete [] readok;
    delete [] TS;
    delete [] buffer;
    cron.stop();
    std::cerr << "Tempo per la lettura dei dati: " << cron.time() << "s\n";
    std::cerr << "Numero di dati letti dal file '"<<filename<<"': "<<data.size()/data_size<<"\n";
    if (data.size()==0) {
        std::cerr << "Attenzione: non sono riuscito a leggere alcun dato. (l'intestazione ha lo stesso numero di colonne dei dati?)\n";
    }

}

template <class TFLOAT> unsigned int ReadLog<TFLOAT>::get_calc_j_index(std::string header) {

    for (unsigned int i=0;i<q_current_type.size();i++) {
        if (q_current_type[i].first==header) {
            return i;
        }
    }
    std::cerr << "Errore: impossibile trovare l'header "<< header<<" nella lista '";
    for (unsigned int i=0;i<q_current_type.size();i++) std::cerr << " '" << q_current_type[i].first<<"'";
    std::cerr << "\n";
    abort();

    return q_current_type.size();
}
/*
template <class TFLOAT> void ReadLog<TFLOAT>::write_currents(std::string filename) {
    if (data_size_from_binary<=0) {
        std::cerr << "Attenzione: non posso scrivere le correnti calcolate quando non ci sono!\n";
        return;
    }

    std::ofstream out(filename);
    for (unsigned int i=0;i<timesteps.size();i++) {
        out << timesteps[i];
        for (unsigned int j=0;j<data_size_from_binary;j++)
            out <<" "<< data[i*data_size+data_size-data_size_from_binary+j];
        out << "\n";
    }
}
*/

//questo analizza la stringa speciale "#traj:JZ N q1 ... qN" e ritorna le cariche
template <class TFLOAT> std::pair<std::string,std::vector<TFLOAT> > ReadLog<TFLOAT>::qs(std::string header) {
    if (header.size()==0 || header[0] != '#')
        return std::pair<std::string,std::vector<TFLOAT> >();
    std::stringstream iss(header);
    std::vector<std::string> t{std::istream_iterator<std::string>{iss},std::istream_iterator<std::string>{}};
    if (t.size()<3){
        std::cerr <<"Errore nel formato dell'header '"<<header<<"'. Deve essere nella forma '#traj:JZ N q1 ... qN'.\n";
        abort();
    }
    int n_charges=std::strtol (t.at(1).c_str(),NULL,0);
    if (n_charges<=0 || errno == ERANGE || t.size() != n_charges + 2) {
        std::cerr << "Errore nel formato dell'header '"<<header<<"'. Deve essere nella forma '#traj:JZ N q1 ... qN'.\n";
        abort();
    }
    std::vector <TFLOAT> c(n_charges);
    for (unsigned int i=2;i<n_charges+2;i++){
        c[i-2]=std::strtold(t.at(i).c_str(),NULL);
        if (errno == ERANGE){
            std::cerr << "Errore nel formato dell'header '"<<header<<"'. Deve essere nella forma '#traj:JZ N q1 ... qN'.\n";
            abort();
        }

    }

    return std::pair<std::string,std::vector<TFLOAT> > (header,c);

}

template <class TFLOAT> int ReadLog<TFLOAT>::need_binary(std::vector<std::string> headers) {
    int res=0;
    for (auto it=headers.begin();it!=headers.end();++it) {
        if (qs(*it).first!="")
            res++;
    }

    return res;

}

template <class TFLOAT> void ReadLog<TFLOAT>::calc_currents(Traiettoria * t,unsigned int n_b){
    traiettoria=t;
    //calcola e legge la corrente partendo dal file binario (nella classe traiettoria sono già presenti le velocità dei centri di massa)
    unsigned int timesteps_tot=data.size()/data_size;
    unsigned int n_data_b=timesteps_tot/n_b;
    unsigned int ntypes=t->get_ntypes();
    //controlla che il numero di correnti corrisponda al numero di coefficienti forniti
    for (auto it = q_current_type.begin();it!=q_current_type.end();++it) {
        if ((*it).second.size()!=ntypes) {
            std::cerr << "Errore: il numero di coefficienti specificati deve corrispondere al numero di tipi di atomi nella traiettoria!\n";
            abort();
        }
    }
    traiettoria->imposta_dimensione_finestra_accesso(n_data_b);
    for (unsigned int ib=0;ib<n_b;ib++){
        unsigned int ultimo=(ib+1)*n_data_b;
        if (ib==n_b-1){
            ultimo=timesteps_tot;
            traiettoria->imposta_dimensione_finestra_accesso(ultimo-ib*n_data_b);
        }
        traiettoria->imposta_inizio_accesso(n_data_b*ib);
        for (unsigned int ts=ib*n_data_b;ts<ultimo;ts++){
            //calcola le varie correnti utilizzando i dati presenti negli header, e copia nello spazio lasciato libero durante la lettura. Poi sono a posto e il resto del codice non cambia
            for (unsigned int i=0;i<q_current_type.size();i++) {
                double * v_cm=t->velocita_cm(ts,0);
                for (unsigned int icoord=0;icoord<3;icoord++)
                    data[ts*data_size+(data_size-data_size_from_binary)+i*3+icoord]=0.0;
                for (unsigned int icm=0;icm<ntypes;icm++) {
                    for (unsigned int icoord=0;icoord<3;icoord++) {
                        data[ts*data_size+(data_size-data_size_from_binary)+i*3+icoord]+=
                                v_cm[icm*3+icoord]*q_current_type[i].second[icm];
                    }
                }
            }
        }
    }
}

template <class TFLOAT> std::pair<unsigned int, bool> ReadLog<TFLOAT>::get_index_of(std::string header) {
    //qui devo controllare la sintassi. Se "#traj:JZ N q1 ... qN", allora calcolo le correnti
    //partendo dalla traiettoria binaria

    unsigned int idx=0;
    if (header.size()==0) {
        std::cerr << "Errore: header di lunghezza nulla! (" __FILE__ <<":" <<__LINE__ <<")\n";
    }
    //TODO: questo va cambiato: l'header delle correnti calcolate è uno solo!
    if (header[0]!='#'){
        for (unsigned int i=0;i<headers.size();i++){
            if (headers.at(i)!="Step"){
                if (headers.at(i)==header){
                    return std::pair<unsigned int,bool>(idx,true);
                } else {
                    idx++;
                }
            }
        }
    }else { // devo andare a vedere la traiettoria binaria

        idx=get_calc_j_index(header);
        return std::pair<unsigned int, bool>(data_size-data_size_from_binary+idx*3,true);

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
