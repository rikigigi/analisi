#include "readlog.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <cerrno>
#include <cstdlib>
#include <iterator>
#include <vector>
#include <thread>

#include "chargefluxts.h"

template <class TFLOAT> ReadLog<TFLOAT>::ReadLog(std::string filename, Traiettoria *t, unsigned int skip, unsigned int nthreads, std::vector<std::string> req_headers):
    traiettoria(t),skip(skip),nthreads(nthreads)
{
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
                 q_current_type.push_back(Qs.second);
             }
         }
         data_size_from_binary=n_columns_from_binary*3;
     }

    unsigned int cont=0;
    data_size= (step_index==std::numeric_limits< unsigned int>::max())?(headers.size()+data_size_from_binary):(headers.size()-1+data_size_from_binary);



    const unsigned int nbatch=50;

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
#ifdef DEBUG
    std::cerr << "Numero di dati letti dal file '"<<filename<<"': "<<data.size()/data_size<<"\n";
#endif

}

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

    return std::pair<std::string,std::vector<TFLOAT> > (t.at(0),c);

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
        if ((*it).size()!=ntypes) {
            std::cerr << "Errore: il numero di coefficienti specificati deve corrispondere al numero di tipi di atomi nella traiettoria!\n";
            abort();
        }
    }
    for (unsigned int ib=0;ib<n_b;ib++){
        for (unsigned int ts=ib*n_data_b;ts<(ib+1)*n_data_b;ts++){
            //calcola le varie correnti utilizzando i dati presenti negli header, e copia nello spazio lasciato libero durante la lettura. Poi sono a posto e il resto del codice non cambia
            for (unsigned int i=0;i<q_current_type.size();i++) {
                double * v_cm=t->velocita_cm(ts,0);
                for (unsigned int icoord=0;icoord<3;icoord++)
                    data[ts*data_size+(data_size-data_size_from_binary)+i*3+icoord]=0.0;
                for (unsigned int icm=0;icm<ntypes;icm++) {
                    for (unsigned int icoord=0;icoord<3;icoord++) {
                        data[ts*data_size+(data_size-data_size_from_binary)+i*3+icoord]+=
                                v_cm[icm*3+icoord]*q_current_type[i][icm];
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
        std::cerr << "Errore: header di lunghezza nulla! ("__FILE__<<":" <<__LINE__ <<")\n";
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
