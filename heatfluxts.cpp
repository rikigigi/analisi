#include "heatfluxts.h"
#include <fstream>

//TODO: anche questo dovrebbe leggere il file a blocchi!


HeatFluxTs::HeatFluxTs(std::string filename, Traiettoria *t,unsigned int skip)
    : traiettoria(t),heatflux(0),skip(skip)
{

    heatflux = new double [traiettoria->get_ntimesteps()*3/skip];
    T = new double [traiettoria->get_ntimesteps()];
    size=traiettoria->get_ntimesteps();

    std::ifstream log(filename);

    //trova l'inizio dei dati
    std::string head_heat =
    "Step Time PotEng TotEng Lx Press Temp c_flusso[1] c_flusso[2] c_flusso[3] ";
    std::string tmp;

    while (log.good()) {
        std::getline(log,tmp);
        if (tmp==head_heat)
            break;
    }

    if (!log.good()) {
        std::cerr << "Errore: impossibile trovare nel file di log \""<<filename<<"\" l'intestazione \""<<head_heat<<"\"\n";
        abort();
    }

    unsigned int cur_ts=0;

    while (log.good() && cur_ts*skip < traiettoria->get_ntimesteps()) {
        int64_t step;
        double time=.0, poteng=.0,toteng=.0,lx=.0,press=.0,temp=.0, j[3]={0.0,0.0,0.0};
        log >> step >> time >> poteng >> toteng >> lx >> press >> temp >>
                j[0] >> j[1] >> j[2];

        if (traiettoria->get_timestep_lammps(cur_ts*skip)>step){
            continue;
        } else if (traiettoria->get_timestep_lammps(cur_ts*skip)<step) {
            std::cerr << "Errore: i frame dei file non corrispondono! (forse il file di log ha una risoluzione temporale minore della traiettoria?)\n";
            abort();
            break;
        }
        //se sono qui vuol dire che i frame corrispondono, deco leggere il dato, che è quello giusto.

        for (unsigned int i=0;i<3;i++)
            heatflux[cur_ts*3+i]=j[i];
        T[cur_ts]=temp;
        L=lx;
        cur_ts++;
    }
}

double * HeatFluxTs::flux(unsigned int ts) {
    if (ts%skip !=0){
        std::cerr << "Errore: richiesto un timestep che non è stato letto (salto ogni "<<skip<<" timesteps)\n";
        abort();
    }

    if (ts<size)
        return &heatflux[ts*3/skip];
    else {
        std::cerr << "Errore: richiesto un flusso di calore fuori dai limiti caricati!\n";
        abort();
        return 0;
    }
}

double * HeatFluxTs::temp(unsigned int ts) {
    if (ts%skip !=0){
        std::cerr << "Errore: richiesto un timestep che non è stato letto (salto ogni "<<skip<<" timesteps)\n";
        abort();
    }

    if (ts<size)
        return &T[ts/skip];
    else {
        std::cerr << "Errore: richiesta una temperatura fuori dai limiti caricati!\n";
        abort();
        return 0;
    }
}
HeatFluxTs::~HeatFluxTs() {
    delete [] heatflux;
    delete [] T;
}
