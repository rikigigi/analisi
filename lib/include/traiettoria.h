/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#ifndef TRAIETTORIA_H
#define TRAIETTORIA_H

#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>

#include "traiettoriabase.h"


#include "lammps_struct.h"

/**
 * Questa classe garantisce un accesso veloce ai timestep caricati con imposta_inizio_accesso
 * , che carica tutti i timestepa partendo da quello specificato in un numero pari a quello richiesto
 * con la funzione imposta_dimensione_finestra_accesso in precedenza. Caricare nuovi pezzi costa nuove
 * allocazioni (e deallocazioni) di memoria e tempo cpu. Il file viene letto tramite la chiamata di sistema mmap.
 * Nel programma vengono allocati con new [] solo degli array dove vengono immagazzinati i dati
 * della finestra (quindi consumando meno memoria che nel caricamento del file in un unico colpo)
**/

class Traiettoria : public TraiettoriaBase<Traiettoria>
{
public:
    Traiettoria(std::string filename);
    ~Traiettoria();

    template <bool SAFE=true, bool ATOM=false>
    double * get_array(double * base, const int & timestep, const int & atomo, const int & stride1, const int &stride2) {
        if constexpr (SAFE) {
            if constexpr (ATOM) {
            if (atomo<0 || atomo >= natoms) {
                std::stringstream ss;
                ss << "Requested atom index ("<< atomo << ") is not in the range 0-"<<natoms<<"\n";
                throw std::runtime_error(ss.str());
            }
            }

            int t=timestep-timestep_corrente;

            if (dati_caricati && t < loaded_timesteps && t>=0) { // vuol dire che ho già caricato i dati

                return base+t*stride1+atomo*stride2;

            } else { // non ho caricato i dati, li carico prima (questo potrebbe essere inefficiente se dopo devo satare di nuovo indietro!
#ifdef DEBUG
                std::cerr << "Warning: loading timesteps not requested with the loading routines!\n";
#endif
                if(imposta_inizio_accesso(timestep)) {
                    t=timestep-timestep_corrente;
                    if (t>0 && t< loaded_timesteps)
                        return base + t*stride1+atomo*stride2;
                    else
                        throw std::runtime_error("requested timestep is out of range");
                } else {
                    throw  std::runtime_error("Error loading the file\n");
                }
            }
        } else {
            return base + (timestep-timestep_corrente)*stride1+atomo*stride2;
        }
    }

    template<bool SAFE=true>
    double * posizioni (const int &timestep, const int &atomo){
        return get_array<SAFE,true>(buffer_posizioni,timestep,atomo,3*natoms,3);
    }
    template<bool SAFE=true>
    double * velocita(const int &timestep, const int &atomo) {
        return get_array<SAFE,true>(buffer_velocita,timestep,atomo,3*natoms,3);
    }
    template<bool SAFE=true>
    double * scatola(const int &timestep) {
        return get_array<SAFE,false>(buffer_scatola,timestep,0,6,0);
    }
    template<bool SAFE=true>
    double * posizioni_cm(const int &timestep, const int &tipo){
        return get_array(buffer_posizioni_cm,timestep,tipo,3*ntypes,3);
    }

    template<bool SAFE=true>
    double * velocita_cm(const int &timestep, const int &tipo){
        return get_array(buffer_velocita_cm,timestep,tipo,3*ntypes,3);
    }
    double *scatola_last(){
        if (loaded_timesteps>0) {
            return &buffer_scatola[(unsigned int) 6*(loaded_timesteps-1)];
        } else {
            throw std::runtime_error("No data is loaded!\n");
        }
    }

    using TraiettoriaBase<Traiettoria>::Errori;
    Traiettoria::Errori imposta_dimensione_finestra_accesso(const int & timesteps);
    Traiettoria::Errori imposta_inizio_accesso(const int & timesteps);
    int64_t get_timestep_lammps(int timestep);
    void index_all();
    //    void set_calculate_center_of_mass(bool);
    //    bool get_calculate_center_of_mass();
private:
    std::map<int,unsigned int>id_map;
    size_t * timesteps; // puntatori (offset rispetto all'inizio) all'inizio di ogni timesteps
    int64_t * timesteps_lammps; // timesteps secondo lammps
    void allunga_timesteps(int nuova_dimensione);
    size_t leggi_pezzo(const size_t & partenza, TimestepManager &timestep, Chunk * &chunk);
    size_t leggi_pezzo_intestazione(const size_t & partenza, TimestepManager &timestep);
    void init_buffer_tipi();
    int fd;
    size_t fsize,tstep_size;
    char * file;
    int timestep_corrente,timestep_indicizzato;
    bool ok,dati_caricati,indexed_all;
    bool triclinic;
    bool calculate_center_of_mass;
    long pagesize;
    size_t allinea_offset(const long &offset, size_t & differenza);

};

#endif // TRAIETTORIA_H
