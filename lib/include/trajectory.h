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

#include "basetrajectory.h"


#include "lammps_struct.h"

/**
 * Questa classe garantisce un accesso veloce ai timestep caricati con set_access_at
 * , che carica tutti i timestepa partendo da quello specificato in un numero pari a quello richiesto
 * con la funzione set_data_access_block_size in precedenza. Caricare nuovi pezzi costa nuove
 * allocazioni (e deallocazioni) di memoria e tempo cpu. Il file viene letto tramite la chiamata di sistema mmap.
 * Nel programma vengono allocati con new [] solo degli array dove vengono immagazzinati i dati
 * della finestra (quindi consumando meno memoria che nel caricamento del file in un unico colpo)
**/

class Trajectory : public BaseTrajectory<Trajectory>
{
public:
    Trajectory(std::string filename);
    ~Trajectory();

    template <bool SAFE=true, bool ATOM=false>
    double * get_array(double * base, const size_t & timestep, const size_t & atomo, const size_t & stride1, const size_t &stride2) {
        if constexpr (SAFE) {
            if constexpr (ATOM) {
                if (atomo >= natoms) {
                    std::stringstream ss;
                    ss << "Requested atom index ("<< atomo << ") is not in the range [0, "<<natoms-1<<"]\n";
                    throw std::runtime_error(ss.str());
                }
            }

            size_t t=timestep-current_timestep; // note that this can be garbage if timestep<current_timestep

            if (dati_caricati && t < loaded_timesteps && timestep>=current_timestep) { // vuol dire che ho già caricato i dati
                return base+t*stride1+atomo*stride2;
            } else { // non ho caricato i dati, li carico prima (questo potrebbe essere inefficiente se dopo devo satare di nuovo indietro!
#ifdef DEBUG
                std::cerr << "Warning: loading timesteps starting from "<<timestep<<" not requested with the loading routines!\n";
#endif
                if(set_access_at(timestep)) {
                    t=timestep-current_timestep;
                    if (timestep>=current_timestep && t< loaded_timesteps){
                        return base + t*stride1+atomo*stride2;
                    }else{
                        throw std::runtime_error("requested timestep is out of range");
                    }
                } else {
                    throw  std::runtime_error("Error loading the file\n");
                }
            }
        } else {
            return base + (timestep-current_timestep)*stride1+atomo*stride2;
        }
    }

    template<bool SAFE=true>
    double * positions (const size_t &timestep, const size_t &atomo){
        return get_array<SAFE,true>(buffer_positions,timestep,atomo,3*natoms,3);
    }
    template<bool SAFE=true>
    double * velocity(const size_t &timestep, const size_t &atomo) {
        return get_array<SAFE,true>(buffer_velocity,timestep,atomo,3*natoms,3);
    }
    template<bool SAFE=true>
    double * box(const size_t &timestep) {
        return get_array<SAFE,false>(buffer_boxes,timestep,0,buffer_boxes_stride,0);
    }
    template<bool SAFE=true>
    double * positions_cm(const size_t &timestep, const size_t &tipo){
        return get_array(buffer_positions_cm,timestep,tipo,3*ntypes,3);
    }

    template<bool SAFE=true>
    double * velocity_cm(const size_t &timestep, const size_t &tipo){
        return get_array(buffer_velocity_cm,timestep,tipo,3*ntypes,3);
    }
    double *box_last(){
        if (loaded_timesteps>0) {
            return &buffer_boxes[ buffer_boxes_stride*(loaded_timesteps-1)];
        } else {
            throw std::runtime_error("No data is loaded!\n");
        }
    }

    using BaseTrajectory<Trajectory>::Errori;
    Trajectory::Errori set_data_access_block_size(const size_t & timesteps);
    Trajectory::Errori set_access_at(const size_t & timesteps);
    int64_t get_timestep_lammps(size_t timestep);
    void index_all();
    int * get_lammps_id();
    int *get_lammps_type();
private:
    std::map<int,unsigned int>id_map;
    size_t * timesteps; // puntatori (offset rispetto all'inizio) all'inizio di ogni timesteps
    int64_t * timesteps_lammps; // timesteps secondo lammps
    void allunga_timesteps(size_t nuova_dimensione);
    template <bool set_chunk=true>
    size_t leggi_pezzo(const size_t & partenza, TimestepManager &timestep, Chunk * &chunk);
    size_t leggi_pezzo_intestazione(const size_t & partenza, TimestepManager &timestep);
    void init_buffer_tipi();
    int fd;
    size_t fsize,tstep_size;
    char * file;
    size_t timestep_indicizzato;
    bool ok,dati_caricati,indexed_all;
    bool calculate_center_of_mass;
    long pagesize;
    size_t allinea_offset(const long &offset, size_t & differenza);

};

#endif // TRAIETTORIA_H
