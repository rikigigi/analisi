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
    double * posizioni (const int &timestep, const int &atomo){

        if (atomo<0 || atomo > natoms) {
            std::cerr << "Atomo richiesto ("<< atomo << ") non è compreso nel range di atomi 0-"<<natoms<<"\n";
            return 0;
        }

        int t=timestep-timestep_corrente;

        if (dati_caricati && t < loaded_timesteps && t>=0) { // vuol dire che ho già caricato i dati

            return &buffer_posizioni[t*3*natoms+atomo*3];

        } else { // non ho caricato i dati, li carico prima (questo potrebbe essere inefficiente se dopo devo satare di nuovo indietro!
    #ifdef DEBUG
            std::cerr << "Attenzione: sto caricando dei timestep non richiesti in precedenza!\n";
    #endif
            if(imposta_inizio_accesso(timestep)) {
                t=timestep-timestep_corrente;
                if (t>0 && t< loaded_timesteps)
                    return &buffer_posizioni[t*3*natoms+atomo*3];
                else
                    abort();
            } else {
                std::cerr << "Errore nel caricamento del file.\n";
                return 0;
            }
        }

    }
    double * velocita(const int &timestep, const int &atomo) {
        if (atomo<0 || atomo > natoms) {
            std::cerr << "Atomo richiesto ("<< atomo << ") non è compreso nel range di atomi 0-"<<natoms<<"\n";
            return 0;
        }


        int t=timestep-timestep_corrente;

        if (dati_caricati && t < loaded_timesteps && t>=0) { // vuol dire che ho già caricato i dati

            return &buffer_velocita[t*3*natoms+atomo*3];

        } else { // non ho caricato i dati, li carico prima (questo potrebbe essere inefficiente se dopo devo satare di nuovo indietro!

    #ifdef DEBUG
            std::cerr << "Attenzione: sto caricando dei timestep non richiesti in precedenza!\n";
    #endif
            if(imposta_inizio_accesso(timestep)){
                t=timestep-timestep_corrente;
                if (t>0 && t< loaded_timesteps)
                    return &buffer_velocita[t*3*natoms+atomo*3];
                else
                    abort();
            } else {
                std::cerr << "Errore nel caricamento del file.\n";
                abort();
                return 0;
            }
        }

    }
    double * scatola(const int &timestep) {
        int t=timestep-timestep_corrente;

        if (dati_caricati && t < loaded_timesteps && t>=0) { // vuol dire che ho già caricato i dati

            return &buffer_scatola[t*6];

        } else { // non ho caricato i dati, li carico prima (questo potrebbe essere inefficiente se dopo devo satare di nuovo indietro!
    #ifdef DEBUG
            std::cerr << "Attenzione: sto caricando dei timestep non richiesti in precedenza!\n";
    #endif
            if(imposta_inizio_accesso(timestep)){
                t=timestep-timestep_corrente;
                if (t>0 && t< loaded_timesteps)
                    return &buffer_scatola[t*6];
                else
                    abort();
            } else {
                std::cerr << "Errore nel caricamento del file.\n";
                abort();
                return 0;
            }
        }
    }
    double * posizioni_cm(const int &timestep, const int &tipo){

        if (tipo<0 || tipo > ntypes) {
            std::cerr << "Tipo richiesto ("<< tipo << ") non è compreso nel range di atomi 0-"<<ntypes<<"\n";
            return 0;
        }

        int t=timestep-timestep_corrente;

        if (dati_caricati && t < loaded_timesteps && t>=0) { // vuol dire che ho già caricato i dati

            return &buffer_posizioni_cm[t*3*ntypes+tipo*3];

        } else { // non ho caricato i dati, li carico prima (questo potrebbe essere inefficiente se dopo devo satare di nuovo indietro!
    #ifdef DEBUG
            std::cerr << "Attenzione: sto caricando dei timestep non richiesti in precedenza!\n";
    #endif
            if(imposta_inizio_accesso(timestep)) {
                t=timestep-timestep_corrente;
                if (t>0 && t< loaded_timesteps)
                    return &buffer_posizioni_cm[t*3*ntypes+tipo*3];
                else
                    abort();
            } else {
                std::cerr << "Errore nel caricamento del file.\n";
                return 0;
            }
        }

    }

    double * velocita_cm(const int &timestep, const int &tipo){

        if (tipo<0 || tipo > ntypes) {
            std::cerr << "Tipo richiesto ("<< tipo << ") non è compreso nel range di atomi 0-"<<ntypes<<"\n";
            return 0;
        }

        int t=timestep-timestep_corrente;

        if (dati_caricati && t < loaded_timesteps && t>=0) { // vuol dire che ho già caricato i dati

            return &buffer_velocita_cm[t*3*ntypes+tipo*3];

        } else { // non ho caricato i dati, li carico prima (questo potrebbe essere inefficiente se dopo devo satare di nuovo indietro!
    #ifdef DEBUG
            std::cerr << "Attenzione: sto caricando dei timestep non richiesti in precedenza!\n";
    #endif
            if(imposta_inizio_accesso(timestep)) {
                t=timestep-timestep_corrente;
                if (t>0 && t< loaded_timesteps)
                    return &buffer_velocita_cm[t*3*ntypes+tipo*3];
                else
                    abort();
            } else {
                std::cerr << "Errore nel caricamento del file.\n";
                return 0;
            }
        }

    }
    double *scatola_last();

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
