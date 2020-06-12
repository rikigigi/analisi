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


struct Intestazione_timestep;
struct Chunk;


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
    double * posizioni (const int & timestep, const int & atomo);
    double * velocita (const int & timestep, const int & atomo);
    double * scatola (const int & timestep);
    double * posizioni_cm(const int & timestep, const int & tipo);
    double * velocita_cm(const int & timestep, const int & tipo);
    double *scatola_last();

    using TraiettoriaBase<Traiettoria>::Errori;
    Traiettoria::Errori imposta_dimensione_finestra_accesso(const int & timesteps);
    Traiettoria::Errori imposta_inizio_accesso(const int & timesteps);
    int64_t get_timestep_lammps(unsigned int timestep);
    void index_all();
//    void set_calculate_center_of_mass(bool);
//    bool get_calculate_center_of_mass();
private:
    std::map<int,unsigned int>id_map;
    size_t * timesteps; // puntatori (offset rispetto all'inizio) all'inizio di ogni timesteps
    int64_t * timesteps_lammps; // timesteps secondo lammps
    void allunga_timesteps(unsigned int nuova_dimensione);
    size_t leggi_pezzo(const size_t & partenza,Intestazione_timestep * &timestep,Chunk * &chunk);
    size_t leggi_pezzo_intestazione(const size_t & partenza, Intestazione_timestep * &timestep);
    void init_buffer_tipi();
    int fd;
    size_t fsize,tstep_size;
    char * file;
    int timestep_corrente,timestep_finestra,timestep_indicizzato;
    bool ok,dati_caricati,indexed_all;
    bool triclinic;
    bool calculate_center_of_mass;
    long pagesize;
    size_t allinea_offset(const size_t & offset, size_t & differenza);

};

#endif // TRAIETTORIA_H
