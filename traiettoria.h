/**
  *
  * (c) Riccardo Bertossa, 2017
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy to receive a copy
  *   of the good modified code, with comments, at
  *    riccardo dot bertossa at gmail dot com
  *
**/



#ifndef TRAIETTORIA_H
#define TRAIETTORIA_H

#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <map>


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

class Traiettoria
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
    double * posizioni_inizio(){return buffer_posizioni;}
    double * velocita_inizio(){return buffer_velocita;}
    int get_ntypes();
    std::vector<unsigned int> get_types();
    unsigned int get_type(const unsigned int &atomo);
    int get_type_min() {return min_type;}
    int get_type_max() {return max_type;}
    enum Errori {non_inizializzato=0,oltre_fine_file=2,Ok=1};
    Traiettoria::Errori imposta_dimensione_finestra_accesso(const int & timesteps);
    Traiettoria::Errori imposta_inizio_accesso(const int & timesteps);
    int get_natoms(){return natoms;}
    int get_ntimesteps(){return n_timesteps;}
    int64_t get_timestep_lammps(unsigned int timestep);
    double get_mass(unsigned int i) {if (i<get_ntypes()) return masse[i]; std::cerr<< "Errore: non posso ritornare una massa per un tipo che non esiste!\n";abort(); return 0.0;}
    void set_mass(unsigned int i,double m) {if (i<get_ntypes()) masse[i]=m;}
    void set_charge(unsigned int i, double c){if (i<get_ntypes()) cariche[i]=c;}
    double get_charge(unsigned int  i){if (i<get_ntypes()) return cariche[i]; std::cerr<< "Errore: non posso ritornare una carica per un tipo che non esiste!\n";abort(); return 0.0;}
    void index_all();
    void set_pbc_wrap(bool);
    double d2_minImage(unsigned int i,unsigned int j, unsigned int itimestep,double *l);
    double d2_minImage(unsigned int i,unsigned int j, unsigned int itimestep,unsigned int jtimestep,double *l);
//    void set_calculate_center_of_mass(bool);
//    bool get_calculate_center_of_mass();
private:
    std::vector<unsigned int> types;
    std::map<unsigned int,unsigned int>type_map;
    double * buffer_posizioni; //velocita' e posizioni copiate dal file caricato con mmap, in ordine (nela traiettoria di LAMMPS sono disordinate)
    double * buffer_velocita;
    double * masse;
    double * cariche;
    double * buffer_scatola; //dimensioni della simulazione ad ogni timestep
    double * buffer_posizioni_cm; // posizioni del centro di massa
    double * buffer_velocita_cm; // velocità del centro di massa

    int * buffer_tipi,*buffer_tipi_id;
    size_t * timesteps; // puntatori (offset rispetto all'inizio) all'inizio di ogni timesteps
    int64_t * timesteps_lammps; // timesteps secondo lammps
    void allunga_timesteps(unsigned int nuova_dimensione);
    int n_timesteps;
    size_t leggi_pezzo(const size_t & partenza,Intestazione_timestep * &timestep,Chunk * &chunk);
    size_t leggi_pezzo_intestazione(const size_t & partenza, Intestazione_timestep * &timestep);
    void init_buffer_tipi();
    int fd;
    size_t fsize,tstep_size;
    char * file;
    int timestep_corrente,timestep_finestra,timestep_indicizzato;
    bool ok,dati_caricati,indexed_all;
    bool triclinic,wrap_pbc;
    bool calculate_center_of_mass;
    int natoms,ntypes,min_type,max_type;
    long pagesize;
    size_t allinea_offset(const size_t & offset, size_t & differenza);

};

#endif // TRAIETTORIA_H
