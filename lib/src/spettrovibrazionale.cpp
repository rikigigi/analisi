/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "spettrovibrazionale.h"
#include <fstream>
#include <stdexcept>

template <class T>
SpettroVibrazionale<T>::SpettroVibrazionale(T * t, bool dump_):VectorOp<SpettroVibrazionale>()
{
    traiettoria=t;
    dump=dump_;
    size=0;
    trasformata=0;
    trasformata_size=0;
    fplan_natoms=0;
    fplan_size=0;
    tipi_atomi=0;
}

template <class T>
SpettroVibrazionale<T>::~SpettroVibrazionale(){

    fftw_free(trasformata);
}

template <class T>
 fftw_plan SpettroVibrazionale<T>::fftw3;
template <class T>
 unsigned int SpettroVibrazionale<T>::fplan_natoms;
template <class T>
 unsigned int SpettroVibrazionale<T>::fplan_size;

template <class T>
void SpettroVibrazionale<T>::deallocate_plan(){
    fftw_destroy_plan(fftw3);
    fplan_natoms=0;
    fplan_size=0;
}

template <class T>
unsigned int SpettroVibrazionale<T>::nExtraTimesteps(unsigned int n_b){
    return 0;
}

template <class T>
void SpettroVibrazionale<T>::reset(const unsigned int numeroTimestepsPerBlocco) {
//inizializzo la memoria per i moduli quadri e basta, se necessario!
// size è la lunghezza della trasformata. La vdata che contiene i moduli quadri sarà size/2+1 (trasformata reale).
    if (numeroTimestepsPerBlocco!=size) {
        size=numeroTimestepsPerBlocco;
        delete [] vdata;
        tipi_atomi=traiettoria->get_ntypes();
        if (tipi_atomi<=0) {
	    throw std::runtime_error("Number of types in the trajectory must be > 0");
            tipi_atomi=1;
        }
        data_length=(size/2+1)*3*tipi_atomi; // uno per ogni direzione dello spazio, per testare l'isotropia, e per tipo di atomo
        vdata = new double[data_length];

    }

}

//prima di chiamare questa la traiettoria deve essere impostata correttamente sulla finestra giusta (la funzione non sa qual'è lo timestep corrente.
template <class T>
void SpettroVibrazionale<T>::calculate(unsigned int primo  ///ignorato: prendo l'inizio di quello che c'è in memoria (attenzione a caricare bene i dati!)
                                  ) {
    //alloca se necessario con fftw_malloc (che alloca la memoria in modo che sia allineata correttaemente per poter sfruttare le istruzioni SIMD del processore
    // e prepara il piano della trasformata. Trasformata_size è la dimensione della trasformata. Deve essere uguale a size, la dimensione dell'array dei moduli quadri.
    // output della trasformata di un array reale: n/2+1 numeri complessi; ho 3*natomi trasformate da fare

    if (fplan_natoms!=traiettoria->get_natoms() || fplan_size != size || trasformata_size!=size) {
        fftw_free(trasformata);
        trasformata = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*(size/2+1)*3*traiettoria->get_natoms());
        trasformata_size=size;
        double * dummy= (double*) fftw_malloc(sizeof(double)*3*traiettoria->get_natoms()*size);
    fftw3 = fftw_plan_many_dft_r2c(1, // rango della trasformata (1D)
                                   (int*)&size, // lunghezza di ciascuna trasformata
                                   3*traiettoria->get_natoms(), // numero di trasformate
                                   /*  ****** input ******  */
                                   dummy,//traiettoria->velocita_inizio(), // puntatore ai dati
                                   NULL, // i dati sono tutti compatti, non fanno parte di array più grandi
                                   3*traiettoria->get_natoms(), // la distanza fra due dati successivi
                                   1, // la distanza fra due serie di dati adiacenti
                                   /*  ****** output ******  */
                                   trasformata, // puntatore all'array di output
                                   NULL,
                                   3*traiettoria->get_natoms(), // la distanza fra due dati successivi
                                   1, // la distanza fra due serie di dati adiacenti
                                   FFTW_PRESERVE_INPUT
                                   );
    fplan_natoms=traiettoria->get_natoms();
    fplan_size=size;
    fftw_free( dummy);
    //devo fare la trasformata della velocità per ogni atomo
    }

    fftw_execute_dft_r2c(fftw3,traiettoria->velocita_inizio(),trasformata);
    // adesso bisogna fare la media del modulo quadro, una per ogni tipo di atomo

    double * media = new double[tipi_atomi];
    unsigned int * cont = new unsigned int [tipi_atomi];

    for (unsigned int iomega=0;iomega<data_length/(tipi_atomi*3);iomega++) {
        //una media per dimensione
        for (unsigned int idim=0;idim<3;idim++){
            //media sugli atomi
            for (unsigned int itipi=0;itipi<tipi_atomi;itipi++) {
                media[itipi]=0.0;
                cont[itipi]=0;
            }
            //qui dentro dovrò distinguere i vari tipi di atomi
            for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
                //parti reale ed immaginaria
                double x=trasformata[iatom*3 // ho x y z , x y z, ... , natoms volte
                        +idim // seleziona x/y/z
                        +iomega*3*traiettoria->get_natoms()
                        ][0];
                double y=trasformata[iatom*3
                        +idim
                        +iomega*3*traiettoria->get_natoms()
                        ][1];
                double modulo2=x*x+y*y;
                int itipo=traiettoria->get_type(iatom);
                if (itipo < 0 ) {
		    throw std::runtime_error("Type index of the atom out of range (!?)");
                    itipo=0;
                } else if (itipo >= tipi_atomi) {
		    throw std::runtime_error("Type index of the atom out of range (non consecutive indexes?!)");
                    itipo=tipi_atomi-1;
                }
                double delta = modulo2-media[itipo];
                media[itipo]+=delta/(++(cont[itipo]));


            }
            for (unsigned int itipo=0;itipo < tipi_atomi;itipo++){
                vdata[itipo*data_length/tipi_atomi+ iomega*3+idim]=media[itipo];
            }
        }
    }
    delete [] media;
    delete [] cont;

    if (dump) {
       std::ofstream out("analisi_vibr.debug",std::fstream::app);
       for (unsigned int iomega=0;iomega<data_length/(tipi_atomi*3);iomega++){
           out << iomega;
           for (unsigned int itype=0;itype<tipi_atomi;itype++) {
               for (unsigned int idim=0;idim<3;idim++)
                    out << " " << spettro(iomega,idim,itype);
           }
           out << "\n";
       }
       out << "\n\n";
    }
}

template <class T>
double SpettroVibrazionale<T>::spettro(unsigned int frequenza, unsigned int dim,unsigned int tipo_atomo) {
    if (frequenza<data_length/(3*tipi_atomi) && dim<3 && tipo_atomo < tipi_atomi) {
        return vdata[tipo_atomo*data_length/tipi_atomi+ frequenza*3+dim];
    } else {
        throw std::runtime_error("Error: index out of range!\n");
    }
}

template <class T>
SpettroVibrazionale<T> & SpettroVibrazionale<T>::operator = (const SpettroVibrazionale<T> & destra) {
    VectorOp<SpettroVibrazionale<T>>::operator =(destra);

    fftw_free(trasformata);
//    deallocate_plan();
    trasformata=0;
    return *this;
}

#ifdef BUILD_MMAP
#include "traiettoria.h"
template class SpettroVibrazionale<Traiettoria>;
#endif

#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class SpettroVibrazionale<Traiettoria_numpy>;
#endif
