#include "testtraiettoria.h"
#include "traiettoria.h"
#include "traiettoria_numpy.h"
#include "config.h"

template <class T>
void TraiettoriaBase<T>::dump_lammps_bin_traj(const std::string &fname, int start_ts, int stop_ts){
    if (start_ts<0 || start_ts>=n_timesteps){
        throw std::runtime_error("You must provide a starting timestep between 0 and the number of timesteps!");
    }
    if (stop_ts<=0)
        stop_ts=n_timesteps;
    std::ofstream out(fname,std::ofstream::binary);
    for (int t=start_ts;t<stop_ts;++t){
        Intestazione_timestep_triclinic head;
        head.natoms=natoms;
        for (unsigned int i=0;i<6;++i)
            head.scatola[i]=scatola(t)[i];
        internal_to_lammps(head.scatola);
        head.timestep=t;
        head.triclinic=triclinic;
        if (triclinic) {
            for (unsigned int i=0;i<3;++i)
                head.xy_xz_yz[i]=scatola(t)[6+i];
        }
        head.condizioni_al_contorno[0]=0;
        head.condizioni_al_contorno[1]=0;
        head.condizioni_al_contorno[2]=0;
        head.condizioni_al_contorno[3]=0;
        head.condizioni_al_contorno[4]=0;
        head.condizioni_al_contorno[5]=0;
        head.dimensioni_riga_output=NDOUBLE_ATOMO;
        head.nchunk=1;
        int n_data=head.natoms*head.dimensioni_riga_output;
        //write timestep header
        head.write(out);
        out.write((char*) &n_data,sizeof(int));

        for (unsigned int iatom=0;iatom<head.natoms;++iatom) {
            double data[NDOUBLE_ATOMO];
            data[0]=iatom;
            data[1]=get_type(iatom);
            for (int i=0;i<3;++i){
                data[2+i]=posizioni(t,iatom)[i];
                data[5+i]=velocita(t,iatom)[i];
            }
            static_assert (NDOUBLE_ATOMO==8, "You have to change the file writing (what do I have to write?) this if you change NDOUBLE_ATOMO" );
            out.write((char*) data,sizeof(double)*NDOUBLE_ATOMO);
        }
    }
}

template <class T>
ssize_t TraiettoriaBase<T>::get_ntypes (){
    if (ntypes==0) {
        types.clear();
        min_type=buffer_tipi[0];
        max_type=buffer_tipi[0];
        bool *duplicati = new bool[natoms];
        for (size_t i=0;i<natoms;i++)
            duplicati[i]=false;
        for (size_t i=0;i<natoms;i++) {
            if (!duplicati[i]) {
                if (buffer_tipi[i]>max_type)
                    max_type=buffer_tipi[i];
                if (buffer_tipi[i]<min_type)
                    min_type=buffer_tipi[i];
                for (size_t j=i+1;j<natoms;j++){
                    if (buffer_tipi[j]==buffer_tipi[i]){
                        duplicati[j]=true;
                    }
                }
                types.push_back(buffer_tipi[i]);
                ntypes++;
            }
        }
        std::sort(types.begin(),types.end());
        type_map.clear();
        for (unsigned int i=0;i<types.size(); i++){
            type_map[types[i]]=i;
        }
        delete [] duplicati;
        masse = new double [ntypes];
        cariche = new double [ntypes];

        for (size_t i=0;i<natoms;i++) {
            buffer_tipi_id[i]=type_map.at(buffer_tipi[i]);
        }
    }
    return ntypes;
}

#ifdef PYTHON_SUPPORT
template class TraiettoriaBase<Traiettoria_numpy>;
#endif

#ifdef BUILD_MMAP
template class TraiettoriaBase<Traiettoria>;
#endif
