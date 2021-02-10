/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "convertibinario.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <utility>
#include "lammps_struct.h"
#include "config.h"

//#define XDR_FILE

#ifdef XDR_FILE

typedef struct XDRFILE XDRFILE;
extern "C" XDRFILE * xdrfile_open(const char *,const char*);
extern "C" int read_trr_natoms(const char *fn,int *natoms);
#define DIM 3
typedef float matrix[DIM][DIM];
typedef float rvec[DIM];
extern "C" int read_trr(XDRFILE *xd,int natoms,int *step,float *t,float *lambda,matrix box,rvec *x,rvec *v,rvec *f);
#endif //XDR_FILE


ConvertiBinario::ConvertiBinario(const std::string filein, const std::string fileout, Type tipo, const std::string typefile)
{
    bool first_message=true;

    if (tipo!=natoms_box_xyz_vxvyvz && tipo!=gromax_trr) {
        std::cerr << "Non implementato!\n";
        abort();
    }


    std::ofstream out(fileout,std::ofstream::binary);
#ifdef XDR_FILE
    XDRFILE * xd;
    int natoms_xdr=0;
    int xdr_res=0;
    matrix box;
    rvec * xx, * vv, * ff;
    int step;
    float time,lambda;
    std::vector<std::pair<int,int>> types;
#endif
    std::ifstream in;
    if (tipo==natoms_box_xyz_vxvyvz){
        in.open(filein);
    } else if (tipo==gromax_trr){
#ifdef XDR_FILE
        if (typefile=="") {
            std::cerr << "Per i file di gromacs è necessario un file aggiuntivo con su ogni riga id e tipo dell'atomo corrispondente!\n";
            abort();
        } else {
            in.open(typefile);
            int type=0,id=0;
            in >> id >> type;
            while (in.good()) {
                types.push_back(std::pair<int,int>(id,type));
                in >> id >> type;
            }
            in.close();
        }
        xdr_res=read_trr_natoms(filein.c_str() , &natoms_xdr);
        xd = xdrfile_open(filein.c_str(),"r");

        if (natoms_xdr!=types.size()) {
            std::cerr << "Errore nella lettura del file dei tipi: il numero di dati ("<<types.size()<<") non corrisponde al numero di atomi letto dal file trr("<<natoms_xdr<<").\n";
            abort();
        }

        xx=new rvec[natoms_xdr];
        vv=new rvec[natoms_xdr];
        ff=0;

#else
        std::cerr << "Non implementato!\n";
        abort();
#endif
    }

    bigint itim=0;
    bool good=true;
    while (good) {
        ++itim;
        //per ogni timestep, legge tutti i dati necessari nella struttura dati
        Intestazione_timestep head;


        if (tipo==natoms_box_xyz_vxvyvz){
            in >> head.natoms >> head.scatola[0] >> head.scatola[1] >> head.scatola[2] >> head.scatola[3] >> head.scatola[4] >> head.scatola[5];
            head.timestep=itim;
            if (!in.good()){
                if (in.eof()) {
                    std::cerr << "Fine del file raggiunta al timestep "<<itim-1<<"\n";
                } else {
                    std::cerr << "Errore nella lettura del file al timestep (timestep "<<itim<<")!\n";
                }
                break;
            }
        } else if (tipo==gromax_trr) {

#ifdef XDR_FILE
            xdr_res=read_trr(xd, natoms_xdr, &step, &time, &lambda, box, xx, vv, ff);
            if (xdr_res!=0) {
                std::cerr << "Fine del file raggiunta dopo "<<itim-1<<" frames.\n";
                break;
            }
            //controlla che la scatola sia un parallelepipedo
            for (unsigned int i=0;i<3;i++)
                for (unsigned int j=0;j<3;j++)
                    if (first_message && i!=j && box[i][j]!=0) {
                        first_message=false;
                        std::cerr << "!!Attenzione!! NON HO IMPLEMENTATO LA LETTURA DI SCATOLE DI SIMULAZIONI NON A FORMA DI PARALLELEPIPEDO: le dimensioni della scatola saranno sbagliate, assicurarsi che i comandi successivi non utilizzino questo dato.\n";
                    }
            head.natoms=natoms_xdr;
            head.scatola[0]=0.0;
            head.scatola[2]=0.0;
            head.scatola[4]=0.0;
            head.scatola[1]=box[0][0];
            head.scatola[3]=box[1][1];
            head.scatola[5]=box[2][2];
            head.timestep=step;


#else
            std::cerr << "Non implementato!\n";
            abort();
#endif



        }

        head.triclinic=false;
        head.condizioni_al_contorno[0]=0;
        head.condizioni_al_contorno[1]=0;
        head.condizioni_al_contorno[2]=0;
        head.condizioni_al_contorno[3]=0;
        head.condizioni_al_contorno[4]=0;
        head.condizioni_al_contorno[5]=0;
        head.dimensioni_riga_output=NDOUBLE_ATOMO;
        head.nchunk=1;
        int n_data=head.natoms*NDOUBLE_ATOMO;
        //qui ho impostato l'intestazione, la scrivo
        out.write((char*) &head,sizeof(Intestazione_timestep));
        out.write((char*) &n_data,sizeof(int));
        //adesso devo scrivere tutte le righe, una dopo l'altra


        for (unsigned int i=0;i<head.natoms;++i) {
            double data[NDOUBLE_ATOMO];
            if (tipo==natoms_box_xyz_vxvyvz){
                for (int j=0; j<NDOUBLE_ATOMO;++j)
                    in >> data[j];
                if (!in.good()){
                    std::cerr << "Errore nella lettura del file (timestep "<<itim<<")!\nATTENZIONE: file di output non completo!\n";
                    break;
                }
                good=in.good();
            } else if (tipo==gromax_trr) {
#ifdef XDR_FILE
                static_assert (NDOUBLE_ATOMO==8, "I don't know what to do with xdr files if NDOUBLE_ATOMO is not 8" );
                data[0]=types[i].first;
                data[1]=types[i].second;
                data[2]=xx[i][0];
                data[3]=xx[i][1];
                data[4]=xx[i][2];
                data[6]=vv[i][0];
                data[5]=vv[i][1];
                data[7]=vv[i][2];
#else
                std::cerr << "Non implementato!\n";
                abort();
#endif
            }
            out.write((char*) data,sizeof(double)*NDOUBLE_ATOMO);
        }
    }
#ifdef XDR_FILE

    if(tipo==gromax_trr) {
        delete [] xx;
        delete [] vv;
    }
#endif
}
