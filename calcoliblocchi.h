#ifndef CALCOLIBLOCCHI_H
#define CALCOLIBLOCCHI_H

#include <vector>
#include <iostream>
#include <assert.h>


template <class T> class MediaVar{
public:
    MediaVar(T*Tmedio, T*Tvar,T*delta,T*tmp) : Tmedio(Tmedio), Tvar(Tvar),delta(delta),tmp(tmp) {}

    void calcola_begin(unsigned int s,T * calc){
        Tmedio->reset(s);
        Tmedio->azzera();
        Tvar->reset(s);
        Tvar->azzera();
        delta->reset(s);
        tmp->reset(s);
        iblock=0;
    }

    void calcola(T* calcolo) {
#ifdef DEBUG2
        std::cerr << "*delta = *calcolo - *Tmedio;\n";
#endif
        //algoritmo media e varianza online
        *delta = *calcolo;
        *delta -= *Tmedio;
        //*delta = *calcolo - *Tmedio
#ifdef DEBUG2
        std::cerr << "*delta = *calcolo - *Tmedio;\n";
#endif
        *tmp = *delta;
        *tmp /= (iblock+1);
        *Tmedio += *tmp;
        //*Tmedio += (*delta)/(iblock+1); // calcolo / double , calcolo+=calcolo
#ifdef DEBUG2
        std::cerr << "*Tvar += (*delta)*(*calcolo-*Tmedio);\n";
#endif
        *tmp = *calcolo;
        *tmp -= *Tmedio;
        *tmp *= *delta;
        *Tvar += *tmp;
        iblock++;
    }
    void calcola_end(unsigned int n_b){
        *Tvar/=((n_b-1)*n_b);
    }
private:
    T *Tmedio,*Tvar,*calcolo,*delta,*tmp;
    unsigned int iblock;
};


/**
  * Calcola, oltre la media e la varianza, anche la covarianza
**/

template <class T> class MediaVarCovar{
public:
    MediaVarCovar(unsigned int n_columns,
                  std::vector< std::pair<unsigned int,unsigned int> > Cvar_list
            ) : n_columns(n_columns),Cvar_list(Cvar_list),dx(0),dx2(0) {

        for (unsigned int i=0;i<Cvar_list.size();i++){
            if(Cvar_list.at(i).first>=n_columns || Cvar_list.at(i).second>=n_columns ){
                std::cerr << "Errore: indice della covarianza fuori dal range!\n";
                abort();
            }
        }

    }

    ~MediaVarCovar(){
        dealloc();
    }

    void calcola_begin(unsigned int si_,T*c){
        dealloc();
        unsigned int size_blockT=c->lunghezza();
        assert(size_blockT%n_columns==0);
        size_block=size_blockT/n_columns;
        iblock=0;
        for (unsigned int icol=0;icol<n_columns;icol++){
            medie.push_back(new double[size_block]);
            varianze.push_back(new double[size_block]);
            for (unsigned int i=0;i<size_block;i++){
                varianze.back()[i]=0;
                medie.back()[i]=0;
            }
        }
        for (unsigned int icvar=0;icvar<Cvar_list.size();icvar++){
            covarianze.push_back(new double[size_block]);
            for (unsigned int i=0;i<size_block;i++){
                covarianze.back()[i]=0;
            }
        }
        dx=new double[n_columns];
        dx2=new double[n_columns];
    }

    void dealloc(){
        for(unsigned int i=0;i<medie.size();i++){
            delete [] medie.at(i);
            delete [] varianze.at(i);
        }
        for(unsigned int i=0;i<covarianze.size();i++){
            delete [] covarianze.at(i);
        }
        medie.clear();
        varianze.clear();
        covarianze.clear();
        delete [] dx;
        delete [] dx2;
    }

    void calcola(T* calcolo) {
        assert(dx2!=0 && dx!=0);
        assert(calcolo->lunghezza()%n_columns==0);
        assert(size_block==calcolo->lunghezza()/n_columns);
        iblock++;


        for (unsigned int ib=0;ib<size_block;ib++){
            for (unsigned int icol=0;icol<n_columns;icol++){

                //calcola dx=x-media
                dx[icol]=calcolo->elemento(ib*n_columns+icol)-medie[icol][ib];
                //aggiorna media+=dx/iblock
                medie[icol][ib]+=dx[icol]/iblock;
                //calcola dx2=x-media
                dx2[icol]=calcolo->elemento(ib*n_columns+icol)-medie[icol][ib];
                //calcola varianze
                varianze[icol][ib]+=dx[icol]*dx2[icol];
            }
            //calcola covarianze
            for (unsigned int icvar=0;icvar<Cvar_list.size();icvar++){
                covarianze[icvar][ib]+=dx[Cvar_list[icvar].first]*dx2[Cvar_list[icvar].second];
            }
        }


    }
    void calcola_end(unsigned int n_b){
        assert(iblock>1);
        assert(n_b==iblock);
        for (unsigned int ib=0;ib<size_block;ib++){
            for (unsigned int icol=0;icol<n_columns;icol++){
                varianze[icol][ib]/=(iblock*(iblock-1));
            }
            for (unsigned int icvar=0;icvar<Cvar_list.size();icvar++){
                covarianze[icvar][ib]/=(iblock*(iblock-1));
            }
        }

        delete [] dx;
        delete [] dx2;
        dx=0;dx2=0;

    }

    double * media(unsigned int col_idx){
        return medie.at(col_idx);
    }
    double * varianza(unsigned int col_idx){
        return varianze.at(col_idx);
    }
    double * covarianza(unsigned int cvar_idx){
        return covarianze.at(cvar_idx);
    }

    unsigned int size(){return size_block;}
    unsigned int n_column(){return medie.size();}
    unsigned int n_cvar(){return covarianze.size();}

private:
    std::vector<double*>medie,varianze,covarianze;
    double *dx,*dx2;
    std::vector<unsigned int> idx;
    std::vector< std::pair<unsigned int,unsigned int> > Cvar_list;
    unsigned int n_columns,size_block;
    unsigned int iblock;
};

#endif // CALCOLIBLOCCHI_H
