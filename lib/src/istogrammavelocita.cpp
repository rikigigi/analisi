/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "istogrammavelocita.h"
#include "mediablocchi.h"


void istg(double a,/// estremo inferiore dell'intervallo da istogrammare
           double b,///estremo superiore dell'intervallo da istogrammare
           unsigned int z,/// numero di intervalli dell'istogramma
           double *ist,/// array che conterrà dopo l'istogramma calcolato. E' responsabilità dell'utente la sua inizializzazione (a zero o a quello che si vuole)
           double val ///elemento da inserire nell'istogramma
           ){
double k;
    k=(val-a)/(b-a)*z;
     if (k<0) {k=0;} else if (k>z-1) {k=z-1;}
    ist[(unsigned int) k]++;

}

IstogrammaVelocita::IstogrammaVelocita(Traiettoria *t, unsigned int nbins,double vminmax_) : OperazioniSuLista<IstogrammaVelocita> ()
{

    traiettoria=t;

    lunghezza_lista=traiettoria->get_ntypes()*3*nbins;
    lista = new double [lunghezza_lista];
    bins=nbins;
    vminmax=vminmax_;


}

unsigned int IstogrammaVelocita::numeroTimestepsOltreFineBlocco(unsigned int n_b) {
    return 0;
}

IstogrammaVelocita & IstogrammaVelocita::operator = (const IstogrammaVelocita & destra) {
    OperazioniSuLista<IstogrammaVelocita>::operator =(destra);
    return *this;
}

void IstogrammaVelocita::reset(const unsigned int numeroTimestepsPerBlocco) {

    for (unsigned int i=0;i<lunghezza_lista;i++)
        lista[i]=0.0;

    ntimestep=numeroTimestepsPerBlocco;

}

void IstogrammaVelocita::calcola(unsigned int primo) {

    for (unsigned int itimestep=primo;itimestep<primo+ntimestep;itimestep++)
        for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
            for (unsigned int icoord=0;icoord<3;icoord++)
                istg(-vminmax,+vminmax,bins,&lista[(traiettoria->get_type(iatom)*3+icoord)*bins],traiettoria->velocita(itimestep,iatom)[icoord]);
        }



}




//template class MediaBlocchi<IstogrammaVelocita,unsigned int,double>;
