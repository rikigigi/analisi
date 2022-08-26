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

IstogrammaVelocita::IstogrammaVelocita(Trajectory *t, unsigned int nbins,double vminmax_) : VectorOp<IstogrammaVelocita> ()
{

    traiettoria=t;

    data_length=traiettoria->get_ntypes()*3*nbins;
    vdata = new double [data_length];
    bins=nbins;
    vminmax=vminmax_;


}

unsigned int IstogrammaVelocita::nExtraTimesteps(unsigned int n_b) {
    return 0;
}

IstogrammaVelocita & IstogrammaVelocita::operator = (const IstogrammaVelocita & destra) {
    VectorOp<IstogrammaVelocita>::operator =(destra);
    return *this;
}

void IstogrammaVelocita::reset(const unsigned int numeroTimestepsPerBlocco) {

    for (unsigned int i=0;i<data_length;i++)
        vdata[i]=0.0;

    ntimestep=numeroTimestepsPerBlocco;

}

void IstogrammaVelocita::calculate(unsigned int primo) {

    for (unsigned int itimestep=primo;itimestep<primo+ntimestep;itimestep++)
        for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
            for (unsigned int icoord=0;icoord<3;icoord++)
                istg(-vminmax,+vminmax,bins,&vdata[(traiettoria->get_type(iatom)*3+icoord)*bins],traiettoria->velocity(itimestep,iatom)[icoord]);
        }



}




//template class BlockAverage<IstogrammaVelocita,unsigned int,double>;
