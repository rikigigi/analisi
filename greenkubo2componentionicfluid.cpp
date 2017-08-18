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

#include "greenkubo2componentionicfluid.h"
#include "chargefluxts.h"
#include "heatfluxts.h"
#include <string>
#include <fstream>

GreenKubo2ComponentIonicFluid::GreenKubo2ComponentIonicFluid(Traiettoria *t, std::string log, double * cariche, unsigned int skip, bool dump,unsigned int lunghezza_funzione_max) : OperazioniSuLista<GreenKubo2ComponentIonicFluid>(),
    traiettoria (t), log(log), ntimesteps(0), skip(skip),je(0),jz(0),scrivi_file(dump),lmax(lunghezza_funzione_max)
{
// ricordarsi di impostare le cariche delgli atomi!
    traiettoria->set_charge(0,cariche[0]);
    traiettoria->set_charge(1,cariche[1]);
    traiettoria->index_all();
    je =new HeatFluxTs(log,traiettoria,skip);
    jz = new ChargeFluxTs(traiettoria);

}

GreenKubo2ComponentIonicFluid::~GreenKubo2ComponentIonicFluid(){
#ifdef DEBUG2
    std::cerr << "Called delete GreenKubo2ComponentIonicFluid, je="<<je<<", jz="<<jz<<"\n";
#endif
    delete je;
    delete jz;
}

GreenKubo2ComponentIonicFluid & GreenKubo2ComponentIonicFluid::operator =(const GreenKubo2ComponentIonicFluid & destra) {
#ifdef DEBUG2
    std::cerr << "Chiamato GreenKubo2ComponentIonicFluid::operator =\n";
#endif
    OperazioniSuLista<GreenKubo2ComponentIonicFluid>::operator =( destra);
    return *this;
}

unsigned int GreenKubo2ComponentIonicFluid::numeroTimestepsOltreFineBlocco(unsigned int n_b) {
    return (traiettoria->get_ntimesteps()/(n_b+1)+1 < lmax || lmax==0)? traiettoria->get_ntimesteps()/(n_b+1)+1 : lmax;
}

void GreenKubo2ComponentIonicFluid::reset(unsigned int numeroTimestepsPerBlocco) {
    leff=(numeroTimestepsPerBlocco<lmax || lmax==0)? numeroTimestepsPerBlocco : lmax;
    lunghezza_lista=(leff/skip+1)*9; // Jee,Jzz,Jez,Jintee,Jintzz,Jintez,lambda
    ntimesteps=numeroTimestepsPerBlocco;
    delete [] lista;
    lista=new double [lunghezza_lista];
    jz->reset(numeroTimestepsPerBlocco+leff);
}

void GreenKubo2ComponentIonicFluid::calcola(unsigned int primo) {

    jz->calcola(primo);

    unsigned int allinea=primo%skip;

    double intzz=0.0;
    double intee=0.0;
    double intez=0.0;
    double intze=0.0;
    double jeeo=0.0,jzzo=0.0,jezo=0.0,jzeo=0.0;
    for (unsigned int itimestep=0;itimestep<leff;itimestep+=skip) {
        //fa la media sulla traiettoria dei vari prodotti,
        //con una differenza di timesteps fissata "itimestep"
        double jee=0.0,jzz=0.0,jez=0.0,jze=0.0;
        unsigned int cont=0;
        for (unsigned int jmedia=allinea;jmedia<ntimesteps;jmedia+=skip) {
            //prodotto JzJz
            double deltazz=(jz->J_z(primo+jmedia)[0]*jz->J_z(primo+jmedia+itimestep)[0]+
                            jz->J_z(primo+jmedia)[1]*jz->J_z(primo+jmedia+itimestep)[1]+
                            jz->J_z(primo+jmedia)[2]*jz->J_z(primo+jmedia+itimestep)[2])/3.0
                                - jzz;
            double deltaez=(je->flux(primo+jmedia)[0]*jz->J_z(primo+jmedia+itimestep)[0]+
                            je->flux(primo+jmedia)[1]*jz->J_z(primo+jmedia+itimestep)[1]+
                            je->flux(primo+jmedia)[2]*jz->J_z(primo+jmedia+itimestep)[2])/3.0
                                - jez;
            double deltaze=(je->flux(primo+jmedia+itimestep)[0]*jz->J_z(primo+jmedia)[0]+
                            je->flux(primo+jmedia+itimestep)[1]*jz->J_z(primo+jmedia)[1]+
                            je->flux(primo+jmedia+itimestep)[2]*jz->J_z(primo+jmedia)[2])/3.0
                                - jze;
            double deltaee=(je->flux(primo+jmedia)[0]*je->flux(primo+jmedia+itimestep)[0]+
                            je->flux(primo+jmedia)[1]*je->flux(primo+jmedia+itimestep)[1]+
                            je->flux(primo+jmedia)[2]*je->flux(primo+jmedia+itimestep)[2])/3.0
                                - jee;
            jzz+=deltazz/(++cont);
            jee+=deltaee/(cont);
            jez+=deltaez/(cont);
            jze+=deltaze/(cont);
        }
        lista[(itimestep/skip)*9+0]=jee;
        lista[(itimestep/skip)*9+1]=jzz;
        lista[(itimestep/skip)*9+2]=jez;
        lista[(itimestep/skip)*9+7]=jze;
        //integrale con il metodo dei trapezi

            intee+=jeeo;
            intez+=jezo;
            intze+=jzeo;
            intzz+=jzzo;

        lista[(itimestep/skip)*9+3]=intee+jee/2.0;
        lista[(itimestep/skip)*9+4]=intzz+jzz/2.0;
        lista[(itimestep/skip)*9+5]=intez+jez/2.0;
        lista[(itimestep/skip)*9+6]=intee+jee/2.0-(intez+jez/2.0)*(intez+jez/2.0)/(intzz+jzz/2.0);
        lista[(itimestep/skip)*9+8]=intze+jze/2.0;

        jeeo=jee;
        jezo=jez;
        jzeo=jze;
        jzzo=jzz;
    }

    if (scrivi_file) {
        std::ofstream outfile(log+".greekdump",std::ios::app);
        for (unsigned int itimestep=0;itimestep<leff;itimestep+=skip) {
            for (unsigned int j=0;j<9;j++){
                outfile << lista[(itimestep/skip)*9+j] << " ";
            }
            outfile << "\n";
        }
        outfile << "\n\n";
    }

}
