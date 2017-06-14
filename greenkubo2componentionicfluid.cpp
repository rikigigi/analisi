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

GreenKubo2ComponentIonicFluid::GreenKubo2ComponentIonicFluid(Traiettoria *t, std::string log, double * cariche, unsigned int skip) : OperazioniSuLista<GreenKubo2ComponentIonicFluid>(),
    traiettoria (t), log(log), ntimesteps(0), skip(skip),je(0),jz(0)
{
// ricordarsi di impostare le cariche delgli atomi!
    traiettoria->set_charge(0,cariche[0]);
    traiettoria->set_charge(1,cariche[1]);
    traiettoria->index_all();
    je =new HeatFluxTs(log,traiettoria,skip);
    jz = new ChargeFluxTs(traiettoria);

}

GreenKubo2ComponentIonicFluid::~GreenKubo2ComponentIonicFluid(){
#ifdef DEBUG
    std::cerr << "Called delete GreenKubo2ComponentIonicFluid\n";
#endif
    delete je;
    delete jz;
}

unsigned int GreenKubo2ComponentIonicFluid::numeroTimestepsOltreFineBlocco(unsigned int n_b) {
    return traiettoria->get_ntimesteps()/(n_b+1)+1;
}

void GreenKubo2ComponentIonicFluid::reset(unsigned int numeroTimestepsPerBlocco) {
    lunghezza_lista=(numeroTimestepsPerBlocco/skip+1)*9; // Jee,Jzz,Jez,Jintee,Jintzz,Jintez,lambda
    ntimesteps=numeroTimestepsPerBlocco;
    delete [] lista;
    lista=new double [lunghezza_lista];
    jz->reset(numeroTimestepsPerBlocco*2);
}

void GreenKubo2ComponentIonicFluid::calcola(unsigned int primo) {

    jz->calcola(primo);

    unsigned int allinea=primo%skip;

    double intzz=0.0;
    double intee=0.0;
    double intez=0.0;
    double intze=0.0;
    for (unsigned int itimestep=0;itimestep<ntimesteps;itimestep+=skip) {
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
        //integrale con il metodo dei trapezi
        if (itimestep==0 || itimestep>=ntimesteps-skip){
            intee+=jee/2.0;
            intez+=jez/2.0;
            intze+=jze/2.0;
            intzz+=jzz/2.0;
        } else {
            intee+=jee;
            intez+=jez;
            intze+=jze;
            intzz+=jzz;
        }
        lista[(itimestep/skip)*9+3]=intee;
        lista[(itimestep/skip)*9+4]=intzz;
        lista[(itimestep/skip)*9+5]=intez;
        lista[(itimestep/skip)*9+6]=intee-intez*intez/intzz;
        lista[(itimestep/skip)*9+7]=jze;
        lista[(itimestep/skip)*9+8]=intze;
    }

}
