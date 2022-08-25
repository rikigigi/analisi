/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include <cmath>
#include<thread>
#include <array>

#include "correlatorespaziale.h"
#include "traiettoria.h"
#include "operazionisulista.h"
#include "config.h"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628

CorrelatoreSpaziale::CorrelatoreSpaziale(Traiettoria *t,
                                                 std::vector< std::array<double,3> >  ks,
                                                 double sigma2,
                                                 unsigned int nthreads,
                                                 unsigned int skip,
                                                 bool debug) :
    t{t},sfac(nullptr), sigma2(sigma2),debug(debug),tipi_atomi(0),ntimesteps(0), CalcolaMultiThread<CorrelatoreSpaziale> (nthreads,skip)
{
    nk=ks.size();
    klist=ks;
    size_k=0;
    size_sfac=0;
}


void CorrelatoreSpaziale::reset(const unsigned int numeroTimestepsPerBlocco) {
    if (ntimesteps>=numeroTimestepsPerBlocco && t->get_ntypes()>=tipi_atomi){
        ntimesteps=numeroTimestepsPerBlocco;
        return;
    }
    delete [] sfac;
    delete [] vdata;
    tipi_atomi=t->get_ntypes();
    ntimesteps=numeroTimestepsPerBlocco;
    size_k=3*tipi_atomi*2;// (complex number)
    size_sfac=size_k*nk;
    data_length=ntimesteps*size_sfac/skip;
    sfac = new double[size_sfac*nthreads];
    vdata = new double[data_length];
}


std::vector<ssize_t> CorrelatoreSpaziale::get_shape() const{
    return {ntimesteps/skip,nk,tipi_atomi,3,2};
}
std::vector<ssize_t> CorrelatoreSpaziale::get_stride() const{
    return {static_cast<long>(sizeof(double)*size_sfac),
                static_cast<long>(sizeof(double)*size_k),
                sizeof(double)*2*3,
                sizeof(double)*2,
                sizeof(double)};
}


void CorrelatoreSpaziale::s_fac_k(const double  k[3],
                                  const unsigned int i_t,
                                  double *out1 ///size 3*type {(x,y,z)_1, ..., (x,y,z)_ntype}
                                  ) const {
    for (unsigned int i_at=0;i_at<t->get_natoms();i_at++){
        const double * pos = t->posizioni(i_t,i_at);
        const double * vel = t->velocita(i_t,i_at);
        unsigned int type=t->get_type(i_at);
        double arg=k[0]*pos[0]+k[1]*pos[1]+k[2]*pos[2];
        double s=sin(arg),c=cos(arg);
        for (unsigned int icoord=0;icoord<3;icoord++){
            out1[3*type*2+icoord*2+0]+=c*vel[icoord];
            out1[3*type*2+icoord*2+1]+=s*vel[icoord];
        }
    }
}

void CorrelatoreSpaziale::calc_single_th(const unsigned int &start, const unsigned int &stop, const unsigned int &primo, const unsigned int &ith) noexcept {

    double * sfact=&sfac[size_sfac*ith];
    int ilista=(start-primo)/skip;
    for (unsigned int itimestep=start;itimestep<stop;itimestep+=skip){
        //calculate at itimestep in vdata ad itimestep+1 in sfac. Then put in itimestep the difference
        //loop over k
        int ik=0;

        for (unsigned int i=0;i<size_sfac;i++) {
            sfact[i]=0.0;
        }
        for (unsigned int i=0;i<size_sfac;i++) {
            vdata[size_sfac*ilista+i]=0.0;
        }
        for (auto &k : klist){
            double kmod=sqrt(k[0]*k[0]+k[1]*k[1]+k[2]*k[2]);
            if (kmod == 0) kmod=1.0;
            s_fac_k(k.data(),itimestep,&vdata[size_sfac*ilista+ik*size_k]);
            s_fac_k(k.data(),itimestep+1,&sfact[ik*size_k]);
            //put the difference inside vdata
            for (unsigned int i=0;i<size_k;i++){
                vdata[size_sfac*ilista+ik*size_k+i]=(sfact[ik*size_k+i]-vdata[size_sfac*ilista+ik*size_k+i])/kmod;
            }
            ++ik;
        }
        ilista++;
    }
} // end calcola

void CorrelatoreSpaziale::print(std::ostream &out){
    //print stuff
    for (unsigned int i=0;i<ntimesteps/skip;++i){
        for (unsigned int j=0;j<size_sfac;++j){
            out << vdata[i*size_sfac+j]<<" ";
        }
        out << std::endl;
    }
}


CorrelatoreSpaziale::~CorrelatoreSpaziale(){
    delete [] sfac;
}


