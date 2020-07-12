#include "sphericalcorrelations.h"
#include "specialfunctions.h"
#include <vector>
#include <thread>
#include <sstream>
#include "config.h"
#include "calc_buffer.h"

template <int l, class TFLOAT, class T>
SphericalCorrelations<l,TFLOAT,T>::SphericalCorrelations(T *t,
                                                         TFLOAT rmin,
                                                         TFLOAT rmax,
                                                         unsigned int nbin,
                                                         unsigned int tmax,
                                                         unsigned int nthreads,
                                                         unsigned int skip,
                                                         bool debug) :
t{*t},rmin{rmin},rmax{rmax},nbin{nbin}, skip{skip}, tmax{tmax}, nthreads{nthreads}, debug{debug}{

}


template <int l, class TFLOAT, class T>
SphericalCorrelations<l,TFLOAT,T>::~SphericalCorrelations() {

}


template <int l, class TFLOAT, class T>
unsigned int SphericalCorrelations<l,TFLOAT,T>::numeroTimestepsOltreFineBlocco(unsigned int n_b) {
    return (t.get_ntimesteps()/(n_b+1)+1 < tmax || tmax==0)? t.get_ntimesteps()/(n_b+1)+1 : tmax;
}

template <int l, class TFLOAT, class T>
void SphericalCorrelations<l,TFLOAT,T>::reset(const unsigned int numeroTimestepsPerBlocco) {


    //description of output
    std::stringstream descr;

    descr << "#TODO"<<std::endl;
    c_descr=descr.str();


    //quanti timestep è lunga la funzione di correlazione
    leff =(numeroTimestepsPerBlocco<tmax || tmax==0)? numeroTimestepsPerBlocco : tmax;
    //numero di timestep su cui fare la media
    ntimesteps=numeroTimestepsPerBlocco;
    natoms=t.get_natoms();
    ntypes=t.get_ntypes();
    lunghezza_lista=leff*ntypes*ntypes*nbin*(l+1);

    dr=(rmax-rmin)/nbin;

    delete [] lista;
    lista=new TFLOAT [lunghezza_lista];
}


template <int lmax, class TFLOAT, class T>
void SphericalCorrelations<lmax,TFLOAT,T>::calc(int timestep, TFLOAT *result, TFLOAT *workspace, TFLOAT * cheby, double *l) {
    //zero result
    for (int i=0;i<(lmax+1)*(lmax+1)*natoms*nbin*ntypes;++i) {
        result[i]=0;
    }
    for (unsigned int iatom=0;iatom<natoms;iatom++) {
        //other atom loop
        for (unsigned int jatom=0;jatom<natoms;jatom++) {
            unsigned int jtype=t.get_type(jatom);
            if (iatom==jatom)
                continue;

            //minimum image distance
            double x[3];
            double d=sqrt(t.d2_minImage(iatom,jatom,timestep,timestep,l,x));
            //bin index
            int idx=(int)floorf((d-rmin)/dr);

            if (idx<nbin && idx >= 0){
                //calculate sin and cos
                //calculate spherical harmonics and add to the correct average
                SpecialFunctions::SphericalHarmonics<lmax,TFLOAT,true,true> sh(x[0],x[1],x[2],cheby,workspace);
                sh.calc();
                //now in workspace you have all the spherical harmonics components of the density of the current jatom around iatom
                //add to the sh density of the current iatom of the current bin of the type jtype
                for (int ll=0;ll<(lmax+1)*(lmax+1);++ll) {
                    result[index_wrk(iatom,jtype,idx)+ll]+=workspace[ll];
                }
            }
        }
    }
}

template <int lmax, class TFLOAT, class T>
void SphericalCorrelations<lmax,TFLOAT,T>::calcola(unsigned int primo) {


    if (nthreads<=1){
        std::cerr << "Attenzione: sto usando un solo thread.\n";
        nthreads=1;
    }

    unsigned int npassith=leff/nthreads;
    std::vector<std::thread> threads;



    for (unsigned int  ith=0;ith<nthreads;ith++) {
        threads.push_back(std::thread([&,ith](){
            unsigned int ultimo= (ith != nthreads-1 )?npassith*(ith+1):leff;
            unsigned int itimestep = primo;
            double l[3]={t.scatola(itimestep)[1]-t.scatola(itimestep)[0],
                         t.scatola(itimestep)[3]-t.scatola(itimestep)[2],
                         t.scatola(itimestep)[5]-t.scatola(itimestep)[4]};

            //working space for spherical harmonics
            TFLOAT workspace[(lmax+1)*(lmax+1)];
            TFLOAT cheby[2*(lmax+1)];
            //working space for averages of correlations over atomic types
            //for every atom I have to allocate space for the (l,m) matrix of coefficients (used to describe the environment around it)
            //then I calculate separately the (l,m) matrix for each type in the environment and for different values of the radial distance (bins)
            int sh_single_type_size=(lmax+1)*(lmax+1)*nbin*ntypes;
            //this is for all atoms. I will compute correlations of objects that are such big
            int sh_snap_size=sh_single_type_size*natoms;
            //at the end I will first sum over all m: the (l,m) matrix will become an (l) vector
            //then I do the average over the central atoms of the same type: the size of the result (for each time lag) is:
            int sh_final_size=sh_single_type_size*ntypes/(lmax+1); // remember: we sum over m before averaging on central atoms

            //allocate space
            TFLOAT *aveWork=new TFLOAT[sh_snap_size +sh_final_size];
            TFLOAT *aveWork1=aveWork;
            TFLOAT *aveTypes=aveWork+  sh_snap_size;
            int *avecont=new int[ntypes];

            //buffer for few sh calculations
           CalcBuffer<TFLOAT> buffer(30,sh_snap_size);

            //loop over time differences -- eheh, did you forgot about time lags?
            for (unsigned int dt=npassith*ith;dt<ultimo;dt++){

                //set to zero the result for this time lag
                azzera(index(dt,0,0),index(dt,ntypes-1,ntypes-1));

                //average over starting timestep
                for (unsigned int imedia=0;imedia<ntimesteps;imedia+=skip){
                    //center atom loop for the snapshot at imedia
                    TFLOAT * sh1=
                    buffer.buffer_calc(*this,imedia+primo,workspace,cheby,l);

                    calc(imedia+primo,aveWork1,workspace,cheby,l);
                    for (int i=0;i<sh_snap_size;++i){
                        if (sh1[i]!=aveWork1[i]){
                            std::cerr << "Error in comparison: it is a bug\n";
                            abort();
                        }
                    }
                    //center atom loop for the snapshot at imedia+dt
                    TFLOAT * sh2=
                    buffer.buffer_calc(*this,imedia+primo+dt,workspace,cheby,l);

                    calc(imedia+primo+dt,aveWork1,workspace,cheby,l);
                    for (int i=0;i<sh_snap_size;++i){
                        if (sh2[i]!=aveWork1[i]){
                            std::cerr << "Error in comparison: it is a bug\n";
                            abort();
                        }
                    }
                    //compute correlations -- this is simple!
                    for (int i=0;i<sh_snap_size;++i) {
                        aveWork1[i]=sh1[i]*sh2[i];
                    }

                    //here we have to average first over m... (otherwise we calculate an average of something that is not invariant under rotation, so we probably get something that always decays to zero, unless we are in a crystal
                    //sum over m: sum everything and put it where the m=0 was placed before
                    //use a nice recurrent structure to keep track of the location of the data
                    using MultiVal = SpecialFunctions::MultiValDynamic::MultiVal<lmax,TFLOAT>;
                    MultiVal v;
                    for (int iatom=0;iatom<natoms;++iatom) {
                        for (int jtype=0;jtype<ntypes;++jtype){
                            for (int ibin=0;ibin<nbin;++ibin){
                                v.init(aveWork1+index_wrk(iatom,jtype,ibin));
                                v.sum_in_m_zero();
                            }
                        }
                    }

                    //do averages over atomic types
                    for (int i=0;i<ntypes;++i){
                        avecont[i]=0;
                    }
                    for (int i=0;i<sh_final_size;++i) {
                        aveTypes[i]=0.0;
                    }
                    for (int iatom=0;iatom<natoms;++iatom) { // loop over center atoms and average correlation functions according to their type
                        int itype=t.get_type(iatom);
                        avecont[itype]++;
                        //average only the sum over m value, that is placed in the memory location corresponding to m=0.
                        // Loop over chunks of (l,m) coefficients ( (l+1)^2 numbers)
                        for (int jtype=0;jtype<ntypes;++jtype)
                            for (int ibin=0;ibin<nbin;++ibin) {
                                v.init(aveWork1+index_wrk(iatom,jtype,ibin)); //set the nice recursive structure that gives some order to the coefficients
                                //the following does: average += (value-average)/counter.
                                //A nice recursive function is used to avoid error prone indexes of arrays with different strides.
                                // (it use metaprogramming tecniques: it should unroll with compiler optimization on!)
                                v.running_average_m_zero_vector(aveTypes+index(0,itype,jtype,ibin),TFLOAT(avecont[itype]));
                            }
//                            for (int ll=0;ll<lmax;++ll){
//                                aveTypes[index(0,itype,jtype,ibin)+ll]+=(aveWork1[index_wrk(iatom,jtype,ibin)+__LL__]-aveTypes[index(0,itype,jtype,ibin)+ll])/TFLOAT(avecont[itype]);
//                            }
//                        for (int ll=0;ll<sh_single_type_size;++ll){
//                            int idx=itype*sh_single_type_size+ll;
//                            aveTypes[idx]+=(aveWork1[sh_single_type_size*iatom+ll]-aveTypes[idx])/TFLOAT(avecont[itype]);
//                        }
                    }
                    //finally add to big average over starting timestep (imedia loop) -- I know, there are a lot of averages and sums and stuff
                    for (int ll=0;ll<sh_final_size;++ll){
                        lista[index(dt,0,0)+ll]+=(aveTypes[ll]-lista[index(dt,0,0)+ll])/TFLOAT((imedia/skip)+1);
                    }

                }

            }

            delete [] aveWork;
            delete [] avecont;

        }));
    }
    //Wasn't it multithreaded?
    for (unsigned int  ith=0;ith<nthreads;ith++)
        threads[ith].join();
    threads.clear();


}

template class SphericalCorrelations<10,double,Traiettoria>;

#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class SphericalCorrelations<10,double,Traiettoria_numpy>;
#endif
