#include "sphericalcorrelations.h"
#include "specialfunctions.h"
#include <vector>
#include <thread>
#include <sstream>
#include "config.h"
#include "calc_buffer.h"
#include "twoloopsplit.h"

template <int l, class TFLOAT, class T>
SphericalCorrelations<l,TFLOAT,T>::SphericalCorrelations(T *t,
                                                         const rminmax_t rminmax,
                                                         size_t nbin,
                                                         size_t tmax,
                                                         size_t nthreads,
                                                         size_t skip,
                                                         size_t buffer_size,
                                                         bool debug, const NeighListSpec neighList) :
t{*t},nbin{nbin}, skip{skip}, tmax{tmax}, nthreads{nthreads}, debug{debug},buffer_size{buffer_size}, neighList{neighList},
SPB{t,nbin,rminmax},
ntypes{static_cast<size_t>(t->get_ntypes())},natoms{t->get_natoms()}{

}


template <int l, class TFLOAT, class T>
SphericalCorrelations<l,TFLOAT,T>::~SphericalCorrelations() {

}


template <int l, class TFLOAT, class T>
unsigned int SphericalCorrelations<l,TFLOAT,T>::nExtraTimesteps(unsigned int n_b) {
    return (t.get_ntimesteps()/(n_b+1)+1 < tmax || tmax==0)? t.get_ntimesteps()/(n_b+1)+1 : tmax;
}

template <int l, class TFLOAT, class T>
void SphericalCorrelations<l,TFLOAT,T>::reset(const unsigned int numeroTimestepsPerBlocco) {


    //description of output
    std::stringstream descr;

    descr << "#TODO"<<std::endl;
    c_descr=descr.str();

    ntimesteps=numeroTimestepsPerBlocco;

    //quanti timestep è lunga la funzione di correlazione
    leff =(numeroTimestepsPerBlocco<tmax || tmax==0)? numeroTimestepsPerBlocco : tmax;
    //numero di timestep su cui fare la media
    data_length=leff*ntypes*ntypes*nbin*(l+1);


    delete [] vdata;
    vdata=new TFLOAT [data_length];
}




template <int lmax, class TFLOAT, class T>
void SphericalCorrelations<lmax,TFLOAT,T>::calculate(unsigned int primo) {

    if (leff+ntimesteps+primo > t.get_ntimesteps()){
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the average or the lenght of the correlation function");
    }

    if (nthreads<=1){
        std::cerr << "Warning: I'm using a single thread.\n";
        nthreads=1;
    }

    if (buffer_size<2){
        std::cerr << "Warning: buffer size cannot be less than 2, setting it to 2"<<std::endl;
        buffer_size=2;
    }


    unsigned int npassith=leff/nthreads;
    std::vector<std::thread> threads;

    //wow! advanced work splitting procedure
    unsigned int block_t=leff/nthreads;
    if (block_t==0) {
        throw std::runtime_error("Don't have enough timesteps to parallelize the calculation. Try again with a smaller number of threads. " AT);
    }
    //buffer_size should be block+1 for maximum efficiency
    if (block_t>buffer_size){
        block_t=buffer_size-1;
    } else {
        buffer_size=block_t+1;
    }
    if (buffer_size<3) {
        buffer_size=3;
        std::cerr << "Warning: setting buffer_size to 3 (not optimal, because there are too few timesteps) " AT;
    }
    TwoLoopSplit<size_t> task_distributer(nthreads,ntimesteps,skip,skip*10,leff,1,block_t);

    //allocate space for per thread averages
    TFLOAT * lista_th = new TFLOAT[data_length*(nthreads-1)];
    unsigned int * lista_th_counters = new unsigned int [leff*nthreads];
    size_t * hit = new size_t[nthreads];
    size_t * miss = new size_t[nthreads];


    for (unsigned int  ith=0;ith<nthreads;ith++) {
        threads.push_back(std::thread([&,ith](){
            Neighbours_T * nns=nullptr;
            if (neighList.size()>0) {
                nns = new Neighbours_T{&t,neighList};
            }
            TFLOAT * lista_th_= ith >0 ? lista_th+data_length*(ith-1) : vdata;
            unsigned int * lista_th_counters_ = lista_th_counters+leff*ith;

            for (unsigned int i=0;i<data_length;++i) lista_th_[i]=0.0;
            for (unsigned int i=0;i<leff;++i) lista_th_counters_[i]=0;


            //working space for spherical harmonics
            TFLOAT workspace[(lmax+1)*(lmax+1)];
            TFLOAT cheby[2*(lmax+1)];
            //working space for averages of correlations over atomic types
            //for every atom I have to allocate space for the (l,m) matrix of coefficients (used to describe the environment around it)
            //then I calculate separately the (l,m) matrix for each type in the environment and for different values of the radial distance (bins)

            //this is for all atoms. I will compute correlations of objects that are such big
            int sh_snap_size=get_snap_size();
            //at the end I will first sum over all m: the (l,m) matrix will become an (l) vector
            //then I do the average over the central atoms of the same type: the size of the result (for each time lag) is:
            int sh_final_size=get_final_snap_size(); // remember: we sum over m before averaging on central atoms

            //allocate space
            TFLOAT *aveWork1=new TFLOAT[sh_snap_size];
            TFLOAT *aveTypes=new TFLOAT[sh_final_size];
            int *avecont=new int[ntypes];
            //buffer for few sh calculations -- note that only the last 2 requests must be in the returned pointer at the same time
            CalcBuffer<TFLOAT,unsigned int > buffer(buffer_size,sh_snap_size);


            //loop over data -- splitted in an efficient (?) way see class TwoLoopSplit. it is the equivalent of the following:
            /*
            for (unsigned int dt=npassith*ith;dt<ultimo;dt++)
                for (unsigned int imedia=0;imedia<ntimesteps;imedia+=skip)
                    */
            //center atom loop for the snapshot at imedia
            bool finished=false;
            size_t t1,t2,dt,t1_old=0;
            task_distributer.get_withoud_advancing(ith,t1_old,t2);
            while(!finished){
                task_distributer.get_next_idx_pair(ith,t1,t2,finished);
                if (t1 > t1_old) {
                    buffer.discard(t1_old+primo);
                } else if (t1<t1_old) {
                    buffer.discard();
                }
                t1_old=t1;
                TFLOAT * sh1=
                        buffer.buffer_calc(* static_cast<SPB*>(this),t1+primo,workspace,cheby,nullptr,nns);

                //center atom loop for the snapshot at imedia+dt
                TFLOAT * sh2=
                        buffer.buffer_calc(* static_cast<SPB*>(this),t2+primo,workspace,cheby,nullptr,nns);

                corr_sh_calc(sh1,sh2,aveTypes,aveWork1, sh_snap_size, sh_final_size, avecont);

                dt=t2-t1;
                lista_th_counters_[dt]++;
                //finally add to big average over starting timestep (imedia loop) -- I know, there are a lot of averages and sums and stuff
                for (int ll=0;ll<sh_final_size;++ll){
                    lista_th_[index(dt,0,0)+ll]+=(aveTypes[ll]-lista_th_[index(dt,0,0)+ll])/TFLOAT(lista_th_counters_[dt]);
                }



            }

            delete [] aveWork1;
            delete [] aveTypes;
            delete [] avecont;
            delete nns;
            hit[ith]=buffer.get_hit();
            miss[ith]=buffer.get_miss();

        }));
    }
    //Wasn't it multithreaded?
    for (unsigned int  ith=0;ith<nthreads;ith++){
        threads[ith].join();
        if (ith>0) {
            hit[0]+=hit[ith];
            miss[0]+=miss[ith];
        }
    }
    threads.clear();
    std::cerr << "Buffer miss/hit: "<<miss[0]<<"/"<<hit[0]<<std::endl;
    delete [] hit;
    delete [] miss;
    //sum up threads averages
    //the data of the first thread is in place
    for (unsigned int dt=0;dt<leff;++dt){
        for (unsigned int ith=1;ith<nthreads;++ith){
            TFLOAT * lista_th_=lista_th+data_length*(ith-1);
            unsigned int * lista_th_counters_ = lista_th_counters+leff*ith;
            if (lista_th_counters_[dt]==0) continue;
            //sum everything in vdata, use lista_th_counters[0] as counter for that accumulator
            if (lista_th_counters[dt]==0) { // just copy
                for (unsigned int i=0;i<get_final_snap_size();++i){
                    vdata[index(dt,0,0)+i]=lista_th_[index(dt,0,0)+i];
                }
            } else if (fabs(double(lista_th_counters[dt])/double(lista_th_counters_[dt])-1.0)<0.01) { // use traditional algorithm to update the mean
                for (unsigned int i=0;i<get_final_snap_size();++i){
                    vdata[index(dt,0,0)+i]=(vdata[index(dt,0,0)+i]*lista_th_counters[dt]+ lista_th_[index(dt,0,0)+i]*lista_th_counters_[dt])/double(lista_th_counters_[dt]+lista_th_counters[dt]);
                }
            } else { // use delta algorithm
                for (unsigned int i=0;i<get_final_snap_size();++i){
                    vdata[index(dt,0,0)+i] += (lista_th_[index(dt,0,0)+i]-vdata[index(dt,0,0)+i])*lista_th_counters_[dt]/double(lista_th_counters_[dt]+lista_th_counters[dt]);
                }
            }
            lista_th_counters[dt]+=lista_th_counters_[dt];
        }
    }


    delete[] lista_th;
    delete[] lista_th_counters;


}


template <int lmax, class TFLOAT, class T>
void SphericalCorrelations<lmax,TFLOAT,T>::corr_sh_calc(const TFLOAT * sh1, const TFLOAT *sh2, TFLOAT * aveTypes, TFLOAT * aveWork1, int sh_snap_size, int sh_final_size, int *avecont ) const noexcept{

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
    }
}

#ifdef BUILD_MMAP
#include "trajectory.h"
template class SphericalCorrelations<10,double,Trajectory>;
template class SphericalCorrelations<6,double,Trajectory>;
#endif

#ifdef PYTHON_SUPPORT
#include "trajectory_numpy.h"
template class SphericalCorrelations<10,double,Trajectory_numpy>;
template class SphericalCorrelations<6,double,Trajectory_numpy>;
#endif
