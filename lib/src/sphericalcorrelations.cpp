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
                                                         unsigned int nbin,
                                                         unsigned int tmax,
                                                         unsigned int nthreads,
                                                         unsigned int skip,
                                                         unsigned int buffer_size,
                                                         bool debug, const NeighListSpec neighList) :
t{*t},rminmax{rminmax},nbin{nbin}, skip{skip}, tmax{tmax}, nthreads{nthreads}, debug{debug},buffer_size{buffer_size}, neighList{neighList} {

   dr.reserve(rminmax.size());
   for (auto & range : rminmax) {
     dr.push_back((range.second-range.first)/nbin);
   }

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

    ntimesteps=numeroTimestepsPerBlocco;

    check_rminmax_size();

    //quanti timestep è lunga la funzione di correlazione
    leff =(numeroTimestepsPerBlocco<tmax || tmax==0)? numeroTimestepsPerBlocco : tmax;
    //numero di timestep su cui fare la media
    lunghezza_lista=leff*ntypes*ntypes*nbin*(l+1);


    delete [] lista;
    lista=new TFLOAT [lunghezza_lista];
}



template <int lmax, class TFLOAT, class T>
void SphericalCorrelations<lmax,TFLOAT,T>::calc(int timestep,
                                                TFLOAT *result,
                                                TFLOAT *workspace,//workspace array
                                                TFLOAT * cheby,//workspace array,
                                                int * counter // counter of #of atoms in the bin, for every atom
                                                ) const {
    //zero result
    for (int i=0;i<(lmax+1)*(lmax+1)*natoms*nbin*ntypes;++i) {
        result[i]=0;
    }
    if (counter){
        for (int i=0;i<natoms*ntypes*nbin;++i){
            counter[i]=0;
        }
    }

    if (neighList.size()==0){
        for (unsigned int iatom=0;iatom<natoms;iatom++) {
            //other atom loop
            const unsigned int itype=t.get_type(iatom);
            for (unsigned int jatom=0;jatom<natoms;jatom++) {
                const unsigned int jtype=t.get_type(jatom);
                const unsigned int pairidx=ntypes*itype+jtype;
                if (iatom==jatom)
                    continue;

                //minimum image distance
                double x[3];
                double d=sqrt(t.d2_minImage(iatom,jatom,timestep,timestep,x));
                //bin index
                int idx=(int)floorf((d-rminmax[pairidx].first)/dr[pairidx]);

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
                    if (counter){
                        counter[index_wrk_counter(iatom,jtype,idx)]++;
                    }
                }
            }
        }
    } else {
        auto nnl = Neighbours<T,double>{&t,neighList};
        nnl.update_neigh(timestep,true); //sorted list of neighbours, to use sann algorithm
        for (unsigned int iatom=0;iatom<natoms;iatom++) {
            for (unsigned int jtype=0;jtype<ntypes;jtype++) {
                //other atom loop
                auto position_iterator = nnl.get_sann_r(iatom,jtype);
                const unsigned int itype=t.get_type(iatom);
                const unsigned int pairidx=ntypes*itype+jtype;
                for (const auto & x : position_iterator){

                    //minimum image distance is already present in the neighbour list:
                    // x is: {r, r_x, r_y, r_z}

                    //calculate sin and cos
                    //calculate spherical harmonics and add to the correct average
                    SpecialFunctions::SphericalHarmonics<lmax,TFLOAT,true,true> sh(x[1],x[2],x[3],cheby,workspace);
                    sh.calc();
                    //now in workspace you have all the spherical harmonics components of the density of the current jatom around iatom
                    //add to the sh density of the current iatom of the current bin of the type jtype
                    for (int ll=0;ll<(lmax+1)*(lmax+1);++ll) {
                        result[index_wrk(iatom,jtype,0)+ll]+=workspace[ll];
                    }
                    if (counter){
                        counter[index_wrk_counter(iatom,jtype,0)]++;
                    }

                }
            }
        }
    }
}

template <int lmax, class TFLOAT, class T>
void SphericalCorrelations<lmax,TFLOAT,T>::calcola(unsigned int primo) {

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
    TwoLoopSplit<unsigned int> task_distributer(nthreads,ntimesteps,skip,skip*10,leff,1,block_t);

    //allocate space for per thread averages
    TFLOAT * lista_th = new TFLOAT[lunghezza_lista*(nthreads-1)];
    unsigned int * lista_th_counters = new unsigned int [leff*nthreads];
    size_t * hit = new size_t[nthreads];
    size_t * miss = new size_t[nthreads];


    for (unsigned int  ith=0;ith<nthreads;ith++) {
        threads.push_back(std::thread([&,ith](){
            TFLOAT * lista_th_= ith >0 ? lista_th+lunghezza_lista*(ith-1) : lista;
            unsigned int * lista_th_counters_ = lista_th_counters+leff*ith;

            for (unsigned int i=0;i<lunghezza_lista;++i) lista_th_[i]=0.0;
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
            unsigned int t1,t2,dt,t1_old=0;
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
                        buffer.buffer_calc(*this,t1+primo,workspace,cheby);

                //center atom loop for the snapshot at imedia+dt
                TFLOAT * sh2=
                        buffer.buffer_calc(*this,t2+primo,workspace,cheby);

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
            TFLOAT * lista_th_=lista_th+lunghezza_lista*(ith-1);
            unsigned int * lista_th_counters_ = lista_th_counters+leff*ith;
            if (lista_th_counters_[dt]==0) continue;
            //sum everything in lista, use lista_th_counters[0] as counter for that accumulator
            if (lista_th_counters[dt]==0) { // just copy
                for (unsigned int i=0;i<get_final_snap_size();++i){
                    lista[index(dt,0,0)+i]=lista_th_[index(dt,0,0)+i];
                }
            } else if (fabs(double(lista_th_counters[dt])/double(lista_th_counters_[dt])-1.0)<0.01) { // use traditional algorithm to update the mean
                for (unsigned int i=0;i<get_final_snap_size();++i){
                    lista[index(dt,0,0)+i]=(lista[index(dt,0,0)+i]*lista_th_counters[dt]+ lista_th_[index(dt,0,0)+i]*lista_th_counters_[dt])/double(lista_th_counters_[dt]+lista_th_counters[dt]);
                }
            } else { // use delta algorithm
                for (unsigned int i=0;i<get_final_snap_size();++i){
                    lista[index(dt,0,0)+i] += (lista_th_[index(dt,0,0)+i]-lista[index(dt,0,0)+i])*lista_th_counters_[dt]/double(lista_th_counters_[dt]+lista_th_counters[dt]);
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
template class SphericalCorrelations<10,double,Traiettoria>;
template class SphericalCorrelations<6,double,Traiettoria>;
#endif

#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class SphericalCorrelations<10,double,Traiettoria_numpy>;
template class SphericalCorrelations<6,double,Traiettoria_numpy>;
#endif
