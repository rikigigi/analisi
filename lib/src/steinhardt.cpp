#include "steinhardt.h"
#include "specialfunctions.h"

#include "config.h"

template <int l, class TFLOAT, class T>
Steinhardt<l,TFLOAT,T>::Steinhardt(T *t,
                                   const rminmax_t rminmax,
                                   unsigned int nbin,
                                   unsigned int nbin_steinhardt,
                                   std::vector<unsigned int> steinhardt_l_histogram,
                                   unsigned int nthreads,
                                   unsigned int skip,
                                   bool debug) :
    SPHC::SphericalCorrelations(t,rminmax,nbin,1,nthreads,skip,2),
    CMT::CalcolaMultiThread(nthreads,skip),
    steinhardt_histogram_size{0},
    steinhardt_l_histogram{steinhardt_l_histogram},
    nbin_steinhardt{nbin_steinhardt}
{

    steinhardt_histogram_size=1;
    for (unsigned i=0;i<steinhardt_l_histogram.size();++i) {
        steinhardt_histogram_size*=nbin_steinhardt;
    }
    lunghezza_lista=0;

}

template <int l, class TFLOAT, class T>
void Steinhardt<l,TFLOAT,T>::reset(const unsigned int numeroTimestepsPerBlocco) {


    SPHC::check_rminmax_size();
    compute_stride();
    //description of output
    std::stringstream descr;

    descr << "# first columns are indexes of the histogram; then:"<<std::endl;
    int cur_col=steinhardt_l_histogram.size();
    for (int ih=0;ih<nbin;++ih){
        descr << "#RADIAL BIN "<<ih<<std::endl;
        for (int i1=0;i1<ntypes;++i1){
            for (int i2=0;i2<ntypes;++i2){
                cur_col +=1;
                descr << "# types("<<i1<<","<<i2<<"): " << cur_col<<std::endl;
                cur_col +=1; //variance
            }
        }
    }
    c_descr=descr.str();

    ntimesteps=numeroTimestepsPerBlocco;
    CMT::ntimesteps=ntimesteps; //TODO: mess here

    size_t new_lungh_lista=steinhardt_histogram_size*nbin*ntypes*ntypes;
    if (new_lungh_lista != lunghezza_lista) {
        delete [] lista;
        lunghezza_lista=new_lungh_lista;
        lista = new TFLOAT [lunghezza_lista];
    }
}

template <int l, class TFLOAT, class T>
void Steinhardt<l,TFLOAT,T>::calc_init(int primo) {
    //check that everything is ok

    //allocate the space for the threads work

    threadResults = new TFLOAT [lunghezza_lista*(nthreads-1)];

    if (ntimesteps/CMT::skip>0) incr=1.0/int(ntimesteps/CMT::skip);
    else                   incr=1;
}
template <int l, class TFLOAT, class T>
void Steinhardt<l,TFLOAT,T>::calc_single_th(int istart,//average index, begin
                                            int istop,//average index, end
                                            int primo,//first index of the block
                                            int ith//thread index
                                            ) {
    //working space for spherical harmonics
    TFLOAT workspace[(l+1)*(l+1)];
    TFLOAT cheby[2*(l+1)];
    TFLOAT *result = new TFLOAT [(l+1)*(l+1)*natoms*nbin*ntypes];
    TFLOAT * threadResult = threadResults + lunghezza_lista*ith;
    int * counter = new int[natoms*ntypes*nbin];
    if (ith==nthreads-1) { // last thread writes directly on the final result memory
        threadResult = lista;
    }

    for (int i = istart; i<istop;++i) {
        SPHC::calc(i,result,workspace,cheby,counter);
        //calculate square, sum all m
        const int sh_size=SPHC::get_snap_size();
        for (int is=0;is<sh_size;++is) {
            result[is]=result[is]*result[is]; //square everything (SH are real here)
        }
        //sum over m for each l:
        using MultiVal = SpecialFunctions::MultiValDynamic::MultiVal<l,TFLOAT>;
        MultiVal v_atomic;
        for (int iatom=0;iatom<natoms;++iatom) {
            int itype=t.get_type(iatom);
            for (int jtype=0;jtype<ntypes;++jtype){
                for (int ibin=0;ibin<nbin;++ibin){
                    v_atomic.init(result+SPHC::index_wrk(iatom,jtype,ibin));
                    v_atomic.sum_in_m_zero();

                    //calculate order parameter and histogram index
                    size_t hist_idx=0;
                    TFLOAT n_atoms=counter[SPHC::index_wrk_counter(iatom,jtype,ibin)];
                    if (n_atoms>0) {
                        for (int jidx=0;jidx<steinhardt_l_histogram.size();++jidx) {
                            int j=steinhardt_l_histogram[jidx];
                            TFLOAT pre_sqrt=v_atomic.get_l_m0(j)/n_atoms/n_atoms*4*PI/(2*j+1);
                            TFLOAT m_steinhardt = sqrt(pre_sqrt);
                            size_t jhidx=floorf(m_steinhardt*nbin_steinhardt);
                            if (jhidx>=nbin_steinhardt) {
                                jhidx=nbin_steinhardt-1;
                            }
                            hist_idx+=jhidx*stride[3+jidx];
                        }
                        //update histogram (threadResult) of this radial bin...
                        size_t idx_all=get_index(ibin,itype,jtype)+hist_idx;
                        threadResult[idx_all]+=incr;
                    }
                }
            }
        }
    }
    delete [] counter;
    delete [] result;
}

template <int l, class TFLOAT, class T>
void Steinhardt<l,TFLOAT,T>::calc_end() {
    //sum all the threads result
    for (int ith=0;ith<nthreads-1;ith++) {
        for (int i=0;i<lunghezza_lista;++i) {
            lista[i]+=threadResults[ith*lunghezza_lista+i];
        }
    }
    delete [] threadResults;
    threadResults=nullptr;
}

#ifdef BUILD_MMAP
template class Steinhardt<6,double,Traiettoria>;
#endif
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class Steinhardt<6,double,Traiettoria_numpy>;
#endif
