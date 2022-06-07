#include "steinhardt.h"
#include "specialfunctions.h"

#include "config.h"

template <int l, class TFLOAT, class T>
Steinhardt<l,TFLOAT,T>::Steinhardt(T *t,
                                   const Rminmax_t rminmax,
                                   unsigned int nbin,
                                   unsigned int nbin_steinhardt,
                                   std::vector<unsigned int> steinhardt_l_histogram,
                                   unsigned int nthreads,
                                   unsigned int skip,
                                   bool do_histogram,
                                   const NeighListSpec nls,
                                   const bool averaged_order) :
    SPB::SphericalBase{t,nbin,rminmax},
    CMT::CalcolaMultiThread{nthreads,skip},
    steinhardt_histogram_size{0},
    steinhardt_l_histogram{steinhardt_l_histogram},
    nbin_steinhardt{nbin_steinhardt},
    natoms{t->get_natoms()},
    ntypes{static_cast<size_t>(t->get_ntypes())},
    nbin{nbin},
    neighListSpec{nls},
    t{*t},
    do_histogram{do_histogram},
    averaged_order{averaged_order}
{

    if (neighListSpec.size()==0 && averaged_order) {
        throw std::runtime_error("You must use a neighbour list when computing a averaged steinhardt order parameter! " AT);
    }
    if (neighListSpec.size()>0 && nbin >1) {
        throw std::runtime_error("When using a neighbour list to compute steinhardt please set nbin to 1. " AT);
    }

    steinhardt_histogram_size=1;
    for (unsigned i=0;i<steinhardt_l_histogram.size();++i) {
        steinhardt_histogram_size*=nbin_steinhardt;
    }
    lunghezza_lista=0;

}

template <int l, class TFLOAT, class T>
void Steinhardt<l,TFLOAT,T>::reset(const unsigned int numeroTimestepsPerBlocco) {

    ntimesteps=numeroTimestepsPerBlocco;
    compute_stride();
    //description of output
    std::stringstream descr;
    size_t new_lungh_lista;
    if (do_histogram){
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
        new_lungh_lista=steinhardt_histogram_size*nbin*ntypes*ntypes;
    } else {
        new_lungh_lista = stride[0]*nbin;
    }


    c_descr=descr.str();
    if (new_lungh_lista != lunghezza_lista) {
        delete [] lista;
        lunghezza_lista=new_lungh_lista;
        lista = new TFLOAT [lunghezza_lista];
    }
}

template <int l, class TFLOAT, class T>
void Steinhardt<l,TFLOAT,T>::calc_init(int primo) {
    //check that everything is ok
    if (ntimesteps+primo > t.get_ntimesteps()){
        throw std::runtime_error("trajectory is too short for this kind of calculation. Select a different starting timestep or lower the size of the average " AT);
    }

    //allocate the space for the threads work
    if (do_histogram){
        threadResults = new TFLOAT [lunghezza_lista*(CMT::nthreads-1)];
    } else {
        threadResults = nullptr;
    }

    if (ntimesteps/CMT::skip>0) incr=1.0/int(ntimesteps/CMT::skip);
    else                   incr=1;
}

template <int l, class TFLOAT, class T>
void Steinhardt<l,TFLOAT,T>::calc_single_th(int istart,//timestep index, begin
                                            int istop,//timestep index, end
                                            int primo,//first index of the block
                                            int ith//thread index
                                            ) {
    //working space for spherical harmonics
    TFLOAT workspace[(l+1)*(l+1)];
    TFLOAT cheby[2*(l+1)];
    TFLOAT *result = new TFLOAT [(l+1)*(l+1)*natoms*nbin*ntypes];
    TFLOAT *result_averaged = nullptr;
    if (averaged_order) result_averaged = new TFLOAT [(l+1)*(l+1)*natoms*nbin*ntypes];
    TFLOAT * threadResult = threadResults + lunghezza_lista*ith;
    int * counter = new int[natoms*ntypes*nbin];
    if (do_histogram){
        if (ith==CMT::nthreads-1) { // last thread writes directly on the final result memory
            threadResult = lista;
        }
        for (size_t i=0;i<lunghezza_lista;++i) {
            threadResult[i]=0;
        }
    } else {
        threadResult = lista; //every thread will write its own timesteps
    }

    Neighbours_T * nns = nullptr;
    if (neighListSpec.size()>0) {
        nns = new Neighbours_T{&t,neighListSpec};
    }

    for (int i = istart; i<istop;i+=CMT::skip) {
        calc(primo+i,result,workspace,cheby,counter,nns);
        //calculate square, sum all m
        const int sh_size=SPB::get_result_size();
        if (result_averaged) {
            //a new loop over neighbours, sum everything again
            for (int iatom=0;iatom<natoms;++iatom ) {
                int itype=t.get_type(iatom);
                for (int jtype=0;jtype<ntypes;++jtype){
                    for (int ibin=0;ibin<nbin;++ibin){

                        TFLOAT *to_sum=result+SPB::index_wrk(iatom,jtype,ibin);
                        TFLOAT *destination = result_averaged+SPB::index_wrk(iatom,jtype,ibin);
                        //size of stuff is (l+1)*(l+1)
                        //copy value of atom itself
                        if (itype==jtype){
                            for (int im=0;im<(l+1)*(l+1);++im){
                                destination[im]=to_sum[im]/counter[SPB::index_wrk_counter(iatom,jtype,ibin)];
                            }
                        }else {
                            for (int im=0;im<(l+1)*(l+1);++im){
                                destination[im]=0;
                            }
                        }
                        //do a loop over neighbours and sum
                        if (nns) {
                            auto index_iterator = nns->get_sann(iatom,jtype);
                            TFLOAT n_atoms_tot=index_iterator.size();
                            if (n_atoms_tot>0) {
                                for (auto nidx : index_iterator){
                                    TFLOAT n_atoms_neigh = counter[SPB::index_wrk_counter(nidx,jtype,ibin)];
                                    to_sum = result+SPB::index_wrk(nidx,jtype,ibin);
                                    for (int im=0;im<(l+1)*(l+1);++im){
                                        destination[im] += to_sum[im]/n_atoms_neigh; // add q_lm(nidx) inn \bar{q}_lm(iatom) accumulator
                                    }
                                }
                            }
                        } else {
                            throw std::runtime_error("Not implemented" AT);
                        }
                    }
                }
            }
        }
        //square everything (SH are real here)
        if (averaged_order) {
            for (int is=0;is<sh_size;++is) {
                result[is]=result_averaged[is]*result_averaged[is]; //use average q_lm. Only differente with not averaged
            }

        } else {
            for (int is=0;is<sh_size;++is) {
                result[is]=result[is]*result[is]; //square everything (SH are real here)
            }
        }


        using MultiVal = SpecialFunctions::MultiValDynamic::MultiVal<l,TFLOAT>;
        MultiVal v_atomic;
        //sum over m for each l:
        for (int iatom=0;iatom<natoms;++iatom) {
            int itype=t.get_type(iatom);
            for (int jtype=0;jtype<ntypes;++jtype){
                for (int ibin=0;ibin<nbin;++ibin){
                    v_atomic.init(result+SPB::index_wrk(iatom,jtype,ibin));
                    v_atomic.sum_in_m_zero();

                    if (do_histogram){
                        //calculate order parameter and histogram index
                        size_t hist_idx=0;
                        TFLOAT n_atoms=counter[SPB::index_wrk_counter(iatom,jtype,ibin)];
                        if (n_atoms>0) {
                            for (int jidx=0;jidx<steinhardt_l_histogram.size();++jidx) {
                                int j=steinhardt_l_histogram[jidx];
                                TFLOAT m_steinhardt=sqrt(v_atomic.get_l_m0(j)/n_atoms/n_atoms*4*PI/(2*j+1));
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
                    }else{
                        //output all steinhardt order parameters in a trajectory-like object
                        TFLOAT n_atoms=counter[SPB::index_wrk_counter(iatom,jtype,ibin)];
                        if (n_atoms>0){
                            for (size_t lidx=0;lidx<l;++lidx){
                                //avoid l=0 that is always 1
                                TFLOAT m_steinhardt=sqrt(v_atomic.get_l_m0(lidx+1)/n_atoms/n_atoms*4*PI/(2*(lidx+1)+1));
                                threadResult[get_index(ibin,i/CMT::skip,iatom,jtype,lidx)]=m_steinhardt;
                            }
                        }else {
                            for (size_t lidx=0;lidx<l;++lidx){
                                threadResult[get_index(ibin,i/CMT::skip,iatom,jtype,lidx)]=0;
                            }
                        }
                    }
                }
            }
        }
    }
    delete [] counter;
    delete [] result;
    delete [] result_averaged;
    delete  nns;
    std::cerr << "Thread " << ith << " istart , istop = "<<istart<<" , "<<istop<< " finished" <<std::endl;
}

template <int l, class TFLOAT, class T>
void Steinhardt<l,TFLOAT,T>::join_data() {
    if (do_histogram){
        //sum all the threads result
        for (int ith=0;ith<CMT::nthreads-1;ith++) {
            for (int i=0;i<lunghezza_lista;++i) {
                lista[i]+=threadResults[ith*lunghezza_lista+i];
            }
        }
    }
    delete [] threadResults;
    threadResults=nullptr;
}

#ifdef BUILD_MMAP
template class Steinhardt<6,double,Traiettoria>;
template class Steinhardt<8,double,Traiettoria>;
template class Steinhardt<10,double,Traiettoria>;
#endif
#ifdef PYTHON_SUPPORT
#include "traiettoria_numpy.h"
template class Steinhardt<6,double,Traiettoria_numpy>;
template class Steinhardt<8,double,Traiettoria_numpy>;
template class Steinhardt<10,double,Traiettoria_numpy>;
#endif
