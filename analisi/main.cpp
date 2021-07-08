/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/


#pragma STDC FENV_ACCESS ON

#include <iostream>
#include <array>
#include <fstream>
#include <boost/program_options.hpp>
#include "testtraiettoria.h"
#include "spettrovibrazionale.h"
#include "modivibrazionali.h"
#include "boost/version.hpp"
#include <boost/lexical_cast.hpp>
#include "config.h"
#include "istogrammavelocita.h"
#include "greenkubo2componentionicfluid.h"
#include "greenkuboNcomponentionicfluid.h"
#include "convolution.h"
#include "heatfluxts.h"
#include "msd.h"
#include "istogrammaatomiraggio.h"
#include "gofrt.h"
#include <functional>
#include "correlatorespaziale.h"
#include "sphericalcorrelations.h"
#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif
#include "convertibinario.h"
#ifdef DEBUG
#include "readlog.h"
#endif
#ifndef NDEBUG
#include <fenv.h>
#endif
#ifdef USE_MPI
#include "mp.h"
#endif
#include <cstdlib>

namespace std{

template<typename A, typename B>
std::ostream & operator << (std::ostream & out, const std::pair<A,B> & p) {
    out << p.first << " "<<p.second;
    return out;
}


/*
template <typename A,typename B>
std::istream & operator >> (std::istream & in, std::pair<A,B> & p) {
    in >> p.first >> p.second;
    return in;
}
*/

template<typename A,std::size_t N>
std::ostream & operator << (std::ostream &out, const std::array<A,N> &c) {
    for (unsigned int i=0;i<N;++i) {
        out << c[i] << " ";
    }
    return out;
}

template<typename A,std::size_t N >
std::ostream & operator << (std::ostream &out, const boost::array<A,N> &c) {
    for (unsigned int i=0;i<N;++i) {
        out << c[i] << " ";
    }
    return out;
}

/*
template<typename A, std::size_t N>
std::istream & operator >> (std::istream &in, std::array<A,N> &c) {
    for (unsigned int i=0;i<N;++i) {
        in >> c[i];
    }
    return in;
}
*/
template<typename A, std::size_t N>
void validate (boost::any& v,
               const std::vector<std::string>& values,
               std::array<A,N> *,int) {
    if (values.size()!=N) {
        throw std::runtime_error("std::array<A,N> needs N arguments");
    }
    std::array<A,N> res;
    for (unsigned int i=0;i<N;++i) {
        res[i]=boost::lexical_cast<A>(values[i]);
    }
    v=res;
}

template<typename A, typename B>
void validate (boost::any& v,
               const std::vector<std::string>& values,
               std::pair<A,B>*, int) {
    if (values.size()!=2){
        throw std::runtime_error("pair<A,B> needs two arguments");
    }
    std::pair<A,B> res;
    res.first = boost::lexical_cast<A>(values[0]);
    res.second = boost::lexical_cast<B>(values[1]);
    v=res;
}
}

int main(int argc, char ** argv)
{

#ifdef USE_MPI
    Mp::mpi(&argc,&argv);
#endif

#ifndef NDEBUG
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif

    std::cerr << _info_msg<<  std::endl;

    std::cerr << "arguments (enclosed by '') were:";
    for (int i=0;i<argc;++i) {
        std::cerr <<" '"<< argv[i]<<"'";
    }
    std::cerr << std::endl;

    //try to read OMP_NUM_THREADS env variable
    int omp_num_threads = 2;
    const char * omp_num_threads_c = getenv("OMP_NUM_THREADS");
    if (omp_num_threads_c != nullptr) {
        int nth_try = atoi(omp_num_threads_c);
        if (nth_try > 0) {
            omp_num_threads = nth_try;
            std::cerr << "Setting default number of threads from OMP_NUM_THREADS="<<omp_num_threads<<std::endl;
        } else {
            std::cerr << "Invalid OMP_NUM_THREADS='"<<omp_num_threads_c<<"' env variable, using a default value of "<<omp_num_threads<<std::endl;
        }
    } else {
        std::cerr << "OMP_NUM_THREADS env variable not found, using a default of "<<omp_num_threads<<std::endl;
    }



    boost::program_options::options_description options("Program to analyze of molecular dynamics trajectories, with multithread and MPI block averages.\n\nAllowed options:");
    std::string input,log_input,corr_out,ris_append_out,ifcfile,fononefile,output_conversion;
    int sub_mean_start=0,numero_frame=0,blocksize=0,elast=0,blocknumber=0,numero_thread,nbins_vel,skip=1,conv_n=20,final=60,stop_acf=0;
    unsigned int n_seg=0,gofrt=0,read_lines_thread=200,sph=0,buffer_size=10;
    bool sub_mean=false,test=false,spettro_vibraz=false,velocity_h=false,heat_coeff=false,debug=false,debug2=false,dumpGK=false,msd=false,msd_cm=false,msd_self=false,bench=false,fpe=false;
    double vmax_h=0,cariche[2],dt=5e-3,vicini_r=0.0;
    std::pair<int,double> nk;
    std::array<double,3> kdir;
    std::vector<unsigned int > cvar_list,kk_l;
    std::vector<double> factors_input;
    std::vector<std::string> headers,output_conversion_gro;
    std::vector< std::pair <unsigned int,unsigned int > > cvar;

    options.add_options()
            ("input,i",boost::program_options::value<std::string>(&input)->default_value(""), "input file in binary LAMMPS format: id type xu yu zu vx vy vz")
            ("loginput,l",boost::program_options::value<std::string>(&log_input),"column formatted file with headers. At the beginning of the file you can have free informations. When the program finds a line with only numbers, it assumes that data starts here, and that the line before contains the headers of the columns.")
            ("help,h", "help message")
            ("thread,N",boost::program_options::value<int>(&numero_thread)->default_value(omp_num_threads),"number of threads to use where supported (defaults to OMP_NUM_THREADS if present)")
            ("blocknumber,B",boost::program_options::value<int>(&blocknumber)->default_value(20),"number of blocks to use to calculate averages and variances and to split the reading of the trajectory")
#ifdef EXPERIMENTAL
            ("ifc,I",boost::program_options::value<std::string>(&ifcfile)->default_value(""),"file with the force constant matrix for phonon coordinate trasformation")
            ("energia-fononi,E",boost::program_options::value<std::string>(&fononefile)->default_value(""),"output file with kinetic and potential energy for each phonon")
#endif
            ("vibrational-spectrum,V",boost::program_options::bool_switch(&spettro_vibraz)->default_value(false),"calculates vibrational spectrum for each atomic type")
            ("velocity-histogram,v",boost::program_options::value<int>(&nbins_vel)->default_value(0),"calculates velocity histogram with the specified number of bins")
            ("histogram-minmax,M",boost::program_options::value<double>(&vmax_h)->default_value(500),"specify a number M. The bins of the histogram will be between -M and +M")
            ("heat-transport-coefficient,H",boost::program_options::bool_switch(&heat_coeff)->default_value(false),"perform green-kubo integral of the provided time series file (specified with the -l option)")
            ("headers,a",boost::program_options::value<std::vector<std::string> > (&headers)->multitoken(),"specify the name of the columns to use in the calculation of the multicomponent green-kubo integrals")
            ("dt,D",boost::program_options::value<double>(&dt)->default_value(5e-3),"timestep (used only in a particular case)")
            ("skip,s",boost::program_options::value<int>(&skip)->default_value(1),"when an average over the trajectory is performed, consecutive steps have a distance specified with this option")
#ifdef EXPERIMENTAL
            ("charge1",boost::program_options::value<double>(&cariche[0])->default_value(1.0),"charge of type 1 (used in a particular case)")
            ("charge2",boost::program_options::value<double>(&cariche[1])->default_value(-1.0),"charge of type 2 (used in a particular case)")
            ("conv_n,C",boost::program_options::value<int>(&conv_n)->default_value(10),"sigma of the gaussian that may be used to compute a convolution with the green-kubo integrals, to smooth out some noise. Number of points units.")
            ("final,f",boost::program_options::value<int>(&final)->default_value(60),"number of points to use to extract the final value of gk integral (used in a particular case)")
#endif
            ("dump-block,d",boost::program_options::bool_switch(&dumpGK)->default_value(false),"dump the data of each block on a file (where implemented)")
            ("stop,S",boost::program_options::value<int>(&stop_acf)->default_value(0),"maximum number of timestep of the function that the particular calculation is going to calculate. If 0, uses the maximum number of timesteps possible.")
            ("covariance,z",boost::program_options::value<std::vector<unsigned int > >(&cvar_list)->multitoken(),"in the block analysis, perform covariance estimation between the pair of column specified here. You must provide an even number of integers. Only for green-kubo calculation")
            ("mean-square-displacement,q",boost::program_options::bool_switch(&msd)->default_value(false),"compute atomic mean square displacement")
            ("mean-square-displacement-cm,Q",boost::program_options::bool_switch(&msd_cm)->default_value(false),"compute atomic mean square displacement and center of mass square displacement for each atomic type")
            ("mean-square-displacement-self",boost::program_options::bool_switch(&msd_self)->default_value(false),"the atomic mean square displacement is calculated in the reference system of the center of mass of the atomic system of the type of the considered atom")
            ("factors,F",boost::program_options::value<std::vector<double> >(&factors_input)->multitoken(),"used to set options particular to the selected calculation. For green-kubo, it select the factors, in order: 1) to be used for both the autocorrelation function and for its integral, and  2) to multiplies a second time only the autocorrelation function integral, this last number for example can be the integration factor. If a single factor is present and a green-kubo analysis is selected, triggers a slightly different calculation (see source code). For the g(r) calculation, specifies the full interval of distances where the istogram is computed")
            ("subtract-mean",boost::program_options::bool_switch(&sub_mean)->default_value(false),"substract from the autocorrelation function a number calculated by performing a mean starting from the timestep index specified with the -u option")
            ("subtract-mean-start,u",boost::program_options::value<int>(&sub_mean_start)->default_value(0),"if --subtract-mean is selected, specifies the starting timestep index to perform the average of the autocorrelation function")
            ("subBlock,k",boost::program_options::value<unsigned int>(&n_seg)->default_value(1),"optimization option for the green-kubo calculation. It is the number of sub-blocks when calculating averages inside the big block used to calculate the variance. The performance depends on the particular system used.")
            ("kk",boost::program_options::bool_switch(&bench)->default_value(false),"small benchmark to find a good value of -k")
            ("kk-range",boost::program_options::value<std::vector<unsigned int > >(&kk_l)->multitoken(),"range where the -k benchmark has to be performed")
            ("binary-convert",boost::program_options::value<std::string>(&output_conversion),"perform the conversion from txt to binary LAMMPS trajectory file. The name of the output binary file is specified here")
#ifdef XDR_FILE	    
            ("binary-convert-gromacs",boost::program_options::value<std::vector<std::string>>(&output_conversion_gro)->multitoken(),"\
	     perform the conversion from gromacs format to the LAMMPS binary. Here you have to specify the output LAMMPS binary and the input type file with format:\n id0 type0\n ...\nidN typeN\nwith atomic types in the same order of the trr gromacs file. The trr gromacs file is specified with -i.")
#endif
            ("neighbor",boost::program_options::value<double>(&vicini_r)->default_value(0.0),"calculate the histogram of the neighbours up to the specified distance")
            ("gofrt,g",boost::program_options::value<unsigned int>(&gofrt)->default_value(0),"calculate the distinctive and non distinctive part of the van Hove correlation function -- a dynamical structure factor. Here you set the size in number of bins of every calculated histogram. Note that it calculate an histogram for each value of t. If you put -S 1, you get only the traditional g(r). You have to specify the minimum and maximum value of the distance in the histogram with -F option. Note that if you don't set a lower bound of 0 on the distance the non distinctive (the self) part will show up only after a many timesteps...")
            ("spherical-harmonics-correlation,Y",boost::program_options::value<unsigned int>(&sph)->default_value(0),"perform the calculation of the correlation function of the atomic density expanded in spherical harmonics. Note that this is a very heavy computation, first do a small test (for example with a very high -S or a low -s value). Here you have to specify the number of different bins of the considered radial distances, specified with -F. The code will calculate a correlation function for each bin.")
            ("buffer-size",boost::program_options::value<unsigned int>(&buffer_size)->default_value(30),"Buffer size for sh frame values. This is machine dependend and can heavily influence the performance of the code")
            ("lt",boost::program_options::value<unsigned int> (&read_lines_thread)->default_value(200),"parameter to read faster the time series column formatted txt file. It specifies the number of lines to read in one after the other for each thread")
#ifdef EXPERIMENTAL
            ("spatial-correlator,A",boost::program_options::value(&nk)->default_value({0,0.0})->multitoken(),"Numero di punti della griglia ...")
            ("spatial-correlator-dir",boost::program_options::value(&kdir)->default_value({1.0,0.0,0.0})->multitoken(),"Direzione di k")
#endif
            ("write-mass-currents",boost::program_options::bool_switch(&debug2)->default_value(false),"write in the output the columns (specified with -a) of the input log file ( specified with -l). Then eventually write an additional current calculated by summing all the velocities from the input trajectory file and multiplying each velocity by the corresponding number specified in the -a string as the following: '#custom_name N q1 ... qN' where # is mandatory, custom_name is a arbitrary string, and q1 ... qN are N floating point numbers that multiplies the velocity of the corresponding atomic type ")
            ("fpe", boost::program_options::bool_switch(&fpe)->default_value(false), "look for floating point exceptions, where implemented")
        #ifdef DEBUG
            ("test-debug",boost::program_options::bool_switch(&debug)->default_value(false),"small test of the trajectory reading routine")
        #endif
            ;

    //            ("nc,non-copressed-input",boost::program_options::bool_switch(&ncompress)->default_value(false),"legge il file di input come file non compresso con gz");

    boost::program_options::variables_map vm;

    try {
        boost::program_options::store(
                    boost::program_options::parse_command_line(argc,argv,options)
                    ,vm
                    );


        boost::program_options::notify(vm);

        if (argc<=1 || ( (output_conversion!="" || output_conversion_gro.size()>0 ) && input=="") ||vm.count("help")|| (vm.count("loginput")==0 && ( debug2 || heat_coeff ) ) || skip<=0 || stop_acf<0 || final<0 || (!sub_mean && (sub_mean_start!=0) ) || sub_mean_start<0 || !(kk_l.size()==0 || kk_l.size()==2)){
            std::cout << options << "\n";
            return 1;
        }

        if (cvar_list.size()%2 != 0) {
            std::cout << "Error: covariance indexes list must contain an even number of elements\n";
            std::cout << options << "\n";
            return 1;
        } else {
            for (unsigned int i=0;i<cvar_list.size()/2;i++) {
                cvar.push_back( std::pair<unsigned int, unsigned int> (cvar_list.at(i*2),cvar_list.at(i*2+1)));
            }
        }



    }

    catch (boost::program_options::error const &e) {
        std::cerr << e.what() << '\n';
        std::cerr << options;
        return 1;
    }


    try {

        if (output_conversion!="") {

            ConvertiBinario conv(input,output_conversion);
            return 0;
        }
        if (output_conversion_gro.size()>0) {
            if (output_conversion_gro.size() !=2) {
                std::cerr << "Error: you must specify 2 files for the gromacs conversion task!\n";
                std::cerr << options;
                return -1;
            } else {
                ConvertiBinario conv(input,output_conversion_gro[0],ConvertiBinario::gromax_trr,output_conversion_gro[1]);
                return 0;
            }
        }



#ifdef FFTW3_THREADS
        fftw_init_threads();
        fftw_plan_with_nthreads(numero_thread);
        std::cerr << "I'm using " << numero_thread << " threads\n";
#ifdef DEBUG
        if (debug) {
            TestTraiettoria ttest(input);
            return 0;
        }
#endif
#endif
#ifdef DEBUG
        if (debug2){

            ReadLog<> test(log_input,0,1,numero_thread,read_lines_thread,headers);
            Traiettoria * binary_traj=NULL;
            if (test.need_binary(headers)>0) {
                binary_traj=new Traiettoria(input);
                test.calc_currents(binary_traj,blocknumber);
            }
            for (unsigned int i=0;i<test.n_timestep();i++){
                std::cout << i << " "<<test.timestep(i)<< " ";
                for (unsigned int j=0;j<test.n_data();j++){
                    std::cout << test.line(i)[j] << " ";
                }
                std::cout << "\n";
            }
            delete binary_traj;

        } else
#endif // DEBUG
            if (heat_coeff) {
                std::cerr << "Green-Kubo heat transport coefficient calculation is beginning...\n";
                ReadLog<> test(log_input,0,1,numero_thread,read_lines_thread,headers);
                Traiettoria * binary_traj=NULL;
                //qui devo aggiungere la traiettoria binaria a ReadLog, qualora ReadLog ne constati la necessità
                if (test.need_binary(headers)>0) {
                    binary_traj=new Traiettoria(input);
                    test.calc_currents(binary_traj,blocknumber);
                }

                double factor_conv=1.0;
                double factor_conv2=1.0;
                double factor_intToCorr=1.0;

                //calcola velocemente la media a blocchi per la temperatura, se non sono specificati manualmente i fattori
                const double const_charge=1.6021765,const_boltzmann=1.38064852;

                unsigned int idx_lx=0;
                double media_=0.0;
                double var_=0.0;
                if (factors_input.size()==1){

                    std::pair<unsigned int,bool> res=test.get_index_of("Temp");
                    idx_lx=res.first;
                    if(!res.second){
                        std::cerr << "Non riesco a trovare la colonna 'Temp' nel file di log '"<<log_input<<"'\n";
                        abort();
                    }
                    unsigned int idx_T=res.first;
                    res=test.get_index_of("Lx");
                    if(!res.second){
                        std::cerr << "Non riesco a trovare la colonna 'Lx' (lato della cella cubica) nel file di log '"<<log_input<<"'\n";
                        abort();
                    }

                    unsigned int cont=0;
                    unsigned int block_size=test.n_timestep()/blocknumber;
                    for (unsigned int iblock=0;iblock<blocknumber;iblock++){
                        unsigned int cont2=0;
                        double media2_=0.0;
                        for (unsigned int i=block_size*iblock;i<block_size*(iblock+1);i+=skip){ // media sul blocco
                            double delta2= test.line(i)[idx_T] - media2_;
                            media2_ = media2_ + delta2/(++cont2);
                        }
                        double delta=media2_-media_;
                        media_=media_+delta/(++cont);
                        var_ = var_ + delta*(media2_ - media_);
                    }
                    var_=var_/(cont*(cont-1));
                    factor_conv=const_charge*const_charge*dt*1e7 / ((pow(test.line(0)[idx_lx],3) )*const_boltzmann*media_*media_);
                    factor_conv2=const_charge*const_charge*dt*1e7 / ((pow(test.line(0)[idx_lx],3) )*const_boltzmann*media_);
                    factor_intToCorr=1.0/(1e-12*dt);
                } else{
                    if (factors_input.size()>2){
                        std::cerr << "Errore: specificare 2 fattori, uno per le funzioni di correlazione e uno per il suo integrale.\n";
                        abort();
                    } else if (factors_input.size()==2){

                        factor_conv=factor_conv2=factors_input[0];
                        factor_intToCorr=factors_input[1];
                    } else {
                        factor_conv=1.0;
                        factor_intToCorr=1.0;
                    }
                    std::cout << "#Factors: "<<factor_conv<<", "<<factor_conv*factor_intToCorr<<"\n";

                }
                //calcola le costanti di conversione da METAL di LAMMPS al sistema internazionale del coefficiente di conducibilità




                if (headers.size()==0){

                    MediaBlocchiG<ReadLog<>,GreenKubo2ComponentIonicFluid,std::string,double*,unsigned int,bool,unsigned int,unsigned int>
                            greenK_c(&test,blocknumber);

                    MediaVarCovar<GreenKubo2ComponentIonicFluid> greenK(GreenKubo2ComponentIonicFluid::narr,cvar);

                    greenK_c.calcola_custom(&greenK,log_input,cariche,skip,dumpGK,stop_acf,numero_thread);

                    double factors[GreenKubo2ComponentIonicFluid::narr]={
                        factor_conv*factor_intToCorr, //Jee
                        factor_conv2*factor_intToCorr, //Jzz
                        factor_conv*factor_intToCorr, //Jez
                        factor_conv, //intee
                        factor_conv2, //intzz
                        factor_conv, //intez
                        factor_conv, //lambda
                        factor_conv*factor_intToCorr, //Jze
                        factor_conv, //intze
                        factor_conv, //int_ein_ee
                        factor_conv, //int_ein_ez
                        factor_conv, //int_ein_ze
                        factor_conv2, //int_ein_zz
                        factor_conv //lambda_einst
                    };
                    double *lambda_conv=new double[greenK.size()];
                    double *lambda_conv_var=new double[greenK.size()];

                    if (final>=greenK.size()) final=greenK.size()-1;
                    if (final <0) final=0;

                    Convolution<double> convoluzione( std::function<double (const double  &)> ([&conv_n](const double & x)->double{
                        return exp(-x*x/(2*conv_n*conv_n));
                    }),(conv_n*6+1),-3*conv_n,3*conv_n,3*conv_n);
                    convoluzione.calcola(greenK.media(6),lambda_conv,greenK.size(),1);
                    convoluzione.calcola(greenK.varianza(6),lambda_conv_var,greenK.size(),1);
                    if (factors_input.size()==0)
                        std::cout << "# T T1sigma  factor1 factor2\n#"
                              <<media_ << " " << sqrt(var_) << " "
                              << factor_conv<<" "<<factor_intToCorr<< "\n";
                    std::cout <<"#valore di kappa a "<<final<< " frame: "<<lambda_conv[final]*factor_conv << " "<< sqrt(lambda_conv_var[final])*factor_conv<<"\n";

                    std::cout << "#Jee,Jzz,Jez,Jintee,Jintzz,Jintez,lambda,jze,Jintze,einst_ee,einst_ez,einst_ze,einst_zz,lambda_einst,lambda_conv,lambda' [,covarianze indicate...]; ciascuno seguito dalla sua varianza\n";
                    for (unsigned int i=0;i<greenK.size();i++) {
                        for (unsigned int j=0;j<GreenKubo2ComponentIonicFluid::narr;j++) {
                            std::cout << greenK.media(j)[i]*factors[j] << " "
                                      << greenK.varianza(j)[i]*factors[j]*factors[j] << " ";
                        }

                        std::cout << lambda_conv[i]*factor_conv<< " "<<lambda_conv_var[i]*factor_conv*factor_conv<<" "
                                  << factor_conv*(greenK.media(3)[i]-pow(greenK.media(5)[i],2)/greenK.media(4)[i])<<" ";
                        for (unsigned int j=0;j<greenK.n_cvar();j++){
                            std::cout << greenK.covarianza(j)[i]*factors[cvar[j].first]*factors[cvar[j].second] << " ";
                        }
                        std::cout  << "\n";

                    }
                    delete [] lambda_conv;
                    delete [] lambda_conv_var;

                } else {

                    MediaBlocchiG<ReadLog<>,GreenKuboNComponentIonicFluid<ReadLog<> >,
                            std::string,
                            unsigned int,
                            std::vector<std::string>,
                            bool,
                            unsigned int,
                            unsigned int,
                            bool,
                            unsigned int,
                            unsigned int,
                            bool,
                            unsigned int,
                            unsigned int>
                            greenK_c(&test,blocknumber);
                    unsigned int narr=headers.size()*headers.size()*3+2;
                    MediaVarCovar<GreenKuboNComponentIonicFluid<ReadLog<> > > greenK(narr,cvar);

                    if(kk_l.size()==0){
                        kk_l.push_back(10);
                        kk_l.push_back(100);
                    }
                    greenK_c.calcola_custom(&greenK,
                                            log_input,
                                            skip,
                                            headers,
                                            dumpGK,
                                            stop_acf,
                                            numero_thread,
                                            sub_mean,
                                            sub_mean_start,
                                            n_seg,
                                            bench,
                                            kk_l[0],kk_l[1]);

                    double *factors= new double [narr];

                    for (unsigned int j=0;j<headers.size()*headers.size();j++)
                        factors[j]=factor_conv*factor_intToCorr;
                    for (unsigned int j=headers.size()*headers.size();j<headers.size()*headers.size()*3+2;j++)
                        factors[j]=factor_conv;

                    if (final>=greenK.size()) final=greenK.size()-1;
                    if (final <0) final=0;

                    std::cout << "# factor integral, factor correlation\n#"
                              << factor_conv << " " << factor_intToCorr <<  "\n";
                    std::cout << greenK_c.puntatoreCalcolo()->get_columns_description();

                    for (unsigned int i=0;i<greenK.size();i++) {
                        for (unsigned int j=0;j<narr;j++) {
                            std::cout << greenK.media(j)[i]*factors[j] << " "
                                      << greenK.varianza(j)[i]*factors[j]*factors[j] << " ";
                        }

                        for (unsigned int j=0;j<greenK.n_cvar();j++){
                            std::cout << greenK.covarianza(j)[i]*factors[cvar[j].first]*factors[cvar[j].second] << " ";
                        }
                        std::cout  << "\n";

                    }
                    delete [] factors;
                }
                delete binary_traj;

            }else if (msd || msd_cm || msd_self){
                std::cerr << "Mean square displacement calculation ";
                unsigned int f_cm=1;
                if (msd_cm) {
                     std::cerr << "of the center of mass and of the atoms is beginning...\n";
                     f_cm=2;
                }else{
                     std::cerr << " of the atoms is beginning...\n";
                     if (msd_self) {
                         std::cerr << "In the reference system of each atomic type center of mass...\n";
                     } else {
                         std::cerr << "In the cell coordinate system...\n";
                     }
                }

                Traiettoria test(input);
#define MSD_(fpe_)\
                {using MSD=MSD<Traiettoria,fpe_>;\
                MediaBlocchi<MSD,unsigned int,unsigned int,unsigned int,bool,bool,bool> Msd(&test,blocknumber);\
                Msd.calcola(skip,stop_acf,numero_thread,msd_cm,msd_self,dumpGK);\
                for (unsigned int i=0;i<Msd.media()->lunghezza()/test.get_ntypes()/f_cm;i++) {\
                    for (unsigned int j=0;j<test.get_ntypes()*f_cm;j++)\
                        std::cout << Msd.media()->elemento(i*test.get_ntypes()*f_cm+j) << " " <<\
                                     Msd.varianza()->elemento(i*test.get_ntypes()*f_cm+j) << " ";\
                    std::cout << "\n";\
                }}
                if (fpe){
                    MSD_(true) //with fpe
                } else {
                    MSD_(false)
                }

            }else if (gofrt>0) {
                if (factors_input.size()!=2){
                    throw std::runtime_error("You have to specify the distance range with the option -F.\n");
                }
                std::cerr << "Calculation of g(r,t) -- distinctive and non distinctive part of the van Hove function...\n";
                Traiettoria tr(input);
                tr.set_pbc_wrap(true); //è necessario impostare le pbc per far funzionare correttamente la distanza delle minime immagini

                MediaBlocchi<Gofrt<double,Traiettoria>,double,double,unsigned int,unsigned int,unsigned int, unsigned int,bool>
                        gofr(&tr,blocknumber);
                gofr.calcola(factors_input[0],factors_input[1],gofrt,stop_acf,numero_thread,skip,dumpGK);

                unsigned int ntyp=tr.get_ntypes()*(tr.get_ntypes()+1);
                unsigned int tmax=gofr.media()->lunghezza()/gofrt/ntyp;

                std::cout << gofr.puntatoreCalcolo()->get_columns_description();
                for (unsigned int t=0;t<tmax;t++) {
                    for (unsigned int r=0;r<gofrt;r++) {
                        std::cout << t << " " << r;
                        for (unsigned int itype=0;itype<ntyp;itype++) {
                            std::cout << " "<< gofr.media()->elemento(
                                             t*ntyp*gofrt+
                                             gofrt*itype+
                                             r)
                                      << " "<< gofr.varianza()->elemento(
                                             t*ntyp*gofrt+
                                             gofrt*itype+
                                             r);
                        }
                        std::cout << "\n";
                    }
                    std::cout << "\n\n";
                }


            }else if (sph>0) {
                if (factors_input.size()!=2){
                    throw std::runtime_error("You have to specify the distance range with the option -F.\n");
                }
                std::cerr << "Calculation of spherical harmonic density correlation function is beginning (this can take a lot of time)...\n";
                Traiettoria tr(input);
                tr.set_pbc_wrap(false); //è necessario impostare le pbc per far funzionare correttamente la distanza delle minime immagini

                MediaBlocchi<SphericalCorrelations<10,double,Traiettoria>,double,double,unsigned int,unsigned int,unsigned int, unsigned int,unsigned int,bool>
                        sh(&tr,blocknumber);
                sh.calcola(factors_input[0],factors_input[1],sph,stop_acf,numero_thread,skip,buffer_size,dumpGK);

                auto shape= sh.media()->get_shape();

                std::cout << sh.puntatoreCalcolo()->get_columns_description();
                auto line_size=shape[1]*shape[2]*shape[3]*shape[4];
                for (unsigned int t=0;t<shape[0];t++) {
                    for (unsigned int r=0;r<line_size;r++) {
                        std::cout << sh.media()->elemento(t*line_size+r) << " "<<
                                     sh.varianza()->elemento(t*line_size+r) << " ";
                    }
                    std::cout << std::endl;
                }


            }else if (vicini_r>0){
                std::cerr << "Beginning of calculation of neighbour histogram\n";
                Traiettoria test(input);
                test.set_pbc_wrap(true);

                IstogrammaAtomiRaggio h(&test,vicini_r,skip,numero_thread);
                unsigned int nt=test.get_ntimesteps(),
                        s=(nt-1)/blocknumber;
                h.reset(s);
                test.imposta_dimensione_finestra_accesso(s);
                for (unsigned int i=0;i<blocknumber;i++){
                    unsigned int t=s*i;
                    test.imposta_inizio_accesso(t);
                    h.calcola(t);
                }
                std::map<unsigned int, unsigned int> *hist=h.get_hist();
                for (unsigned int i=0;i<test.get_ntypes();i++){
                    std::cout << "\""<<i<<"\"\n";
                    for (auto it=hist[i].begin();it!=hist[i].end();++it){
                        std::cout << it->first << " " << it->second << "\n";
                    }
                    std::cout << "\n\n";
                }


            }else if (nbins_vel>0) {
                std::cerr << "Velocity histogram...\n";
                Traiettoria test(input);
                MediaBlocchi<IstogrammaVelocita,unsigned int,double> istogramma_vel(&test,blocknumber);
                istogramma_vel.calcola(nbins_vel,vmax_h);

                //stampa i risultati
                unsigned int lungh_sing_h=istogramma_vel.media()->lunghezza()/(3*test.get_ntypes());
                for (unsigned int i=0;i<lungh_sing_h;i++){
                    std::cout << -vmax_h + 2*vmax_h*i/lungh_sing_h << " ";
                    for (unsigned int it=0;it<3*test.get_ntypes();it++)
                        std::cout << istogramma_vel.media()->elemento(lungh_sing_h*it+i) << " "
                                  << istogramma_vel.varianza()->elemento(lungh_sing_h*it+i) << " ";
                    std::cout << "\n";
                }

            } else if (fononefile != "") {
                std::cerr << "Normal mode coordinate transformation...\n";
                Traiettoria test(input);
		blocksize = test.get_ntimesteps()/blocknumber;
                ModiVibrazionali test_normali(&test,ifcfile,fononefile,numero_thread,blocksize);
                test_normali.reset(stop_acf);
                test_normali.calcola(0);

                return 1;
            } else if (spettro_vibraz) {
                std::cerr << "Vibrational spectrum...\n";
                Traiettoria test(input);
                //            SpettroVibrazionale test_spettro(&test);
                MediaBlocchi<SpettroVibrazionale<Traiettoria>,bool> test_spettro_blocchi(&test,blocknumber);

                test_spettro_blocchi.calcola(dumpGK);
                for (unsigned int i=0;i<test_spettro_blocchi.media()->lunghezza()/(3*test.get_ntypes());i++) {
                    std::cout << i << " ";
                    for (unsigned int itipo=0;itipo<test.get_ntypes();itipo++)
                        std::cout << test_spettro_blocchi.media()->spettro(i,0,itipo) << " " << test_spettro_blocchi.varianza()->spettro(i,0,itipo) << " "
                                  << test_spettro_blocchi.media()->spettro(i,1,itipo) << " " << test_spettro_blocchi.varianza()->spettro(i,1,itipo) << " "
                                  << test_spettro_blocchi.media()->spettro(i,2,itipo) << " " << test_spettro_blocchi.varianza()->spettro(i,2,itipo) << " ";
                    std::cout << "\n";
                }

            } else if (nk.first>0) {
                std::cerr << "Spatial correlations of velocities...\n";
                //generate k-list
                std::vector<std::array<double,3> > klist;
                klist.reserve(nk.first);
                for (unsigned int i=0;i<nk.first;++i){
                    auto k=kdir;
                    for (auto & ki : k) {
                        ki=ki*nk.second*double(i)/double(nk.first);
                    }
                    klist.push_back(k);
                }
                Traiettoria t(input);
                t.set_pbc_wrap(true);

                CorrelatoreSpaziale corr(&t,klist,0.0,numero_thread,skip,false);
                unsigned int nt=t.get_ntimesteps(),
                        s=(nt-1)/blocknumber;
                corr.reset(s);
                t.imposta_dimensione_finestra_accesso(s+corr.numeroTimestepsOltreFineBlocco(blocknumber));
                //print stuff
                std::cout << "#kdir= "<<kdir[0]<<" "<<kdir[1]<<" "<<kdir[2]<<std::endl
                          << "#dk= "<<nk.second<<std::endl
                          << "#nk= "<<nk.first<<std::endl;
                for (unsigned int i=0;i<blocknumber;++i) {
                    unsigned int ts=s*i;
                    t.imposta_inizio_accesso(ts);
                    corr.calcola(ts);
                    corr.print(std::cout);
                }
            }
    }

    catch (const std::exception& e) {
        std::cerr << e.what()<<"\n";
        return 1;
    }




    return 0;
}
