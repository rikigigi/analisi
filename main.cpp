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



#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include "interfaccia.h"
#include "testtraiettoria.h"
#include "spettrovibrazionale.h"
#include "modivibrazionali.h"
#include "boost/version.hpp"
#include "config.h"
#include "istogrammavelocita.h"
#include "greenkubo2componentionicfluid.h"
#include "convolution.h"
#include "heatfluxts.h"
#include <functional>
#include <fftw3.h>


int main(int argc, char ** argv)
{
    boost::program_options::options_description options("Opzioni consentite");
    std::string input,log_input,corr_out,ris_append_out,ifcfile,fononefile;
    int numero_frame=0,blocksize=0,elast=0,blocknumber=0,numero_thread,nbins,skip=1,conv_n=20,final=60,stop_acf=0;
    bool test=false,spettro_vibraz=false,velocity_h=false,heat_coeff=false,debug=false,dumpGK=false;
    double vmax_h=0,cariche[2];
    options.add_options()
#if BOOST_VERSION >= 105600
            ("input,i",boost::program_options::value<std::string>(&input)->default_value(""), "file di input nel formato binario di LAMMPS: id type xu yu zu vx vy vz")
            ("loginput,l",boost::program_options::value<std::string>(&log_input)->required(),"file di input con il log di LAMMPS e i dati termodinamici nel formato: Step Time PotEng TotEng Lx Press Temp c_flusso[1] c_flusso[2] c_flusso[3]")
#else
            ("input,i",boost::program_options::value<std::string>(&input)->default_value(""), "file di input nel formato binario di LAMMPS: id type xu yu zu vx vy vz")
            ("loginput,l",boost::program_options::value<std::string>(&log_input),"file di input con il log di LAMMPS e i dati termodinamici nel formato: Step Time PotEng TotEng Lx Press Temp c_flusso[1] c_flusso[2] c_flusso[3]")
#endif
            ("help,h", "messaggio di aiuto")
            ("thread,N",boost::program_options::value<int>(&numero_thread)->default_value(2),"numero di thread da usare (dove supportato)")
            ("timesteps,n",boost::program_options::value<int>(&numero_frame)->default_value(50000), "numero di timestep che il programma legge dal file")
            ("blocksize,b",boost::program_options::value<int>(&blocksize)->default_value(1000),"numero di timestep per blocco (per calcolare l'errore con l'analisi a blocchi, o per quei conti che devono essere fatti pezzo per pezzo)")
            ("blocknumber,B",boost::program_options::value<int>(&blocknumber)->default_value(20),"numero di blocchi")
            ("errorblocklast,e",boost::program_options::value<int>(&elast)->default_value(1000),"numero di timestep finali di ogni blocco per calcolare l'errore sull'integrale della funzione di autocorrelazione (il programma fa una media della varianza dell'integrale punto per punto)")
            ("corr_out,c",boost::program_options::value<std::string>(&corr_out)->default_value(""),"file dove scrivere la funzione di autocorrelazione e il suo integrale nel formato: t acf acf_var int_acf int_acf_var")
            ("output,o",boost::program_options::value<std::string>(&ris_append_out)->default_value(""),"file dove scrivere i risultati finali (aggiungendoli alla fine")
            ("test,t",boost::program_options::bool_switch(&test)->default_value(false),"esegue le procedure di test del programma")
            ("ifc,I",boost::program_options::value<std::string>(&ifcfile)->default_value(""),"file della matrice con la derivata della forza nelle varie dimensioni spaziali (per l'analisi dei modi normali)")
            ("energia-fononi,E",boost::program_options::value<std::string>(&fononefile)->default_value(""),"nome del file da scrivere con dentro le energie cinetica e potenziale per singolo fonone, riga per riga al passare dei timestep")
            ("vibrational-spectrum,V",boost::program_options::bool_switch(&spettro_vibraz)->default_value(false),"esegue l'analisi a blocchi, secondo i parametri, dello spettro vibrazionale, per ogni tipo di atomo della simulazione. Scrive sullo stdout media e varianza per ciascuna delle tre direzioni spaziali.")
            ("velocity-histogram,v",boost::program_options::bool_switch(&velocity_h)->default_value(false),"calcola l'istogramma delle velocità con una media a blocchi per calcolare l'errore")
            ("bin-number,m",boost::program_options::value<int>(&nbins)->default_value(50),"numero di intervalli da usare nell'istogramma delle velocità")
            ("histogram-minmax,M",boost::program_options::value<double>(&vmax_h)->default_value(500),"l'istogramma viene calcolato nell'intervallo compreso fra -<valore> e +<valore>, dove <valore> è quello qui specificato")
            ("heat-transport-coefficient,H",boost::program_options::bool_switch(&heat_coeff)->default_value(false),"calcola il coefficiente di trasporto termico per sale fuso a due componenti")
            ("heat-transport-skip,s",boost::program_options::value<int>(&skip)->default_value(1),"incremento minimo di timesteps da usare nel calcolo delle funzioni di correlazione")
            ("charge1",boost::program_options::value<double>(&cariche[0])->default_value(1.0),"carica in unità elementari del tipo 1")
            ("charge2",boost::program_options::value<double>(&cariche[1])->default_value(-1.0),"carica in unità elementari del tipo 2")
            ("conv_n,C",boost::program_options::value<int>(&conv_n)->default_value(10),"sigma della gaussiana con cui viene fatta la convoluzione del coefficiente di trasporto termico al variare del tempo di integrazione (in numero di frame)")
            ("final,f",boost::program_options::value<int>(&final)->default_value(60),"numero di frame a cui estrarre il risultato finale")
            ("dump-block-H,d",boost::program_options::bool_switch(&dumpGK)->default_value(false),"scrive in dei file (aggiungendo ogni volta alla fine) i calcoli di ogni signolo blocco per il calcolo del coefficiente di trasporto termico in un sale fuso")
            ("stop-acf,S",boost::program_options::value<int>(&stop_acf)->default_value(0),"lunghezza massima, in frame, di tutte le funzioni di correlazione e relativi integrali. Se posto a zero è pari alle dimensioni del blocco")
#ifdef DEBUG
            ("test-debug",boost::program_options::bool_switch(&debug)->default_value(false),"test vari")
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

        if (vm.count("help")||vm.count("input")==0 || vm.count("loginput")==0 || skip<=0 || stop_acf<0 || final<0){
            std::cout << "COMPILED AT " __DATE__ " " __TIME__ " by " CMAKE_CXX_COMPILER " whith flags " CMAKE_CXX_FLAGS  " on a " CMAKE_SYSTEM " whith processor " CMAKE_SYSTEM_PROCESSOR ".\n";
            std::cout << options << "\n";
            return 1;
        }
    }

    catch (boost::program_options::error const &e) {
        std::cerr << e.what() << '\n';
        std::cout << "COMPILED AT " __DATE__ " " __TIME__ " by " CMAKE_CXX_COMPILER " whith flags " CMAKE_CXX_FLAGS  " on a " CMAKE_SYSTEM " whith processor " CMAKE_SYSTEM_PROCESSOR ".\n";
        std::cerr << options;
        return 1;
    }


    try {

#ifdef FFTW3_THREADS
        fftw_init_threads();
        fftw_plan_with_nthreads(numero_thread);
#ifdef DEBUG
        std::cerr << "fftw_plan_with_nthreads("<<numero_thread<<")\n";
        if (debug) {
            TestTraiettoria ttest(input);
            return 0;
        }
#endif
#endif

        if (heat_coeff) {
            std::cerr << "Inizio del calcolo del coefficiente di trasporto termico per un sale a due componenti...\n";
            Traiettoria test(input);
            test.imposta_dimensione_finestra_accesso(1);
            test.imposta_inizio_accesso(0);
            MediaBlocchi<GreenKubo2ComponentIonicFluid,std::string,double*,unsigned int,bool,unsigned int> greenK(&test,blocknumber,log_input,cariche,skip,dumpGK,stop_acf);
            greenK.calcola();
            greenK.puntatoreCalcolo()->puntatoreHeatFluxTs()->temp(0);
            //calcola velocemente la media a blocchi per la temperatura

            double media_=0.0;
            double var_=0.0;
            unsigned int cont=0;
            unsigned int block_size=test.get_ntimesteps()/skip/blocknumber;
            for (unsigned int iblock=0;iblock<blocknumber;iblock++){
                unsigned int cont2=0;
                double media2_=0.0;
                for (unsigned int i=block_size*iblock;i<block_size*(iblock+1);i++){ // media sul blocco
                    double delta2= greenK.puntatoreCalcolo()->puntatoreHeatFluxTs()->temp(0)[i] - media2_;
                    media2_ = media2_ + delta2/(++cont2);
                }
                double delta=media2_-media_;
                media_=media_+delta/(++cont);
                var_ = var_ + delta*(media2_ - media_);
            }
            var_=var_/(cont*(cont-1));

            double factor_conv=1.6022*1.6022*5*skip / ((pow(greenK.puntatoreCalcolo()->puntatoreHeatFluxTs()->get_L(),3) )*1.38064852e-4*media_*media_);
            double factor_conv2=1.6022*1.6022*5*skip / ((pow(greenK.puntatoreCalcolo()->puntatoreHeatFluxTs()->get_L(),3) )*1.38064852e-4*media_);
            double factor_intToCorr=1.0/(1e-15*5*skip);  //0.01 ps è l'intervallo di integrazione
            double factors[9]={
                factor_conv*factor_intToCorr, //Jee
                factor_conv2*factor_intToCorr, //Jzz
                factor_conv*factor_intToCorr, //Jez
                factor_conv, //intee
                factor_conv2, //intzz
                factor_conv, //intez
                factor_conv, //lambda
                factor_conv*factor_intToCorr, //Jze
                factor_conv //intze
            };
            double *lambda_conv=new double[greenK.media()->lunghezza()/9];
            double *lambda_conv_var=new double[greenK.media()->lunghezza()/9];

            if (final>=greenK.media()->lunghezza()/9) final=greenK.media()->lunghezza()/9-1;
            if (final <0) final=0;

            Convolution<double> convoluzione( std::function<double (const double  &)> ([&conv_n](const double & x)->double{
                return exp(-x*x/(2*conv_n*conv_n));
            }),(conv_n*6+1),-3*conv_n,3*conv_n,3*conv_n);
            convoluzione.calcola(&greenK.media()->accesso_lista()[6],lambda_conv,greenK.media()->lunghezza()/9,9);
            convoluzione.calcola(&greenK.varianza()->accesso_lista()[6],lambda_conv_var,greenK.media()->lunghezza()/9,9);
            std::cout << "# T T1sigma  atomi/volume densitNaCL\n#"
                      <<media_ << " " << sqrt(var_) << " "
                      << test.get_natoms()/pow(greenK.puntatoreCalcolo()->puntatoreHeatFluxTs()->get_L(),3)<< " "<<
                         test.get_natoms()/pow(greenK.puntatoreCalcolo()->puntatoreHeatFluxTs()->get_L(),3)*(22.990+35.453)/2.0*1.66054<<"\n"
                      <<"#valore di kappa a "<<final<< " frame: "<<lambda_conv[final]*factor_conv << " "<< sqrt(lambda_conv_var[final])*factor_conv<<"\n";

            std::cout << "#Jee,Jzz,Jez,Jintee,Jintzz,Jintez,lambda,jze,Jintze,lambda_conv; ciascuno seguito dalla sua varianza\n";
            for (unsigned int i=0;i<greenK.media()->lunghezza()/9;i++) {
                for (unsigned int j=0;j<9;j++) {
                    std::cout << greenK.media()->elemento(i*9+j)*factors[j] << " "
                              << greenK.varianza()->elemento(i*9+j)*factors[j]*factors[j] << " ";
                }

                std::cout << lambda_conv[i]*factor_conv<< " "<<lambda_conv_var[i]*factor_conv*factor_conv<<" "
                          << factor_conv*(greenK.media()->elemento(i*9+3)-pow(greenK.media()->elemento(i*9+5),2)/greenK.media()->elemento(i*9+4)) << "\n";

            }

        }else if (velocity_h) {
            std::cerr << "Inizio del calcolo dell'istogramma della velocità...\n";
            Traiettoria test(input);
            test.imposta_dimensione_finestra_accesso(1);
            test.imposta_inizio_accesso(0);
            MediaBlocchi<IstogrammaVelocita,unsigned int,double> istogramma_vel(&test,blocknumber,nbins,vmax_h);
            istogramma_vel.calcola();

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
            std::cerr << "Inizio analisi dei modi normale di vibrazione del cristallo...\n";
            Traiettoria test(input);
            test.imposta_dimensione_finestra_accesso(1);
            test.imposta_inizio_accesso(0);
            ModiVibrazionali test_normali(&test,ifcfile,fononefile,numero_thread,blocksize);
            test_normali.reset(numero_frame);
            test_normali.calcola(0);

            return 1;
        } else if (spettro_vibraz) {
            std::cerr << "Inizio del calcolo dello spettro vibrazionale...\n";
            Traiettoria test(input);
            test.imposta_dimensione_finestra_accesso(1);
            test.imposta_inizio_accesso(0);
//            SpettroVibrazionale test_spettro(&test);
            MediaBlocchi<SpettroVibrazionale> test_spettro_blocchi(&test,blocknumber);

            test_spettro_blocchi.calcola();
            for (unsigned int i=0;i<test_spettro_blocchi.media()->lunghezza()/(3*test.get_ntypes());i++) {
                std::cout << i << " ";
                for (unsigned int itipo=0;itipo<test.get_ntypes();itipo++)
                    std::cout << test_spettro_blocchi.media()->spettro(i,0,itipo) << " " << test_spettro_blocchi.varianza()->spettro(i,0,itipo) << " "
                              << test_spettro_blocchi.media()->spettro(i,1,itipo) << " " << test_spettro_blocchi.varianza()->spettro(i,1,itipo) << " "
                              << test_spettro_blocchi.media()->spettro(i,2,itipo) << " " << test_spettro_blocchi.varianza()->spettro(i,2,itipo) << " ";
                std::cout << "\n";
            }

        } else {
            std::cerr << "Inizio del calcolo del coefficiente di trasporto termico per un fluido a una componente...\n";
            Analisi traiettoria(input,log_input,numero_frame);
            //traiettoria.traiettoria.allocate_J();
            traiettoria.traiettoria.calc_J_autocorrelation(blocksize,elast);

            if (corr_out!="") {
                std::ofstream out;
                out.open(corr_out.c_str(),std::ofstream::app);
                if (out.is_open()) {
                    out << "\n\n#block size "<< blocksize <<"\n";
                    for (int i=0;i<traiettoria.traiettoria.get_n_J_autocorr();i++) {
                        out<<traiettoria.time_history[i]-traiettoria.time_history[0]<< " " << traiettoria.traiettoria.J_correlation(i) << " "<< traiettoria.traiettoria.J_correlation_var(i)<< " "<<traiettoria.traiettoria.J_correlation_int(i) <<" " << traiettoria.traiettoria.J_correlation_int_var(i) <<"\n";
                    }
                    out.flush();
                }
            }
            if (ris_append_out!="") {
                std::ofstream out;
                out.open(ris_append_out.c_str(),std::ofstream::app);
                if (out.is_open()) {
                    double T,T_var,P,P_var;
                    traiettoria.T_history.media_var_calc_block(T,T_var,blocksize);
                    traiettoria.pressure_history.media_var_calc_block(P,P_var,blocksize);
                    out << "\n\n#block size "<< blocksize <<"\n";
                    out << T << " " << T_var << " " << P << " " << P_var<< " ";
                    out.flush();
                }
            }
        }

    }

    catch (const std::exception& e) {
        std::cerr << e.what()<<"\n";
        return 1;
    }




    return 0;
}
