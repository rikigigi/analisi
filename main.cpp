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
//#include "interfaccia.h"
#include "testtraiettoria.h"
#include "spettrovibrazionale.h"
#include "modivibrazionali.h"
#include "boost/version.hpp"
#include "config.h"
#include "istogrammavelocita.h"
#include "greenkubo2componentionicfluid.h"
#include "greenkuboNcomponentionicfluid.h"
#include "convolution.h"
#include "heatfluxts.h"
#include "msd.h"
#include "istogrammaatomiraggio.h"
#include <functional>
#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif
#include "convertibinario.h"
#ifdef DEBUG
#include "readlog.h"
#endif


int main(int argc, char ** argv)
{
    boost::program_options::options_description options("Riccardo Bertossa, 2017\nProgramma per l'analisi di traiettorie di LAMMPS, principalmente finalizzato al calcolo del coefficiente di conducibilità termica tramite l'analisi a blocchi.\n\nOpzioni consentite");
    std::string input,log_input,corr_out,ris_append_out,ifcfile,fononefile,output_conversion;
    int sub_mean_start=0,numero_frame=0,blocksize=0,elast=0,blocknumber=0,numero_thread,nbins,skip=1,conv_n=20,final=60,stop_acf=0;
    unsigned int n_seg=0;
    bool sub_mean=false,test=false,spettro_vibraz=false,velocity_h=false,heat_coeff=false,debug=false,debug2=false,dumpGK=false,msd=false,bench=false;
    double vmax_h=0,cariche[2],dt=5e-3,vicini_r=0.0;
    std::vector<unsigned int > cvar_list,kk_l;
    std::vector<double> factors_input;
    std::vector<std::string> headers;
    std::vector< std::pair <unsigned int,unsigned int > > cvar;
    options.add_options()
        #if BOOST_VERSION >= 105600
            ("input,i",boost::program_options::value<std::string>(&input)->default_value(""), "file di input nel formato binario di LAMMPS: id type xu yu zu vx vy vz")
            ("loginput,l",boost::program_options::value<std::string>(&log_input),"file di input con il log di LAMMPS e i dati termodinamici nel formato: Step Time PotEng TotEng Lx Press Temp c_flusso[1] c_flusso[2] c_flusso[3]")
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
            ("headers,a",boost::program_options::value<std::vector<std::string> > (&headers)->multitoken(),"Specifica le colonne del log di LAMMPS da utilizzare per il calcolo del coefficiente. La prima colonna specificata deve essere quella della corrente di energia. Utilizza il codice generale con l'inversione della matrice.")
            ("dt,D",boost::program_options::value<double>(&dt)->default_value(5e-3),"Intervallo di dump in picosecondi delle correnti nel file di log di LAMMPS")
            ("skip,s",boost::program_options::value<int>(&skip)->default_value(1),"incremento minimo di timesteps da usare nel calcolo delle funzioni di correlazione o nello spostamento quadratico medio")
            ("charge1",boost::program_options::value<double>(&cariche[0])->default_value(1.0),"carica in unità elementari del tipo 1")
            ("charge2",boost::program_options::value<double>(&cariche[1])->default_value(-1.0),"carica in unità elementari del tipo 2")
            ("conv_n,C",boost::program_options::value<int>(&conv_n)->default_value(10),"sigma della gaussiana con cui viene fatta la convoluzione del coefficiente di trasporto termico al variare del tempo di integrazione (in numero di frame)")
            ("final,f",boost::program_options::value<int>(&final)->default_value(60),"numero di frame a cui estrarre il risultato finale")
            ("dump-block-H,d",boost::program_options::bool_switch(&dumpGK)->default_value(false),"scrive in dei file (aggiungendo ogni volta alla fine) i calcoli di ogni signolo blocco per il calcolo del coefficiente di trasporto termico in un sale fuso")
            ("stop,S",boost::program_options::value<int>(&stop_acf)->default_value(0),"lunghezza massima, in frame, di tutte le funzioni di correlazione e relativi integrali o dello spostamento quadratico medio. Se posto a zero è pari alle dimensioni del blocco")
            ("covarianze,z",boost::program_options::value<std::vector<unsigned int > >(&cvar_list)->multitoken(),"nel calcolo del coefficiente di conducibilità, oltre alla media e alla varianza di tutte le variabili calcola anche la covarianza della coppia di quantità calcolate indicate. Deve essere un numero pari di numeri")
            ("mean-square-displacement,q",boost::program_options::bool_switch(&msd)->default_value(false),"calcola e stampa nell'output lo spostamento quadratico medio per ogni specie atomica")
            ("factors,F",boost::program_options::value<std::vector<double> >(&factors_input)->multitoken(),"imposta i fattori dall'esterno. (in ordine: fattore, fattore di integrazione). Le funzioni di autocorrelazione vengono moltiplicate per il fattore, e gli integrali per fattore*fattore di integrazione. Legge solo le colonne delle correnti.")
            ("subtract-mean",boost::program_options::bool_switch(&sub_mean)->default_value(false),"sottrae la media dalla funzione di correlazione, calcolata a partire dal timestep specificato con -u")
            ("subtract-mean-start,u",boost::program_options::value<int>(&sub_mean_start)->default_value(0),"timestep nella funzione di correlazione a partire dal quale iniziare a calcolare la media")
            ("subBlock,k",boost::program_options::value<unsigned int>(&n_seg)->default_value(1),"opzione di ottimizzazione del calcolo della funzione di correlazione. Indica in quanti blocchetti suddividere il calcolo delle medie(influenza l'efficienza della cache della CPU)")
            ("kk",boost::program_options::bool_switch(&bench)->default_value(false),"esegue un benchmark per il valore ottimale di k")
            ("kk-range",boost::program_options::value<std::vector<unsigned int > >(&kk_l)->multitoken(),"valore minimo e massimo del range in cui testare k")
            ("binary-convert",boost::program_options::value<std::string>(&output_conversion),"esegui la conversione nel formato binario di lammps del file specificato come input nel file qui specificato")
            ("neighbor",boost::program_options::value<double>(&vicini_r)->default_value(0.0),"Se impostato calcola l'istogramma del numero di vicini per tipo entro il raggio specificato.")
        #ifdef DEBUG
            ("test-debug",boost::program_options::bool_switch(&debug)->default_value(false),"test vari")
            ("test-debug2",boost::program_options::bool_switch(&debug2)->default_value(false),"test vari 2")
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

        if ((output_conversion!="" && input=="") ||vm.count("help")|| (vm.count("loginput")==0 && (output_conversion=="" && !velocity_h) ) || skip<=0 || stop_acf<0 || final<0 || (!sub_mean && (sub_mean_start!=0) ) || sub_mean_start<0 || !(kk_l.size()==0 || kk_l.size()==2)){
            std::cout << "COMPILED AT " __DATE__ " " __TIME__ " by " CMAKE_CXX_COMPILER " whith flags " CMAKE_CXX_FLAGS  " on a " CMAKE_SYSTEM " whith processor " CMAKE_SYSTEM_PROCESSOR ".\n";
            std::cout << options << "\n";
            return 1;
        }

        if (cvar_list.size()%2 != 0) {
            std::cout << "Errore: la lista degli indici delle covarianze deve contenere un numero pari di numeri\n";
            std::cout << "COMPILED AT " __DATE__ " " __TIME__ " by " CMAKE_CXX_COMPILER " whith flags " CMAKE_CXX_FLAGS  " on a " CMAKE_SYSTEM " whith processor " CMAKE_SYSTEM_PROCESSOR ".\n";
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
        std::cout << "COMPILED AT " __DATE__ " " __TIME__ " by " CMAKE_CXX_COMPILER " whith flags " CMAKE_CXX_FLAGS  " on a " CMAKE_SYSTEM " whith processor " CMAKE_SYSTEM_PROCESSOR ".\n";
        std::cerr << options;
        return 1;
    }


    try {

        if (output_conversion!="") {

            ConvertiBinario conv(input,output_conversion);
            return 0;
        }


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
#ifdef DEBUG
        if (debug2){

            ReadLog<> test(log_input);
            for (unsigned int i=0;i<test.n_timestep();i++){
                std::cerr << i << " "<<test.timestep(i)<< " ";
                for (unsigned int j=0;j<test.n_data();j++){
                    std::cerr << test.line(i)[j] << " ";
                }
                std::cerr << "\n";
            }

        } else
#endif // DEBUG
            if (heat_coeff) {
                std::cerr << "Inizio del calcolo del coefficiente di trasporto termico...\n";
                ReadLog<> test(log_input);

                //qui devo aggiungere la traiettoria binaria a ReadLog, qualora ReadLog ne constati la necessità -- NOT IMPLEMENTED
                if (test.need_binary(headers)) {
                    //fai cose -- NOT IMPLEMENTED
                }

                double factor_conv=1.0;
                double factor_conv2=1.0;
                double factor_intToCorr=1.0;

                //calcola velocemente la media a blocchi per la temperatura, se non sono specificati manualmente i fattori
                const double const_charge=1.6021765,const_boltzmann=1.38064852;

                unsigned int idx_lx=0;
                double media_=0.0;
                double var_=0.0;
                if (factors_input.size()==0){

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
                    if (factors_input.size()!=2){
                        std::cerr << "Errore: specificare 2 fattori, uno per le funzioni di correlazione e uno per il suo integrale.\n";
                        abort();
                    }
                    factor_conv=factor_conv2=factors_input[0];
                    factor_intToCorr=factors_input[1];
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
                        std::cout << "# T T1sigma  atomi/volume factor1 factor2\n#"
                              <<media_ << " " << sqrt(var_) << " "
                             << test.get_natoms()/pow(test.line(0)[idx_lx] ,3)<< " "<< factor_conv<<" "<<factor_intToCorr<< "\n";
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

                    MediaBlocchiG<ReadLog<>,GreenKuboNComponentIonicFluid<>,
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
                    MediaVarCovar<GreenKuboNComponentIonicFluid<> > greenK(narr,cvar);

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

            }else if (msd){

                std::cerr << "Inizio del calcolo dello spostamento quadratico medio...\n";
                Traiettoria test(input);
                test.imposta_dimensione_finestra_accesso(1);
                test.imposta_inizio_accesso(0);

                MediaBlocchi<MSD,unsigned int,unsigned int,unsigned int> Msd(&test,blocknumber);
                Msd.calcola(skip,stop_acf,numero_thread);
                for (unsigned int i=0;i<Msd.media()->lunghezza()/test.get_ntypes();i++) {
                    for (unsigned int j=0;j<test.get_ntypes();j++)
                        std::cout << Msd.media()->elemento(i*test.get_ntypes()+j) << " " <<
                                     Msd.varianza()->elemento(i*test.get_ntypes()+j) << " ";
                    std::cout << "\n";
                }

            }else if (vicini_r>0){
                std::cerr << "Inizio del calcolo dell'istogramma del numero di vicini per tipo di tutti gli atomi\n";
                Traiettoria test(input);
                test.set_pbc_wrap(true);
                test.imposta_dimensione_finestra_accesso(1);
                test.imposta_inizio_accesso(0);

                IstogrammaAtomiRaggio h(&test,vicini_r,skip,numero_thread);
                unsigned int nt=test.get_ntimesteps(),
                        s=(nt-1)/blocknumber;
                h.reset(s);
                test.imposta_dimensione_finestra_accesso(s);
                for (unsigned int i=0;i<blocknumber;i++){
                    unsigned int t=blocknumber*i;
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


            }else if (velocity_h) {
                std::cerr << "Inizio del calcolo dell'istogramma della velocità...\n";
                Traiettoria test(input);
                test.imposta_dimensione_finestra_accesso(1);
                test.imposta_inizio_accesso(0);
                MediaBlocchi<IstogrammaVelocita,unsigned int,double> istogramma_vel(&test,blocknumber);
                istogramma_vel.calcola(nbins,vmax_h);

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

            }

    }

    catch (const std::exception& e) {
        std::cerr << e.what()<<"\n";
        return 1;
    }




    return 0;
}
