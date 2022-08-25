/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "modivibrazionali.h"
#include <fstream>
#include <complex>
#include <vector>
#include <thread>
#include <mutex>
#include <iostream>


const double evjoules  = 1.60217733e-19, // from eV to Joules
massunit  = 1.6605402e-27,  //# from atomic mass unit to kilogramm
angsmeter = 1.0e-20;       // # from Angstrom^2 to meter^2
const double conv=evjoules/(angsmeter*massunit*1e24); //time in picoseconds


void istg(double a,/// estremo inferiore dell'intervallo da istogrammare
           double b,///estremo superiore dell'intervallo da istogrammare
           unsigned int z,/// numero di intervalli dell'istogramma
           double *ist,/// array che conterrà dopo l'istogramma calcolato. E' responsabilità dell'utente la sua inizializzazione (a zero o a quello che si vuole)
           double val ///elemento da inserire nell'istogramma
           );



//unità immaginaria
const std::complex<double> I(0.0,1.0);

ModiVibrazionali::ModiVibrazionali(Trajectory * tr, std::string ifcfile, std::string fononefile, unsigned int n_threads,unsigned int timestep_blocco) : VectorOp<ModiVibrazionali>()
{
    if (n_threads==0)
        numero_threads=1;
    else
        numero_threads=n_threads;
    traiettoria=tr;
    fileFononi=fononefile;
    posizioni_equilibrio=new PosizioniEquilibrio(tr,timestep_blocco);
    numero_timesteps=0;
    //il numero di modi normali totali è N, dove N è il numero di atomi totali del cristallo
    read_force_file(ifcfile);
    timestepBlocco=timestep_blocco;
}

void ModiVibrazionali::read_force_file(std::string f) {
    //leggi il file ifc prodotto dal programma in python
    ifc.resize(3*traiettoria->get_natoms(),3*traiettoria->get_natoms());
    std::ifstream ifcfile(f.c_str());
    double t;
    if (!ifcfile.good()) //controlla che il file non sia arrivato alla fine
    {
        std::cerr << "File '"<<f<<"' non valido!\n";
        abort();
    }
    //scarta le prime 6 righe
    std::string line;
    for (unsigned int i=0;i<6;i++)
        std::getline(ifcfile,line);
    int nat_f=0,nty_f=0;
    ifcfile >> nty_f;
    ifcfile >> nat_f;
    if (nty_f != traiettoria->get_ntypes() || nat_f != traiettoria->get_natoms()) {
        std::cerr << "Errore: il numero di atomi (o di tipi) non corrisponde fra la traiettoria fornita ed il file ifc!\n";
        abort();
    }

    // legge anche le masse degli atomi
    for (unsigned int i=0;i<nat_f;i++) {
        int tipo;
        double m,x,y,z;
        ifcfile >> tipo >> m >> x >> y >> z;
        if (tipo<=0 || tipo > traiettoria->get_ntypes()) {
            std::cerr << "Errore: numero del tipo di atomo non compreso fra 1 e "<< traiettoria->get_ntypes()<<".\n";
        } else {
            traiettoria->set_mass(tipo-1,m);
        }
    }

    //legge la mega matrice delle derivate delle forze
    for (unsigned int i=0;i<3*nat_f;i++)
        for (unsigned int j=0;j<3*nat_f;j++) {
            double t;
            if (!ifcfile.good()) //controlla che il file non sia arrivato alla fine
            {
                std::cerr << "Errore di I/O sul file '"<<f<<"' !\n";
                abort();
            }
            ifcfile >> t;
            ifc(i,j)=t;
        }

}

unsigned int ModiVibrazionali::nExtraTimesteps(unsigned int n_b){
    return 0;
}

void ModiVibrazionali::reset(unsigned int s) {
    numero_timesteps=s;
    data_length=traiettoria->get_natoms()*3;
    delete [] vdata;
    vdata=new double[data_length];
    posizioni_equilibrio->reset(s);
}

void ModiVibrazionali::azzera() {
    for (unsigned int i=0;i<data_length;i++) {
        vdata[i]=0;
    }
}

void ModiVibrazionali::calculate(unsigned int primo) {
    posizioni_equilibrio->calculate(primo);

    /*calcola, per ogni atomo e per ogni modo normale, il prodotto scalare
     * fra lo spostamento dell'atomo dalla sua posizione d'equilibrio e
     * l'autovettore di quel modo normale. L'autovettore avrà tre componenti
     * (x y z) per ciascun atomo della base del reticolo cristallino.
     * Per esempio NaCl ha una base fatta da Na e Cl (2 atomi). La matrice
     * dinamica, che dipende dal modo normale scelto, sarà una matrice 6x6.
     * Tale matrice avrà al massimo 6 autovalori e 6 autovettori. Questi
     * autovettori rappresentano l'ampiezza dello spostamento per ciascun
     * atomo della base per ciascuna direzione. Questa poi andrà moltiplicata
     * per la fase e^{ik\dot r}, dove k è lo stesso. L'autovettore può avere
     * delle componenti complesse: questo indica uno sfasamento rispetto alla
     * fase data dal vettore k.
     * NOTA: assumo che il numero di atomi della base è uguale al numero di
     * tipi di atomi della traiettoria.
    */

    /*adesso devo calcolare il passo reticolare, e quindi i vari vettori k
     * che vado a considerare
    */

    posizioni_equilibrio->fit_nacl();
    posizioni_equilibrio->calculate_new_pos();

    /*
     * ho un numero di modi normali pari al numero di celle che considero.
     * Per esempio se ho 1728 atomi in una cella di simulazione cubica di NaCl,
     * ho 216 celle elementari di 8 atomi (considero la cella con la stessa simmetria
     * della simulazione). Avrò 216 vettori d'onda nella prima zona di Brillouin.
     *
     * Comunque alla fine il numero di autovettori e autovalori che posso ottenere sarà
     * 3*atomi_per_cella*celle=3*numero_totale_atomi
    */

    // il *2 è per sicurezza
    autovettori.resize(3*posizioni_equilibrio->get_atoms_cell(),traiettoria->get_natoms()*3*2);
    autovalori.resize(3*traiettoria->get_natoms()*2);
    vettori_onda.resize(3,posizioni_equilibrio->get_number_cells()*2);

    /*
     * adesso ho il passo reticolare,
     * per ogni k all'interno della zona di brillouin calcolo il modo di vibrazione
    */

    std::ofstream eigenfile((fileFononi+".eigenvalues").c_str());

    unsigned int mcont=0;
    double kmin[3],kmax[3];
    posizioni_equilibrio->get_brillouin_limit(kmax,kmin);
    double nmin[3]={kmin[0]*posizioni_equilibrio->get_simulation_size(0)/(2*std::acos(-1))-1,
               kmin[1]*posizioni_equilibrio->get_simulation_size(1)/(2*std::acos(-1))-1,
               kmin[2]*posizioni_equilibrio->get_simulation_size(2)/(2*std::acos(-1))-1
                };
    double nmax[3]={kmax[0]*posizioni_equilibrio->get_simulation_size(0)/(2*std::acos(-1))+1,
                 kmax[1]*posizioni_equilibrio->get_simulation_size(1)/(2*std::acos(-1))+1,
                 kmax[2]*posizioni_equilibrio->get_simulation_size(2)/(2*std::acos(-1))+1
                };
    for (int nx=nmin[0];nx<=nmax[0];nx++)
        for (int ny=nmin[1];ny<=nmax[1];ny++)
            for (int nz=nmin[2];nz<=nmax[2];nz++) {
                double k[3]={2*std::acos(-1)*double(nx)/posizioni_equilibrio->get_simulation_size(0),
                             2*std::acos(-1)*double(ny)/posizioni_equilibrio->get_simulation_size(1),
                             2*std::acos(-1)*double(nz)/posizioni_equilibrio->get_simulation_size(2)};
                if (posizioni_equilibrio->zona_brillouin(k)) {
#ifdef DEBUG
                    std::cerr << "k[3]= {" << k[0]<< ", "<<k[1]<<", "<<k[2]<<")\n";
#endif
                    //calcola il modo di vibrazione, iniziando dalla matrice dinamica
                    Eigen::MatrixXcd dynMat,dynMat_h;
                    dynMat.resize(3*posizioni_equilibrio->get_atoms_cell(),3*posizioni_equilibrio->get_atoms_cell());
                    dynMat_h.resize(3*posizioni_equilibrio->get_atoms_cell(),3*posizioni_equilibrio->get_atoms_cell());


/**
  *****************************************************************************************************************
  *****************************************************************************************************************
        Begin calculation of dynamical matrix
  *****************************************************************************************************************
  *****************************************************************************************************************
**/

                    for (unsigned int a=0;a<3*posizioni_equilibrio->get_atoms_cell();a++)
                        for (unsigned int b=0;b<3*posizioni_equilibrio->get_atoms_cell();b++)
                            dynMat(a,b)=0;


                    for (unsigned int iatom_cell=0;iatom_cell<posizioni_equilibrio->get_atoms_cell();iatom_cell++)
                        for (unsigned int icoord=0;icoord<3;icoord++)
                            for (unsigned int jcoord=0;jcoord<3;jcoord++)
                                // la matrice dinamica è fatta di blocchi 3x3 dove si considera l'atomo della cella i e j
                                for (unsigned int jatom_all=0;jatom_all<traiettoria->get_natoms();jatom_all++) {
                                    std::complex<double> element=0.0;
                                    for (unsigned int jimage=0;jimage<posizioni_equilibrio->get_atom_nearest_image_translation(jatom_all).size();jimage++){
                                        //calcola r-r'
                                        double rk=0.0; //= k*(r_i-r_0)
                                        for (unsigned int ii=0;ii<3;ii++){ //(loop over coordinates)
                                            rk+=(posizioni_equilibrio->get_atom_position_origin_cell(iatom_cell)[ii]-
                                                  (posizioni_equilibrio->get_fitted_pos(jatom_all)[ii] + posizioni_equilibrio->get_atom_nearest_image_translation(jatom_all).at(jimage)[ii] )
                                                 )   *   k[ii];
                                            // this r_i-r_0 is the relative position of the cell of the j atom
                                        }
                                        element +=
                                                ifc(posizioni_equilibrio->get_atom_index_origin_cell(iatom_cell)*3+icoord,
                                                    //          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                                    //               atom index of the 0th cell
                                                    jatom_all*3+jcoord)*exp(I*rk);
                                        //                                  ^^^^
                                        //                                 i k*(r_i-r_0)
                                    }
                                    element/=posizioni_equilibrio->get_atom_nearest_image_translation(jatom_all).size();
                                    dynMat(icoord+3*iatom_cell,
                                           jcoord+3*posizioni_equilibrio->get_atom_base_index(jatom_all)) +=element;
                                }


                    // la radice delle masse è già inserica nel file ifc
/*
                    //moltiplica per il coefficiente 1/sqrt(m_i m_j) le comoponenti della matrice
                    for (unsigned int iatom=0;iatom<posizioni_equilibrio->get_atoms_cell();iatom++)
                        for (unsigned int jatom=0;jatom<posizioni_equilibrio->get_atoms_cell();jatom++)
                        {
                            double factor=1.0/std::sqrt(traiettoria->get_mass(  posizioni_equilibrio->get_type_base(iatom) )
                                                        *traiettoria->get_mass( posizioni_equilibrio->get_type_base(jatom) ));
                            for (unsigned int icoord=0;icoord<3;icoord++)
                                for (unsigned int jcoord=0;jcoord<3;jcoord++)
                                    dynMat(icoord+iatom*3,jcoord+jatom*3) *= factor;
                        }
*/
                    //voglio che la matrice sia autoaggiunta

                    dynMat_h= (dynMat + dynMat.adjoint())*0.5;

                    /**
                      *****************************************************************************************************************
                      *****************************************************************************************************************
                            End calculation of dynamical matrix,
                            begin of diagonalization
                      *****************************************************************************************************************
                      *****************************************************************************************************************
                    **/


                    //diagonalizza la matrice e trova autovettori e autovalori

                    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> diagonalizzatore;
                    diagonalizzatore.compute(dynMat_h);

                    // adesso decompone gli spostamenti e le velocità, per calcolare l'energia per modo

                    if (mcont >=3*traiettoria->get_natoms()) {
                        std::cerr << "Raggiunto il numero massimo di modi normali (e ancora ne trovo?!)\n";
                        goto fine_calcolo_fononi;
                    }
                    //write in sub block of mega-eigenvalue column vector the eigenvalues
                    /*  ___
                        |e|
                        |i|
                        |g|
                        |1|
                        ---
                        |e|
                        |i|
                        |g|
                        |2|
                        ---
                        |e|
                        |i|
                        |g|
                        |3|
                        ---
                        |.|
                        |.|
                        |.|
                        | |
                        ---
                         .
                         .
                         .
                    */
                    autovalori.block(mcont *3*posizioni_equilibrio->get_atoms_cell(),0,
                                     3*posizioni_equilibrio->get_atoms_cell(),1)=diagonalizzatore.eigenvalues();

                    //write in sub block of mega-eigenvector matrix eigenvectors
                    //  |      |      | ... |  ...
                    //  | eig1 | eig2 | ... |  ...
                    //  |      |      | ... |  ...
                    autovettori.block(0,3*posizioni_equilibrio->get_atoms_cell()*mcont,//indici da cui iniziare a prendere la sottomatrice
                                      3*posizioni_equilibrio->get_atoms_cell(),3*posizioni_equilibrio->get_atoms_cell())=
                            diagonalizzatore.eigenvectors();
                    vettori_onda.col(mcont)=Eigen::Map<Eigen::Vector3d>(k);

#ifdef DEBUG
                    /*
                    std::cerr << "Test autovettori: ";
                    for (unsigned int iaut=0;iaut<posizioni_equilibrio->get_atoms_cell()*3 ;iaut++)
//                        std::cerr << dynMat_h*diagonalizzatore.eigenvectors().col(iaut)-diagonalizzatore.eigenvalues()(iaut)*diagonalizzatore.eigenvectors().col(iaut)<< "\n";
                        std::cerr << dynMat_h*autovettori.col(3*posizioni_equilibrio->get_atoms_cell()*mcont+iaut)-
                                     autovalori(3*posizioni_equilibrio->get_atoms_cell()*mcont+iaut)*
                                     autovettori.col(3*posizioni_equilibrio->get_atoms_cell()*mcont+iaut)<< "\n";
                    std::cerr << "\n";
*/
                    eigenfile << k[0]<< " "<< k[1] << " " << k[2];
                    for (unsigned int ieigen=0;ieigen<3*posizioni_equilibrio->get_atoms_cell();ieigen++)
                            eigenfile << " "<< diagonalizzatore.eigenvalues()(ieigen);
                    eigenfile << "\n";
                    eigenfile.flush();
#endif

                    mcont++;
                }
            }
fine_calcolo_fononi:



#ifdef DEBUG
    std::cerr << "Numero vettori d'onda calcolati: "<<mcont<<"\n";
    std::ofstream spostamenti((fileFononi+".spostamenti_normCoord").c_str());
#ifdef DEBUG2
    for (unsigned int imode=0;imode<mcont;imode++)
    for (unsigned int ivector=0;ivector<3*posizioni_equilibrio->get_atoms_cell();ivector++)
        for (unsigned int jvector=0;jvector<3*posizioni_equilibrio->get_atoms_cell();jvector++){
//                        std::cerr << dynMat_h*diagonalizzatore.eigenvectors().col(iaut)-diagonalizzatore.eigenvalues()(iaut)*diagonalizzatore.eigenvectors().col(iaut)<< "\n";
            std::cerr << autovettori.col(ivector+imode*3*posizioni_equilibrio->get_atoms_cell()).dot (
                             autovettori.col(jvector+imode*3*posizioni_equilibrio->get_atoms_cell())
                             ) << "\n";


            std::cerr <<
            autovettori.block(0,imode*3*posizioni_equilibrio->get_atoms_cell(),
                     3*posizioni_equilibrio->get_atoms_cell(),3*posizioni_equilibrio->get_atoms_cell()).row(ivector).dot(
                        autovettori.block(0,imode*3*posizioni_equilibrio->get_atoms_cell(),
                                          3*posizioni_equilibrio->get_atoms_cell(),3*posizioni_equilibrio->get_atoms_cell()).row(jvector)
                        ) <<  "\n";
        }
    std::cerr << "\n";
#endif

#endif

    //adesso devo scomporre le vibrazioni usando gli autovettori calcolati, per ogni timestep

    std::ofstream fononi(fileFononi.c_str());
    std::ofstream fononi_etot((fileFononi+".TVtot").c_str());


    double *TV = new double [numero_threads*mcont*3*posizioni_equilibrio->get_atoms_cell()*2];
    double *TVEmedio = new double [mcont*3*posizioni_equilibrio->get_atoms_cell()*3];
    double *TVEvar = new double [mcont*3*posizioni_equilibrio->get_atoms_cell()*3];
    for (unsigned int i=0;i<mcont*3*posizioni_equilibrio->get_atoms_cell()*3;i++)
        TVEmedio[i]=TVEvar[i]=0;
    unsigned int nmean=0;


#ifdef DEBUG
    std::complex<double> *Q = new std::complex<double> [numero_threads*traiettoria->get_natoms()*3];
    std::complex<double> *Qdot = new std::complex<double> [numero_threads*traiettoria->get_natoms()*3];
    double *U = new double [numero_threads*traiettoria->get_natoms()*3];
#endif

    std::mutex mut;

    // se ho una traiettoria troppo lunga, che non sta nella ram in un colpo solo, devo caricare il file a pezzi (compito svolto da traiettoria)
    unsigned int ntimestep_load=timestepBlocco + numero_threads - timestepBlocco%numero_threads;
    if (timestepBlocco!=0)
        traiettoria->set_data_access_block_size(ntimestep_load); //quanti dati carico in memoria

    for (unsigned int istep=0;istep<numero_timesteps;istep+=numero_threads) {
        if (timestepBlocco!=0 && istep%ntimestep_load==0){
            traiettoria->set_access_at(istep+primo); // carica il prossimo blocco di dati
        }


        //calcola le coordinate dei modi normali
        std::vector<std::thread> threads;
        for (unsigned int ith=istep;ith-istep<numero_threads && ith<numero_timesteps;ith++) {
            threads.push_back(std::thread([&,ith](){
#ifdef DEBUG
                mut.lock();
                std::cerr << "Passo ith="<<ith<<"\n";
                mut.unlock();
                // mostra gli spostamenti
                for (unsigned int i=0;i<traiettoria->get_natoms();i++){
                    double ut[3];
                    posizioni_equilibrio->get_displacement(i,ith+primo,ut);
                    for (unsigned int c=0;c<3;c++)
                        U[(ith-istep)*3*traiettoria->get_natoms()+i*3+c]=ut[c];
//                        spostamenti << ut[c]<<" ";
                }
#endif
                //loop over normal modes
                for (unsigned int imodo=0;imodo<mcont*3*posizioni_equilibrio->get_atoms_cell();imodo++){
                    //normal mode coordinates, initialize to zero
                    std::complex<double> q(0,0),q_punto(0,0);
                    //loop over atoms
                    for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++){
                        //get displacements over equilibrium positions
                        double u[3];
#ifdef DEBUG
#endif
                        posizioni_equilibrio->get_displacement(iatom,ith+primo,u); // calculate displacement of iatom, current timestep

                        q += sqrt(traiettoria->get_mass( traiettoria->get_type(iatom) )) //square root of mass
                                // I is imaginary unit
                                *exp(I*vettori_onda.col(imodo/(3*posizioni_equilibrio->get_atoms_cell())).dot(      Eigen::Map<Eigen::Vector3d>(posizioni_equilibrio->get_fitted_pos(iatom)) - Eigen::Map<Eigen::Vector3d>(posizioni_equilibrio->get_atom_position_origin_cell(posizioni_equilibrio->get_atom_base_index(iatom))) ) )*
                                //      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                //       k vector of the normal mode                                  scalar product                    reticula position of iatom
                                (autovettori.block<3,1>(0+3*posizioni_equilibrio->get_atom_base_index(iatom),imodo).dot(
                            //   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^scalar product
                            //      eigenvector of the normal mode corrisponding to the base index of iatom
                                     Eigen::Map<Eigen::Vector3d>(u)) );
                            //       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ displacement
                        q_punto += sqrt(traiettoria->get_mass( traiettoria->get_type(iatom) ))
                                *exp(I*vettori_onda.col(imodo/(3*posizioni_equilibrio->get_atoms_cell())).dot(       Eigen::Map<Eigen::Vector3d>(posizioni_equilibrio->get_fitted_pos(iatom)) - Eigen::Map<Eigen::Vector3d>(posizioni_equilibrio->get_atom_position_origin_cell(posizioni_equilibrio->get_atom_base_index(iatom)))   )    )*
                                (autovettori.block<3,1>(0+3*posizioni_equilibrio->get_atom_base_index(iatom),imodo).dot(
                                 Eigen::Map<Eigen::Vector3d>(traiettoria->velocita(ith+primo,iatom))) );
                        //        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                        //                           velocity of iatom, current timestep
                    }
//                    q /= sqrt(traiettoria->get_natoms());
//                    q_punto /= sqrt(traiettoria->get_natoms());
#ifdef DEBUG
                    Q[(ith-istep)*traiettoria->get_natoms()*3+imodo]=q;
                    Qdot[(ith-istep)*traiettoria->get_natoms()*3+imodo]=q_punto;
#endif
                    //            double T,V;
                    TV[(ith-istep)*mcont*3*posizioni_equilibrio->get_atoms_cell()*2+imodo*2] = 0.5*conv*autovalori(imodo)*(q.real()*q.real()+q.imag()*q.imag());
                    TV[(ith-istep)*mcont*3*posizioni_equilibrio->get_atoms_cell()*2+imodo*2+1] =0.5*(q_punto.real()*q_punto.real()+q_punto.imag()*q_punto.imag());
                    //            fononi << T<<" " <<V <<" "<< T+V<<" ";


                }
#ifdef DEBUG

#endif

            }));
        }
        for (unsigned int ithread=0;ithread<threads.size();ithread++){
            threads[ithread].join();
#ifdef DEBUG
            std::cerr << "Thread ithread="<<ithread<<" finito\n";
#endif
        }
        for (unsigned int itvv=0;itvv<threads.size();itvv++){
            double Ttot=0.0,Vtot=0.0;
            for (unsigned int itv=0;itv<mcont*3*posizioni_equilibrio->get_atoms_cell();itv++){
                fononi << TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2]<<" " <<TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2+1] <<" "<<TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2] + TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2+1]<<" ";
                Ttot+=TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2];
                Vtot+=TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2+1];

                double delta,delta2;


                delta=TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2]-TVEmedio[itv*3+0];
                TVEmedio[itv*3+0]+=delta/(++nmean);
                delta2=TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2]-TVEmedio[itv*3+0];
                TVEvar[itv*3+0]+=delta*delta2;


                delta=TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2 +1]-TVEmedio[itv*3+1];
                TVEmedio[itv*3+1]+=delta/(nmean);
                delta2=TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2 +1]-TVEmedio[itv*3+1];
                TVEvar[itv*3+1]+=delta*delta2;


                delta=TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2]+TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2 +1]
                        -TVEmedio[itv*3+2];
                TVEmedio[itv*3+2]+=delta/(nmean);
                delta2=TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2]+TV[itvv*mcont*3*posizioni_equilibrio->get_atoms_cell()*2 +itv*2 +1]
                        -TVEmedio[itv*3+2];
                TVEvar[itv*3+2]+=delta*delta2;
            }

            fononi_etot << Ttot << " " << Vtot << " "<< Ttot+Vtot<< "\n";
            fononi << "\n";
#ifdef DEBUG
            for (unsigned int iq=0;iq<traiettoria->get_natoms();iq++)
                spostamenti <<Q[itvv*traiettoria->get_natoms()*3+iq].real()<< " "<<  Q[itvv*traiettoria->get_natoms()*3+iq].imag() << " "
                            << Qdot[itvv*traiettoria->get_natoms()*3+iq].real() << " "<< Qdot[itvv*traiettoria->get_natoms()*3+iq].imag() << " ";
            spostamenti << "  ";
            for (unsigned int iu=0;iu<3*traiettoria->get_natoms();iu++)
                spostamenti << " " << U[3*traiettoria->get_natoms()*itvv+iu];
            spostamenti << "\n";
#endif
        }
        fononi.flush();
        fononi_etot.flush();
#ifdef DEBUG
        spostamenti.flush();
#endif
        threads.clear();
    }
    fononi.close();
    fononi_etot.close();

    std::ofstream energiaMedia((fileFononi+".energia_media").c_str());
    energiaMedia << "#kx ky kz  T Tvar V Vvar E Evar\n";
    for (unsigned int imodo=0;imodo<mcont;imodo++){
        for (unsigned int icell=0;icell<posizioni_equilibrio->get_atoms_cell()*3;icell++){
            for (unsigned int i=0;i<3;i++)
                energiaMedia << vettori_onda.col(imodo)(i) << " ";
            for (unsigned int i=0;i<3;i++)
                energiaMedia << TVEmedio[ 3*posizioni_equilibrio->get_atoms_cell()*imodo *3
                                         +icell*3
                                         + i] << " "
                            << TVEvar[3*posizioni_equilibrio->get_atoms_cell()*imodo *3 +icell*3 + i]/((nmean-1)*nmean) << " ";
            energiaMedia << "\n";
        }
    }




    delete [] TV;
    delete [] TVEmedio;
    delete [] TVEvar;

    delete [] Q;
    delete [] Qdot;



}











/* questo probabilmente non serve, basta fare l'analisi a blocchi del risultato finale, senza includere anche le posizioni medie
ModiVibrazionali & ModiVibrazionali::operator =(const ModiVibrazionali &d){
    VectorOp<ModiVibrazionali> operator= (d);
    posizioni_equilibrio
}

ModiVibrazionali & ModiVibrazionali::operator+= (const ModiVibrazionali&d){

}

ModiVibrazionali & ModiVibrazionali::operator-= (const ModiVibrazionali&d){

}

ModiVibrazionali & ModiVibrazionali::operator*= (const ModiVibrazionali&d){

}

ModiVibrazionali & ModiVibrazionali::operator/= (const ModiVibrazionali&d){

}
*/
