/**
  *
  * (c) Riccardo Bertossa, 2019
  *
  *   Use at  your own risk.
  *
  *   If you modified the code, I could be happy if you contribute on github!
  *
**/



#include "posizioniequilibrio.h"
#include <limits>
#include <cmath>
#include <cfenv>
#include <map>
#include "doubleround.h"

PosizioniEquilibrio::PosizioniEquilibrio(Trajectory * tr,unsigned int timesteps_sottoblocco)
{
    traiettoria=tr;
    data_length=traiettoria->get_natoms()*3;
    vdata = new double [data_length];
    posizioni_fittate_reticolo=new double[data_length];
    for (unsigned int i=0;i<data_length;i++){
        posizioni_fittate_reticolo[i]=0.0;
        vdata[i]=0.0;
    }
    traslation= new std::vector< std::array<double,3> > [traiettoria->get_natoms()];
    calcolato=false;
    lunghezza_media=traiettoria->get_ntimesteps();
    origin_cell_idx=0;
    coordinate_reticolo_atomi.resize(traiettoria->get_natoms());
    simulation_size[0]=0.0;
    simulation_size[1]=0.0;
    simulation_size[2]=0.0;
    base_tipi=0;
    atomi_per_cella=0;
    timestepSottoblocco=timesteps_sottoblocco;

}

/**
 * calcola le traslazioni delle immagini più vicine di tutti gli atomi che interagiscono con gli atomi della cella centrale.
 * In pratica costruisce una specie di sfera intorno alla sfera centrale. In questo modo, quando calcolo la matrice dinamica usando
 * un vettore d'onda che non è fra quelli che rendono il reticolo periodico, non vengono introdotte fasi che non centrano niente.
 * (nella formula della matrice la fase è equivalente alla fase dell'atomo spostato di una quantità pari alle dimensioni della
 * cella solo per i vettori d'onda compatibili con la cella stessa). Usando questa funzione, anche se sbaglio i vettori d'onda,
 * sono sicuro che dopo gli autovettori e gli autovalori sono comunque corretti (a parte il vettore d'onda).
**/
void PosizioniEquilibrio::calculate_new_pos() {

    for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
        //cerca l'immagine dell'atomo più vicina all'atomo considerato. Se ne trova più di una alla stessa distanza, conta anche quante sono.
        double norm_min=std::numeric_limits<double>::max();
        double new_pos[3]={0,0,0};
        for (int im1=-1;im1<=1;im1++)
         for (int im2=-1;im2<=1;im2++)
          for (int im3=-1;im3<=1;im3++) {
               double tmp_vec[3]={get_fitted_pos(iatom)[0]+get_simulation_size(0)*im1-get_atom_position_origin_cell(get_atom_base_index(iatom))[0],
                                  get_fitted_pos(iatom)[1]+get_simulation_size(1)*im2-get_atom_position_origin_cell(get_atom_base_index(iatom))[1],
                                  get_fitted_pos(iatom)[2]+get_simulation_size(2)*im3-get_atom_position_origin_cell(get_atom_base_index(iatom))[2],
                                  };
               double tmp_norm=sqrt(tmp_vec[0]*tmp_vec[0]+tmp_vec[1]*tmp_vec[1]+tmp_vec[2]*tmp_vec[2]);
               DoubleRound::round(&tmp_norm); //arrotonda gli ultimi bit del numero a virgola mobile
               if (tmp_norm<norm_min) {
                   norm_min=tmp_norm;
                   traslation[iatom].clear();
                   for (unsigned int j=0;j<3;j++)
                    new_pos[j]=tmp_vec[j];
                   traslation[iatom].push_back(std::array<double,3>{{get_simulation_size(0)*im1,get_simulation_size(1)*im2,get_simulation_size(2)*im3}});


               } else if (tmp_norm==norm_min) {
                   traslation[iatom].push_back(std::array<double,3>{{get_simulation_size(0)*im1,get_simulation_size(1)*im2,get_simulation_size(2)*im3}});
               }
          }
    }

#ifdef DEBUG2
    for (unsigned int i=0;i<traiettoria->get_natoms();i++) {
        for (unsigned int j=0;j<3;j++)
            std::cout << get_fitted_pos(i)[j]+traslation[i].at(0)[j]<<" ";
        std::cout << "\n";
    }
#endif
}

std::vector< std::array <double,3> > & PosizioniEquilibrio::get_atom_nearest_image_translation(unsigned int iatom){
    if (iatom<traiettoria->get_natoms()) {
        return traslation[iatom];
    } else {
        std::cerr << "Errore: richiesta la traslazione più vicina di un atomo fuori dal range ("<<iatom<<")!\n";
        abort();
    }
}

PosizioniEquilibrio::~PosizioniEquilibrio(){
    delete [] origin_cell_idx;
    delete [] posizioni_fittate_reticolo;
    delete [] traslation;
}

unsigned int PosizioniEquilibrio::nExtraTimesteps(unsigned int n_b){
    return 0;
}

void PosizioniEquilibrio::reset(const unsigned int numeroTimestepsPerBlocco){
    lunghezza_media=numeroTimestepsPerBlocco;
}

void PosizioniEquilibrio::calculate(unsigned int primo) {


    azzera();

    if (timestepSottoblocco>0)
        traiettoria->set_data_access_block_size(timestepSottoblocco);

    for (unsigned int i=0;i<lunghezza_media;i++){
        if (timestepSottoblocco>0 && i%timestepSottoblocco==0){
            traiettoria->set_access_at(primo+i);
        }
        for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++){
            for (unsigned int icoord=0;icoord<3;icoord++){
                double delta = traiettoria->positions(primo+i,iatom)[icoord]
                               - vdata[iatom*3+icoord];
                vdata[iatom*3+icoord] += delta/(i+1);
            }
        }
    }
    calcolato=true;
}

double * PosizioniEquilibrio::get_atom_position(unsigned int iatom){
    if (iatom >= traiettoria->get_natoms()){
        std::cerr << "Errore: indice dell'atomo fuori dal range!\n";
        abort();
        return 0;
    }

    return & vdata[iatom*3];
}

void PosizioniEquilibrio::reticolo_inizializza(double *base_r, double *base,unsigned int *base_type, const unsigned int &nbase) {
    reticolo << base_r[0], base_r[3], base_r[6],
                base_r[1], base_r[4], base_r[7],
                base_r[2], base_r[5], base_r[8];

    reticolo_inv=reticolo.inverse();
    reticolo_reciproco=reticolo_inv.transpose()*(2*std::acos(-1));
    base_reticolo.resize(3,nbase);
    delete [] base_tipi;
    base_tipi=new unsigned int[nbase];
    delete [] origin_cell_idx;
    origin_cell_idx=new unsigned int [nbase];
    atomi_per_cella=nbase;
    for (unsigned int i=0;i<nbase;i++){
        if (base_type[i]>=traiettoria->get_ntypes()) {
            std::cerr << "Errore: hai specificato un tipo di atomo nella base che non è compreso fra 0 e "<< traiettoria->get_ntypes()<<"!\n";
            abort();
        }
        base_tipi[i]=base_type[i];
        for (unsigned int ii=0;ii<3;ii++)
            base_reticolo(ii,i)=base[i*3+ii];
    }
}

//sicuramente questi k non appartengono alla zona di brillouin
void PosizioniEquilibrio::get_brillouin_limit(double *kmax, double *kmin){
    int idx[3];
    for (unsigned int i=0;i<3;i++) {
        kmax[i]=-std::numeric_limits<double>::max();
        kmin[i]=std::numeric_limits<double>::max();
    }
    for (idx[0]=-1;idx[0]<=1;idx[0]++)
        for (idx[1]=-1;idx[1]<=1;idx[1]++)
            for (idx[2]=-1;idx[2]<=1;idx[2]++) {
                double idxd[3]={(double)idx[0],(double)idx[1],(double)idx[2]};
                Eigen::Vector3d k=reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd);
                for (unsigned int i=0;i<3;i++) {
                    if (k(i)<kmin[i])
                        kmin[i]=k(i);
                    if (k(i)>kmax[i])
                        kmax[i]=k(i);
                }

            }
}

bool PosizioniEquilibrio::zona_brillouin(double *k_test){
    // qui, date due facce opposte parallele, devo prendere i vettori su una e togliere i vettori sull'altra.
    // devo stare attento anche al fatto che zero, con i numeri a virgola mobile, non è sempre zero...
    int idx[3];
    double k_test_n=sqrt(k_test[0]*k_test[0]+k_test[1]*k_test[1]+k_test[2]*k_test[2]);

/*
    for (idx[0]=-1;idx[0]<=1;idx[0]++)
        for (idx[1]=-1;idx[1]<=1;idx[1]++)
            for (idx[2]=-1;idx[2]<=1;idx[2]++) {
                if (idx[0]==0 && idx[1]==0 && idx[2]==0) continue;
                //
                //  l'equazione del piano perpendicolare al vettore "a" e passante per il punto "x0":
                //  a . ( x - x0 ) = 0
                //
                double idxd[3]={(double)idx[0],(double)idx[1],(double)idx[2]};
                double idxd2[3]={(double)idx[0]/2.0,(double)idx[1]/2.0,(double)idx[2]/2.0};
                double res=((reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd)).dot
                         (  Eigen::Map<Eigen::Vector3d> (k_test) - reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd2) )
                            );
                if (res>std::numeric_limits<double>::epsilon()*k_test_n*10 ) return false;
            }
*/

#define DO_IT   double idxd[3]={(double)idx[0],(double)idx[1],(double)idx[2]}; \
    double idxd2[3]={(double)idx[0]/2.0,(double)idx[1]/2.0,(double)idx[2]/2.0}; \
    double res=((reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd)).dot(  Eigen::Map<Eigen::Vector3d> (k_test) - reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd2) ) ); \
    if (res> std::numeric_limits<double>::epsilon()*k_test_n*10 ) return false;\
    double res1=-((reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd)).dot(  Eigen::Map<Eigen::Vector3d> (k_test) + reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd2) ) ); \
    if (res1> -std::numeric_limits<double>::epsilon()*k_test_n*10 ) return false;

    idx[0]=1;
       for (idx[1]=-1;idx[1]<=1;idx[1]++)
           for (idx[2]=-1;idx[2]<=1;idx[2]++) {
               DO_IT
           }
    idx[0]=0;
       idx[1]=1;
           for (idx[2]=-1;idx[2]<=1;idx[2]++) {
               DO_IT
           }
    idx[0]=0;
       idx[1]=0;
          idx[2]=1;{
              DO_IT
          }

    return true;
}

bool PosizioniEquilibrio::zona_brillouin_simple_cubic(double *k_test){
    // qui, date due facce opposte parallele, devo prendere i vettori su una e togliere i vettori sull'altra.
    // devo stare attento anche al fatto che zero, con i numeri a virgola mobile, non è sempre zero...
    double k_test_n=sqrt(k_test[0]*k_test[0]+k_test[1]*k_test[1]+k_test[2]*k_test[2]);

                /*
                 * l'equazione del piano perpendicolare al vettore "a" e passante per il punto "x0":
                 * a . ( x - x0 ) = 0
                */
 #define DO_IT_S double idxd2[3]={idxd[0]/2.0,idxd[1]/2.0,idxd[2]/2.0};\
                double res=((reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd)).dot\
                (  Eigen::Map<Eigen::Vector3d> (k_test) - reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd2) )\
                );\
                double res1= -((reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd)).dot\
                (  Eigen::Map<Eigen::Vector3d> (k_test) + reticolo_reciproco * Eigen::Map<Eigen::Vector3d>(idxd2) )\
                );\
                if (res>std::numeric_limits<double>::epsilon()*k_test_n ) return false; \
                if (res1>-std::numeric_limits<double>::epsilon()*k_test_n ) return false;
    {
        double idxd[3]={1.0,0.0,0.0};\
        DO_IT_S
    }
    {
        double idxd[3]={0.0,1.0,0.0};\
        DO_IT_S
    }
    {
        double idxd[3]={0.0,0.0,1.0};\
        DO_IT_S
    }


    return true;
}

void PosizioniEquilibrio::coord_reticolo(double *xyz, double *uvw_min, double *uvw_max) {
    Eigen::Vector3d uvw=reticolo_inv*Eigen::Map<Eigen::Vector3d>(xyz);
    for(unsigned int i=0;i<3;i++) {
        if (uvw(i)<uvw_min[i]) {
            uvw_min[i]=uvw(i);
        }
        if (uvw(i)>uvw_max[i]) {
            uvw_max[i]=uvw(i);
        }
    }
}



double PosizioniEquilibrio::d2_reticolo_spostamento_medio(double * min, double *max, /// limiti della scatola
                                        double *spostamento ) {

    //obiettivo: voglio assegnare ad ogni atomo il suo indice di cella,
    //e calcolare esattamente le positions degli atomi secondo il reticolo a temperatura zero

    //assumo che il reticolo sia fisso. devo solo trovare come traslarlo per centrare nel migliore dei modi gli atomi, che non si trovano esattamente sulle positions di temperatura zero (la temperatura della simulazione è finita)

    //sceglie un atomo di uno spigolo e lo considera come origine. devo controllare di poter costruire una cella con questo atomo
    //gli atomi dello spigolo sono quelli che hanno meno primi vicini.

    //identifica gli atomi ai bordi guardando il numero di primi vicini. Inizierà a fittare la cella usando un atomo al borda come base di partenza.

    std::multimap<double,unsigned int> distanze;
    std::multimap<unsigned int,std::vector<unsigned int> > numero_primi_vicini;

    for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++){
        distanze.clear();
        for (unsigned int jatom=0;jatom<traiettoria->get_natoms();jatom++) {
            if (iatom!=jatom) {
                double d2=0.0;
                for (unsigned int i=0;i<3;i++)
                    d2+=(vdata[iatom*3+i]-vdata[jatom*3+i])*(vdata[iatom*3+i]-vdata[jatom*3+i]);
                distanze.insert(std::pair<double,unsigned int>(d2,jatom));
            }
        }
        //guarda i primi vicini
        double first=distanze.begin()!=distanze.end() ? distanze.begin()->first : 0.0;
        unsigned int vicinicont=0;
        for (auto it= distanze.begin();it!=distanze.end();it++) {
            if (fabs(it->first-first)/first < 0.05){
                vicinicont++;
            } else {
                break;
            }
        }
        auto find = numero_primi_vicini.find(vicinicont);
        if (find==numero_primi_vicini.end())
            numero_primi_vicini.insert(std::pair<unsigned int, std::vector<unsigned int> >(vicinicont,std::vector<unsigned int> (1,iatom)));
        else
            find->second.push_back(iatom);
    }

    //guardo gli atomi che hanno meno primi vicini -- std::map ordina gli elementi secondo la chiave
    bool trovato=false;
    auto iterator=numero_primi_vicini.begin()->second.begin();
    for (;iterator!=numero_primi_vicini.begin()->second.end();iterator++) {
        //controllo che aggiungendo a questo atomo i vettori della base finisco ragionevolmente vicino ad altri atomi del tipo corretto
        for (unsigned int ibase=0;ibase<atomi_per_cella;ibase++) {
            //cerca fra tutti gli atomi quello più vicino alla posizione data
            distanze.clear();
            for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
                double d2=0.0;
                for (unsigned int i=0;i<3;i++ )
                    d2+=(vdata[iatom*3+i]-(vdata[*iterator *3+i]+base_reticolo(i,ibase)) )*
                        (vdata[iatom*3+i]-(vdata[*iterator *3+i]+base_reticolo(i,ibase)) );
                distanze.insert(std::pair<double,unsigned int>(d2,iatom));
            }
            //prendi il più vicino. se la distanza è maggiore di una cerca soglia, scarta questo atomo dell'angolo
            auto it=distanze.begin();
            for (;it!=distanze.end();it++) {
                if (it->first >base_reticolo.col(ibase).dot(base_reticolo.col(ibase))*(0.05*0.05))
                    goto prossimo_angolo;
                if (traiettoria->get_type(it->second)==base_tipi[ibase]){
                    origin_cell_idx[ibase]=it->second; //ricorda gli indici della base per un uso successivo
                    break;
                }
            }
        }
        // se sono qui tutti gli atomi sono a posto
        trovato=true;

        break;
        prossimo_angolo:
        1;
    }

    if (!trovato) {
        std::cerr << "Impossibile identificare una cella elementare sugli angoli della simulazione! Impostare il primo vettore di base a 0,0,0 o provare a cambiare il tipo di atomo alla posizione 0,0,0, e scegliere un atomo che sta in un angolo della simulazione.\n";
        abort();
    }

    double origin[3]={vdata[3* *iterator+0],vdata[3* *iterator+1],vdata[3* *iterator+2]};

    //stima molto rozza dei limiti degli indici che coprono completamente la cella di simulazione
    double coord[]={max[0]-origin[0],max[1]-origin[1],max[2]-origin[2]},
            uvw_max[]={-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max(),-std::numeric_limits<double>::max()},
            uvw_min[]={std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),std::numeric_limits<double>::max()};
    coord_reticolo(coord,uvw_min,uvw_max); // m,m,m
    coord[0]=min[0]-origin[0];
    coord_reticolo(coord,uvw_min,uvw_max); // M,m,m
    coord[1]=min[1]-origin[1];
    coord_reticolo(coord,uvw_min,uvw_max); // M,M,m
    coord[2]=min[2]-origin[2];
    coord_reticolo(coord,uvw_min,uvw_max); // M,M,M
    coord[1]=max[1]-origin[1];
    coord_reticolo(coord,uvw_min,uvw_max); // M,m,M
    coord[2]=max[2]-origin[2];
    coord[1]=min[1]-origin[1];
    coord[0]=max[0]-origin[0];
    coord_reticolo(coord,uvw_min,uvw_max); // m,M,m
    coord[1]=max[1]-origin[1];
    coord[2]=min[2]-origin[2];
    coord_reticolo(coord,uvw_min,uvw_max); // m,m,M
    coord[1]=min[1]-origin[1];
    coord_reticolo(coord,uvw_min,uvw_max); // m,M,M

    /*adesso ho i limiti nelle nuove coordinate, genero i punti del reticolo
     * nella scatola e poi calcolo la distanza da quelli della traiettoria
    */

    spostamento[0]=0.0;
    spostamento[1]=0.0;
    spostamento[2]=0.0;
    unsigned int cont=0;
    reticolo_xyz.clear();
    std::fesetround(FE_TONEAREST);
    for (int u=std::lrint( uvw_min[0]);u<=std::lrint(uvw_max[0]);u++){
        for (int v=std::lrint( uvw_min[1]);v<=std::lrint(uvw_max[1]);v++){
            for (int w=std::lrint( uvw_min[2]);w<=std::lrint(uvw_max[2]);w++){
                //controlla che l'atomo sia dentro la scatola, se si lo aggiunge
                //per ogni atomo della base...
                Eigen::Vector3d uvw;
                uvw << double(u),double(v),double(w);
                Eigen::Vector3d xyzR=reticolo*uvw;
                for (unsigned int i=0;i<base_reticolo.cols();i++){
                    Eigen::Vector3d xyz=xyzR+base_reticolo.col(i)+Eigen::Map<Eigen::Vector3d>(origin);
                    if (xyz(0)>=min[0] && xyz(0)<=max[0] &&
                        xyz(1)>=min[1] && xyz(1)<=max[1] &&
                        xyz(2)>=min[2] && xyz(2)<=max[2]) {
                        //aggiungi l'atomo
                        reticolo_xyz.insert(std::pair<std::array<int,4>,std::pair<Eigen::Vector3d,unsigned int>>
                                            ( std::array<int,4>{{u,v,w,(int)i}},std::pair <Eigen::Vector3d,unsigned int> (xyz,traiettoria->get_natoms()) )
                                            );
                    }
                }
            }
        }
    }

    if (reticolo_xyz.size()< traiettoria->get_natoms()) {
        std::cerr << "La scatola non è abbastanza grande per contenere il reticolo intero ("<<reticolo_xyz.size() <<" atomi contro "<< traiettoria->get_natoms()<<" atomi della traiettoria)!\n";
    }

    //guardo la distanza dei punti da quelli della media (da calcolare in precedenza)
    double d2=0.0,d2t,d2m;
    for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
        std::array<int,4> uvwi=reticolo_xyz.begin()->first;
            // cerco l'atomo più vicino dello stesso tipo, per restituire comunque un numero sensato
            auto iteratore = reticolo_xyz.begin();
            Eigen::Vector3d d = iteratore->second.first - Eigen::Map<Eigen::Vector3d>(&vdata[iatom*3]);
            Eigen::Vector3d dmin;
            d2m= d.dot(d)*4;
            for (iteratore++;iteratore != reticolo_xyz.end();iteratore++) {
                if (base_tipi[iteratore->first[3]]==traiettoria->get_type(iatom) && iteratore->second.second == traiettoria->get_natoms()){
                    d = iteratore->second.first - Eigen::Map<Eigen::Vector3d>(&vdata[iatom*3]);
                    d2t= d.dot(d);
                    if (d2t<d2m) {
                        d2m=d2t;
                        uvwi=iteratore->first;
                        dmin=d;
                    }
                }
            }
            if (reticolo_xyz[uvwi].second!=traiettoria->get_natoms()) {
                std::cerr << "Errore: questo punto del reticolo ha già un atomo assegnato (provare ad aumentare le dimensioni della scatola dove si genera il reticolo?)\n";
                abort();
            }
            reticolo_xyz[uvwi].second=iatom;
            d2+=d2m;
            //aggiorna lo spostamento medio dal reticolo
            double delta;
            delta=dmin(0)-spostamento[0];
            spostamento[0] += delta/++cont;
            delta=dmin(1)-spostamento[1];
            spostamento[1] += delta/cont;
            delta=dmin(2)-spostamento[2];
            spostamento[2] += delta/cont;
            coordinate_reticolo_atomi[iatom]=uvwi;
            posizioni_fittate_reticolo[iatom*3+0]=reticolo_xyz[uvwi].first(0);
            posizioni_fittate_reticolo[iatom*3+1]=reticolo_xyz[uvwi].first(1);
            posizioni_fittate_reticolo[iatom*3+2]=reticolo_xyz[uvwi].first(2);

    }

    //trasla tutto dello spostamento indicato per centrare correttamente il reticolo

    for (unsigned int iatom=0;iatom<traiettoria->get_natoms();iatom++) {
        for (unsigned int icord=0;icord<3;icord++) {
            posizioni_fittate_reticolo[iatom*3+icord] += spostamento[icord];
        }
    }

    return d2;

}

double * PosizioniEquilibrio::get_fitted_pos(unsigned int iatom) {
    if (iatom < traiettoria->get_natoms()) {
        return &posizioni_fittate_reticolo[3*iatom];
    } else {
        return 0;
    }
}

void PosizioniEquilibrio::get_displacement(unsigned int iatom, unsigned int tstep,double * displ) {

    for (unsigned int i=0;i<3;i++) {
        displ[i]=traiettoria->positions(tstep,iatom)[i]-get_fitted_pos(iatom)[i];
    }
}


std::map<std::array<int,4>,std::pair<Eigen::Vector3d,unsigned int> > & PosizioniEquilibrio::get_reticolo() {
    return reticolo_xyz;
}


double * PosizioniEquilibrio::get_atom_position_origin_cell(unsigned int icell){
    if (icell>=atomi_per_cella){
        std::cerr << "Errore: richiesta l'indice della cella origine per un indice della cella fuori range ("<<icell<<", il valore massimo è "<< atomi_per_cella-1<<" atomi per cella)\n";
        abort();
        return 0;
    }
    return &posizioni_fittate_reticolo[origin_cell_idx[icell]*3];
}

unsigned int PosizioniEquilibrio::get_type_base(unsigned int ibase) {
    if (ibase<atomi_per_cella)
    return base_tipi[ibase];
    else {
        std::cerr << "Errore: richiesto il tipo di atomo della base con un indice fuori dai valori consentiti epr la base.\n";
        abort();
        return 0;
    }
}

unsigned int PosizioniEquilibrio::get_atom_index_origin_cell(unsigned int icell){
    if (icell>=atomi_per_cella) {
        std::cerr << "Errore: richiesta l'indice della cella origine per un indice della cella fuori range ("<<icell<<", il valore massimo è "<< atomi_per_cella-1<<" atomi per cella)\n";
        abort();
        return 0;
    }
    return origin_cell_idx[icell];
}

unsigned int PosizioniEquilibrio::get_atom_base_index(unsigned int iatom) {
    if (iatom<traiettoria->get_natoms()) {
        return coordinate_reticolo_atomi[iatom][3];
    } else {
        std::cerr << "Errore: richiesto un atomo fuori dal range!\n";
        abort();
        return 0;
    }
}

unsigned int PosizioniEquilibrio::get_number_cells(){
    if (traiettoria->get_natoms()%atomi_per_cella == 0) {
        return traiettoria->get_natoms()/atomi_per_cella;
    } else {
        std::cerr << "Errore: il numero di atomi totale non è divisibile per il numero di atomi per cella!\n";
        abort();
        return traiettoria->get_natoms()/atomi_per_cella;
    }
}

double PosizioniEquilibrio::get_simulation_size(unsigned int icoord) {
    if (icoord <3) {
        double *b=traiettoria->scatola_last();
        return b[icoord*2+1]-b[icoord*2];
    } else {
        std::cerr << "Errore: richiesta una dimensione che non esiste in get_simulation_size()";
        abort();
        return 0;
    }
}

void PosizioniEquilibrio::fit_nacl(){


    if (calcolato){
        /*
     * TODO: calcola la distanza dei primi vicini media,
     * con questa fai il passo reticolare,
     * poi per calcolare lo spostamento rispetto all'origine,
     * prima imposta un atomo del tipo giusto all'origine,
     * poi calcola lo spostamento medio delgi atomi rispetto al reticolo.
     * Lo spostamento di prima più questo nuovo sarà lo spostamento effettivo del reticolo rispetto all'origine.
    */

         //trova le dimensioni della scatola, leggendole dalla traiettoria
        double max[3]={traiettoria->scatola_last()[1],traiettoria->scatola_last()[3],traiettoria->scatola_last()[5]},
               min[3]={traiettoria->scatola_last()[0],traiettoria->scatola_last()[2],traiettoria->scatola_last()[4]};


#ifdef DEBUG
        std::cerr << "Limiti della scatola (minxyz),(maxxyz) ("<< min[0]<< ", " << min[1] << ", " << min[2] << "), (" << max[0] << ", " << max[1] << ", " << max[2]<<")\n";
#endif

        /*    ______
         *   /____ /|
         *  |     | | natoms = prod (max-min)/dminM sulle coordinate
         *  |     | |
         *  |_____|/
         *
         */

        //assumo scatola quadrata
        double media_primi_vicini= cbrt( (max[0]-min[0])*(max[1]-min[1])*(max[2]-min[2])/(traiettoria->get_natoms()/8))/2.0;

        for (unsigned int i=0;i<3;i++) {
            max[i]=max[i]+media_primi_vicini;
            min[i]=min[i]-media_primi_vicini;
        }

#ifdef DEBUG
        std::cerr << "Distanza primi vicini: "<< media_primi_vicini<<"\n";
#endif

        /*

        //devo usare un reticolo cubico semplice, che è la simmetria della cella di simulazione
        double brav[] = {media_primi_vicini*2,0,0,
                         0,media_primi_vicini*2,0,
                         0,0,media_primi_vicini*2};
        double bas[] = {0,0,0, //atomi di tipo A
                        media_primi_vicini,media_primi_vicini,0,
                        0, media_primi_vicini, media_primi_vicini,
                        media_primi_vicini,0,media_primi_vicini,
                        //atomi di tipo B
                        media_primi_vicini,0,0,
                        0,media_primi_vicini,0,
                        0,0,media_primi_vicini,
                        media_primi_vicini,media_primi_vicini,media_primi_vicini
                       };
        unsigned int base_t [] = {0,0,0,0,1,1,1,1};

        reticolo_inizializza(brav,bas,base_t,8);
        */

        double brav[] = {media_primi_vicini,media_primi_vicini,0,
                         media_primi_vicini,0,media_primi_vicini,
                         0,media_primi_vicini,media_primi_vicini};
        double bas[] = {0,0,0,
                        media_primi_vicini,media_primi_vicini,media_primi_vicini
                       };
        unsigned int base_t [] = {0,1};
        reticolo_inizializza(brav,bas,base_t,2);

        double d2=d2_reticolo_spostamento_medio(min,max,origine_reticolo);

        std::cerr <<"Distanza del primo vicino: "<<media_primi_vicini<<
              "; d2: "<<d2<<"\n";


    }
}
