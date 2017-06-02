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



#include <stdint.h>
#include <iostream>
#include "rnd.h"
#include "interfaccia.h"
#include "cronometro.h"
#include "fit.h"
#include <omp.h>


double zero_=0.0;

///
/// Resetta i contatori interni
template <typename T,typename T2> void History<T,T2>::reset(){
        cur=0;
        full=false;
}

///
/// (Re)inizializza le variabili che sono usate per memorizzare i dati
template <typename T,typename T2> void History<T,T2>::reset_size(unsigned int s, ///< numero massimo di elementi
                                                                 T * ini 		 ///< eventuale puntatore all'elemento a cui devono essere inizializzati gli elementi (per esempio se voglio che inizialmente siano tutti zero). Se non c'è gli elementi non vengono inizializzati a nessun valore.
                                                                 ) {
        if (history!=NULL) delete [] history;
        full=false;
        history=NULL;
        history=new T [s];
        if (history==NULL) std::cerr << "Errore nella reinizializzazione di History ("<<s<<") "<<this<<"\n";
        size=s;
        cur=0;
        if (ini!=NULL)
                for (unsigned int i=0;i<size;i++) history[i]=*ini;
        if (s_y!=NULL) delete [] s_y;
        if (s_x!=NULL) delete [] s_x;
        s_xy_size=0;
        s_x=NULL;
        s_y=NULL;
}
///
/// Stessa funzione di \ref reset_size()
template <typename T,typename T2> History<T,T2>::History(unsigned int s, T * ini) {
        full=false;
        history=NULL;
        history=new T [s];
        if (history==NULL) std::cerr << "Errore nell'inizializzazione di History ("<<s<<") "<<this<<"\n";
        size=s;
        cur=0;
        if (ini!=NULL)
                for (unsigned int i=0;i<size;i++) history[i]=*ini;
        s_xy_size=0;
        s_x=NULL;
        s_y=NULL;

}

///
/// Usando questo costruttore dopo bisogna usare \ref reset_size() per inizializzare correttamente l'oggetto!
template <typename T,typename T2> History<T,T2>::History() {
        full=false;
        history=NULL;
        size=0;
        cur=0;
        s_xy_size=0;
        s_x=NULL;
        s_y=NULL;

}

///
/// libera lo spazio
template <typename T,typename T2> History<T,T2>::~History(){
        delete [] history ;
}

///
/// ritorna vero se ho esaurito lo spazio (e quindi sto scrivendo sopra i dati più vecchi)
template <typename T,typename T2> bool History<T,T2>::is_full(){
        return full;
}


///
/// numero di dati scritti
template <typename T,typename T2> unsigned int History<T,T2>::get_n(){
        if (full){
                return size;
        }else {
                return cur;
        }
}

///
/// aggiunge un nuovo elemento in coda. Se ho finito gli spazi disponibili scrivo sopra quello più vecchio.
template <typename T,typename T2> void History<T,T2>::add(T nuovo) {
        if (size>0){
                history[cur%size]=nuovo;
        }
        cur++;
        if (!full && cur >= size) {
                full=true;
        }
        media_var_fatto=false;
}

///
/// Restituisce un riferimento a un elemento.
template <typename T,typename T2> T & History<T,T2>::operator [] (unsigned int n ///< 0 è l'elemento più vecchio, in ordine di entrata
                                                                                                                                  ) {
        if (size==0) std::cerr << "Richiesto un elemento a un oggetto History senza elementi ("<<this<<")\n";
        if (full) {
                return history[(cur+n)%size];
        } else {
                return history[n%size];
        }
}

/** calcola, se non sono già state calcolate, media e varianza del campione memorizzato
 * e le scrivono nelle variabili fornite
* */
template <typename T,typename T2> void History<T,T2>::media_var(T &media_, ///< qui scrivo la media
                                                                T2 & var_, ///< qui scrivo la varianza
                                                                unsigned int ini_, ///< primo elemento da cui inizio a calcolare media e varianza.
                                                                unsigned int fine_ ///< calcolo la media e la varianza fino al campione precedente a questo (l'elemento etichettato con questo valore non viene considerato). Se vale zero considero tutti gli elementi inseriti (tiene anche conto del caso in cui l'oggetto non è completamente pieno).
                                                                ){
        if (!media_var_fatto || fine_!=media_fine || ini_!=media_ini){
                media_var_calc(media,var,ini_,fine_,true);
                media_ini=ini_;
                media_fine=fine_;
                media_var_fatto=true;
        }
        media_=media;
        var_=var;
}

template <typename T,typename T2> void History<T,T2>::media_var_calc(T &media_, ///< qui scrivo la media
                                                                                                                                         T2 & var_, ///< qui scrivo la varianza
                                                                                                                                         unsigned int ini, ///< primo elemento da cui inizio a calcolare media e varianza.
                                                                     unsigned int fine_, ///< calcolo la media e la varianza fino al campione precedente a questo (l'elemento etichettato con questo valore non viene considerato). Se vale zero considero tutti gli elementi inseriti (tiene anche conto del caso in cui l'oggetto non è completamente pieno).
                                                                     bool calc_var ///< se impostato a true, calcola anche la varianza
                                                                     ){
        unsigned int fine,cont=0;
    if (fine_!=0 && fine_ <= size) {
        std::cerr << "Attenzione: richiesta una media di un oggetto History che non finisce alla fine dello spazio allocato per i dati. ("<<ini<<"-"<<fine_<<" oggetto: "<<this<<")\n";
    }
        if (fine_==0) {
                if (full) fine=size; else fine=cur;
        } else if (fine_>size) {
                fine=size;
        } else {
                fine=fine_;
        }

    if (ini>=fine) {
        std::cerr << "Richiesta una media di un oggetto History a partire da un elemento dopo l'ultimo ("<<ini<<"-"<<fine<<" oggetto: "<<this<<"). Dovresti preoccuparti seriamente di questo messaggio!\n";
    }

        media_=0.0;
        var_=0.0;
        if (calc_var){
                for (unsigned int i=ini;i<fine;i++){
                        cont++;
                        T delta= operator [] (i) - media_;
                        media_ = media_ + delta/cont;
                        var_ = var_ + delta*(operator [] (i) - media_);
                }
                var_ = var_/(cont-1);
        }else {
                for (unsigned int i=ini;i<fine;i++){
                        cont++;
                        T delta= operator [] (i) - media_;
                        media_ = media_ + delta/cont;
                }
        }
}



template <typename T,typename T2> void History<T,T2>::media_var_calc_block(T &media_, ///< qui scrivo la media
                                                                     T2 & var_, ///< qui scrivo la varianza
                                                                     const unsigned int & block_size ///< dimensione di un blocco in timesteps
                                                                     ){
    unsigned int fine,cont=0;

    if (full) fine=size; else fine=cur;
    unsigned int n_b=fine/block_size;

    media_=0.0;
    var_=0.0;
    for (unsigned int iblock=0;iblock<n_b;iblock++){
        unsigned int cont2=0;
        T media2_=0.0;
        for (unsigned int i=block_size*iblock;i<block_size*(iblock+1);i++){ // media sul blocco
            T delta2= operator [] (i) - media2_;
            media2_ = media2_ + delta2/(++cont2);
        }
        T delta=media2_-media_;
        media_=media_+delta/(++cont);
        var_ = var_ + delta*(media2_ - media_);
    }
    var_=var_/(cont*(cont-1));
}

/** effettua la block analysis per il calcolo della varianza
*/
template <typename T,typename T2> void History<T,T2>::stat_ineff_calc(unsigned int N_punti){

        //alloca lo spazio richiesto
        if (s_xy_size!=N_punti){
                if (s_y!=NULL) delete [] s_y;
                if (s_x!=NULL) delete [] s_x;
                s_x=NULL;
                s_y=NULL;
                s_x=new unsigned int [N_punti];
                s_y=new T2 [N_punti];
                s_xy_size=N_punti;
        }
        unsigned int fine;
        if (full) fine=size; else fine=cur;
        unsigned int N_max=fine/10; //numero massimo di misure per blocco
        unsigned int N_min=2;       //numero minimo di misure per blocco
        T media_,media_b,var_;
        T var_fittizia;
        // calcola la media e la varianza del campione
                if (!media_var_fatto || 0!=media_fine || 0!=media_ini){
                media_var_calc(media,var,0,0,true);
                media_ini=0;
                media_fine=0;
                media_var_fatto=true;
        }
        for (unsigned int pcont=0;pcont<N_punti;pcont++){
                unsigned int N=(float)(N_max-N_min)*(float)(pcont)/(float)N_punti + N_min;
                unsigned int n_b=fine/N;
                var_=0.0;
                media_=0.0;
                //stima della varianza delle medie dei blocchi
                for (unsigned int j=0;j<n_b;j++){
                        // per ogni blocco calcola la media
                        media_var_calc(media_b,var_fittizia,j*N,(j+1)*N,false);
                        T delta= media_b - media_;
                        media_ = media_ + delta/(j+1);
                        var_ = var_ + delta*(media_b - media_);
                }
                var_=var_/(n_b-1);
                s_x[pcont]=N;
                s_y[pcont]=(T)N*var_/var;
        }
}

/** fa la media degli ultimi n punti per stimare l'inefficenza statistica del campione.
 * Nota: è una stima molto cruda, niente mi garantisce che la stima sia corretta. Sicuramente è una stima per difetto.
*/
template <typename T,typename T2> T2 History<T,T2>::stat_ineff_media_ultimi(unsigned int N_punti){
    if (N_punti>=s_xy_size/2 || N_punti==0) {
        std::cerr << "Questa media per la stima di 's' non può funzionare. Gli errori sono sbagliati, per difetto (oggetto History"<<this<<")\n";
    }
    if (s_xy_size==0) std::cerr << "Non ho ancora fatto l'analisi a blocchi! ( richiesta media per la stima di 's', oggetto History "<<this<<")\n";
    T2 res=0.0;
    for (unsigned int i=s_xy_size-N_punti;i<s_xy_size;i++){
        res+=s_y[i];
    }
    res=res/(T2)N_punti;
    return res;
}

Position_history::Position_history(){
        fitted_f=NULL;
        fit1=0;
        fit2=0;
        J=NULL;
        types=NULL;
        media_traiettoria_skip=1;
        media_sezioni_scatole=true;
                cur_n=0;
                size=0;
                natoms=0;
                pos=NULL;
                vel=NULL;
                t=NULL;
                d=NULL;
        d_var=NULL;
                full=false;
                medie_cont=0;
        idx_A=NULL;
                A=NULL;
        A=new double[FIT_LIN_NPAR];
        A_err=NULL;
        A_err=new double[FIT_LIN_NPAR];
                tau=0.0;
                n_at=NULL;
        n_at_var=NULL;
                max_d=0;
                max_d_sc=0;
                D=0.0;
        m=0.0;
                corr_v=NULL;
        corr_v_var=NULL;
        corr_J_int=NULL;
        corr_J=NULL;
        corr_J_var=NULL;
        corr_v_int=NULL;
        corr_v_int_dif=NULL;
        corr_v_int_ecc=NULL;
}

/** effettua il fit del profilo multiesponenziale ottenuto in teoria (vedi appunti).
 * Prende solo il primo terzo dei dati.
 * */
void Position_history::scatola_exp_fit(){

}

/** copia i risultati del fit A[0] è il termine costante, A[1,2,...,7] sono i termini dell'esponenziale con tau_1 e tau_2
 * */
void Position_history::get_param_scatola_exp_fit(double *A_,double &tau_,double * A__, double & tau__){
    for (unsigned int i=0;i<FIT_LIN_NPAR;i++) {
        A_[i]=A[i];
        A__[i]=A_err[i];
    }
        tau_=tau;
    tau__=tau_err;
}

void Position_history::reset_types(int *types_) {
    if (types!=NULL)
    for (unsigned int i=0;i<natoms;i++) {
        types[i]=types_[i];
    }
    else
        std::cerr << "Errore: richiesto di impostare i tipi degli atomi senza aver prima allocato l'array che li contiene.\n";
}

void Position_history::reset_natoms_size_cell_types(unsigned int n,unsigned int s,double cell_[3],int * types_) {
    reset_natoms_size_cell(n,s,cell_);
    reset_types (types_);

}

void Position_history::reset_natoms_size_cell(unsigned int n,unsigned int s,double cell_[3]){
    fit1=0;
    fit2=0;
        cur_n=0;
        medie_cont=0;
        max_d=0;
        n_autoc_v=0;
    n_autoc_J=0;
        full=false;
        if (n!=natoms || s!=size){
                size=s;
                natoms=n;
                types=new int [natoms];
                pos=new double * [size];
                vel=new double * [size];
                idx_A=new unsigned int * [size];
                d=new double [size];
        d_var=new double [size];
                n_at=new double [size];
        n_at_var=new double [size];
                corr_v=new double [size];
        corr_v_var=new double[size];
        corr_v_int=new double[size];
        corr_v_int_var=new double[size];
        corr_v_int_ecc=new double[size];
        corr_v_int_dif=new double[size];
                t = new double [size];
                for (unsigned int i=0;i<size;i++){
                        pos[i]=new double[natoms*3];
                        vel[i]=new double[natoms*3];
                        idx_A[i]=new unsigned int [natoms];
                }
        }
        for (unsigned int i=0;i<3;i++){
                cell[i]=cell_[i];
        }
}

double Position_history::operator [](unsigned int i){
        return d[i];
}

void Position_history::reset(){
        cur_n=0;
        medie_cont=0;
    full=false;
}

void Position_history::add(double t_, double *J_){
    t[cur_n]=t_;
    for (unsigned int i=0;i<3;i++){
        J[cur_n*3+i]=J_[i];
    }
    cur_n=(cur_n+1)%size;
    if (cur_n==0){
            full=true;
            //media_temporale_pos2();
    }
 }

void Position_history::add(double t_, double * pos_, double * vel_, double *J_){
        for (unsigned int i=0;i<natoms*3;i++){
                pos[cur_n][i]=pos_[i];
                vel[cur_n][i]=vel_[i];
        }
        for (unsigned int i=0;i<natoms;i++){
                idx_A[cur_n][i]=find_atom(i,pos_);
        }
        t[cur_n]=t_;
    if (J_!=NULL) {
        if (J==NULL) std::cerr <<"J==NULL! Adesso probabilmente SIGSEGV. Ciao ciao!\n";
        for (unsigned int i=0;i<3;i++){
            J[cur_n*3+i]=J_[i];
        }
    }
        cur_n=(cur_n+1)%size;
        if (cur_n==0){
                full=true;
                //media_temporale_pos2();
        }
}

double Position_history::rimasti(unsigned int * idxA, unsigned int * idxB){
        unsigned int ris=0;
    if (media_sezioni_scatole){
        for (unsigned int i=0;i<natoms;i++){
            if (idxA[i]==idxB[i]) ris++;
        }
        // diviso 2 per avere una media (ho considerato entrambi gli scompartimenti)
        return (double)ris/2.0;
    } else {
        for (unsigned int i=0;i<natoms;i++){
            if (idxA[i]==1 && 1==idxB[i]) ris++;
        }
        return (double) ris;
    }
}

unsigned int Position_history::destra(unsigned int *idxA){
    unsigned int cont=0;
    for (unsigned int iatom=0;iatom<natoms;iatom++){
        if (idxA[iatom]==1) cont++;
    }
    return cont;
}

void Position_history::hist_natoms_dx(unsigned int *& ist, double &media, double &var, unsigned int &min,unsigned int &max, unsigned int skip_tstep){
    unsigned int s=0;
    max=0;
    min=natoms;
    if (full) s=size; else s=cur_n;
    unsigned int n_s=s/skip_tstep,cont=0;
    unsigned int *n_dx= new unsigned int[n_s];
    media=0.0;
    var=0.0;
    for (unsigned int i=0;i<n_s;i++){
        n_dx[i]=destra(idx_A[skip_tstep*i]);
        if (n_dx[i]<min) min=n_dx[i];
        if (n_dx[i]>max) max=n_dx[i];
        double delta=(double)n_dx[i] - media;
        media +=delta/(++cont);
        var+=delta*((double)n_dx[i] - media);
    }
    var=var/((double)(cont*(cont-1)));

    if (ist!=NULL) delete [] ist;
    ist=new unsigned int [max-min+1];
    for (unsigned int i=0;i<max-min+1;i++) ist[i]=0;
    for (unsigned int i=0;i<n_s;i++) {
        ist[n_dx[i]-min]++;
    }

    delete [] n_dx;
}

/**Trova dov'è l'atomo. Restituisce l'indice del settore corrispondente all'atomo richiesto.
 * */
unsigned int Position_history::find_atom(unsigned int i, ///< indice dell'atomo
                                                                                 double * pos_ ///< posizioni degli atomi
                                     ){
        unsigned int n_parti=2;
        double pb[3];
        for (unsigned int jj=0;jj<3;jj++){
            pb[jj]=pos_[i*3+jj]-round(pos_[i*3+jj]/cell[jj])*cell[jj];
        }
        // la cella è centrata in zero! (per questo aggiungo metà cella per avere solo coordinate
        pb[0]+=cell[0]/2.0;
        for (unsigned int j=0;j<n_parti;j++){
                // guardo solo la prima coordinata per dividere positive)

                if (pb[0] >= cell[0]*(double)j/(double)n_parti
                 && pb[0] <  cell[0]*(double)(j+1)/(double)n_parti){
                         return j;
                }
        }
        std::cerr << "Manca un atomo all'appello: "<<i<<" posizionato in "<<pos_[i*3]<<"(Position_history: "<<this<<")\n";
        return 0;
}


double * Position_history::element(unsigned int j){
        if (!full)
                return pos[j%size];
        else
                return pos[(cur_n+j)%size];
}

unsigned int Position_history::max_(unsigned int d_max){
        unsigned int s;
        if (d_max==0){
                if (full) s=size; else s=cur_n;
        } else {
                if (full) {
                        if (d_max<size)
                                s=d_max;
                        else
                                s=size;
                } else {
                        if (d_max>cur_n)
                                s=cur_n;
                        else
                                s=d_max;
                }
        }
        return s;
}

void Position_history::media_temporale_pos2_(unsigned int d_max){
        double media;
        d[0]=0.0;
        unsigned int ss,s=max_(d_max);
        if (full) ss=size; else ss=cur_n;
        max_d=s;
    if ((t[s]-t[0])<=0.25) {
                std::cerr << "Non dovrei (ma ci provo) fare una media su un tempo così corto:\
 la funzione qui non è ancora lineare! (oggetto Position_history "<<this<<")\n";
        }
        unsigned int n_025=0;
        for (unsigned int j=1;j<s;j++){
                d[j]=0.0;
        if (n_025==0 && (t[j]-t[0]) >=0.25) n_025=j;
        //media sulla traccia
                for (unsigned int i=0;i<ss-j;i++) {
                        media = 0.0;
            //media sugli atomi
                        for (unsigned int at=0;at<natoms*3;at++){
                                media+=(element(i) [at]-element(i+j)[at])*(element(i) [at]-element(i+j)[at]);
                        }
                        media=media/(natoms-1);
                        d[j]=d[j]+media;
                }
                d[j]=d[j]/(double)(ss-j);
        }
        // scarta i primi 0.25/tstep passi
        ols (s-n_025, &t[n_025], &d[n_025], &m, &q, &em, &eq, &chi2);

}

void Position_history::set_media_traiettoria_skip(unsigned int n_){
    media_traiettoria_skip=n_;
    if (n_>=size) {
        std::cerr << "Attenzione: non sto facendo la media sulla traiettoria.\n";
    }
}

void Position_history::set_scatola_media_sezioni(bool n_){
    media_sezioni_scatole = n_;
}

void Position_history::media_temporale_pos2(unsigned int d_max){
    d[0]=0.0;
    unsigned int ss,s=max_(d_max);
    if (full) ss=size; else ss=cur_n;
    max_d=s;
    //trasla l'origine dei tempi
    double *x = new double [max_d];
    unsigned int i0=0;
    if (full) i0=cur_n;
    unsigned int n_025=0;
    for (unsigned int i=0;i<max_d;i++){
        x[i]=t[i]-t[i0];
        if (n_025==0 && x[i] >=1.5) n_025=i;
    }
    if ((x[s-1])<=1.5) {
        std::cerr << "Non dovrei (ma ci provo) fare una media su un tempo così corto (x["<<s-1<<"] = "<<x[s-1]<<"): la funzione qui non è ancora lineare! (oggetto Position_history "<<this<<")\n";
    }
    unsigned int n_b=(ss-s)/s;
    if (n_b==0) std::cerr <<"Attenzione: non posso fare l'analisi a blocchi della lunghezza richiesta!\n";
    double *d_t=new double[n_b*s];
    double *M=new double[n_b],*M_var=new double[n_b],*Q=new double[n_b],*Q_var=new double[n_b],*C=new double[n_b];
    cronometro cron;
    cron.set_print_interval(10.0);
    cron.set_expected(1.0/((double)n_b*s));
    std::cerr << "Inizio del calcolo del MSD (Einstein) con analisi a blocchi...\n";
    cron.start();
#pragma omp parallel for
    for (unsigned int iblock=0;iblock<n_b;iblock++){
        for (unsigned int t=0;t<s;t++){
            unsigned int cont=0;
            d_t[iblock*s+t]=0.0;
            //1. media sulla traiettoria
            for (unsigned int i=iblock*s;i<(iblock+1)*s;i=i+media_traiettoria_skip){
                //2. media sugli atomi,
                for (unsigned int iat=0;iat<natoms*3;iat++){
                    double delta1=  (element(i) [iat]-element(i+t)[iat])*(element(i) [iat]-element(i+t)[iat])- d_t[iblock*s+t];
                    d_t[iblock*s+t] = d_t[iblock*s+t] + delta1/(++cont);
                }
            }
            cron.stop(omp_get_thread_num());
            cron.print_interval(omp_get_thread_num());
        }

    }
    std::cerr << "Tempo di calcolo: " << cron.time()<<"s\nEseguo i fit lineari, calcolo le medie e le varianze...\n";
    cron.reset();
    cron.start();
    unsigned int fit_ini=0,fit_end=max_d-1;
    if (fit2 !=0){ // quindi l'utente ha richiesto uno specifico intervallo
        if (fit1<max_d && fit1!=fit2){
            fit_ini=fit1;
        } else {
            std::cerr << "Il range di fit richiesto e' fuori dal blocco (inizio a "<<fit1<<" richiesto, l'ultimo timestep e' "<<max_d-1<<"). Ignorata la richiesta.\n";
        }

        if (fit1<max_d && fit1!=fit2){
            fit_end=fit2;
        }else {
            std::cerr << "Il range di fit richiesto e' fuori dal blocco (fine a "<<fit2<<" richiesta, l'ultimo timestep e' "<<max_d-1<<"). Ignorata la richiesta.\n";
        }
    }
    if (fit_end < fit_ini){
        std::cerr << "Il range di fit richiesto e' al contrario. (da "<<fit1<<" a "<<fit2<<") Adesso lo inverto.\n";
        unsigned int fit_t_=fit_ini;
        fit_ini=fit_end;
        fit_end=fit_t_;
    }
    for (unsigned int t=1;t<s;t++){
        unsigned int cont=0;
        d[t]=0.0;
        d_var[t]=0.0;
        for (unsigned int iblock=0;iblock<n_b;iblock++){
            double delta1= d_t[iblock*s+t] - d[t];
            d[t]+= delta1/(++cont);
            d_var[t]+=delta1*(d_t[iblock*s+t] - d[t]);
        }
        d_var[t]/=((cont-1)*cont);
    }
    m=0.0;
    em=0.0;
    q=0.0;
    eq=0.0;
    chi2=0.0;
    unsigned int cont=0;
    for (unsigned int iblock=0;iblock<n_b;iblock++){
        ols (fit_end-fit_ini, &x[fit_ini], &d_t[iblock*s+fit_ini], &M[iblock], &Q[iblock], &M_var[iblock], &Q_var[iblock], &C[iblock]);
        double delta1= M[iblock] - m;
        m+=delta1/(++cont);
        em+=delta1*(M[iblock] - m);
        double delta2= Q[iblock] - q;
        q+=delta2/(cont);
        eq+=delta2*(Q[iblock] - q);
    }
    em=sqrt(em/(double)((cont-1)*cont));
    eq=sqrt(eq/(double)((cont-1)*cont));
    std::cerr << "Einstein fatto.\n";
    delete [] d_t;
    delete [] M;
    delete [] M_var;
    delete [] Q;
    delete [] Q_var;
    delete [] C;
    delete [] x;
}

void Position_history::media_temporale_scatole(unsigned int d_max){
    unsigned int ss,s=max_(d_max);
    max_d_sc=s;
    if (full) ss=size; else ss=cur_n;
    unsigned int n_b=(ss-s)/s;
    if (n_b==0) std::cerr <<"Attenzione: non posso fare l'analisi a blocchi della lunghezza richiesta!\n";
    double *n_at_b=new double[s*n_b];
    double *x = new double [max_d_sc];
    unsigned int i0=0;
    //trasla l'origine dei tempi
    if (full) i0=cur_n;
    for (unsigned int i=0;i<max_d_sc;i++){
        x[i]=t[i]-t[i0];
    }
    if (s<40) std::cerr << "Attenzione: numero di punti troppo basso ("<<s<<")\n";
    cronometro cron;
    cron.set_print_interval(10.0);
    cron.set_expected(1.0/((double)n_b*s));
    std::cerr << "Inizio del calcolo della diffusione degli atomi fra le parti destra e sinistra della simulazione con analisi a blocchi...\n";
    cron.start();
#pragma omp parallel for
    for (unsigned int iblock=0;iblock<n_b;iblock++){
        n_at_b[iblock*s]=(double) natoms/2.0;
        for (unsigned int t=1;t<s;t++){
            unsigned int cont=0;
            n_at_b[iblock*s+t]=0.0;
            //1. media sulla traiettoria
            for (unsigned int i=iblock*s;i<(iblock+1)*s;i=i+media_traiettoria_skip){
                double delta=rimasti(idx_A[i],idx_A[i+t])-n_at_b[iblock*s+t];
                n_at_b[iblock*s+t]+=delta/(++cont);
            }
        }
        cron.stop(omp_get_thread_num());
        cron.print_interval(omp_get_thread_num());
    }
    std::cerr << "Tempo di calcolo: " << cron.time()<<"s\nEseguo i fit e calcolo le medie e le varianze...\n";
    cron.reset();
    cron.set_expected(1.0/(float)n_b);
    cron.start();
    n_at[0]=(double)natoms/2.0;
    n_at_var[0]=0.0;
    for (unsigned int t=1;t<s;t++){
        unsigned int cont=0;
        n_at[t]=0.0;
        n_at_var[t]=0.0;
        for (unsigned int iblock=0;iblock<n_b;iblock++){
            double delta1= n_at_b[iblock*s+t] - n_at[t];
            n_at[t]+= delta1/(++cont);
            n_at_var[t]+=delta1*(n_at_b[iblock*s+t] - n_at[t]);
        }
        n_at_var[t]/=((cont-1)*cont);
    }
    double *A_b=new double[n_b*FIT_LIN_NPAR];
    double *tau_b=new double[n_b];
    for (unsigned int j=0;j<FIT_LIN_NPAR;j++) {
        A[j]=0.0;
        A_err[j]=0.0;
    }
    tau=0.0;
    tau_err=0.0;
    unsigned int cont=0;
    unsigned int fit_ini=0,fit_end=max_d_sc-1;
    if (fit2 !=0){ // quindi l'utente ha richiesto uno specifico intervallo
        if (fit1<max_d_sc && fit1!=fit2){
            fit_ini=fit1;
        } else {
            std::cerr << "Il range di fit richiesto e' fuori dal blocco (inizio a "<<fit1<<" richiesto, l'ultimo timestep e' "<<max_d_sc-1<<"). Ignorata la richiesta.\n";
        }

        if (fit1<max_d_sc && fit1!=fit2){
            fit_end=fit2;
        }else {
            std::cerr << "Il range di fit richiesto e' fuori dal blocco (fine a "<<fit2<<" richiesta, l'ultimo timestep e' "<<max_d_sc-1<<"). Ignorata la richiesta.\n";
        }
    }
    if (fit_end < fit_ini){
        std::cerr << "Il range di fit richiesto e' al contrario. (da "<<fit1<<" a "<<fit2<<") Adesso lo inverto.\n";
        unsigned int fit_t_=fit_ini;
        fit_ini=fit_end;
        fit_end=fit_t_;
    }
#pragma omp parallel for
    for (unsigned int iblock=0;iblock<n_b;iblock++){
        if (m==0.0) m=0.066;
        tau_b[iblock]=(cell[0]/(2*PI))*(cell[0]/(2*PI))/(m/2.0);
        fit<FIT_LIN_NPAR>::fit_lin_exp(fit_end-fit_ini,&x[fit_ini],&n_at_b[iblock*s+fit_ini],&A_b[iblock*FIT_LIN_NPAR],tau_b[iblock],true,natoms);
        double delta2=tau_b[iblock]-tau;
        tau+=delta2/(++cont);
        tau_err+=delta2*(tau_b[iblock]-tau);
        for (unsigned int j=0;j<FIT_LIN_NPAR;j++){
            double delta1=A_b[iblock*FIT_LIN_NPAR+j]-A[j];
            A[j]+=delta1/cont;
            A_err[j]+=delta1*(A_b[iblock*FIT_LIN_NPAR+j]-A[j]);
        }
        cron.stop(omp_get_thread_num());
        cron.print_interval(omp_get_thread_num());
    }
    tau_err=sqrt(tau_err/((cont-1)*cont));
    delete [] fitted_f;
    fitted_f=new double [max_d_sc];
    fit<FIT_LIN_NPAR>::get_fit_func(max_d_sc,x,fitted_f,A,tau);
    for (unsigned int j=0;j<FIT_LIN_NPAR;j++) A_err[j]=sqrt(A_err[j]/((cont-1)*cont));
    std::cerr<<"Tempo trascorso: "<<cron.time()<<"s. Scatole fatte.\n";
    delete [] x;
    delete [] A_b;
    delete [] tau_b;
    delete [] n_at_b;
}

void Position_history::media_temporale_scatole_set_fit_range(unsigned int f1, unsigned int f2){
    if (f1>=f2){
        std::cerr <<"Specificato un range di fit non valido! (da "<<f1<<" a "<<f2<<") Ignoro la richiesta\n";
    } else {
        fit1=f1;
        fit2=f2;
    }
}

void Position_history::media_temporale_scatole2(unsigned int d_max){
        unsigned int ss,s=max_(d_max);
        max_d_sc=s;
        if (full) ss=size; else ss=cur_n;
        n_at[0]=(double) natoms / 2.0;
        for (unsigned int j=1;j<s;j++){
                n_at[j]=0.0;
                for (unsigned int i=0;i<ss-j;i++) {
                        n_at[j]=n_at[j]+rimasti(idx_A[i],idx_A[i+j]);
                }
                n_at[j]=n_at[j]/(double)(ss-j);
        }

}

/** Calcola la funzione di autocorrelazione delle velocità e il suo integrale, tutto con i relativi errori
 *  Computazionalmente pesante! Da snellire un po'.
*/
void Position_history::calc_v_autocorrelation_old(unsigned int d_max) {
        unsigned int ss,s=max_(d_max);
        n_autoc_v=s;
        if (full) ss=size; else ss=cur_n;
    // per stimare la varianza
    History<double,double> block_an(ss,&zero_);
        D=0.0;
    double D_dif=0.0,D_ecc=0.0;
    double d1,d2,non_usato,var;
    // integrale della media sulla traiettoria di <v(i)*v(i+j)>
        for (unsigned int j=0;j<s;j++){
        d2=0.0;
        // media sulla traiettoria di <v(i)*v(i+j)>
        if (j%10==0) block_an.reset();
                for (unsigned int i=0;i<ss-s;i++){
            d1=0.0;
                        //media sugli atomi del prodotto v(i)*v(i+j)
            for (unsigned int iat=0;iat<natoms*3;iat++){
                double delta1= vel[i][iat]*vel[i+j][iat] - d1;
                                d1 = d1 + delta1/(iat+1);
                        }
            if (j%10==0) block_an.add(d1);
            double delta2=d1-d2;
                        d2 = d2 + delta2/(i+1);
                }
        /** Ogni tanto fa la block analysis per stimare l'errore sulla funzione di autocorrelazione, e alla fine l'errore sull'integrale.
        * Quest'ultimo può essere ricavato integrando la funzione di correlazione con tutti gli errori statistici
        * prima sommati e poi sottratti (così ho una valore 'massimo' e uno 'minimo' dell'integrale).
        * Non posso sommare le varianze dei vari termini della somma che approssima l'integrale perchè due termini
        * successivi sono chiaramente fortemente correlati.
        */
        if (j%10==0) {
            block_an.stat_ineff_calc(50);
            block_an.media_var(non_usato,var);
            corr_v_var[j]=block_an.stat_ineff_media_ultimi(5)*var/(double)block_an.get_n();
        }else {
            corr_v_var[j]=corr_v_var[j-1];
        }
                corr_v[j]=d2;
        D+=d2;
        D_dif+=d2-sqrt(corr_v_var[j]);
        D_ecc+=d2+sqrt(corr_v_var[j]);
        corr_v_int_dif[j]=D_dif;
        corr_v_int_ecc[j]=D_ecc;
        corr_v_int[j]=D;
    }
}

void Position_history::calc_v_autocorrelation(unsigned int d_max, unsigned int err_last){
    unsigned int ss,s=max_(d_max);
    n_autoc_v=s;
    if (full) ss=size; else ss=cur_n;
    unsigned int n_b=(ss-s)/s; // per fare in modo che tutti i punti hanno la stessa statistica tolgo la parte finale della traccia
    //divido la traccia in n_b blocchi di lunghezza s, e su ognuno calcolo la funzione di autocorrelazione e il suo integrale
    double *c_t = new double[s*n_b];
    double *c_t_int_T = new double[n_b*s];
    if (n_b==0) std::cerr <<"Attenzione: non posso fare l'analisi a blocchi della lunghezza richiesta!\n";
    cronometro cron;
    cron.set_print_interval(10.0);
    cron.set_expected(1.0/((double)n_b*s));
    std::cerr << "Inizio del calcolo della funzione di autocorrelazione delle velocità e relativo integrale con analisi a blocchi...\n";
    cron.start();
#pragma omp parallel for
    for (unsigned int iblock=0;iblock<n_b;iblock++){
        // integrale di c(t)
        c_t_int_T[(iblock+1)*s-1]=0.0; // l'ultimo valore è l'integrale totale
        for (unsigned int t=0;t<s;t++){
            //calcola prima c(t):
            c_t[iblock*s+t]=0.0;
            unsigned int cont=0;
            //1. media sulla traiettoria
            for (unsigned int i=iblock*s;i<(iblock+1)*s;i=i+media_traiettoria_skip){
                //2. media sugli atomi, c(t)
                for (unsigned int iat=0;iat<natoms*3;iat++){
                    double delta1= vel[i][iat]*vel[i+t][iat] - c_t[iblock*s+t];
                    c_t[iblock*s+t] = c_t[iblock*s+t] + delta1/(++cont);
                }
            }
            // integrale con i trapezi
            if (t==0 || t==(s-1)){
                c_t_int_T[(iblock+1)*s-1]+=c_t[iblock*s+t]/2.0;
            }else {
                c_t_int_T[(iblock+1)*s-1]+=c_t[iblock*s+t];
            }
            c_t_int_T[iblock*s+t]=c_t_int_T[(iblock+1)*s-1];
            cron.stop(omp_get_thread_num());
            cron.print_interval(omp_get_thread_num());
        }
    }
    std::cerr << "Tempo di calcolo: " << cron.time()<<"s\nCalcolo le medie e le varianze...";
    cron.reset();
    cron.set_expected(1.0/(double)s);
    cron.start();
    //adesso ho varie funzioni di correlazione e integrali, posso calcolare media e varianza per ogni punto
    for (unsigned int t=0;t<s;t++){
        corr_v[t]=0.0;
        corr_v_var[t]=0.0;
        corr_v_int[t]=0.0;
        corr_v_int_var[t]=0.0;
        unsigned int cont=0;
        for (unsigned int iblock=0;iblock<n_b;iblock++){
            //c(t)
            double delta=c_t[iblock*s+t] - corr_v[t];
            corr_v[t]+=delta/(++cont);
            corr_v_var[t]+=delta*(c_t[iblock*s+t]-corr_v[t]);
            //l'integrale
            double delta1=c_t_int_T[iblock*s+t] - corr_v_int[t];
            corr_v_int[t]+=delta1/(cont);
            corr_v_int_var[t]+=delta1*(c_t_int_T[iblock*s+t]-corr_v_int[t]);
        }
        corr_v_var[t]/=((n_b-1)*n_b);
        corr_v_int_var[t]/=((n_b-1)*n_b);
        corr_v_int_ecc[t]=corr_v_int[t]+sqrt(corr_v_int_var[t]);
        corr_v_int_dif[t]=corr_v_int[t]-sqrt(corr_v_int_var[t]);
        cron.stop();
        cron.print_interval();
    }
    D=corr_v_int[s-1];
    D_err=0.0;
    unsigned int cont=0;
    for (unsigned int i=s-err_last;i<s;i++){
        double delta=corr_v_int_var[i]-D_err;
        D_err+=delta/(++cont);
    }
    D_err=sqrt(D_err);
    std::cerr << "Tempo di calcolo: " << cron.time()<<"s\nGreen-Kubo velocità fatto!\n";
    delete [] c_t;
    delete [] c_t_int_T;
}

/**
  * calcola la funzione di autocorrelazione di J (nel tempo) e il relativo integrale usando l'analisi a blocchi.
*/
void Position_history::calc_J_autocorrelation(unsigned int d_max /// lunghezza dei blocchi che uso per l'analisi
                                              ,unsigned int err_last /// per calcolare l'errore faccio la media dell'errore per la coda dell'integrale della funzione di autocorrelazione. Questo numero dice quanti timestep prima della fine del blocco inizio a farela media per l'errore finale.
                                              ){
    if (J==NULL) {
        std::cerr << "J==NULL true; SIGSEGV a breve. Ciao ciao! (forse era meglio prima chiamare allocate_J(), e caricare J con un adeguato file log di lammps)\n";
    }
    unsigned int ss,s=max_(d_max);
    n_autoc_J=s;
    if (full) ss=size; else ss=cur_n;
    unsigned int n_b=(ss-s)/s; // per fare in modo che tutti i punti hanno la stessa statistica tolgo la parte finale della traccia
    //divido la traccia in n_b blocchi di lunghezza s, e su ognuno calcolo la funzione di autocorrelazione e il suo integrale
    double *c_t = new double[s*n_b];
    double *c_t_int_T = new double[n_b*s];
    if (n_b==0) std::cerr <<"Attenzione: non posso fare l'analisi a blocchi della lunghezza richiesta!\n";
    cronometro cron;
    cron.set_print_interval(10.0);
    cron.set_expected(1.0/((double)n_b*s));
    std::cerr << "Inizio del calcolo della funzione di autocorrelazione di J e relativo integrale con analisi a blocchi...\n";
    cron.start();
//#pragma omp parallel for
    for (unsigned int iblock=0;iblock<n_b;iblock++){
        // integrale di c(t)
        c_t_int_T[(iblock+1)*s-1]=0.0; // l'ultimo valore è l'integrale totale
        for (unsigned int t=0;t<s;t++){
            //calcola prima c(t):
            c_t[iblock*s+t]=0.0;
            unsigned int cont=0;
            //1. media sulla traiettoria
            for (unsigned int i=iblock*s;i<(iblock+1)*s;i=i+media_traiettoria_skip){
                for (unsigned int icoord=0;icoord<3;icoord++){
                    double delta1= J[i*3+icoord]*J[(i+t)*3+icoord] - c_t[iblock*s+t];
                    c_t[iblock*s+t] = c_t[iblock*s+t] + delta1/(++cont);
                }
            }
            // integrale con i trapezi
            if (t==0 || t==(s-1)){
                c_t_int_T[(iblock+1)*s-1]+=c_t[iblock*s+t]/2.0;
            }else {
                c_t_int_T[(iblock+1)*s-1]+=c_t[iblock*s+t];
            }
            c_t_int_T[iblock*s+t]=c_t_int_T[(iblock+1)*s-1];
            cron.stop(omp_get_thread_num());
            cron.print_interval(omp_get_thread_num());
        }
    }
    std::cerr << "Tempo di calcolo: " << cron.time()<<"s\nCalcolo le medie e le varianze...";
    cron.reset();
    cron.set_expected(1.0/(double)s);
    cron.start();
    //adesso ho varie funzioni di correlazione e integrali, posso calcolare media e varianza per ogni punto
    for (unsigned int t=0;t<s;t++){
        corr_J[t]=0.0;
        corr_J_var[t]=0.0;
        corr_J_int[t]=0.0;
        corr_J_int_var[t]=0.0;
        unsigned int cont=0;
        for (unsigned int iblock=0;iblock<n_b;iblock++){
            //c(t)
            double delta=c_t[iblock*s+t] - corr_J[t];
            corr_J[t]+=delta/(++cont);
            corr_J_var[t]+=delta*(c_t[iblock*s+t]-corr_J[t]);
            //l'integrale
            double delta1=c_t_int_T[iblock*s+t] - corr_J_int[t];
            corr_J_int[t]+=delta1/(cont);
            corr_J_int_var[t]+=delta1*(c_t_int_T[iblock*s+t]-corr_J_int[t]);
        }
        corr_J_var[t]/=((n_b-1)*n_b);
        corr_J_int_var[t]/=((n_b-1)*n_b);
        cron.stop();
        cron.print_interval();
    }
    DJ=corr_J_int[s-1];
    DJ_err=0.0;
    unsigned int cont=0;
    for (unsigned int i=s-err_last;i<s;i++){
        double delta=corr_J_int_var[i]-DJ_err;
        DJ_err+=delta/(++cont);
    }
    DJ_err=sqrt(DJ_err);
    std::cerr << "Tempo di calcolo: " << cron.time()<<"s\nGreen-Kubo - calore fatto!\n";
    delete [] c_t;
    delete [] c_t_int_T;
}
void Position_history::get_m_q_em_eq_chi2(double &M,double &Q,double &EM,double &EQ,double &CHI2){
        M=m;
        Q=q;
        EM=em;
        EQ=eq;
        CHI2=chi2;
}

void Position_history::allocate_J(){
    if (J!=NULL) delete [] J;
    J=new double [3*size];

    corr_J=new double [size];
    corr_J_var=new double[size];
    corr_J_int=new double[size];
    corr_J_int_var=new double[size];
}

void Position_history::allocate_J(unsigned int size_){
    size=size_;
    delete [] t;
    t = new double [size];
    allocate_J();
}


Position_history::~Position_history(){
    if (J!=NULL) delete [] J;
    if (pos!=0){
        for (unsigned int i=0;i<size;i++){
            delete [] pos[i];
        }
        delete [] pos;
    }
    if (vel!=0){
        for (unsigned int i=0;i<size;i++){
            delete [] vel[i];
        }
        delete [] vel;
    }
    if (idx_A!=0){
        for (unsigned int i=0;i<size;i++){
            delete [] idx_A[i];
        }
        delete [] idx_A;
    }
    if (d!= NULL) {
        delete [] d;
        delete [] d_var;
    }
    if (n_at != NULL){
        delete [] n_at;
        delete [] n_at_var;
    }
    if (corr_v!=NULL){
        delete [] corr_v;
        delete [] corr_v_var;
        delete [] corr_v_int;
        delete [] corr_v_int_var;
        delete [] corr_v_int_dif;
        delete [] corr_v_int_ecc;
    }
    delete [] corr_J;
    delete [] corr_J_var;
    delete [] corr_J_int;
    delete [] corr_J_int_var;
    if (t!=NULL) delete [] t;
}

Analisi::Analisi(std::string traiettoria_fname,std::string energia_fname,unsigned int size){
    double *positions=NULL,*velocities=NULL,J[3];
    int *types=NULL;
    FILE *fp = fopen(traiettoria_fname.c_str(),"rb");
    std::string head = "Step Time Temp E_pair E_mol TotEng Press \n";
//    std::string head_heat = "Step Time Temp E_pair E_mol TotEng Press Jx Jy Jz \n";
    std::string head_heat =   "Step Time PotEng TotEng Lx Press Temp c_flusso[1] c_flusso[2] c_flusso[3] \n";
    std::string tmp;
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    FILE *fp2 = fopen(energia_fname.c_str(),"r");
    cronometro cron;
    cron.set_expected(1.0/(double)size);
    cron.set_print_interval(10.0);
    std::cerr << "Inizio lettura della traiettoria (puo' richiedere del tempo)...\n";
    cron.start();
    time_history.reset_size(size,&zero_);
    T_history.reset_size(size,&zero_);
    pressure_history.reset_size(size,&zero_);
    kin_energy_history.reset_size(size,&zero_);
    pot_energy_history.reset_size(size,&zero_);


    if (!fp) {
        std::cerr << "impossibile aprire il file della traiettoria " << traiettoria_fname << "\n";
    }
    if (!fp2) {
        std::cerr << "impossibile aprire il file di log" << energia_fname << "\n";
    }

    bool fail = true;
    while ((read = getline(&line, &len, fp2)) != -1){
        if (feof(fp2)) {
            fclose(fp2);
            break;
        }
        tmp.assign(line);
        if (tmp==head){
            fail=false;
            heat_autocorr=false;
            break;
        } else if (tmp==head_heat){
            fail=false;
            heat_autocorr=true;
            break;
        }// sono arrivato all'inizio dei dati
    }

    if(fail) {
        std::cerr << "Impossibile trovare i dati termodinamici della simulazione con il formato corretto!\n";
    }

    bool first=true;
    bigint nntimestep;
    int error_cont=0;
    size_t st;
    while (1) {
        error_cont=0;
        if(fp && fread(&nntimestep,sizeof(bigint),1,fp)!=1) {std::cerr << "Errore nella lettura del numero di timestep corrente! Termino subito\n";break;} //ntimestep indica il numero del timestep che adesso il programma leggerà
#ifdef DEBUG_
        std::cerr << "Timestep n° "<< ntimestep << "\n";
#endif
        // detect end-of-file

        if ((fp && feof(fp))||feof(fp2)) {
            std::cerr << "Fine dei file raggiunta (timestep n° "<<ntimestep<<")\n";
            fclose(fp);
            fclose(fp2);
            break;
        }
        //"Step Time Temp E_pair E_mol TotEng Press"//
        //std::string head_heat = "  Step Time Temp   E_pair E_mol TotEng Press Jx          Jy          Jz \n";
        //std::string head_heat =   "Step Time PotEng TotEng Lx    Press  Temp  c_flusso[1] c_flusso[2] c_flusso[3] \n";

        int Step;
        double Time, Temp, E_pair, E_mol, TotEng, Press;
        int res=0;
        if (heat_autocorr){
            res=fscanf(fp2,"%i %lf %lf %lf %lf %lf %lf %lf %lf %lf",&Step,&Time, &E_pair, &TotEng, &E_mol, &Press, &Temp, &J[0], &J[1], &J[2]);
        } else {
            res=fscanf(fp2,"%i %lf %lf %lf %lf %lf %lf",&Step,&Time, &Temp, &E_pair, &E_mol, &TotEng, &Press);
        }
        if (Step < nntimestep){ //leggo ancora uno
            if (heat_autocorr){
                res=fscanf(fp2,"%i %lf %lf %lf %lf %lf %lf %lf %lf %lf",&Step,&Time, &E_pair, &TotEng, &E_mol, &Press, &Temp, &J[0], &J[1], &J[2]);
            } else {
                res=fscanf(fp2,"%i %lf %lf %lf %lf %lf %lf",&Step,&Time, &Temp, &E_pair, &E_mol, &TotEng, &Press);
            }
        }
#ifdef DEBUG_
        std::cerr << head <<"\n" << Step << " " << Time << " " <<  Temp << " " <<  E_pair << " " <<  E_mol << " " <<  TotEng << " " <<  Press << "\n";
#endif
        if (feof(fp2)) {
            std::cerr << "Fine del file di log raggiunta\n";
            if(fp) fclose(fp);
            fclose(fp2);
        }
        if (res != 7 && (!heat_autocorr)){
            std::cerr << "Numero sbagliato di dati termodinamici ("<< res<<" al posto di 7) letti dal file " << energia_fname << "\n usare \"thermo_style custom step time temp epair emol etotal press\" in lammps\n";
            if(fp) fclose(fp);
            fclose(fp2);
            break;
        } else if (res != 10 && (heat_autocorr)){
            std::cerr << "Numero sbagliato di dati termodinamici ("<< res<<" al posto di 10) letti dal file " << energia_fname << "\n usare \"thermo_style custom step time temp epair emol etotal press\" in lammps\n";
            if(fp) fclose(fp);
            fclose(fp2);
            break;
        }

        if (fp) while (nntimestep<Step) {
            if(fread(&natoms_,sizeof(bigint),1,fp)!=1) {std::cerr << "Errore nella lettura del numero di atomi! Termino subito\n";break;}
            if(fread(&triclinic_,sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura del tipo di cella! Termino subito\n";break;}
            fseek(fp,sizeof(int)*6+sizeof(double)*6,SEEK_CUR); // intestazione del frame
            if (triclinic_) {
                fseek(fp,sizeof(double)*3,SEEK_CUR);
            }
            if (ferror(fp)!=0) perror("Errore:");
            if(fread(&size_one,sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura del numero di dati per atomo! Termino subito\n";break;}
            if(fread(&nchunk,sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura del numero di pezzi per processore! Termino subito\n";break;}
            if (ferror(fp)!=0) perror("Errore:");
            for (unsigned int i = 0; i < nchunk; i++) {
                int n=0;
                if(fread(&n,sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura del numero di dati per pezzo "<<i<<" ! Termino subito\n";break;}
                if (ferror(fp)!=0) perror("Errore:");
                fseek(fp,sizeof(double)*n,SEEK_CUR);
            }
            if(fread(&nntimestep,sizeof(bigint),1,fp)!=1) {std::cerr << "Errore nella lettura del numero di timestep corrente! Termino subito\n";break;} //ntimestep indica il numero del timestep che adesso il programma leggerà
        }

        if (fp && Step!=nntimestep){
            std::cerr << "Il timestep non corrisponde fra la traiettoria e il file log! ("<<nntimestep <<" e " << Step <<")\n";
            error_cont++;
        }
        T_history.add(Temp);
        time_history.add(Time);
        pressure_history.add(Press);
        pot_energy_history.add(E_pair);
        kin_energy_history.add(TotEng);

        if(fp){
            if(fread(&natoms_,sizeof(bigint),1,fp)!=1) {std::cerr << "Errore nella lettura del numero di atomi! Termino subito\n";break;}
            if(fread(&triclinic_,sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            if(fread(&boundary_[0][0],6*sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            if(fread(&xlo_,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            if(fread(&xhi_,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            if(fread(&ylo_,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            if(fread(&yhi_,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            if(fread(&zlo_,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            if(fread(&zhi_,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            cell_[0]=xhi_-xlo_;
            cell_[1]=yhi_-ylo_;
            cell_[2]=zhi_-zlo_;

#ifdef DEBUG_
            std::cerr <<"natoms: "<< natoms_ << " triclinic: " << triclinic_ << " boundary[3][2] & cell[3]: ";
            for (unsigned int jj=0;jj<3;jj++) for (unsigned int kk=0;kk<2;kk++) std::cerr << boundary_[jj][kk] <<" ";
            for (unsigned int jj=0;jj<3;jj++) std::cerr << cell_[jj]<<" ";
            std::cerr <<"\n";
#endif
            if (triclinic_) {
                if(fread(&xy,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
                if(fread(&xz,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
                if(fread(&yz,sizeof(double),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            }
            if(fread(&size_one,sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
            if(fread(&nchunk,sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
#ifdef DEBUG_
            std::cerr << "size_one: " <<size_one << " nchunk: " << nchunk << "\n";
#endif
            int atomc=0;
            if (first){
                traiettoria.reset_natoms_size_cell(natoms_,size,cell_);
                traiettoria.allocate_J();
                positions=new double[natoms_*3];
                velocities=new double[natoms_*3];
                types=new int[natoms_];
                for (unsigned int i=0;i<natoms_;i++) {
                    types[i]=0;
                }
            }

            // loop over processor chunks in file

            int n=0,maxbuf=0;
            double *buf=NULL;
            for (unsigned int i=0;i<natoms_*3;i++){
                positions[i]=0.0;
                velocities[i]=0.0;
            }
            for (unsigned int i = 0; i < nchunk; i++) {

                if(fread(&n,sizeof(int),1,fp)!=1) {std::cerr << "Errore nella lettura! Termino subito\n";break;}
#ifdef DEBUG_
                std::cerr << "chunk n° " << i <<" n: "<<n<<"\n";
#endif
                // extend buffer to fit chunk size

                if (n > maxbuf) {
                    if (buf){
                        delete [] buf;
                        buf=NULL;
                    }
                    buf = new double[n];
                    if (buf==NULL) {
                        std::cerr << "Allocazione di buf fallita!\n";
                        error_cont++;
                    }
                    maxbuf = n;
                } else if (!buf) {
                    std::cerr << "Qualcosa di strano è accaduto!\nn = " << n <<"; maxbuf = "<<maxbuf<<"\n";
                    error_cont++;
                }

                // read chunk and write as size_one values per line

                if(fread(buf,sizeof(double),n,fp)!=n){std::cerr << "Errore nella lettura del chunk! Termino subito\n";break;}
                unsigned int nn =n/ size_one;
                for (int j = 0; j < nn; j++) { //loop sugli atomi
                    // in buf[j*size_one] ci saranno le posizioni e le velocità da leggere, o qualsiasi altra cosa specificata in dump_custom
                    // per le varie cose che voglio calcolare
                    size_t size_one_t=2+6;
                    if (size_one!=size_one_t){
                        std::cerr << "Le dimensioni dei dati di un atomo non corrispondono ("<<size_one<<" invece di " <<size_one_t<<")\nusare \"dump 1 all custom ${dump_interval} dump.${simname}.bin id type xu yu zu vx vy vz\" in lammps\n";
                        fclose(fp);
                        fclose(fp2);
                        break;
                    }
#ifdef DEBUG_
                    std::cerr << "Atomo "<< buf[j*size_one] << "\n";
#endif
                    int iatom=round(buf[j*size_one])-1;
                    int typeatom=round(buf[j*size_one]);
                    if (!first && typeatom!=types[iatom]) {
                        std::cerr << "L'atomo con id " << iatom << " ha cambiato tipo!\n";
                        error_cont++;
                    }

                    if (iatom>=natoms_ || iatom<0){
                        std::cerr << "Attenzione: atomo fuori dal range 0-"<<natoms_-1<<" ("<< iatom <<"). Qualcosa di brutto accadra'.\n";
                        error_cont++;
                    }
                    //id type xu yu zu vx vy vz
                    types[iatom]=typeatom;
                    for (unsigned int jj=0;jj<3;jj++){
                        positions[(iatom)*3+jj]=buf[j*size_one+2+jj];
                        velocities[(iatom)*3+jj]=buf[j*size_one+5+jj];
#ifdef DEBUG_
                        std::cerr<<" x("<<jj<<") = " << positions[(iatom)*3+jj]<<" v("<<jj<<") = " << velocities[(iatom)*3+jj]<< " ";
#endif
                    }
                    atomc++;
#ifdef DEBUG_
                    std::cerr <<"\n";
#endif
                }

            }
#ifdef DEBUG_
            std::cerr <<"atomc = "<<atomc<<"\n";
#endif
            if (atomc!=natoms_){
                std::cerr << "Il numero di atomi letti non corrisponde a quello dichiatato! ("<<atomc<<" contati contro "<<natoms_<<" dichiarati)\n";
                error_cont++;
            }
            for (unsigned int i=0;i<natoms_;i++){
                if(positions[i*3]==0.0 && velocities[i*3]==0.0&&positions[i*3+1]==0.0 && velocities[i*3+1]==0.0&&positions[i*3+2]==0.0 && velocities[i*3+2]==0.0){
#ifdef DEBUG_
                    std::cerr<<"atomo "<<i<<": "<<positions[i]<<" "<<velocities[i]<<"\n";
#endif
                    error_cont++;
                }
            }
            if (error_cont>0) {
                std::cerr << "Troppi errori ("<<error_cont<<"), termino qui la lettura (ultimo step valido "<< ntimestep<<")\n";
                break;
            }
            ntimestep=nntimestep;
            if (first) {
                traiettoria.reset_types(types);
            }
            if (!heat_autocorr){
                traiettoria.add(Time,positions,velocities);
            } else {
                traiettoria.add(Time,positions,velocities,J);
            }
            for (unsigned ii=0;ii<3;ii++) cell[ii]=cell_[ii];
            natoms=natoms_;
            if (traiettoria.is_full()){
                std::cerr << "Ho riempito la traiettoria fino alla dimensione richiesta. Termino la lettura.\n";
                break;
            }
        } else if (heat_autocorr && first) {
            traiettoria.allocate_J(size);
            traiettoria.add(Time,J);
        } else if (heat_autocorr) {
            traiettoria.add(Time,J);

        }
        if (first) {
            first=false;
        }
        cron.stop();
        cron.print_interval();
    }
    std::cerr << "Tempo per la lettura della traiettoria: " << cron.time()<<"s\n";
    delete [] positions;
    delete [] velocities;
    //if (line!=NULL) delete [] line;
}

Analisi::~Analisi(){
    // per il momento nulla da deallocare.
}

//template instantiation
template class History<double,double>;
