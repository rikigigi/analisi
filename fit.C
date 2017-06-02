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



#include "fit.h"
#include "cronometro.h"
#include "interfaccia.h"
#include <cmath>

#include <iostream>
#include <cstring>
#include <cstdio>


template <unsigned int N_lin> bool fit<N_lin>::output;
template <unsigned int N_lin> unsigned int fit<N_lin>::fitcont=0;
template <unsigned int N_lin> unsigned int fit<N_lin>::fitskip=4;
template <unsigned int N_lin> char * fit<N_lin>::fout;


template <unsigned int N_lin> void fit<N_lin>::set_fitskip(unsigned int Nf){
    fitskip=Nf;
}

/** Implementazione con Eigen del metodo del linear least squares.
 * Siano \f$ (x_i,y_i)\,\,i\in\{1,\dots,m\} \f$ i dati misurati,
 *  e \f$ X_{ij}(x_i)\,\,i\in\{1,\dots,m\}\,\,J\in\{1,\dots,n\}\f$ una funzione generica dei dati.
 * La funzione risolve il problema di rendere minimo l'errore che si ottiene assumendo una relazione del tipo
 * \f$ \sum_{j=1}^n X_{ij}\beta_j=y_i \,\,i\in\{1,\dots,m\}\f$, restituendo i \f$\beta_j$\f$
 * */
template <unsigned int N_lin > void fit<N_lin>::ols_multidim(const Matrix<double,Dynamic,N__lin> &X, ///< matrice degli \f$X_{ij}\f$
																	const Matrix<double,Dynamic,1> &Y, ///< array degli \f$y_i\f$
                                                                    Matrix<double,N__lin,1> &B ///< restituisce i parametri \f$\beta_j\f$
																	){
    Matrix<double,N__lin,Dynamic> Xt=X.transpose();
    Matrix<double,N__lin,N__lin> M=Xt*X;
    Matrix<double,N__lin,1> A=Xt*Y;
	B=M.fullPivHouseholderQr().solve(A);
}

/**calcola la matrice
 * */
template <unsigned int N_lin> double fit<N_lin>::genX(const unsigned int &i,
												 const unsigned int &j,
												 const double &tau,
												 const double *x
												 ){
	if (j==0) return 1.0;
	return exp(-x[i]*(2*j-1)*(2*j-1)/tau);
}

/**La funzione che si sta fittando
 * */
template <unsigned int N_lin> double fit<N_lin>::fit_func(const double &t,
                                                     const Matrix<double,N__lin,1> &B,
													 const double &tau
													 ){
	double ret=0.0;
    for (int i=N_lin-1;i>=0;i--) {
		ret+=genX(0,i,tau,&t)*B(i);
	}
	return ret;
}

template <unsigned int N_lin> double fit<N_lin>::get_fit_func(const unsigned int & N, const double *x, double * y, const double * b, const double & tau){
    Matrix<double,N__lin,1> B(N_lin,1);
    for (unsigned int i=0;i<N_lin;i++) B(i)=b[i];
    for (unsigned int i=0;i<N;i++){
        y[i]=fit_func(x[i],B,tau);
    }
}

/**funzione ausiliaria per il calcolo della derivata rispetto a 1/tau (che è più semplice)
 * */
template <unsigned int N_lin> double fit<N_lin>::aus_func(const double &t,
                                                     const Matrix<double,N__lin,1> &B,
													 const double &tau
													 ){
	double ret=0.0; // il primo è il termine costante rispetto a lambda e la sua derivata fa zero
    for (int i=N_lin-1;i>=1;i--) {
		ret+=genX(0,i,tau,&t)*B(i)*(-(2*i-1)*(2*i-1)*t);
	}
	return ret;
}

template <unsigned int N_lin> double fit<N_lin>::chi2(const unsigned int &N,
											     const double *x,
											     const double *y, 
                                                 const Matrix<double,N__lin,1> &B,
											     const double &tau
											     ){
	double ris=0.0;
    for (unsigned int i=0;i<N;i=i+fitskip){
		double f=fit_func(x[i],B,tau);
		ris+=(f-y[i])*(f-y[i]);
	}
	return ris;
}

/**Derivata del chi2 rispetto a 1/tau
 * */
template <unsigned int N_lin> double fit<N_lin>::chi2_derivat (const unsigned int &N,
														  const double *x,
														  const double *y, 
                                                          const Matrix<double,N__lin,1> &B,
														  const double &tau
														  ){
	double ris=0.0;
    for (unsigned int i=0;i<N;i=i+fitskip){
		double f=fit_func(x[i],B,tau);
		ris+=(f-y[i])*aus_func(x[i],B,tau);
	}
	return ris*2.0;
}


/**Fitta la funzione per ricavare tau e i vari coefficienti. Assume che lo stato iniziale dello sviluppo in serie di fourier ha solo termini dispari.
 * Può essere facilmente modificata per includere anche i
 *termini che vengono dai termini pari dello sviluppo in serie di fourier dello stato iniziale.
 * */
template <unsigned int N_lin> void fit<N_lin>::fit_lin_exp(const unsigned int &N, ///< numero di dati
                                                           const double *x, ///< array con le ascisse (per esempio i tempi)
                                                           const double *y, ///< array con le misure (per esempio l'integrale della densità su metà scatola)
                                                           double *b, ///< array di lunghezza N_lin dove verranno scritti i coefficienti che stanno davanti alle esponenziali
                                                           double &t, ///< usato come valore iniziale di tau, poi qui verrà scritto il valore di tau calcolato
                                                           bool teor_lin, ///< se true non esegue il fit per la parte lineare ma usa i valori calcolati teoricamente
                                                           unsigned int natoms ///< numero di atomi nella scatola
                                                           ){
    Matrix<double,N__lin,1> B(N_lin);
    Matrix<double,Dynamic,1> Y(N,1);
    Matrix<double,Dynamic,N__lin> X(N,N_lin);
    //valori iniziali
    double tau=t,lambda=1.0/tau,tmin=0.0;

    B(0)=(double)natoms/4.0;
    for (unsigned int i = 1;i<N_lin;i++){
        B(i)=(double)natoms*2/((2*i-1)*(2*i-1)*PI*PI);
    }
    //copia le y
    for (unsigned int i=0;i<N;i++){
        Y(i)=y[i];
    }

    std::ofstream F;
    char * fout_n=NULL;
    if (output){
        fout_n= new char[strlen(fout)+5];
        strcpy(fout_n,fout);
        char tt[5];
        sprintf(tt,"%d",++fitcont);
        strcat(fout_n,tt);
        F.open(fout_n,std::ofstream::app);
    }

    unsigned int step=0;
    const unsigned int n_seg=42;
    double *coeff=new double[n_seg];
    for (unsigned int i=0;i<n_seg;i++) coeff[i]=(double)i/((double)(n_seg-1));
    while(1){
      //  double coeff=1.0;
        if (!teor_lin){
            //calcola gli elementi di matrice in base al tau (nuovo tau)
            for (unsigned int i=0;i<N;i++){
                for (unsigned int j=0;j<N_lin;j++){
                    X(i,j)=genX(i,j,tau,x);
                }
            }
            //risolve il sistema dei minimi quadrati per la parte lineare
            ols_multidim(X,Y,B);
        }
        /*
        //steepest discend + termine di momento per 1/tau = lambda
        if (step==0) {
            derivat0= chi2_derivat(N,x,y,B,tau);
            d=-0.001*derivat0/fabs(derivat0);
            derivat0=fabs(derivat0);
        }else{
            // derivat0 viene usato per tenere conto dell'ordine di grandezza della derivata
            d=-0.00001 * chi2_derivat(N,x,y,B,tau)/derivat0;
        }
        // ho dimenticato un meno da qualche parte?
        lambda = lambda-  d;
        tau=1.0/lambda;
#ifdef DEBUG_
        std::cerr << "========\n" <<"passo "<<step+1<<"\n"<<"parametri trovati: \n"<<B<<"\n"<<"delta_lamda = "<<d<<"\nvalore di lamda e tau:" << lambda <<" " <<tau<<"\n";
        //DEBUG_
#endif
        if(++step % 30 == 0){
            // quando l'incremento di lambda è troppo piccolo l'algoritmo si ferma
            if (fabs(d)<1e-14) break;
            cron.stop();
            if (cron.time()>120.0){
                std::cerr << "Il tempo di fit ha superato il ragionevole tempo di 120s ("<<cron.time()<<"). Termino qui il fit. (d="<<d<<")\n";
                break;
            }
           // coeff=0.01+exp(-((double)step)/500.0);
        }*/
        double a=0.0,b=1000.0,b_=b,tcur=0.0;
        double min=1.0/0.0;
        for (unsigned int j=0;j<15;j++){
            tcur=a;
            tmin=a;
            //min=chi2(N,x,y,B,a);
            if (output) F << tcur <<" " <<min <<"\n";
            for (unsigned int i=0;i<n_seg;i++){
                tcur=a+(b-a)*coeff[i];
                double C=chi2(N,x,y,B,tcur);
                if (output) F << tcur <<" " <<C <<" \n";
                if (C<min && C>=0 && std::isfinite(C)) {
                    min=C;
                    tmin=tcur;
                }
            }
            b_=b;
            b=tmin+(b_-a)*coeff[1];
            a=tmin-(b_-a)*coeff[1];
        }

        if (teor_lin) {break;} else {
            if (fabs(tau-tmin)<1e-6) break;
            tau=tmin;
        }
    }
    //scrive i risultati
    for (unsigned int i=0;i<N_lin;i++) b[i]=B(i);
    t=tmin;
    if (output) F.close();
    delete [] coeff;
}

template <unsigned int N_lin> void fit<N_lin>::set_output(char *o){
    fout=o;
    output=true;
}

//template instantation
template class fit<FIT_LIN_NPAR>;
