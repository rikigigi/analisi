// numeri pseudo-casuali
// Riccardo Bertossa, 01/11/2013
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#define UNI (rnd_double()-0.5)*2
extern "C" {
uint32_t lcg(void);
void init_cmwc4096(void); // inizializza cmwc4096 con lcg()
uint32_t cmwc4096(void); //cmwc4096 x1 (32 bit)
uint64_t rnd64(void); //cmwc4096 x2 (64 bit)
double rnd_double(); // double fra 0.5 e 1
float rnd_single(); // single fra 0.5 e 1
void rnd_gauss(double * ris);
void istg(double a, double b, unsigned int z, unsigned int nc,unsigned int N, unsigned int *ist, void (*sorgente)(double *),double *mu, double *var);//genera l'istogramma (intervallo [a,b], numero di intervalli dell'istogramma z, numero di gruppi nc di N numeri da ricevere da *sorgente, array dell'istogramma, funzione generatrice(double[N]) )
void istg_(double a, double b, unsigned int z, unsigned int nc, unsigned int *ist, double (*sorgente)(),double *mu, double *var); // se la funzione non e' void
void istg__(double a, double b, unsigned int z, unsigned int nc, unsigned int *ist, float (*sorgente)(),float *mu, float *var); // float
// istogramma nell'intervallo [a,b], xx[z+1] contiene gli estremi di ogni intervallo,
// z è il numero di intervalli dell'istogramma, nc il numero di campioni contenuti in  *sorgente
// in mu e var vengono scritti valore medio e varianza
// *ist[z] è l'istogramma
void istogramma(float a, float b, float *xx, unsigned int z, unsigned int nc, unsigned int *ist, float *sorgente,float *mu, float *var); // float, array
float normal_gauss();
void ols (unsigned int n, double *x, double *y, double *m, double *q, double *em, double *eq, double *chi2);
void ols2 (unsigned int n, float *x, float *y, float *ey, float *m, float *q, float *em, float *eq,float *chi2);
#ifndef ANALISI
void pearson_C(unsigned int n, float *x, unsigned int *ist,float *chi2,unsigned int *nu, float (*fdistrib)(float,void *),void*);
void pearson(unsigned int n, float *x, unsigned int *ist,float *chi2,unsigned int *nu, float (*fdistrib)(float)); // calcola il chi^2 per fare il test di pearson
double gamma_n2 (unsigned int n); // calcola la gamma di n/2
float chi2_quantile(float p,unsigned int nu);
#endif
void covarianza(unsigned int n,float *x, float *y, float *mu, float *var, float * covar);
}
