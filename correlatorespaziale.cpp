#include "correlatorespaziale.h"
#include "traiettoria.h"
#include "operazionisulista.h"

#include "config.h"
#ifdef HAVEfftw3
#include <fftw3.h>
#else
#include <fftw.h>
#endif
CorrelatoreSpaziale::CorrelatoreSpaziale(Traiettoria *t,
                                          unsigned int n,
                                          TFLOAT sigma2,
                                          unsigned int nthreads=0,
                                          unsigned int skip=1,
                                          bool debug=false) :
    e_of_r(0), n(n), sigma2(sigma2),nthreads(nthreads),debug(debug),skip(skip)
{


}

template  <class TFLOAT>
void CorrelatoreSpaziale<TFLOAT>::reset(const unsigned int numeroTimestepsPerBlocco) {
    if (e_of_r==0)
        e_of_r = (TFLOAT *) fftw_malloc(sizeof(TFLOAT)*n*n*n);
}

CorrelatoreSpaziale::~CorrelatoreSpaziale(){
    fftw_free(e_of_r);
}
