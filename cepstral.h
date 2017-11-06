#ifndef CEPSTRAL_H
#define CEPSTRAL_H
#include <vector>
#include <utility>
#include <array>

struct Serie{
    std::array<double*,2> serie;
    unsigned int length;
    std::array<int,2> fft_idx;
};

class Cepstral
{
public:
    /**
      * nota: allocare gli array usati con fftw_alloc_real(N) o fftw_alloc_complex(N) per assicurare l'allineamnto necessario ad utilizzare le istruzioni SIMD, dove supportate.
    **/
    Cepstral(std::vector<Serie> in_);
    void calcola();
private:
    std::vector<Serie> correnti;
    double * in;
    std::vector<double *> in_src;
    unsigned int size_in;
    std::vector <fftw_complex *> fft,fft_conv;
    fftw_plan fplan;

};

#endif // CEPSTRAL_H
