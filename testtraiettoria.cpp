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



#include "testtraiettoria.h"
#include <iostream>
TestTraiettoria::TestTraiettoria(std::string filename) : Traiettoria(filename)
{

    imposta_dimensione_finestra_accesso(1000);
    imposta_inizio_accesso(0);
    std::cout << "#---0-\n";
    for (int i=0;i<get_natoms();i++) {
        std::cout << posizioni(167,i)[0] << " " << posizioni(167,i)[1] << " " <<  posizioni(167,i)[2] << " " <<
        velocita(167,i)[0] << " " <<velocita(167,i)[1] << " " <<velocita(167,i)[2] << "\n";
    }
    std::cout << "#----887-\n";
    imposta_inizio_accesso(800);
    for (int i=0;i<get_natoms();i++) {
        std::cout << posizioni(887,i)[0] << " " << posizioni(887,i)[1] << " "<<  posizioni(887,i)[2] << " " <<
        velocita(887,i)[0] << " " <<velocita(887,i)[1] << " " <<velocita(887,i)[2] << "\n";
    }
    std::cout << "#----1773-\n";
    imposta_inizio_accesso(1600);
    for (int i=0;i<get_natoms();i++) {
        std::cout << posizioni(1773,i)[0] << " " << posizioni(1773,i)[1] << " "<<  posizioni(1773,i)[2] << " " <<
        velocita(1773,i)[0] << " " <<velocita(1773,i)[1] << " " <<velocita(1773,i)[2] << "\n";
    }
/*
    std::cout << "#----1734573-\n";
    for (int i=0;i<get_natoms();i++) {
        std::cout << posizioni(1734573,i)[0] << " " << posizioni(1734573,i)[1] << " "<<  posizioni(1734573,i)[2] << " " <<
        velocita(1734573,i)[0] << " " <<velocita(1734573,i)[1] << " " <<velocita(1734573,i)[2] << "\n";
    }
*/

    std::cout << "#--------\n";

}
