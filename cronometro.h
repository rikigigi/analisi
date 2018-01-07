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

#ifndef CRONOMETRO_H
#define CRONOMETRO_H

#include <chrono>


class cronometro {
public:
	void start(); // registra l'istante
    void stop(unsigned int nthread=0); // calcola il tempo trascorso dall'ultimo start e lo aggiunge al tempo totale
	void reset(); // resetta il tempo totale
    double time(); // ritorna il tempo totale
    void set_expected(double); // imposta il calcolo della stima del tempo rimanente, il parametro è la frazione di lavoro fatta a ogni chiamata di stop
	void unset_expected();
    double expected(); // è il tempo stimato calcolato all'ultima chiamata di stop
	cronometro();
    void set_print_interval(double); // imposta l'intervallo di tempo minimo affinché print_interval stampi qualcosa
    void print_interval(unsigned int nthread=0); // stampa qualcosa o no guardando la differenza di tempo fra il tempo di stop precedente all'ultima volta che ha stampato e l'ultima chiamata di stop
private:
    double tempo;
    double tempo_exp;
    double p, p_cont,pi,last_print;
    std::chrono::high_resolution_clock::time_point t0;
	bool calc_exp;
};

#endif
