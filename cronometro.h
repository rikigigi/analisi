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



#include <time.h>


class cronometro {
public:
	void start(); // registra l'istante
    void stop(unsigned int nthread=0); // calcola il tempo trascorso dall'ultimo start e lo aggiunge al tempo totale
	void reset(); // resetta il tempo totale
	float time(); // ritorna il tempo totale
	void set_expected(float); // imposta il calcolo della stima del tempo rimanente, il parametro è la frazione di lavoro fatta a ogni chiamata di stop
	void unset_expected();
	float expected(); // è il tempo stimato calcolato all'ultima chiamata di stop
	cronometro();
	void set_print_interval(float); // imposta l'intervallo di tempo minimo affinché print_interval stampi qualcosa
    void print_interval(unsigned int nthread=0); // stampa qualcosa o no guardando la differenza di tempo fra il tempo di stop precedente all'ultima volta che ha stampato e l'ultima chiamata di stop
private:
	float tempo;
	float tempo_exp;
	float p, p_cont,pi,last_print;
	clock_t t0;
	bool calc_exp;
};

