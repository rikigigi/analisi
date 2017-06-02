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



#include "cronometro.h"
#include <iostream>

using namespace std;

cronometro::cronometro(){
	tempo=0.0;
	t0=clock();
	calc_exp=false;
	tempo_exp=42.0;
	pi=5.0;
	last_print=0;
}

void cronometro::set_print_interval(float P){
	pi=P;
}

void cronometro::print_interval(unsigned int nthread){
    if (tempo-last_print>pi && nthread==0){
		cerr << "passati: "<<tempo<<"s stimati rimanenti: "<<tempo_exp<<"s\n";
		last_print=tempo;
	}
}

void cronometro::unset_expected(){
	calc_exp=false;
}

//ppc è la percentuale del lavoro fatta a ogni chiamata di stop
void cronometro::set_expected(float ppc){
	p=ppc;
	p_cont=0.0;
	calc_exp=true;
}

float cronometro::time(){
	return tempo;
}

float cronometro::expected(){
	return tempo_exp;
}

void cronometro::reset(){
	tempo=0.0;
}

void cronometro::start(){
	t0=clock();
}

void cronometro::stop(unsigned int nthread){
    if (nthread==0){
        clock_t t=clock();
        tempo+=((float)t-t0)/CLOCKS_PER_SEC;
        t0=t;
    }
	if (calc_exp){
		p_cont+=p;
		tempo_exp=tempo*(1.0-p_cont)/p_cont;
	}
}

