#ifndef CALCOLIBLOCCHI_H
#define CALCOLIBLOCCHI_H




template <class T> class MediaVar{
public:
    MediaVar(T*Tmedio, T*Tvar,T*delta,T*tmp) : Tmedio(Tmedio), Tvar(Tvar),delta(delta),tmp(tmp) {}

    void calcola_begin(unsigned int s){
        Tmedio->reset(s);
        Tmedio->azzera();
        Tvar->reset(s);
        Tvar->azzera();
        delta->reset(s);
        tmp->reset(s);
        iblock=0;
    }

    void calcola(T* calcolo) {
#ifdef DEBUG2
        std::cerr << "*delta = *calcolo - *Tmedio;\n";
#endif
        //algoritmo media e varianza online
        *delta = *calcolo;
        *delta -= *Tmedio;
        //*delta = *calcolo - *Tmedio
#ifdef DEBUG2
        std::cerr << "*delta = *calcolo - *Tmedio;\n";
#endif
        *tmp = *delta;
        *tmp /= (iblock+1);
        *Tmedio += *tmp;
        //*Tmedio += (*delta)/(iblock+1); // calcolo / double , calcolo+=calcolo
#ifdef DEBUG2
        std::cerr << "*Tvar += (*delta)*(*calcolo-*Tmedio);\n";
#endif
        *tmp = *calcolo;
        *tmp -= *Tmedio;
        *tmp *= *delta;
        *Tvar += *tmp;
        iblock++;
    }
    void calcola_end(unsigned int n_b){
        *Tvar/=((n_b-1)*n_b);
    }
private:
    T *Tmedio,*Tvar,*calcolo,*delta,*tmp;
    unsigned int iblock;
};


template <class T> class MediaVarCovar{
public:
    MediaVarCovar() {

    }

    void calcola_begin(unsigned int s){
        Tmedio->reset(s);
        Tmedio->azzera();
        Tvar->reset(s);
        Tvar->azzera();
        delta->reset(s);
        tmp->reset(s);
        iblock=0;
    }

    void calcola(T* calcolo) {
#ifdef DEBUG2
        std::cerr << "*delta = *calcolo - *Tmedio;\n";
#endif
        //algoritmo media e varianza online
        *delta = *calcolo;
        *delta -= *Tmedio;
        //*delta = *calcolo - *Tmedio
#ifdef DEBUG2
        std::cerr << "*delta = *calcolo - *Tmedio;\n";
#endif
        *tmp = *delta;
        *tmp /= (iblock+1);
        *Tmedio += *tmp;
        //*Tmedio += (*delta)/(iblock+1); // calcolo / double , calcolo+=calcolo
#ifdef DEBUG2
        std::cerr << "*Tvar += (*delta)*(*calcolo-*Tmedio);\n";
#endif
        *tmp = *calcolo;
        *tmp -= *Tmedio;
        *tmp *= *delta;
        *Tvar += *tmp;
        iblock++;
    }
    void calcola_end(unsigned int n_b){
        *Tvar/=((n_b-1)*n_b);
    }
private:
    T *Tmedio,*Tvar,*calcolo,*delta,*tmp;
    unsigned int iblock;
};

#endif // CALCOLIBLOCCHI_H
