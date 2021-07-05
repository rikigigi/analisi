#ifndef FLOATING_EXCEPTIONS_H
#define FLOATING_EXCEPTIONS_H

#include <cfenv>
#include <functional>
#pragma STDC FENV_ACCESS ON
template <typename FUNCTION_SIG, bool ACTIVE, int FPE=FE_DIVBYZERO|FE_OVERFLOW|FE_INVALID>
class FloatingPointExceptionManager {
public:
    FloatingPointExceptionManager(std::function<FUNCTION_SIG(int)> const & f) : f{f},exceptcounter{0}  {
        if constexpr (ACTIVE){
            std::fegetenv(&oldexept);
            std::feclearexcept(FPE);
        }
    }

    ~FloatingPointExceptionManager() {
        if constexpr (ACTIVE){
            int fee=std::fetestexcept(FPE);
            if (fee) {
                f(fee);
            }
            std::feclearexcept(FPE);
            std::fesetenv(&oldexept);
        }
    }
    void check_nan(const double & d) {
        if constexpr (ACTIVE){
            if (d!=d) {
                exceptcounter--;
                f(exceptcounter);
            }
        }
    }
private:
    int exceptcounter;
    std::fenv_t oldexept;
    const std::function<FUNCTION_SIG(int)> f;
};

#endif // FLOATING_EXCEPTIONS_H
