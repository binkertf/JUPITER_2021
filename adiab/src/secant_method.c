#include "jupiter.h"


real secant_method(real e, real rho, real rel_err){
    /*
    energy e in code units
    density rho in code units
    pressure p in code units
    */
    real a, b, c;
    int max_it = 20; int i = 1;
    
    real zero_func(real p){
        return (p - e*(Gamma1(p * MU / rho * TEMP0, rho)-1.));
    }
    
    a = 1.0001*e*(GAMMA -1.);
    b = 0.9999*e*(GAMMA -1.);
    c = b - zero_func(b) * (b - a) / (zero_func(b) - zero_func(a));
    
    while(zero_func(c)/c > rel_err && i<max_it){
        a = b;
        b = c;
        c = b - zero_func(b) * (b - a) / (zero_func(b) - zero_func(a));
        ++i;
    }
    return c;
}