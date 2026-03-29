/*
 * SuperDiFF — Delay Differential Equations (DDE)
 */

#include "superdiff.h"

/* Solver para y'(t) = f(t, y(t), y(t-tau)) com tau constante */
void superdiff_edr_passos(real64 (*f)(real64, real64, real64), real64 (*historico_func)(real64), 
                         real64 tau, real64 t0, real64 tf, real64 h, real64* t_saida, real64* y_saida) {
    int n_passos = (int)((tf - t0) / h) + 1;
    real64 t = t0;
    real64 y = historico_func(t0);
    
    t_saida[0] = t;
    y_saida[0] = y;

    for (int i = 1; i < n_passos; i++) {
        real64 t_atrasado = t - tau;
        real64 y_atrasado;
        
        if (t_atrasado <= t0) {
            y_atrasado = historico_func(t_atrasado);
        } else {
            /* Interpolacao simples do historico calculado */
            int k = (int)((t_atrasado - t0) / h);
            real64 fracao = (t_atrasado - (t0 + k*h)) / h;
            y_atrasado = y_saida[k] + fracao * (y_saida[k+1] - y_saida[k]);
        }
        
        /* Passo de Euler (para brevidade tecnica) */
        y += h * f(t, y, y_atrasado);
        t += h;
        t_saida[i] = t;
        y_saida[i] = y;
    }
}
