/*
 * SuperDiFF — Stochastic Differential Equations (SDE) Solvers
 */

#include "superdiff.h"
#include <time.h>

/* Gerador de numeros aleatorios Gaussiano (Box-Muller) */
real64 superdiff_aleatorio_gaussiano() {
    static int tem_proximo = 0;
    static real64 valor_proximo;
    if (tem_proximo) {
        tem_proximo = 0;
        return valor_proximo;
    }
    real64 u, v, s;
    do {
        u = 2.0 * rand() / (real64)RAND_MAX - 1.0;
        v = 2.0 * rand() / (real64)RAND_MAX - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    real64 multiplicador = sqrt(-2.0 * log(s) / s);
    valor_proximo = v * multiplicador;
    tem_proximo = 1;
    return u * multiplicador;
}

/* Euler-Maruyama: dX = a(X,t)dt + b(X,t)dW */
void superdiff_ede_euler_maruyama(real64 (*a)(real64, real64), real64 (*b)(real64, real64), 
                                  real64 x0, real64 t0, real64 tf, real64 dt, real64* t_saida, real64* x_saida) {
    int n_passos = (int)((tf - t0) / dt);
    real64 t = t0;
    real64 x = x0;
    srand(time(NULL));

    t_saida[0] = t;
    x_saida[0] = x;

    for (int i = 1; i <= n_passos; i++) {
        real64 dw = sqrt(dt) * superdiff_aleatorio_gaussiano();
        x += a(x, t) * dt + b(x, t) * dw;
        t += dt;
        t_saida[i] = t;
        x_saida[i] = x;
    }
}
