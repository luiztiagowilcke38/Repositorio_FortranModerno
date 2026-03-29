/*
 * SuperDiFF — Fractional Differential Equations
 */

#include "superdiff.h"

/* Gamma function approximation (Lanczos) */
real64 superdiff_gamma(real64 x) {
    static const real64 p[] = {
        676.5203681218851, -1259.1392167224028, 771.32342877765313,
        -176.61502916214059, 12.507343278686905, -0.13857109526572012,
        9.9843695780195716e-6, 1.5056327351493116e-7
    };
    if (x < 0.5) return M_PI / (sin(M_PI * x) * superdiff_gamma(1.0 - x));
    x -= 1.0;
    real64 a = 0.99999999999980993;
    real64 t = x + 7.5;
    for (int i = 0; i < 8; i++) a += p[i] / (x + i + 1.0);
    return sqrt(2.0 * M_PI) * pow(t, x + 0.5) * exp(-t) * a;
}

/* Grunwald-Letnikov para derivada fracionaria de ordem alfa */
void superdiff_fracionaria_grunwald(real64 alfa, real64 (*f)(real64), real64 t, real64 h, real64* resultado) {
    int n = (int)(t / h);
    real64 soma = 0.0;
    real64 w = 1.0;
    for (int k = 0; k <= n; k++) {
        soma += w * f(t - k * h);
        w *= (k - alfa) / (k + 1.0);
    }
    *resultado = soma / pow(h, alfa);
}
