#include "superdiff.h"
#include <stdlib.h>
#include <string.h>

t_matriz* superdiff_criar_matriz(int linhas, int colunas) {
    t_matriz* m = (t_matriz*)malloc(sizeof(t_matriz));
    m->n_linhas = linhas;
    m->n_colunas = colunas;
    m->dados = (real64*)calloc(linhas * colunas, sizeof(real64));
    return m;
}

void superdiff_liberar_matriz(t_matriz* m) {
    if (m) {
        free(m->dados);
        free(m);
    }
}

// Multiplicacao de matriz por vetor (Operacao Atomica para PDE/ODE)
void superdiff_matriz_vetor_mult(t_matriz* m, const real64* x, real64* y) {
    for (int i = 0; i < m->n_linhas; i++) {
        y[i] = 0.0;
        for (int j = 0; j < m->n_colunas; j++) {
            y[i] += m->dados[i * m->n_colunas + j] * x[j];
        }
    }
}
