#include "superdiff.h"

//==================================================================================================
// Runge-Kutta 4 Turbo (Otimizacao SIMD-Ready)
//==================================================================================================
t_solucao* superdiff_rk4_turbo(t_problema_edo* problema, real64 h_passo) {
    t_solucao* solucao = (t_solucao*)malloc(sizeof(t_solucao));
    int n_eq = problema->n_equacoes;
    int n_passos = (int)((problema->t_final - problema->t_inicial) / h_passo);
    
    solucao->n_passos = n_passos;
    solucao->t = MALLOC_ALINHADO(n_passos * sizeof(real64));
    solucao->y = (real64**)malloc(n_passos * sizeof(real64*));
    
    real64* restrict y_ativo = MALLOC_ALINHADO(n_eq * sizeof(real64));
    memcpy(y_ativo, problema->y_inicial, n_eq * sizeof(real64));
    
    real64 t_iter = problema->t_inicial;
    real64 *restrict k1 = MALLOC_ALINHADO(n_eq*sizeof(real64));
    real64 *restrict k2 = MALLOC_ALINHADO(n_eq*sizeof(real64));
    real64 *restrict k3 = MALLOC_ALINHADO(n_eq*sizeof(real64));
    real64 *restrict k4 = MALLOC_ALINHADO(n_eq*sizeof(real64));
    real64 *restrict y_tmp = MALLOC_ALINHADO(n_eq*sizeof(real64));

    real64 h2 = 0.5 * h_passo;
    real64 h6 = h_passo / 6.0;

    for (int i = 0; i < n_passos; i++) {
        solucao->t[i] = t_iter;
        solucao->y[i] = MALLOC_ALINHADO(n_eq * sizeof(real64));
        memcpy(solucao->y[i], y_ativo, n_eq * sizeof(real64));

        // Estagio 1
        problema->funcao_f(t_iter, y_ativo, k1, problema->parametros);
        #pragma omp simd
        for(int j=0; j<n_eq; j++) y_tmp[j] = y_ativo[j] + h2*k1[j];
        
        // Estagio 2
        problema->funcao_f(t_iter + h2, y_tmp, k2, problema->parametros);
        #pragma omp simd
        for(int j=0; j<n_eq; j++) y_tmp[j] = y_ativo[j] + h2*k2[j];
        
        // Estagio 3
        problema->funcao_f(t_iter + h2, y_tmp, k3, problema->parametros);
        #pragma omp simd
        for(int j=0; j<n_eq; j++) y_tmp[j] = y_ativo[j] + h_passo*k3[j];
        
        // Estagio 4
        problema->funcao_f(t_iter + h_passo, y_tmp, k4, problema->parametros);

        // Atualizacao Final (Vetorizada)
        #pragma omp simd
        for (int j = 0; j < n_eq; j++) {
            y_ativo[j] += h6 * (k1[j] + 2.0*k2[j] + 2.0*k3[j] + k4[j]);
        }
        t_iter += h_passo;
    }

    free(k1); free(k2); free(k3); free(k4); free(y_tmp); free(y_ativo);
    return solucao;
}

void superdiff_liberar_solucao(t_solucao* sol) {
    if (!sol) return;
    for (int i = 0; i < sol->n_passos; i++) free(sol->y[i]);
    free(sol->y);
    free(sol->t);
    free(sol);
}
