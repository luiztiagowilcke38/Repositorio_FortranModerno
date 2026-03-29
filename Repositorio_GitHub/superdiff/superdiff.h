#ifndef SUPERDIFF_H
#define SUPERDIFF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <immintrin.h> // AVX support

typedef double real64;

// Alinhamento para SIMD (32 bytes para AVX, 64 para AVX-512)
#define ALINHAMENTO 64
#define MALLOC_ALINHADO(tamanho) aligned_alloc(ALINHAMENTO, (tamanho))

//--- Estruturas de Dados Otimizadas ---

typedef struct {
    int n_linhas;
    int n_colunas;
    real64* restrict dados __attribute__((aligned(ALINHAMENTO)));
} t_matriz;

typedef struct {
    void (*funcao_f)(real64, const real64* restrict, real64* restrict, void*);
    real64* restrict y_inicial;
    int n_equacoes;
    real64 t_inicial;
    real64 t_final;
    real64 tol_absoluta;
    real64 tol_relativa;
    void* parametros;
} t_problema_edo;

typedef struct {
    real64* restrict t;
    real64** restrict y;
    int n_passos;
    int n_avaliacoes;
    int status_erro;
} t_solucao;

typedef struct {
    int nx, ny, nz;
    real64 dx, dy, dz;
    real64* restrict dados_grade __attribute__((aligned(ALINHAMENTO)));
} t_grade_edp;

//--- Motores de Alta Performance (SIMD + OpenMP) ---

t_solucao* superdiff_rk4_turbo(t_problema_edo* problema, real64 h_passo);
void superdiff_edp_calor_tiled_2d(t_grade_edp* grade, real64 dt, int passos, real64 difusiv);
void superdiff_poisson_avx_2d(t_grade_edp* grade, real64 f_fuente);

//--- Algebra Linear de Baixo Nivel ---

void superdiff_gemv_turbo(const t_matriz* m, const real64* restrict x, real64* restrict y);
void superdiff_liberar_solucao(t_solucao* sol);

#endif
