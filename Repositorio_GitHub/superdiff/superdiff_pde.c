#include "superdiff.h"
#include <omp.h>

#define TAMANHO_BLOCO 16

//==================================================================================================
// Poisson 2D via SOR com Cache-Blocking (Tiling) - Performance Computacional
//==================================================================================================
void superdiff_poisson_avx_2d(t_grade_edp* grade, real64 f_fuente) {
    int nx = grade->nx;
    int ny = grade->ny;
    real64 dx2 = grade->dx * grade->dx;
    real64* restrict dados = grade->dados_grade;
    
    // Divisao em blocos para maximizar reuso de cache L1/L2
    #pragma omp parallel
    {
        for (int iter = 0; iter < 200; iter++) {
            #pragma omp for collapse(2) schedule(static)
            for (int bi = 1; bi < nx - 1; bi += TAMANHO_BLOCO) {
                for (int bj = 1; bj < ny - 1; bj += TAMANHO_BLOCO) {
                    
                    // Processamento Interno do Bloco
                    int i_max = (bi + TAMANHO_BLOCO < nx - 1) ? bi + TAMANHO_BLOCO : nx - 1;
                    int j_max = (bj + TAMANHO_BLOCO < ny - 1) ? bj + TAMANHO_BLOCO : ny - 1;
                    
                    for (int i = bi; i < i_max; i++) {
                        #pragma omp simd
                        for (int j = bj; j < j_max; j++) {
                            real64 r = (dados[(i+1)*ny + j] + dados[(i-1)*ny + j] +
                                        dados[i*ny + (j+1)] + dados[i*ny + (j-1)] -
                                        4.0 * dados[i*ny + j]) / dx2 - f_fuente;
                            dados[i*ny + j] += 1.5 * r * dx2 / 4.0; // SOR com omega=1.5
                        }
                    }
                }
            }
            #pragma omp barrier
        }
    }
}
