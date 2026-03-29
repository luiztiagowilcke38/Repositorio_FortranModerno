!==================================================================================================
! MODULO: MOD_STATS_TESTS (TESTES ESTATISTICOS DE KOLMOGOROV-SMIRNOV)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Implementacao do teste de bondade de ajuste de Kolmogorov-Smirnov (K-S).
! Compara uma distribuicao amostral com uma teorica ou duas amostras entre si.
!
! REFERENCIA: Numerical Recipes in Fortran 90 / 14.3. Kolmogorov-Smirnov Test.
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE mod_stats_tests
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ks_test_1d, prob_ks

  INTEGER, PARAMETER :: DP = real64

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! TESTE K-S PARA UMA AMOSTRA CONTRA DISTRIBUICAO TEORICA
  !------------------------------------------------------------------------------------------------
  SUBROUTINE ks_test_1d(data, n, func, d_stat, prob)
    REAL(DP), INTENT(IN)  :: data(n)
    INTEGER,  INTENT(IN)  :: n
    REAL(DP), INTENT(OUT) :: d_stat, prob
    
    INTERFACE
       FUNCTION func(x)
          IMPORT :: DP
          REAL(DP), INTENT(IN) :: x
          REAL(DP) :: func
       END FUNCTION func
    END INTERFACE

    REAL(DP) :: data_sorted(n), en, ff, fn, dt
    INTEGER  :: j

    data_sorted = data
    CALL sort(data_sorted, n)
    
    en = REAL(n, DP)
    d_stat = 0.0_DP
    DO j = 1, n
       fn = j / en
       ff = func(data_sorted(j))
       dt = MAX(ABS((j-1)/en - ff), ABS(fn - ff))
       IF (dt > d_stat) d_stat = dt
    END DO
    
    prob = prob_ks(SQRT(en) * d_stat)
  END SUBROUTINE ks_test_1d

  !------------------------------------------------------------------------------------------------
  ! FUNCAO DE PROBABILIDADE K-S (SERIE DE KOLMOGOROV)
  !------------------------------------------------------------------------------------------------
  FUNCTION prob_ks(alam) RESULT(prob)
    REAL(DP), INTENT(IN) :: alam
    REAL(DP) :: prob
    REAL(DP), PARAMETER :: EPS1 = 0.001_DP, EPS2 = 1.0E-8_DP
    REAL(DP) :: a2, fac, term, term_old
    INTEGER  :: j

    a2 = -2.0_DP * alam**2
    fac = 2.0_DP
    prob = 0.0_DP
    term_old = 0.0_DP
    
    DO j = 1, 100
       term = fac * EXP(a2 * REAL(j**2, DP))
       prob = prob + term
       IF (ABS(term) <= EPS1 * term_old .OR. ABS(term) <= EPS2) RETURN
       fac = -fac
       term_old = ABS(term)
    END DO
    prob = 1.0_DP ! Falha na convergencia
  END FUNCTION prob_ks

  !------------------------------------------------------------------------------------------------
  ! UTILITARIO DE ORDENACAO (SHELL SORT)
  !------------------------------------------------------------------------------------------------
  SUBROUTINE sort(arr, n)
    REAL(DP), INTENT(INOUT) :: arr(n)
    INTEGER,  INTENT(IN)    :: n
    INTEGER :: i, j, h
    REAL(DP) :: v
    h = 1
    DO WHILE (h <= n/3) ; h = 3*h + 1 ; END DO
    DO WHILE (h >= 1)
       DO i = h + 1, n
          v = arr(i) ; j = i
          DO WHILE (j > h .AND. arr(j-h) > v)
             arr(j) = arr(j-h) ; j = j - h
          END DO
          arr(j) = v
       END DO
       h = h / 3
    END DO
  END SUBROUTINE sort

END MODULE mod_stats_tests
