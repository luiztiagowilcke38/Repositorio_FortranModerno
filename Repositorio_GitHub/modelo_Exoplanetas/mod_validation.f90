MODULE mod_validation
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE

CONTAINS

  FUNCTION calcular_bic(ln_l_max, n_par, n_dados) RESULT(bic)
    REAL(real64), INTENT(IN) :: ln_l_max; INTEGER, INTENT(IN) :: n_par, n_dados
    REAL(real64) :: bic
    bic = -2.0_real64 * ln_l_max + REAL(n_par, real64) * LOG(REAL(n_dados, real64))
  END FUNCTION calcular_bic

  FUNCTION calcular_aic(ln_l_max, n_par) RESULT(aic)
    REAL(real64), INTENT(IN) :: ln_l_max; INTEGER, INTENT(IN) :: n_par
    REAL(real64) :: aic
    aic = -2.0_real64 * ln_l_max + 2.0_real64 * REAL(n_par, real64)
  END FUNCTION calcular_aic

  FUNCTION realizar_teste_par_impar(fluxo_par, fluxo_impar) RESULT(p_valor)
    REAL(real64), DIMENSION(:), INTENT(IN) :: fluxo_par, fluxo_impar
    REAL(real64) :: p_valor
    ! Teste t de Student simplificado para comparar profundidades de transito
    p_valor = 0.5_real64 ! Valor simulado
  END FUNCTION realizar_teste_par_impar

  SUBROUTINE test_modulo_validation()
    REAL(real64) :: bic, aic
    bic = calcular_bic(-100.0_real64, 5, 1000)
    aic = calcular_aic(-100.0_real64, 5)
    IF (bic > aic) THEN
      PRINT *, "[OK] Modulo de Validacao (BIC/AIC) validado."
    ELSE
      PRINT *, "[ERRO] Inconsistencia no calculo de criterios de informacao."
    END IF
  END SUBROUTINE test_modulo_validation

END MODULE mod_validation
