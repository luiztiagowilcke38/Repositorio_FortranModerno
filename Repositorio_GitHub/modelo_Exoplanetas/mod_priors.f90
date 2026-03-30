MODULE mod_priors
  USE iso_fortran_env, ONLY: real64
  USE mod_constants, ONLY: PI_CONST
  IMPLICIT NONE
  
  ENUM, BIND(C)
    ENUMERATOR :: TIPO_UNIFORME = 1
    ENUMERATOR :: TIPO_LOG_UNIFORME = 2
    ENUMERATOR :: TIPO_GAUSSIANO = 3
    ENUMERATOR :: TIPO_GAUSSIANO_TRUNCADO = 4
    ENUMERATOR :: TIPO_JEFFREYS = 5
    ENUMERATOR :: TIPO_BETA = 6
  END ENUM

  TYPE :: estrutura_prior
    INTEGER :: tipo
    REAL(real64) :: p1, p2 ! Parametros da distribuicao (ex: min/max ou media/sigma)
    REAL(real64) :: limite_min, limite_max
  END TYPE estrutura_prior

CONTAINS

  FUNCTION calcular_log_prior(valor, prior) RESULT(log_p)
    REAL(real64), INTENT(IN) :: valor
    TYPE(estrutura_prior), INTENT(IN) :: prior
    REAL(real64) :: log_p
    
    ! Verificacao de limites basicos
    IF (valor < prior%limite_min .OR. valor > prior%limite_max) THEN
      log_p = -HUGE(1.0_real64)
      RETURN
    END IF
    
    SELECT CASE (prior%tipo)
    CASE (TIPO_UNIFORME)
      log_p = -LOG(prior%p2 - prior%p1)
      
    CASE (TIPO_LOG_UNIFORME)
      log_p = -LOG(valor * LOG(prior%p2 / prior%p1))
      
    CASE (TIPO_GAUSSIANO)
      log_p = -0.5_real64 * (((valor - prior%p1)/prior%p2)**2 + LOG(2.0_real64 * PI_CONST * prior%p2**2))
      
    CASE DEFAULT
      log_p = 0.0_real64
    END SELECT
  END FUNCTION calcular_log_prior

  SUBROUTINE test_modulo_priors()
    TYPE(estrutura_prior) :: p_teste
    REAL(real64) :: res
    p_teste%tipo = TIPO_UNIFORME
    p_teste%p1 = 0.0; p_teste%p2 = 10.0; p_teste%limite_min = 0.0; p_teste%limite_max = 10.0
    res = calcular_log_prior(5.0_real64, p_teste)
    IF (ABS(res + LOG(10.0_real64)) < 1.0E-10_real64) THEN
      PRINT *, "[OK] Modulo de Priors validado."
    ELSE
      PRINT *, "[ERRO] Falha no calculo do log-prior."
    END IF
  END SUBROUTINE test_modulo_priors

END MODULE mod_priors
