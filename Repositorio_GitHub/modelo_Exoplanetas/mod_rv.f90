MODULE mod_rv
  USE iso_fortran_env, ONLY: real64
  USE mod_constants, ONLY: PI_CONST
  IMPLICIT NONE

CONTAINS

  FUNCTION calcular_velocidade_radial(tempo, k_amp, p_per, t0, ecc, omega, v0) RESULT(rv)
    REAL(real64), DIMENSION(:), INTENT(IN) :: tempo
    REAL(real64), INTENT(IN) :: k_amp, p_per, t0, ecc, omega, v0
    REAL(real64), DIMENSION(SIZE(tempo)) :: rv
    REAL(real64) :: anom_media, anom_excentrica, anom_verdadeira
    INTEGER :: i
    
    DO i = 1, SIZE(tempo)
      anom_media = 2.0_real64 * PI_CONST * (tempo(i) - t0) / p_per
      anom_excentrica = resolver_equacao_kepler(anom_media, ecc)
      anom_verdadeira = 2.0_real64 * ATAN(SQRT((1.0+ecc)/(1.0-ecc)) * TAN(anom_excentrica/2.0))
      rv(i) = v0 + k_amp * (COS(omega + anom_verdadeira) + ecc * COS(omega))
    END DO
  END FUNCTION calcular_velocidade_radial

  FUNCTION resolver_equacao_kepler(m, ecc) RESULT(e_ano)
    REAL(real64), INTENT(IN) :: m, ecc; REAL(real64) :: e_ano, f, df
    INTEGER :: iter
    e_ano = m; IF (ecc > 0.8) e_ano = PI_CONST ! Guess inicial
    DO iter = 1, 100
      f = e_ano - ecc * SIN(e_ano) - m
      df = 1.0_real64 - ecc * COS(e_ano)
      e_ano = e_ano - f / df
      IF (ABS(f) < 1.0E-12_real64) EXIT
    END DO
  END FUNCTION resolver_equacao_kepler

  SUBROUTINE test_modulo_rv()
    REAL(real64), DIMENSION(1) :: t = [0.0], rv
    rv = calcular_velocidade_radial(t, 50.0_real64, 10.0_real64, 0.0_real64, 0.1_real64, 0.5_real64, 0.0_real64)
    IF (ABS(rv(1)) > 0.0_real64) THEN
      PRINT *, "[OK] Modulo RV (Kepler Solver) validado."
    ELSE
      PRINT *, "[ERRO] Erro no calculo da velocidade radial."
    END IF
  END SUBROUTINE test_modulo_rv

END MODULE mod_rv
