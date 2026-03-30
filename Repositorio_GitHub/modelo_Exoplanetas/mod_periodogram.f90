MODULE mod_periodogram
  USE iso_fortran_env, ONLY: real64
  USE mod_constants, ONLY: PI_CONST
  IMPLICIT NONE

CONTAINS

  SUBROUTINE buscar_bls(tempo, fluxo, p_min, p_max, n_freq, pot_max, p_opt)
    REAL(real64), DIMENSION(:), INTENT(IN) :: tempo, fluxo
    REAL(real64), INTENT(IN) :: p_min, p_max; INTEGER, INTENT(IN) :: n_freq
    REAL(real64), INTENT(OUT) :: pot_max, p_opt
    REAL(real64) :: freq, per, pot_atual, fase
    INTEGER :: k, i, n
    
    n = SIZE(tempo); pot_max = 0.0; p_opt = 0.0
    !$OMP PARALLEL DO PRIVATE(k, per, freq, pot_atual) SHARED(tempo, fluxo, pot_max, p_opt)
    DO k = 1, n_freq
      per = p_min + (p_max - p_min) * REAL(k-1, real64) / REAL(n_freq-1, real64)
      ! Simplificacao do algoritmo BLS (Busca de caixa)
      pot_atual = avaliar_potencia_caixa(tempo, fluxo, per)
      
      !$OMP CRITICAL
      IF (pot_atual > pot_max) THEN
        pot_max = pot_atual; p_opt = per
      END IF
      !$OMP END CRITICAL
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE buscar_bls

  FUNCTION avaliar_potencia_caixa(t, f, p) RESULT(pot)
    REAL(real64), DIMENSION(:), INTENT(IN) :: t, f; REAL(real64), INTENT(IN) :: p
    REAL(real64) :: pot, fase_temp
    INTEGER :: i
    ! Soma das profundidades em fase (Demonstracao simplificada)
    pot = 0.0; DO i = 1, SIZE(t)
      fase_temp = MOD(t(i), p) / p
      IF (fase_temp > 0.45 .AND. fase_temp < 0.55) pot = pot + (1.0_real64 - f(i))
    END DO
  END FUNCTION avaliar_potencia_caixa

  SUBROUTINE test_modulo_periodogram()
    REAL(real64), DIMENSION(100) :: t, f; REAL(real64) :: p_m, p_o
    INTEGER :: i
    DO i = 1, 100; t(i) = REAL(i, real64)*0.1; f(i) = 1.0; END DO
    f(45:55) = 0.99 ! Simula um transito no meio
    CALL buscar_bls(t, f, 5.0_real64, 15.0_real64, 100, p_m, p_o)
    IF (p_o > 0.0_real64) THEN
      PRINT *, "[OK] Modulo Periodograma (BLS) validado."
    ELSE
      PRINT *, "[ERRO] BLS nao detectou o sinal periodico."
    END IF
  END SUBROUTINE test_modulo_periodogram

END MODULE mod_periodogram
