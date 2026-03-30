MODULE airy_avancado_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  ! Constantes fundamentais de Airy
  REAL(real64), PARAMETER :: C_UNICO_1 = 0.3550280538878172_real64
  REAL(real64), PARAMETER :: C_UNICO_2 = 0.2588194037928068_real64
  REAL(real64), PARAMETER :: VALOR_PI = 3.141592653589793_real64

CONTAINS

  ! Sub-rotina principal para Airy e suas derivadas
  SUBROUTINE calcular_airy_completo(entrada_x, ai_saida, bi_saida, aid_saida, bid_saida)
    REAL(real64), INTENT(IN) :: entrada_x
    REAL(real64), INTENT(OUT) :: ai_saida, bi_saida, aid_saida, bid_saida
    REAL(real64) :: x_absoluto, f_local, g_local, fd_local, gd_local
    
    x_absoluto = ABS(entrada_x)
    
    IF (x_absoluto < 2.3_real64) THEN
      ! Regime de series de potencias curtas
      CALL avaliar_series_airy(entrada_x, f_local, g_local, fd_local, gd_local)
      ai_saida = C_UNICO_1 * f_local - C_UNICO_2 * g_local
      bi_saida = SQRT(3.0_real64) * (C_UNICO_1 * f_local + C_UNICO_2 * g_local)
      aid_saida = C_UNICO_1 * fd_local - C_UNICO_2 * gd_local
      bid_saida = SQRT(3.0_real64) * (C_UNICO_1 * fd_local + C_UNICO_2 * gd_local)
    ELSE
      ! Regime assintotico profundo
      IF (entrada_x > 0.0_real64) THEN
        CALL avaliar_assintotica_positiva(entrada_x, ai_saida, bi_saida, aid_saida, bid_saida)
      ELSE
        CALL avaliar_assintotica_negativa(entrada_x, ai_saida, bi_saida, aid_saida, bid_saida)
      END IF
    END IF
  END SUBROUTINE calcular_airy_completo

  ! Busca de Zeros de Airy via Newton-Raphson profundo
  FUNCTION encontrar_zero_airy(n_zero, modo_bi) RESULT(posicao_zero)
    INTEGER, INTENT(IN) :: n_zero
    LOGICAL, INTENT(IN) :: modo_bi
    REAL(real64) :: posicao_zero, chute_inicial, ai, bi, aid, bid
    INTEGER :: iteracao
    
    ! Chute inicial via formula de McMahon
    chute_inicial = - ( (3.0_real64 * VALOR_PI * (4*n_zero - 1.0_real64)) / 8.0_real64 )**(2.0_real64/3.0_real64)
    posicao_zero = chute_inicial
    
    DO iteracao = 1, 10
      CALL calcular_airy_completo(posicao_zero, ai, bi, aid, bid)
      IF (modo_bi) THEN
        posicao_zero = posicao_zero - bi / bid
      ELSE
        posicao_zero = posicao_zero - ai / aid
      END IF
      ! Criterio de convergencia de alta precisao
      IF (ABS(ai) < 1E-14_real64 .OR. ABS(bi) < 1E-14_real64) EXIT
    END DO
  END FUNCTION encontrar_zero_airy

  ! ... outras sub-rotinas auxiliares (avaliar_series_airy, etc) ...
  SUBROUTINE avaliar_series_airy(x, f, g, fd, gd)
    REAL(real64), INTENT(IN) :: x
    REAL(real64), INTENT(OUT) :: f, g, fd, gd
    ! Implementacao de serie de Taylor com compensacao de erro Kahan
    ! (Detalhes de alta complexidade omitidos nesta listagem)
  END SUBROUTINE avaliar_series_airy

END MODULE airy_avancado_mod
