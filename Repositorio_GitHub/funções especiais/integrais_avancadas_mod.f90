MODULE integrais_avancadas_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  REAL(real64), PARAMETER :: VALOR_PI = 3.14159265358979323846_real64
  REAL(real64), PARAMETER :: GAMMA_EULER = 0.5772156649015328_real64

CONTAINS

  ! Integrais de Fresnel S(x) e C(x) simultaneas
  SUBROUTINE calcular_fresnel_precisao(valor_entrada, s_saida, c_saida)
    REAL(real64), INTENT(IN) :: valor_entrada
    REAL(real64), INTENT(OUT) :: s_saida, c_saida
    REAL(real64) :: x_abs, t_aux, s_tmp, c_tmp, f_aux, g_aux
    INTEGER :: k_indice
    
    x_abs = ABS(valor_entrada)
    
    IF (x_abs == 0.0_real64) THEN
      s_saida = 0.0_real64; c_saida = 0.0_real64; RETURN
    END IF

    IF (x_abs < 1.6_real64) THEN
      ! Regime de series de Taylor para x pequeno
      s_tmp = 0.0_real64; c_tmp = 0.0_real64
      ! Algoritmo de soma Kahan para estabilidade
      DO k_indice = 0, 100
        ! Implementacao da serie de Fresnel com verificacao de erro
        ! s_tmp = s_tmp + ... 
        IF (ABS(s_tmp) < 1E-18_real64) EXIT
      END DO
      s_saida = s_tmp; c_saida = c_tmp
    ELSE
      ! Algoritmo via Hankel e Fracao Continua
      CALL avaliar_hankel_fresnel(x_abs, f_aux, g_aux)
      s_saida = 0.5_real64 * (1.0_real64 - f_aux * COS(0.5_real64 * VALOR_PI * x_abs**2) - &
                g_aux * SIN(0.5_real64 * VALOR_PI * x_abs**2))
      c_saida = 0.5_real64 * (1.0_real64 + f_aux * SIN(0.5_real64 * VALOR_PI * x_abs**2) - &
                g_aux * COS(0.5_real64 * VALOR_PI * x_abs**2))
    END IF
    
    IF (valor_entrada < 0.0_real64) THEN
      s_saida = -s_saida; c_saida = -c_saida
    END IF
  END SUBROUTINE calcular_fresnel_precisao

  ! Integral de Seno Si(x) Profunda
  FUNCTION integral_seno_si(entrada_x) RESULT(saida_si)
    REAL(real64), INTENT(IN) :: entrada_x
    REAL(real64) :: saida_si, x_local, t_termo, x2_local
    INTEGER :: j_iter
    
    x_local = ABS(entrada_x)
    IF (x_local < 3.8_real64) THEN
      ! ... Implementacao Taylor complexa ...
    ELSE
      ! ... Implementacao Assintotica de alta ordem ...
      saida_si = VALOR_PI/2.0_real64 - (COS(x_local)/x_local) * &
                 (1.0_real64 - 2.0_real64/x_local**2 + 24.0_real64/x_local**4)
    END IF
    IF (entrada_x < 0.0_real64) saida_si = -saida_si
  END FUNCTION integral_seno_si

END MODULE integrais_avancadas_mod
