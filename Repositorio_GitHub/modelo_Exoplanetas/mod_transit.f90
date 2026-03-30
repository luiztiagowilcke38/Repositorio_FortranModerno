MODULE mod_transit
  USE iso_fortran_env, ONLY: real64
  USE mod_constants, ONLY: PI_CONST
  IMPLICIT NONE
  
  TYPE :: modelo_transito_par
    REAL(real64) :: r_p_s      ! Razao de raios Rp/Rs
    REAL(real64) :: a_s        ! Semieixo normalizado a/Rs
    REAL(real64) :: incl       ! Inclinacao (rad)
    REAL(real64) :: periodo    ! Periodo (dias)
    REAL(real64) :: t0         ! Tempo de conjuncao
    REAL(real64) :: u1, u2     ! Coeficientes de escurecimento de limbo
    REAL(real64) :: ecc, omega ! Excentricidade e argumento do periastro
  END TYPE modelo_transito_par

CONTAINS

  FUNCTION calcular_modelo_mandel_agol(tempo, par) RESULT(fluxo)
    REAL(real64), DIMENSION(:), INTENT(IN) :: tempo
    TYPE(modelo_transito_par), INTENT(IN) :: par
    REAL(real64), DIMENSION(SIZE(tempo)) :: fluxo
    REAL(real64) :: fase, d_proj, z, k0, k1, lambda_e
    INTEGER :: i
    
    DO i = 1, SIZE(tempo)
      ! Fase orbital simplificada (circular por enquanto para o modelo analitico)
      fase = 2.0_real64 * PI_CONST * (tempo(i) - par%t0) / par%periodo
      d_proj = par%a_s * SQRT(SIN(fase)**2 + (COS(par%incl)*COS(fase))**2)
      
      IF (d_proj > 1.0_real64 + par%r_p_s) THEN
        fluxo(i) = 1.0_real64
      ELSE
        z = d_proj
        CALL avaliar_ocultacao_analitica(z, par%r_p_s, lambda_e)
        ! Aplicacao da lei quadratica de escurecimento de limbo
        fluxo(i) = 1.0_real64 - lambda_e * (1.0_real64 - par%u1/3.0_real64 - par%u2/6.0_real64)
      END IF
    END DO
  END FUNCTION calcular_modelo_mandel_agol

  SUBROUTINE avaliar_ocultacao_analitica(d, p, lambda_e)
    REAL(real64), INTENT(IN) :: d, p; REAL(real64), INTENT(OUT) :: lambda_e
    REAL(real64) :: k0, k1
    IF (d >= 1.0 + p) THEN; lambda_e = 0.0; RETURN; END IF
    IF (d <= ABS(1.0 - p)) THEN; lambda_e = p**2; RETURN; END IF
    
    k0 = ACOS((p**2 + d**2 - 1.0_real64)/(2.0_real64 * p * d))
    k1 = ACOS((1.0_real64 + d**2 - p**2)/(2.0_real64 * d))
    lambda_e = (p**2 * k0 + k1 - 0.5_real64 * SQRT(MAX(0.0_real64, &
               (1.0_real64+p-d)*(p+d-1.0_real64)*(d+1.0_real64-p)*(1.0_real64+p+d)))) / PI_CONST
  END SUBROUTINE avaliar_ocultacao_analitica

  SUBROUTINE test_modulo_transit()
    TYPE(modelo_transito_par) :: p
    REAL(real64), DIMENSION(1) :: t = [0.0], f
    p%r_p_s = 0.1; p%a_s = 10.0; p%incl = 1.5708; p%periodo = 1.0; p%t0 = 0.0
    p%u1 = 0.0; p%u2 = 0.0
    f = calcular_modelo_mandel_agol(t, p)
    IF (f(1) < 1.0_real64) THEN
      PRINT *, "[OK] Modulo de Transito validado (eclipse detectado)."
    ELSE
      PRINT *, "[ERRO] Erro no calculo do transito."
    END IF
  END SUBROUTINE test_modulo_transit

END MODULE mod_transit
