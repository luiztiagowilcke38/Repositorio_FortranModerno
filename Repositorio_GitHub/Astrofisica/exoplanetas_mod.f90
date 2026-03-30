MODULE exoplanetas_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  REAL(real64), PARAMETER :: PI = 3.14159265358979323846_real64

CONTAINS

  FUNCTION calcular_fluxo_transito(distancia_r, raio_p, u1, u2) RESULT(fluxo)
    REAL(real64), INTENT(IN) :: distancia_r, raio_p, u1, u2
    REAL(real64) :: fluxo, lambda_e, d, p
    
    d = distancia_r; p = raio_p; fluxo = 1.0_real64
    IF (d > 1.0_real64 + p) RETURN ! Fora do transito
    
    IF (d > ABS(1.0_real64 - p) .AND. d <= 1.0_real64 + p) THEN
      lambda_e = calcular_area_intersecao(d, p) / PI
    ELSE IF (d <= ABS(1.0_real64 - p)) THEN
      lambda_e = p**2
    END IF
    
    ! Aplicacao do escurecimento de limbo (Simplificacao Mandel-Agol 2002)
    fluxo = 1.0_real64 - lambda_e * (1.0_real64 - u1/3.0_real64 - u2/6.0_real64)
  END FUNCTION calcular_fluxo_transito

  FUNCTION calcular_area_intersecao(d, p) RESULT(area)
    REAL(real64), INTENT(IN) :: d, p; REAL(real64) :: area, k0, k1
    IF (d >= 1.0 + p) THEN; area = 0.0; RETURN; END IF
    IF (d <= ABS(1.0 - p)) THEN; area = PI * MIN(1.0_real64, p)**2; RETURN; END IF
    
    k0 = ACOS((p**2 + d**2 - 1.0_real64)/(2.0_real64 * p * d))
    k1 = ACOS((1.0_real64 + d**2 - p**2)/(2.0_real64 * d))
    area = p**2 * k0 + k1 - 0.5_real64 * SQRT(MAX(0.0_real64, (1.0_real64+p-d)*(p+d-1.0_real64)*(d+1.0_real64-p)*(1.0_real64+p+d)))
  END FUNCTION calcular_area_intersecao

  SUBROUTINE gerar_curva_luz_orbita(periodo, r_estrela, r_planeta, inclina, a_semi, n_pontos)
    REAL(real64), INTENT(IN) :: periodo, r_estrela, r_planeta, inclina, a_semi
    INTEGER, INTENT(IN) :: n_pontos
    REAL(real64) :: t, dt, phi, r_dist, p_ratio, u1, u2, f
    INTEGER :: i
    
    dt = periodo / REAL(n_pontos, real64); p_ratio = r_planeta / r_estrela
    u1 = 0.3_real64; u2 = 0.2_real64 ! Coeficientes tipicos
    
    DO i = 0, n_pontos - 1
      t = i * dt; phi = 2.0_real64 * PI * t / periodo
      ! Projecao da distancia no plano do ceu (Kepleriana circular)
      r_dist = a_semi * SQRT(SIN(phi)**2 + (COS(inclina)*COS(phi))**2) / r_estrela
      f = calcular_fluxo_transito(r_dist, p_ratio, u1, u2)
      ! Imprimir resultado (tempo, fluxo)
      IF (MOD(i, 10) == 0) PRINT "(F10.4, F12.8)", t, f
    END DO
  END SUBROUTINE gerar_curva_luz_orbita

END MODULE exoplanetas_mod
