MODULE harmonicos_profissionais_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  REAL(real64), PARAMETER :: VALOR_PI = 3.14159265358979323846_real64

CONTAINS

  SUBROUTINE avaliar_ylm_completo(grau_l, ordem_m, theta_coord, fi_coord, &
                                  y_real, y_imag)
    INTEGER, INTENT(IN) :: grau_l, ordem_m
    REAL(real64), INTENT(IN) :: theta_coord, fi_coord
    REAL(real64), INTENT(OUT) :: y_real, y_imag
    REAL(real64) :: norm_f, p_lm, cosseno_t, m_fi
    
    cosseno_t = COS(theta_coord)
    p_lm = calcular_legendre_estavel(grau_l, ABS(ordem_m), cosseno_t)
    
    norm_f = SQRT(((2*grau_l + 1) / (4.0_real64 * VALOR_PI)) * &
             EXP(log_fatoreal(grau_l - ABS(ordem_m)) - log_fatoreal(grau_l + ABS(ordem_m))))
    
    m_fi = REAL(ordem_m, real64) * fi_coord
    y_real = norm_f * p_lm * COS(m_fi)
    y_imag = norm_f * p_lm * SIN(m_fi)
    
    IF (ordem_m < 0 .AND. MOD(ABS(ordem_m), 2) /= 0) THEN
      y_real = -y_real; y_imag = -y_imag
    END IF
  END SUBROUTINE avaliar_ylm_completo

  FUNCTION calcular_legendre_estavel(l, m, x) RESULT(plm)
    INTEGER, INTENT(IN) :: l, m
    REAL(real64), INTENT(IN) :: x
    REAL(real64) :: plm, p_mm, p_mmp1, somx2, fat
    INTEGER :: j, ll
    
    p_mm = 1.0_real64
    IF (m > 0) THEN
      somx2 = SQRT((1.0_real64 - x) * (1.0_real64 + x))
      fat = 1.0_real64
      DO j = 1, m
        p_mm = -p_mm * fat * somx2
        fat = fat + 2.0_real64
      END DO
    END IF
    
    IF (l == m) THEN
      plm = p_mm
    ELSE
      p_mmp1 = x * (2*m + 1) * p_mm
      IF (l == m + 1) THEN
        plm = p_mmp1
      ELSE
        DO ll = m + 2, l
          plm = (x * (2*ll - 1) * p_mmp1 - (ll + m - 1) * p_mm) / (ll - m)
          p_mm = p_mmp1
          p_mmp1 = plm
        END DO
      END IF
    END IF
  END FUNCTION calcular_legendre_estavel

  FUNCTION log_fatoreal(n_ent) RESULT(log_f)
    INTEGER, INTENT(IN) :: n_ent
    REAL(real64) :: log_f
    INTEGER :: i
    log_f = 0.0_real64
    IF (n_ent <= 1) RETURN
    DO i = 2, n_ent
      log_f = log_f + LOG(REAL(i, real64))
    END DO
  END FUNCTION log_fatoreal

END MODULE harmonicos_profissionais_mod
