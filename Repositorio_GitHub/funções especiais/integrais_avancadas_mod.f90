MODULE integrais_avancadas_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  REAL(real64), PARAMETER :: PI = 3.14159265358979323846_real64
  REAL(real64), PARAMETER :: EPS = 1.0E-15_real64
  REAL(real64), PARAMETER :: FPMIN = 1.0E-30_real64

CONTAINS

  SUBROUTINE calcular_fresnel_total(x, s, c)
    REAL(real64), INTENT(IN) :: x; REAL(real64), INTENT(OUT) :: s, c
    REAL(real64) :: abs_x, x2, t, f, g, ang, soma_s, soma_c, termo
    INTEGER :: k
    abs_x = ABS(x)
    IF (abs_x < 1.5_real64) THEN
      soma_s = 0.0_real64; soma_c = abs_x; x2 = abs_x*abs_x; t = (PI/2.0_real64)*x2; termo = abs_x
      DO k = 1, 100
        termo = -termo * t**2 / REAL((2*k)*(2*k-1), real64)
        soma_c = soma_c + termo / REAL(4*k + 1, real64)
        IF (ABS(termo) < EPS) EXIT
      END DO
      termo = (PI/2.0_real64) * abs_x**3; soma_s = termo / 3.0_real64
      DO k = 1, 100
        termo = -termo * t**2 / REAL((2*k)*(2*k+1), real64)
        soma_s = soma_s + termo / REAL(4*k + 3, real64)
        IF (ABS(termo) < EPS) EXIT
      END DO
      s = soma_s; c = soma_c
    ELSE
      CALL avaliar_aux_lentz(abs_x, f, g)
      ang = 0.5_real64 * PI * abs_x**2
      s = 0.5_real64 * (1.0_real64 - f*COS(ang) - g*SIN(ang))
      c = 0.5_real64 * (1.0_real64 + f*SIN(ang) - g*COS(ang))
    END IF
    IF (x < 0.0_real64) THEN; s = -s; c = -c; END IF
  END SUBROUTINE calcular_fresnel_total

  SUBROUTINE avaliar_aux_lentz(x, f, g)
    REAL(real64), INTENT(IN) :: x; REAL(real64), INTENT(OUT) :: f, g
    COMPLEX(real64) :: h, c, d, delta, a, b, z_inv
    INTEGER :: n
    z_inv = CMPLX(1.0_real64 / (2.0_real64 * PI * x**2), 0.0_real64, real64)
    h = CMPLX(1.0_real64, 0.0_real64, real64)
    IF (ABS(h) < FPMIN) h = CMPLX(FPMIN, 0.0_real64, real64)
    c = h; d = 0.0_real64
      DO n = 1, 100
        a = -CMPLX(REAL((2*n-1)**2, real64), 0.0_real64, real64) * z_inv
        b = CMPLX(1.0_real64, 0.0_real64, real64)
        d = b + a * d
        IF (ABS(d) < FPMIN) d = CMPLX(FPMIN, 0.0_real64, real64)
        c = b + a / c
        IF (ABS(c) < FPMIN) c = CMPLX(FPMIN, 0.0_real64, real64)
        d = 1.0_real64 / d
        delta = c * d; h = h * delta
        IF (ABS(REAL(delta)-1.0_real64) < EPS .AND. ABS(AIMAG(delta)) < EPS) EXIT
      END DO
      h = h * CMPLX(1.0_real64, 1.0_real64, real64) / CMPLX(2.0_real64 * PI * x, 0.0_real64, real64)
      f = REAL(h, real64)
      g = AIMAG(h)
    END SUBROUTINE avaliar_aux_lentz
  
    FUNCTION calcular_si(x) RESULT(res)
      REAL(real64), INTENT(IN) :: x; REAL(real64) :: res, t, t2; INTEGER :: i
      res = x; t = x; t2 = x**2
      IF (ABS(x) < 4.0_real64) THEN
        DO i = 1, 60
          t = -t * t2 / REAL((2*i)*(2*i+1), real64)
          res = res + t
          IF (ABS(t) < EPS) EXIT
        END DO
    ELSE
      res = PI/2.0_real64 - (COS(x)/x)*(1.0_real64 - 2.0_real64/x**2 + 24.0_real64/x**4) &
                          - (SIN(x)/x)*(1.0_real64/x - 6.0_real64/x**3 + 120.0_real64/x**5)
    END IF
  END FUNCTION calcular_si

END MODULE integrais_avancadas_mod
