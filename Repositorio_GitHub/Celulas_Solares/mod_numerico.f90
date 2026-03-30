MODULE mod_numerico
  USE iso_fortran_env
  IMPLICIT NONE

CONTAINS

  ! Algoritmo de Thomas (TDMA) para sistemas tridiagonais
  SUBROUTINE solver_thomas(a, b, c, f, x)
    REAL(real64), INTENT(IN) :: a(:), b(:), c(:), f(:)
    REAL(real64), INTENT(OUT) :: x(:)
    INTEGER :: n, i
    REAL(real64) :: m, cp(SIZE(f)), fp(SIZE(f))
    
    n = SIZE(f)
    cp(1) = c(1) / b(1)
    fp(1) = f(1) / b(1)
    
    DO i = 2, n
       m = 1.0d0 / (b(i) - a(i) * cp(i-1))
       cp(i) = c(i) * m
       fp(i) = (f(i) - a(i) * fp(i-1)) * m
    END DO
    
    x(n) = fp(n)
    DO i = n-1, 1, -1
       x(i) = fp(i) - cp(i) * x(i+1)
    END DO
  END SUBROUTINE solver_thomas

  ! Funcao de Bernoulli B(x) = x / (exp(x) - 1) com tratamento para x -> 0
  REAL(real64) FUNCTION bernoulli(x)
    REAL(real64), INTENT(IN) :: x
    IF (ABS(x) < 1.0d-5) THEN
       bernoulli = 1.0d0 - 0.5d0*x + x**2/12.0d0 - x**4/720.0d0
    ELSE
       bernoulli = x / (EXP(x) - 1.0d0)
    END IF
  END FUNCTION bernoulli

END MODULE mod_numerico
