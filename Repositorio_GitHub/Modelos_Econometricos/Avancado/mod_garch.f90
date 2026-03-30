MODULE mod_garch
  USE iso_fortran_env
  IMPLICIT NONE

CONTAINS

  ! Filtro de Variancia GARCH(1,1)
  SUBROUTINE filtrar_garch(omega, alfa, beta, retornos, sigma2)
    REAL(real64), INTENT(IN) :: omega, alfa, beta, retornos(:)
    REAL(real64), INTENT(OUT) :: sigma2(:)
    INTEGER :: t, n
    
    n = SIZE(retornos)
    sigma2 = 0.0d0
    sigma2(1) = omega / (1.0d0 - alfa - beta) ! Variancia incondicional
    
    DO t = 2, n
       sigma2(t) = omega + alfa * retornos(t-1)**2 + beta * sigma2(t-1)
    END DO
  END SUBROUTINE filtrar_garch

  ! Log-Verossimilhanca GARCH(1,1) com t-Student
  REAL(real64) FUNCTION logL_garch_student(theta, ret)
    REAL(real64), INTENT(IN) :: theta(4) ! [omega, alfa, beta, nu]
    REAL(real64), INTENT(IN) :: ret(:)
    REAL(real64), ALLOCATABLE :: sig2(:)
    REAL(real64) :: nu, val
    INTEGER :: t, n
    
    n = SIZE(ret)
    ALLOCATE(sig2(n))
    CALL filtrar_garch(theta(1), theta(2), theta(3), ret, sig2)
    
    nu = theta(4)
    val = 0.0d0
    DO t = 1, n
       val = val + LOG(GAMMA((nu+1.0d0)/2.0d0)) - LOG(GAMMA(nu/2.0d0)) - 0.5d0*LOG(3.14159d0*(nu-2.0d0)) &
             - 1.0d0/2.0d0*LOG(sig2(t)) - (nu+1.0d0)/2.0d0 * LOG(1.0d0 + ret(t)**2 / ((nu-2.0d0)*sig2(t)))
    END DO
    logL_garch_student = val
  END FUNCTION logL_garch_student

END MODULE mod_garch
