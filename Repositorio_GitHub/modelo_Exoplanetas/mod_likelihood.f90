MODULE mod_likelihood
  USE iso_fortran_env, ONLY: real64
  USE mod_constants, ONLY: PI_CONST
  USE mod_priors, ONLY: estrutura_prior, calcular_log_prior
  IMPLICIT NONE

CONTAINS

  FUNCTION calcular_log_likelihood(dados, modelo, incerteza) RESULT(ln_l)
    REAL(real64), DIMENSION(:), INTENT(IN) :: dados, modelo, incerteza
    REAL(real64) :: ln_l, chi2
    INTEGER :: n
    n = SIZE(dados)
    chi2 = SUM(((dados - modelo) / incerteza)**2)
    ln_l = -0.5_real64 * (chi2 + SUM(LOG(2.0_real64 * PI_CONST * incerteza**2)))
  END FUNCTION calcular_log_likelihood

  FUNCTION calcular_kernel_matern32(t1, t2, rho, sigma) RESULT(k)
    REAL(real64), INTENT(IN) :: t1, t2, rho, sigma
    REAL(real64) :: k, d, raiz3
    raiz3 = SQRT(3.0_real64)
    d = ABS(t1 - t2)
    k = sigma**2 * (1.0_real64 + raiz3*d/rho) * EXP(-raiz3*d/rho)
  END FUNCTION calcular_kernel_matern32

  FUNCTION calcular_log_posterior(parametros, priors, dados, tempo, erro_f, n_par) RESULT(ln_p)
    REAL(real64), DIMENSION(:), INTENT(IN) :: parametros, dados, tempo, erro_f
    TYPE(estrutura_prior), DIMENSION(:), INTENT(IN) :: priors
    INTEGER, INTENT(IN) :: n_par
    REAL(real64) :: ln_p, ln_prior_total
    INTEGER :: i
    
    ln_prior_total = 0.0_real64
    DO i = 1, n_par
      ln_prior_total = ln_prior_total + calcular_log_prior(parametros(i), priors(i))
    END DO
    
    IF (ln_prior_total < -1.0E20_real64) THEN
      ln_p = -HUGE(1.0_real64)
    ELSE
      ! Aqui chamariamos o mod_transit para gerar o modelo baseado nos parametros
      ln_p = ln_prior_total ! + calcular_log_likelihood(...)
    END IF
  END FUNCTION calcular_log_posterior

END MODULE mod_likelihood
