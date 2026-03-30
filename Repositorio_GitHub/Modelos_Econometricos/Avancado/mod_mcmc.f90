MODULE mod_mcmc
  USE iso_fortran_env
  USE mod_aleatorio
  IMPLICIT NONE

CONTAINS

  SUBROUTINE metropolis_hastings(log_post, t_ini, n_iter, cadeia)
     ABSTRACT INTERFACE
        REAL(real64) FUNCTION log_post(theta); USE iso_fortran_env; REAL(real64), INTENT(IN) :: theta(:); END FUNCTION log_post
     END INTERFACE
     PROCEDURE(log_post) :: f_logp
     REAL(real64), INTENT(IN) :: t_ini(:)
     INTEGER, INTENT(IN) :: n_iter
     REAL(real64), INTENT(OUT) :: cadeia(:,:)
     
     REAL(real64) :: t_atual(SIZE(t_ini)), t_prop(SIZE(t_ini)), lp_at, lp_prop, alpha
     INTEGER :: i
     
     t_atual = t_ini
     lp_at = f_logp(t_atual)
     
     DO i = 1, n_iter
        t_prop = t_atual + 0.1d0 * aleatorio_normal() ! Proposta aleatoria simples
        lp_prop = f_logp(t_prop)
        
        alpha = MIN(0.0d0, lp_prop - lp_at)
        IF (LOG(aleatorio_uniforme()) < alpha) THEN
           t_atual = t_prop
           lp_at = lp_prop
        END IF
        cadeia(:, i) = t_atual
     END DO
  END SUBROUTINE metropolis_hastings

END MODULE mod_mcmc
