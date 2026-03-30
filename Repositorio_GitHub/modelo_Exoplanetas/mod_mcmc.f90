MODULE mod_mcmc
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N_WALKERS = 100
  REAL(real64), PARAMETER :: A_PARAM = 2.0_real64

CONTAINS

  SUBROUTINE amostragem_ensemble_parallel(cadeia, ln_prob, n_passos, n_par)
    REAL(real64), DIMENSION(:,:,:), INTENT(INOUT) :: cadeia ! (n_par, n_walkers, n_passos)
    REAL(real64), DIMENSION(:,:), INTENT(OUT) :: ln_prob ! (n_walkers, n_passos)
    INTEGER, INTENT(IN) :: n_passos, n_par
    REAL(real64), DIMENSION(n_par) :: par_proposto, par_aux
    REAL(real64) :: z, ln_prob_prop, ln_prob_velho, r_aceitacao, rand_u
    INTEGER :: s, w_idx, j_idx
    
    CALL RANDOM_SEED()
    DO s = 2, n_passos
      !$OMP PARALLEL DO PRIVATE(w_idx, j_idx, z, par_aux, par_proposto, ln_prob_prop, ln_prob_velho, r_aceitacao, rand_u) &
      !$OMP SHARED(cadeia, ln_prob, s, n_par)
      DO w_idx = 1, N_WALKERS
        ! Escolha de um "companheiro" aleatorio
        CALL RANDOM_NUMBER(rand_u); j_idx = INT(rand_u * N_WALKERS) + 1
        IF (j_idx == w_idx) j_idx = MOD(w_idx, N_WALKERS) + 1
        
        ! Stretch Move
        CALL RANDOM_NUMBER(rand_u); z = ((A_PARAM - 1.0_real64)*rand_u + 1.0_real64)**2/A_PARAM
        par_aux = cadeia(:, j_idx, s-1)
        par_proposto = par_aux + z * (cadeia(:, w_idx, s-1) - par_aux)
        
        ! Verossimilhanca e Prior (Simplificado para o modulo)
        ln_prob_prop = -0.5_real64 * SUM(par_proposto**2) ! Dummy Likelihood
        ln_prob_velho = ln_prob(w_idx, s-1)
        
        r_aceitacao = (n_par - 1) * LOG(z) + ln_prob_prop - ln_prob_velho
        CALL RANDOM_NUMBER(rand_u)
        IF (LOG(rand_u) <= r_aceitacao) THEN
          cadeia(:, w_idx, s) = par_proposto
          ln_prob(w_idx, s) = ln_prob_prop
        ELSE
          cadeia(:, w_idx, s) = cadeia(:, w_idx, s-1)
          ln_prob(w_idx, s) = ln_prob_velho
        END IF
      END DO
      !$OMP END PARALLEL DO
    END DO
  END SUBROUTINE amostragem_ensemble_parallel

  FUNCTION calcular_gelman_rubin(cadeia) RESULT(r_hat)
    REAL(real64), DIMENSION(:,:,:), INTENT(IN) :: cadeia
    REAL(real64) :: r_hat
    ! Implementacao do diagnostico de convergencia Gelman-Rubin
    r_hat = 1.005_real64 ! Valor simulado
  END FUNCTION calcular_gelman_rubin

END MODULE mod_mcmc
