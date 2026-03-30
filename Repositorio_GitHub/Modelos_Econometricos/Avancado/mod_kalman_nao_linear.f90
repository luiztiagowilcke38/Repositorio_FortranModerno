MODULE mod_kalman_nao_linear
  USE iso_fortran_env
  USE mod_numerico ! Reutilizando helper de Cholesky se disponivel
  IMPLICIT NONE

CONTAINS

  ! Filtro de Kalman Unscented (UKF) - Um passo
  SUBROUTINE ukf_passo(f_trans, h_obs, x, P, Q, R, y, n, W_m, W_c, kappa)
     ABSTRACT INTERFACE
        SUBROUTINE f_trans(x, x_novo); USE iso_fortran_env; REAL(real64), INTENT(IN) :: x(:); REAL(real64), INTENT(OUT) :: x_novo(:); END SUBROUTINE f_trans
        SUBROUTINE h_obs(x, y_novo); USE iso_fortran_env; REAL(real64), INTENT(IN) :: x(:); REAL(real64), INTENT(OUT) :: y_novo(:); END SUBROUTINE h_obs
     END INTERFACE
     
     PROCEDURE(f_trans) :: f_func
     PROCEDURE(h_obs) :: h_func
     REAL(real64), INTENT(INOUT) :: x(:), P(:,:)
     REAL(real64), INTENT(IN) :: Q(:,:), R(:,:), y(:), W_m(:), W_c(:)
     INTEGER, INTENT(IN) :: n
     REAL(real64), INTENT(IN) :: kappa
     
     REAL(real64) :: chi(n, 2*n+1), chi_pred(n, 2*n+1), zeta(SIZE(y), 2*n+1)
     REAL(real64) :: x_pred(n), y_pred(SIZE(y)), P_pred(n,n), P_yy(SIZE(y), SIZE(y)), P_xy(n, SIZE(y))
     INTEGER :: i
     
     ! 1. Geracao de Pontos Sigma (via Cholesky de P)
     ! ... (Omitido p/ brevidade no esqueleto, mas implementado no capitulo real) ...
     
     ! 2. Propagacao Nao-Linear
     DO i = 1, 2*n+1
        CALL f_func(chi(:,i), chi_pred(:,i))
     END DO
     
     ! 3. Reconstrucao de Media e Covariancia
     x_pred = 0.0d0
     DO i = 1, 2*n+1
        x_pred = x_pred + W_m(i) * chi_pred(:,i)
     END DO
     ! ... (P_pred follow similar logic) ...
     
  END SUBROUTINE ukf_passo

END MODULE mod_kalman_nao_linear
