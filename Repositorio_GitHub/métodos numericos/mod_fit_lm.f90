!==================================================================================================
! MODULO: MOD_FIT_LM (ALGORITMO DE LEVENBERG-MARQUARDT DE ALTA FIDELIDADE)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Implementacao robusta para ajuste de curvas nao-lineares. Combina o metodo do
! gradiente descendente com o metodo de Gauss-Newton via parametro de amortecimento.
!
! REFERENCIA: Numerical Recipes in Fortran 90 / 22. Levenberg-Marquardt Method.
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE mod_fit_lm
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fitting_lm, lm_control_type

  INTEGER, PARAMETER :: DP = real64

  TYPE :: lm_control_type
     REAL(DP) :: tol = 1.0E-8_DP
     REAL(DP) :: lambda_init = 0.001_DP
     REAL(DP) :: lambda_factor = 10.0_DP
     INTEGER  :: max_iter = 100
  END TYPE lm_control_type

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! INTERFACE PRINCIPAL DE AJUSTE
  !------------------------------------------------------------------------------------------------
  SUBROUTINE fitting_lm(x, y, sigma, n_data, params, n_params, func, control, chi_square)
    REAL(DP), INTENT(IN)    :: x(n_data), y(n_data), sigma(n_data)
    INTEGER,  INTENT(IN)    :: n_data, n_params
    REAL(DP), INTENT(INOUT) :: params(n_params)
    TYPE(lm_control_type), INTENT(IN) :: control
    REAL(DP), INTENT(OUT)   :: chi_square
    
    INTERFACE
       SUBROUTINE func(x_val, p, y_model, dy_dp)
          IMPORT :: DP
          REAL(DP), INTENT(IN)  :: x_val, p(:)
          REAL(DP), INTENT(OUT) :: y_model, dy_dp(SIZE(p))
       END SUBROUTINE func
    END INTERFACE

    REAL(DP) :: alpha(n_params, n_params), beta(n_params)
    REAL(DP) :: covar(n_params, n_params), d_params(n_params), params_trial(n_params)
    REAL(DP) :: lambda, chi_old, chi_new
    INTEGER  :: iter

    ! Inicializacao
    lambda = control%lambda_init
    CALL MRQC_COEFF(x, y, sigma, n_data, params, n_params, func, alpha, beta, chi_old)
    
    DO iter = 1, control%max_iter
       covar = alpha
       ! Aplicacao do fator de Levenberg na diagonal
       DO j = 1, n_params ; covar(j,j) = alpha(j,j) * (1.0_DP + lambda) ; END DO
       
       ! Resolve o sistema linear Covar * d_params = beta
       CALL GAUSS_JORDAN_SOLVER(covar, n_params, beta, d_params)
       
       params_trial = params + d_params
       CALL MRQC_COEFF(x, y, sigma, n_data, params_trial, n_params, func, covar, d_params, chi_new)
       
       IF (chi_new < chi_old) THEN
          ! Aceita a melhoria
          lambda = lambda / control%lambda_factor
          chi_old = chi_new
          alpha = covar
          beta = d_params
          params = params_trial
          IF (ABS(chi_old - chi_new) < control%tol) EXIT
       ELSE
          ! Rejeita e aumenta o amortecimento (mais para o gradiente)
          lambda = lambda * control%lambda_factor
       END IF
    END DO
    chi_square = chi_old
  END SUBROUTINE fitting_lm

  !------------------------------------------------------------------------------------------------
  ! CALCULO DOS COEFICIENTES CURVATURA (ALPHA) E GRADIENTE (BETA)
  !------------------------------------------------------------------------------------------------
  SUBROUTINE MRQC_COEFF(x, y, sigma, n_data, p, m, func, alpha, beta, chisq)
    REAL(DP), INTENT(IN)  :: x(:), y(:), sigma(:)
    INTEGER,  INTENT(IN)  :: n_data, m
    REAL(DP), INTENT(IN)  :: p(m)
    REAL(DP), INTENT(OUT) :: alpha(m,m), beta(m), chisq
    
    INTERFACE
       SUBROUTINE func(x_val, p, y_model, dy_dp)
          IMPORT :: DP
          REAL(DP), INTENT(IN)  :: x_val, p(:)
          REAL(DP), INTENT(OUT) :: y_model, dy_dp(SIZE(p))
       END SUBROUTINE func
    END INTERFACE

    REAL(DP) :: y_mod, dy_dp(m), sig2_inv, dy
    INTEGER  :: i, k

    alpha = 0.0_DP ; beta = 0.0_DP ; chisq = 0.0_DP
    DO i = 1, n_data
       CALL func(x(i), p, y_mod, dy_dp)
       sig2_inv = 1.0_DP / (sigma(i)**2)
       dy = y(i) - y_mod
       DO j = 1, m
          DO k = 1, j
             alpha(j,k) = alpha(j,k) + dy_dp(j) * dy_dp(k) * sig2_inv
          END DO
          beta(j) = beta(j) + dy * dy_dp(j) * sig2_inv
       END DO
       chisq = chisq + dy**2 * sig2_inv
    END DO
    ! Simetria da Hessiana aproximada
    DO j = 2, m ; DO k = 1, j-1 ; alpha(k,j) = alpha(j,k) ; END DO ; END DO
  END SUBROUTINE MRQC_COEFF

  !------------------------------------------------------------------------------------------------
  ! SOLUCIONADOR LINEAR ROBUSTO (GAUSS-JORDAN COM PIVOTEAMENTO)
  !------------------------------------------------------------------------------------------------
  SUBROUTINE GAUSS_JORDAN_SOLVER(A, n, b, x_sol)
    REAL(DP), INTENT(INOUT) :: A(n,n), b(n)
    INTEGER,  INTENT(IN)    :: n
    REAL(DP), INTENT(OUT)   :: x_sol(n)
    INTEGER :: i, j, k, pivot_row
    REAL(DP) :: pivot_val, factor

    DO k = 1, n
       ! Pivoteamento parcial
       pivot_row = k ; pivot_val = ABS(A(k,k))
       DO i = k+1, n
          IF (ABS(A(i,k)) > pivot_val) THEN
             pivot_val = ABS(A(i,k)) ; pivot_row = i
          END IF
       END DO
       IF (pivot_row /= k) THEN
          CALL swap_rows(A(k,:), A(pivot_row,:))
          CALL swap_elements(b(k), b(pivot_row))
       END IF
       
       ! Normalizacao
       factor = A(k,k)
       A(k,:) = A(k,:) / factor
       b(k) = b(k) / factor
       
       ! Eliminacao
       DO i = 1, n
          IF (i /= k) THEN
             factor = A(i,k)
             A(i,:) = A(i,:) - factor * A(k,:)
             b(i) = b(i) - factor * b(k)
          END IF
       END DO
    END DO
    x_sol = b
  END SUBROUTINE GAUSS_JORDAN_SOLVER

  SUBROUTINE swap_rows(r1, r2)
    REAL(DP), INTENT(INOUT) :: r1(:), r2(:)
    REAL(DP) :: tmp(SIZE(r1))
    tmp = r1 ; r1 = r2 ; r2 = tmp
  END SUBROUTINE swap_rows

  SUBROUTINE swap_elements(e1, e2)
    REAL(DP), INTENT(INOUT) :: e1, e2
    REAL(DP) :: tmp
    tmp = e1 ; e1 = e2 ; e2 = tmp
  END SUBROUTINE swap_elements

END MODULE mod_fit_lm
