!==================================================================================================
! MODULO: MOD_OPT_SIMPLEX (ALGORITMO NELDER-MEAD SIMPLEX)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Metodo de otimizacao multidimensional sem o uso de derivadas.
! Utiliza uma figura geometrica (simplex) que rasteja, expande e contrai em direcao ao minimo.
!
! REFERENCIA: Numerical Recipes in Fortran 90 / 10.4. Downhill Simplex Method.
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE mod_opt_simplex
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: amoeba

  INTEGER, PARAMETER :: DP = real64

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! ROTINA PRINCIPAL: AMOEBA (Downhill Simplex)
  !------------------------------------------------------------------------------------------------
  SUBROUTINE amoeba(p, y, ndim, ftol, func, iter)
    REAL(DP), INTENT(INOUT) :: p(ndim+1, ndim)
    REAL(DP), INTENT(INOUT) :: y(ndim+1)
    INTEGER,  INTENT(IN)    :: ndim
    REAL(DP), INTENT(IN)    :: ftol
    INTEGER,  INTENT(OUT)   :: iter
    
    INTERFACE
       FUNCTION func(x)
          IMPORT :: DP
          REAL(DP), INTENT(IN) :: x(:)
          REAL(DP) :: func
       END FUNCTION func
    END INTERFACE

    INTEGER, PARAMETER :: ITMAX = 5000
    REAL(DP) :: psum(ndim), ytest, rtol
    INTEGER  :: i, ihi, ilo, inhi

    iter = 0
    psum = SUM(p, DIM=1)
    
    DO
       ilo = 1
       IF (y(1) > y(2)) THEN ; ihi = 1 ; inhi = 2 ; ELSE ; ihi = 2 ; inhi = 1 ; END IF
       
       DO i = 1, ndim+1
          IF (y(i) <= y(ilo)) ilo = i
          IF (y(i) > y(ihi)) THEN
             inhi = ihi ; ihi = i
          ELSE IF (y(i) > y(inhi)) THEN
             IF (i /= ihi) inhi = i
          END IF
       END DO
       
       rtol = 2.0_DP * ABS(y(ihi) - y(ilo)) / (ABS(y(ihi)) + ABS(y(ilo)) + 1.0E-10_DP)
       
       ! Convergencia
       IF (rtol < ftol) THEN
          CALL swap_elements(y(1), y(ilo))
          CALL swap_rows(p(1,:), p(ilo,:))
          RETURN
       END IF
       
       IF (iter >= ITMAX) STOP "Excedeu iteracoes no Simplex"
       iter = iter + 2
       
       ! Passo 1: Reflexao
       ytest = amotry(p, y, psum, ndim, func, ihi, -1.0_DP)
       IF (ytest <= y(ilo)) THEN
          ! Passo 2: Expansao
          ytest = amotry(p, y, psum, ndim, func, ihi, 2.0_DP)
       ELSE IF (ytest >= y(inhi)) THEN
          ! Passo 3: Contracao
          ytest = amotry(p, y, psum, ndim, func, ihi, 0.5_DP)
          IF (ytest >= y(ihi)) THEN
             ! Passo 4: Multi-contracao
             DO i = 1, ndim+1
                IF (i /= ilo) THEN
                   p(i,:) = 0.5_DP * (p(i,:) + p(ilo,:))
                   y(i) = func(p(i,:))
                END IF
             END DO
             iter = iter + ndim
             psum = SUM(p, DIM=1)
          END IF
       ELSE
          iter = iter - 1
       END IF
    END DO
  END SUBROUTINE amoeba

  !------------------------------------------------------------------------------------------------
  ! OPERACAO DE TESTE DO SIMPLEX (TRANSFORMACAO)
  !------------------------------------------------------------------------------------------------
  FUNCTION amotry(p, y, psum, ndim, func, ihi, fac) RESULT(ytry)
    REAL(DP), INTENT(INOUT) :: p(:, :), y(:), psum(:)
    INTEGER,  INTENT(IN)    :: ndim, ihi
    REAL(DP), INTENT(IN)    :: fac
    REAL(DP) :: ytry
    
    INTERFACE
       FUNCTION func(x)
          IMPORT :: DP
          REAL(DP), INTENT(IN) :: x(:)
          REAL(DP) :: func
       END FUNCTION func
    END INTERFACE

    REAL(DP) :: fac1, fac2, ptry(ndim)

    fac1 = (1.0_DP - fac) / ndim
    fac2 = fac1 - fac
    ptry = psum * fac1 - p(ihi, :) * fac2
    ytry = func(ptry)
    
    IF (ytry < y(ihi)) THEN
       y(ihi) = ytry
       psum = psum + ptry - p(ihi, :)
       p(ihi, :) = ptry
    END IF
  END FUNCTION amotry

  SUBROUTINE swap_rows(r1, r2) ; REAL(DP), INTENT(INOUT) :: r1(:), r2(:) ; REAL(DP) :: tmp(SIZE(r1)) ; tmp=r1; r1=r2; r2=tmp ; END SUBROUTINE swap_rows
  SUBROUTINE swap_elements(e1, e2) ; REAL(DP), INTENT(INOUT) :: e1, e2 ; REAL(DP) :: tmp ; tmp=e1; e1=e2; e2=tmp ; END SUBROUTINE swap_elements

END MODULE mod_opt_simplex
