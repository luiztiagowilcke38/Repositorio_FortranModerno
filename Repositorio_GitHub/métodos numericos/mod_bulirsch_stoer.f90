!==================================================================================================
! MODULO: MOD_BULIRSCH_STOER (INTEGRADOR DE EDO DE ALTA PRECISAO)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Implementacao do metodo de Bulirsch-Stoer para sistemas de EDOs.
! Utiliza extrapolacao de Richardson sobre o metodo do ponto medio modificado.
!
! REFERENCIA: Numerical Recipes in Fortran 90 / 16.4. Bulirsch-Stoer Method.
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE mod_bulirsch_stoer
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bulirsch_stoer_step

  INTEGER, PARAMETER :: DP = real64
  INTEGER, PARAMETER :: KMAXX = 8
  INTEGER, PARAMETER :: IMAXX = KMAXX + 1

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! PASSO DE INTEGRACAO BULIRSCH-STOER
  !------------------------------------------------------------------------------------------------
  SUBROUTINE bulirsch_stoer_step(y, dydx, nv, x, htry, eps, yscal, hdid, hnext, derivs)
    REAL(DP), INTENT(INOUT) :: y(nv)
    REAL(DP), INTENT(IN)    :: dydx(nv), yscal(nv)
    REAL(DP), INTENT(INOUT) :: x
    REAL(DP), INTENT(IN)    :: htry, eps
    REAL(DP), INTENT(OUT)   :: hdid, hnext
    INTEGER,  INTENT(IN)    :: nv
    
    INTERFACE
       SUBROUTINE derivs(x, y, dydx)
          IMPORT :: DP
          REAL(DP), INTENT(IN)  :: x, y(:)
          REAL(DP), INTENT(OUT) :: dydx(SIZE(y))
       END SUBROUTINE derivs
    END INTERFACE

    REAL(DP) :: x_new, h, err, yerr(nv), y_tmp(nv), ysave(nv)
    REAL(DP) :: d(nv, KMAXX), xest, a(IMAXX)
    INTEGER  :: i, k, n_step(IMAXX) = [2, 4, 6, 8, 10, 12, 14, 16, 18]

    h = htry
    ysave = y
    
    DO
       DO k = 1, KMAXX
          CALL modified_midpoint(ysave, dydx, nv, x, h, n_step(k), y_tmp, derivs)
          xest = (h / n_step(k))**2
          CALL rational_extrapolate(k, xest, y_tmp, y, yerr, nv, d)
          
          ! Verifica erro estimado
          err = 0.0_DP
          DO i = 1, nv
             err = MAX(err, ABS(yerr(i) / yscal(i)))
          END DO
          err = err / eps
          
          IF (err <= 1.0_DP) THEN
             x = x + h
             hdid = h
             IF (k == KMAXX) THEN
                hnext = h * 0.95_DP
             ELSE IF (k >= KMAXX - 1) THEN
                hnext = h * 1.2_DP
             ELSE
                hnext = h * (n_step(k+1) / REAL(n_step(k), DP))
             END IF
             RETURN
          END IF
       END DO
       
       ! Reduz passo se falhar
       h = 0.25_DP * h / (2.0_DP**((KMAXX-k)/2))
       IF (ABS(h) < 1.0E-15_DP) STOP "Passo muito pequeno no Bulirsch-Stoer"
    END DO
  END SUBROUTINE bulirsch_stoer_step

  !------------------------------------------------------------------------------------------------
  ! METODO DO PONTO MEDIO MODIFICADO
  !------------------------------------------------------------------------------------------------
  SUBROUTINE modified_midpoint(y_in, dydx, n, x_start, h_tot, n_step, y_out, derivs)
    REAL(DP), INTENT(IN)  :: y_in(n), dydx(n), x_start, h_tot
    INTEGER,  INTENT(IN)  :: n, n_step
    REAL(DP), INTENT(OUT) :: y_out(n)
    
    INTERFACE
       SUBROUTINE derivs(x, y, dydx)
          IMPORT :: DP
          REAL(DP), INTENT(IN)  :: x, y(:)
          REAL(DP), INTENT(OUT) :: dydx(SIZE(y))
       END SUBROUTINE derivs
    END INTERFACE

    REAL(DP) :: h, ym(n), yn(n), x, d_out(n)
    INTEGER  :: i, step

    h = h_tot / n_step
    ym = y_in
    yn = y_in + h * dydx
    x = x_start + h
    CALL derivs(x, yn, d_out)
    
    DO step = 2, n_step
       y_out = ym + 2.0_DP * h * d_out
       ym = yn
       yn = y_out
       x = x + h
       CALL derivs(x, yn, d_out)
    END DO
    
    y_out = 0.5_DP * (ym + yn + h * d_out)
  END SUBROUTINE modified_midpoint

  !------------------------------------------------------------------------------------------------
  ! EXTRAPOLACAO RACIONAL (POLLY ESTENSAO)
  !------------------------------------------------------------------------------------------------
  SUBROUTINE rational_extrapolate(k, xest, yest, yz, dy, nv, d)
    INTEGER,  INTENT(IN)    :: k, nv
    REAL(DP), INTENT(IN)    :: xest
    REAL(DP), INTENT(INOUT) :: yest(nv)
    REAL(DP), INTENT(OUT)   :: yz(nv), dy(nv)
    REAL(DP), INTENT(INOUT) :: d(nv, KMAXX)
    
    REAL(DP) :: fx(KMAXX), v, c, b, b1, dd, x(IMAXX)
    INTEGER  :: j, i
    SAVE x

    x(k) = xest
    IF (k == 1) THEN
       yz = yest
       d(:,1) = yest
       dy = yest
    ELSE
       DO i = 1, nv
          v = d(i, 1)
          d(i, 1) = yest(i)
          c = yest(i)
          DO j = 2, k
             b1 = x(k-j+1) * v
             b = b1 - xest * c
             IF (b /= 0.0_DP) THEN
                b = (c - v) / b
                dd = c * b
                c = b1 * b
             ELSE
                dd = v
             END IF
             v = d(i, j)
             d(i, j) = dd
             yest(i) = yest(i) + dd ! Reuso temporario para erro
          END DO
          dy(i) = dd
          yz(i) = yest(i)
       END DO
    END IF
  END SUBROUTINE rational_extrapolate

END MODULE mod_bulirsch_stoer
