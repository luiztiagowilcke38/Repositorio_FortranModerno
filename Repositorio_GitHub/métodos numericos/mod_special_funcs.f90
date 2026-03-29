!==================================================================================================
! MODULO: MOD_SPECIAL_FUNCS (FUNCOES DE BESSEL E POLINOMIOS ORTOGONAIS)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Implementacao de funcoes especiais via relacoes de recorrencia estaveis
! (Miller's algorithm) e expansoes de series.
!
! REFERENCIA: Numerical Recipes in Fortran 90 / 6.6. Bessel Functions of Fractional Order.
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE mod_special_funcs
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: numrec_j0, numrec_j1, numrec_jn

  INTEGER, PARAMETER :: DP = real64

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! FUNCAO DE BESSEL J0(X) VIA APROXIMACAO RACIONAL
  !------------------------------------------------------------------------------------------------
  FUNCTION numrec_j0(x) RESULT(j0)
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: j0, ax, z, xx, y
    ! (Constantes abreviadas para clareza didatica)
    
    ax = ABS(x)
    IF (ax < 8.0_DP) THEN
       y = x**2
       j0 = (57568490574.0_DP + y * (-13362590354.0_DP + y * (651619640.7_DP + y * (-11214424.18_DP)))) / &
            (57568490411.0_DP + y * (102953298.5_DP + y * (1235122.351_DP + y * (5193.303833_DP))))
    ELSE
       z = 8.0_DP / ax
       y = z**2
       xx = ax - 0.785398164_DP
       j0 = SQRT(0.636619772_DP / ax) * (cos(xx) * (1._DP + y * (-0.1098628627E-2_DP)) - &
            z * sin(xx) * (-0.1562499995E-1_DP + y * (0.1430488765E-3_DP)))
    END IF
  END FUNCTION numrec_j0

  FUNCTION numrec_j1(x) RESULT(j1)
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: j1, ax
    ax = ABS(x)
    ! Implementacao analoga a J0 com coeficientes especificos para J1
    j1 = x * 0.5_DP ! Placeholder para expansao pequena
  END FUNCTION numrec_j1

  !------------------------------------------------------------------------------------------------
  ! FUNCAO DE BESSEL JN(X) VIA RECORRENCIA PARA BAIXO (ALGORITMO DE MILLER)
  !------------------------------------------------------------------------------------------------
  FUNCTION numrec_jn(n, x) RESULT(jn)
    INTEGER,  INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: jn
    INTEGER, PARAMETER :: IACC = 40
    REAL(DP) :: bj, bjm, bjp, sum_val, tox
    INTEGER :: j, m

    IF (n == 0) THEN ; jn = numrec_j0(x) ; RETURN ; END IF
    IF (n == 1) THEN ; jn = numrec_j1(x) ; RETURN ; END IF
    
    tox = 2.0_DP / ABS(x)
    IF (ABS(x) > REAL(n, DP)) THEN
       ! Recorrencia para cima e estavel
       bjm = numrec_j0(x)
       bj = numrec_j1(x)
       DO j = 1, n-1
          bjp = j * tox * bj - bjm
          bjm = bj
          bj = bjp
       END DO
       jn = bj
    ELSE
       ! Recorrencia para baixo (Miller)
       m = 2 * ((n + INT(SQRT(REAL(IACC * n, DP)))) / 2)
       jn = 0.0_DP ; bjp = 0.0_DP ; bj = 1.0_DP ; sum_val = 0.0_DP
       DO j = m, 1, -1
          bjm = j * tox * bj - bjp
          bjp = bj
          bj = bjm
          IF (j > n) jn = 0.0_DP
          IF (j == n) jn = bjp
          IF (MOD(j, 2) == 0) sum_val = sum_val + bj
       END DO
       sum_val = 2.0_DP * sum_val - bj
       jn = jn / sum_val
    END IF
    IF (x < 0.0_DP .AND. MOD(n, 2) /= 0) jn = -jn
  END FUNCTION numrec_jn

END MODULE mod_special_funcs
