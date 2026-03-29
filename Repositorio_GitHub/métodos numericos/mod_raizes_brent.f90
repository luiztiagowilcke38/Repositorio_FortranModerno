!==================================================================================================
! MODULO: MOD_RAIZES_BRENT (METODO DE BRENT PARA BUSCA DE RAIZES)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Algoritmo hibrido que combina bissecao, secante e interpolacao quadratica
! inversa. E considerado um dos mais robustos para busca de raizes em 1D.
!
! REFERENCIA: Numerical Recipes in Fortran 90 / 9.3. Van Wijngaarden-Dekker-Brent Method.
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE mod_raizes_brent
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: zbrent

  INTEGER, PARAMETER :: DP = real64

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! FUNCAO ZBRENT: Encontra a raiz no intervalo [x1, x2]
  !------------------------------------------------------------------------------------------------
  FUNCTION zbrent(func, x1, x2, tol, ierr) RESULT(root)
    REAL(DP), INTENT(IN)  :: x1, x2, tol
    INTEGER,  INTENT(OUT) :: ierr
    REAL(DP) :: root
    
    INTERFACE
       FUNCTION func(x)
          IMPORT :: DP
          REAL(DP), INTENT(IN) :: x
          REAL(DP) :: func
       END FUNCTION func
    END INTERFACE

    INTEGER, PARAMETER :: MAXIT = 100
    REAL(DP), PARAMETER :: EPS = EPSILON(1.0_DP)
    REAL(DP) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
    INTEGER  :: iter

    a = x1 ; b = x2 ; fa = func(a) ; fb = func(b)
    ierr = 0

    ! Verifica se o intervalo contem a raiz (mudanca de sinal)
    IF ((fa > 0.0_DP .AND. fb > 0.0_DP) .OR. (fa < 0.0_DP .AND. fb < 0.0_DP)) THEN
       ierr = 1 ! Erro: Raiz nao cercada
       root = b
       RETURN
    END IF

    c = a ; fc = fa ; d = b - a ; e = d
    DO iter = 1, MAXIT
       IF ((fb > 0.0_DP .AND. fc > 0.0_DP) .OR. (fb < 0.0_DP .AND. fc < 0.0_DP)) THEN
          c = a ; fc = fa ; d = b - a ; e = d
       END IF
       IF (ABS(fc) < ABS(fb)) THEN
          a = b ; b = c ; c = a
          fa = fb ; fb = fc ; fc = fa
       END IF
       
       tol1 = 2.0_DP * EPS * ABS(b) + 0.5_DP * tol
       xm = 0.5_DP * (c - b)
       
       ! Convergencia atingida?
       IF (ABS(xm) <= tol1 .OR. fb == 0.0_DP) THEN
          root = b
          RETURN
       END IF
       
       IF (ABS(e) >= tol1 .AND. ABS(fa) > ABS(fb)) THEN
          s = fb / fa
          IF (a == c) THEN
             ! Tentativa de interpolacao linear (metodo da secante)
             p = 2.0_DP * xm * s
             q = 1.0_DP - s
          ELSE
             ! Tentativa de interpolacao quadratica inversa
             q = fa / fc
             r = fb / fc
             p = s * (2.0_DP * xm * q * (q - r) - (b - a) * (r - 1.0_DP))
             q = (q - 1.0_DP) * (r - 1.0_DP) * (s - 1.0_DP)
          END IF
          
          IF (p > 0.0_DP) q = -q
          p = ABS(p)
          
          ! Verifica se a interpolacao e aceitavel
          IF (2.0_DP * p < MIN(3.0_DP * xm * q - ABS(tol1 * q), ABS(e * q))) THEN
             e = d ; d = p / q
          ELSE
             ! Fallback para bissecao
             d = xm ; e = d
          END IF
       ELSE
          ! Passo de bissecao se o decréscimo for lento
          d = xm ; e = d
       END IF
       
       a = b ; fa = fb
       IF (ABS(d) > tol1) THEN
          b = b + d
       ELSE
          b = b + SIGN(tol1, xm)
       END IF
       fb = func(b)
    END DO
    
    ierr = 2 ! Erro: Nao convergiu no numero maximo de iteracoes
    root = b
  END FUNCTION zbrent

END MODULE mod_raizes_brent
