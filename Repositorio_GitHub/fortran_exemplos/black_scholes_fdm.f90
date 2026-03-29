!==================================================================================================
! Programa: black_scholes_fdm.f90
! Resolve a EDP de Black-Scholes usando o Metodo de Diferencas Finitas (Crank-Nicolson)
! Equacao: dV/dt + r*S*dV/dS + 0.5*sigma^2*S^2*d2V/dS2 - r*V = 0
!==================================================================================================
PROGRAM black_scholes_fdm
  USE iso_fortran_env
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = real64
  REAL(DP), PARAMETER :: S_MAX = 200.0_DP, T_MAX = 1.0_DP
  REAL(DP), PARAMETER :: K = 100.0_DP, R = 0.05_DP, SIGMA = 0.2_DP
  INTEGER, PARAMETER :: NS = 200, NT = 100 ! Grid points
  
  REAL(DP) :: ds, dt, s(0:NS), v(0:NS), v_new(0:NS)
  REAL(DP) :: a(NS-1), b(NS-1), c(NS-1), rhs(NS-1)
  INTEGER :: i, j

  PRINT *, "===================================================="
  PRINT *, " SOLVER BLACK-SCHOLES: METODO CRANK-NICOLSON "
  PRINT *, "===================================================="

  ds = S_MAX / REAL(NS, DP)
  dt = T_MAX / REAL(NT, DP)
  
  ! 1. Inicializacao do Grid e Condicao Final (Payoff no vencimento)
  DO i = 0, NS
    s(i) = REAL(i, DP) * ds
    v(i) = MAX(s(i) - K, 0.0_DP) ! Call Option payoff
  END DO

  ! 2. Loop de Tempo (Retrogrado de T para 0)
  DO j = 1, NT
    ! Montagem do sistema tridiagonal (A*V_new = B*V)
    ! Coeficientes para Crank-Nicolson
    DO i = 1, NS-1
      a(i) = 0.25_DP * dt * (SIGMA**2 * (REAL(i, DP)**2) - R * REAL(i, DP))
      b(i) = -0.5_DP * dt * (SIGMA**2 * (REAL(i, DP)**2) + R)
      c(i) = 0.25_DP * dt * (SIGMA**2 * (REAL(i, DP)**2) + R * REAL(i, DP))
      
      ! Lado direito (Explicit part)
      rhs(i) = v(i) + 0.5_DP * dt * ( &
               0.5_DP * SIGMA**2 * s(i)**2 * (v(i+1) - 2.0_DP*v(i) + v(i-1))/ds**2 + &
               R * s(i) * (v(i+1) - v(i-1))/(2.0_DP*ds) - R*v(i) )
    END DO

    ! Matriz implicita L = (I - 0.5*dt*L_op)
    ! Usamos Thomas Algorithm para resolver: (1-b)*v_new - a*v_new_prev ...
    ! Redefinindo a, b, c para o solver tridiagonal
    rhs(1) = rhs(1) + a(1) * 0.0_DP ! Condicao de contorno S=0
    ! rhs(NS-1) = rhs(NS-1) + c(NS-1) * (S_MAX - K*exp(-r*dt*j)) ! Condicao S_MAX
    
    ! Resolve o sistema tridiagonal
    CALL thomas_solver(NS-1, -a, 1.0_DP-b, -c, rhs, v_new(1:NS-1))
    
    ! Condicoes de Contorno
    v_new(0) = 0.0_DP
    v_new(NS) = S_MAX - K * EXP(-R * REAL(j, DP) * dt)
    
    v = v_new
  END DO

  PRINT *, "Preco da Opcao em S = K (100):", v(NS/2)
  PRINT *, "Solucao Analitica Estimada:", analytical_bs(100.0_DP, K, R, SIGMA, T_MAX)

CONTAINS

  SUBROUTINE thomas_solver(n, a, b, c, r, x)
    INTEGER, INTENT(IN) :: n
    REAL(DP), INTENT(IN) :: a(n), b(n), c(n), r(n)
    REAL(DP), INTENT(OUT) :: x(n)
    REAL(DP) :: cp(n), rp(n), m
    INTEGER :: i
    
    cp(1) = c(1) / b(1)
    rp(1) = r(1) / b(1)
    DO i = 2, n
      m = b(i) - a(i) * cp(i-1)
      cp(i) = c(i) / m
      rp(i) = (r(i) - a(i) * rp(i-1)) / m
    END DO
    x(n) = rp(n)
    DO i = n-1, 1, -1
      x(i) = rp(i) - cp(i) * x(i+1)
    END DO
  END SUBROUTINE

  REAL(DP) FUNCTION analytical_bs(S0, K_val, r_val, sigma_val, T_val)
    REAL(DP), INTENT(IN) :: S0, K_val, r_val, sigma_val, T_val
    REAL(DP) :: d1, d2
    d1 = (LOG(S0/K_val) + (r_val + 0.5_DP*sigma_val**2)*T_val) / (sigma_val*SQRT(T_val))
    d2 = d1 - sigma_val*SQRT(T_val)
    analytical_bs = S0 * normal_cdf(d1) - K_val * EXP(-r_val*T_val) * normal_cdf(d2)
  END FUNCTION

  REAL(DP) FUNCTION normal_cdf(x)
    REAL(DP), INTENT(IN) :: x
    normal_cdf = 0.5_DP * (1.0_DP + erf(x / SQRT(2.0_DP)))
  END FUNCTION

END PROGRAM black_scholes_fdm
