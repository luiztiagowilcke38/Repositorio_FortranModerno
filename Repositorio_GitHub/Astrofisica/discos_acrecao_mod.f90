MODULE discos_acrecao_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  REAL(real64), PARAMETER :: MASSA_CENTRAL = 1.0_real64 ! Em massas solares
  REAL(real64), PARAMETER :: ALPHA_VISC = 0.01_real64
  INTEGER, PARAMETER :: N_RADIAL = 200

CONTAINS

  SUBROUTINE evoluir_disco_viscoso(sigma, r, dt, t_fim)
    REAL(real64), DIMENSION(N_RADIAL), INTENT(INOUT) :: sigma
    REAL(real64), DIMENSION(N_RADIAL), INTENT(IN) :: r
    REAL(real64), INTENT(IN) :: dt, t_fim
    REAL(real64) :: t_atual, dr
    REAL(real64), DIMENSION(N_RADIAL) :: viscosidade, b_vetor, delta_sigma
    REAL(real64), DIMENSION(N_RADIAL, N_RADIAL) :: matriz_dif
    INTEGER :: i
    
    t_atual = 0.0_real64; dr = r(2) - r(1)
    DO i = 1, N_RADIAL; viscosidade(i) = calcular_visc_local(r(i)); END DO
    
    DO WHILE (t_atual < t_fim)
      CALL construir_matriz_difusao(r, viscosidade, dt, dr, matriz_dif)
      CALL calcular_lado_direito(sigma, r, viscosidade, dt, dr, b_vetor)
      
      ! Resolver sistema linear (Matriz Tri-diagonal simplificada)
      CALL resolver_tridiagonal(matriz_dif, b_vetor, sigma)
      t_atual = t_atual + dt
    END DO
  END SUBROUTINE evoluir_disco_viscoso

  FUNCTION calcular_visc_local(r_coord) RESULT(nu)
    REAL(real64), INTENT(IN) :: r_coord; REAL(real64) :: nu, h, cs
    cs = 1.0_real64 / SQRT(r_coord) ! Perfil simples de temperatura
    h = 0.05_real64 * r_coord      ! Aspect ratio constante
    nu = ALPHA_VISC * cs * h
  END FUNCTION calcular_visc_local

  SUBROUTINE construir_matriz_difusao(r, nu, dt, dr, m)
    REAL(real64), DIMENSION(:), INTENT(IN) :: r, nu; REAL(real64), INTENT(IN) :: dt, dr
    REAL(real64), DIMENSION(:,:), INTENT(OUT) :: m; INTEGER :: i, n; REAL(real64) :: diff_c
    n = SIZE(r); m = 0.0
    DO i = 2, n - 1
      diff_c = 3.0_real64 * dt / (r(i) * dr**2)
      m(i, i-1) = -diff_c * SQRT(r(i)) * nu(i-1) * SQRT(r(i-1))
      m(i, i) = 1.0_real64 + 2.0_real64 * diff_c * nu(i) * r(i)
      m(i, i+1) = -diff_c * SQRT(r(i)) * nu(i+1) * SQRT(r(i+1))
    END DO
    m(1,1) = 1.0; m(n,n) = 1.0 ! Condicoes de contorno simples
  END SUBROUTINE construir_matriz_difusao

  SUBROUTINE calcular_lado_direito(sigma, r, nu, dt, dr, b)
    REAL(real64), DIMENSION(:), INTENT(IN) :: sigma, r, nu; REAL(real64), INTENT(IN) :: dt, dr
    REAL(real64), DIMENSION(:), INTENT(OUT) :: b; b = sigma
  END SUBROUTINE calcular_lado_direito

  SUBROUTINE resolver_tridiagonal(a, b, x)
    REAL(real64), DIMENSION(:,:), INTENT(IN) :: a; REAL(real64), DIMENSION(:), INTENT(IN) :: b
    REAL(real64), DIMENSION(:), INTENT(OUT) :: x; INTEGER :: i, n
    REAL(real64), DIMENSION(SIZE(b)) :: cp, dp; n = SIZE(b)
    cp(1) = a(1,2)/a(1,1); dp(1) = b(1)/a(1,1)
    DO i = 2, n-1
      cp(i) = a(i,i+1) / (a(i,i) - a(i,i-1)*cp(i-1))
      dp(i) = (b(i) - a(i,i-1)*dp(i-1)) / (a(i,i) - a(i,i-1)*cp(i-1))
    END DO
    dp(n) = (b(n) - a(n,n-1)*dp(n-1)) / (a(n,n) - a(n,n-1)*cp(n-1))
    x(n) = dp(n); DO i = n-1, 1, -1; x(i) = dp(i) - cp(i)*x(i+1); END DO
  END SUBROUTINE resolver_tridiagonal

END MODULE discos_acrecao_mod
