MODULE mod_poisson
  USE iso_fortran_env
  USE mod_parametros_si
  USE mod_malha
  USE mod_numerico
  IMPLICIT NONE
  
  REAL(real64) :: potencial(n_pontos)

CONTAINS

  ! Resolve Poisson: div(eps * grad(phi)) = -q(p - n + Nd - Na)
  SUBROUTINE resolver_poisson(n, p, temperatura)
    REAL(real64), INTENT(IN) :: n(:), p(:), temperatura
    REAL(real64) :: a(n_pontos), b(n_pontos), c(n_pontos), f(n_pontos)
    REAL(real64) :: delta_phi(n_pontos)
    INTEGER :: i
    REAL(real64) :: vt
    
    vt = (k_boltzmann * temperatura) / q_eletron
    
    ! Montagem da matriz tridiagonal (Linearizacao de Gummel)
    DO i = 2, n_pontos - 1
       a(i) = (eps0 * eps_si) / (dx(i-1) * (dx(i-1) + dx(i)) / 2.0d0)
       c(i) = (eps0 * eps_si) / (dx(i) * (dx(i-1) + dx(i)) / 2.0d0)
       b(i) = -(a(i) + c(i)) - (q_eletron / vt) * (n(i) + p(i))
       f(i) = -(eps0 * eps_si * (potencial(i+1)-2*potencial(i)+potencial(i-1))/dx(i)**2) & ! (Simplificado)
              - q_eletron * (p(i) - n(i) + Nd(i) - Na(i))
    END DO
    
    ! Condicao de Contorno Ôhmica (Dirichlet)
    b(1) = 1.0d0; f(1) = potencial(1)
    b(n_pontos) = 1.0d0; f(n_pontos) = potencial(n_pontos)
    
    CALL solver_thomas(a, b, c, f, delta_phi)
    potencial = potencial + delta_phi
  END SUBROUTINE resolver_poisson

END MODULE mod_poisson
