MODULE mod_transporte
  USE iso_fortran_env
  USE mod_parametros_si
  USE mod_malha
  USE mod_poisson
  USE mod_numerico
  IMPLICIT NONE

CONTAINS

  ! Resolve Equacoes de Continuidade usando Scharfetter-Gummel
  SUBROUTINE resolver_continuidade_eletrons(n, geracao, temperatura)
    REAL(real64), INTENT(INOUT) :: n(:)
    REAL(real64), INTENT(IN) :: geracao(:), temperatura
    REAL(real64) :: a(n_pontos), b(n_pontos), c(n_pontos), f(n_pontos)
    REAL(real64) :: vt, d_phi, mu_n
    INTEGER :: i
    
    vt = (k_boltzmann * temperatura) / q_eletron
    mu_n = mobilidade_eletrons(1.0d17, temperatura) ! Dopagem media p/ simplificar
    
    DO i = 2, n_pontos - 1
       d_phi = (potencial(i+1) - potencial(i)) / vt
       a(i) = (mu_n * vt / dx(i-1)) * bernoulli(d_phi)
       c(i) = (mu_n * vt / dx(i)) * bernoulli(-d_phi)
       b(i) = -(a(i) + c(i))
       f(i) = -geracao(i)
    END DO
    
    ! Contatos
    b(1) = 1.0d0; f(1) = n(1)
    b(n_pontos) = 1.0d0; f(n_pontos) = n(n_pontos)
    
    CALL solver_thomas(a, b, c, f, n)
  END SUBROUTINE resolver_continuidade_eletrons

END MODULE mod_transporte
