!==================================================================================================
! Modulo: libnegf_observables.f90
! Calculo de Grandezas Fisicas (Transmissao, Corrente, DOS)
!==================================================================================================
MODULE libnegf_observables
  USE libnegf_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: transmission, density_of_states, landauer_current, fermi_dirac
  PUBLIC :: transmission_fast, spectral_function, current_spectrum

CONTAINS

  ! Coeficiente de Transmissao: T(E) = Tr[ GammaL * Gr * GammaR * Ga ]
  FUNCTION transmission(gr, gamma_l, gamma_r) RESULT(t)
    COMPLEX(DP), INTENT(IN) :: gr(:,:), gamma_l(:,:), gamma_r(:,:)
    REAL(DP) :: t
    COMPLEX(DP), ALLOCATABLE :: tmp(:,:)
    INTEGER :: n, i
    
    n = SIZE(gr, 1)
    ALLOCATE(tmp(n,n))
    
    ! Calculo otimizado da traca: Tr(A*B*C*D)
    ! T = Tr[ GammaL * Gr * GammaR * Gr^dagger ]
    tmp = MATMUL(gamma_l, MATMUL(gr, MATMUL(gamma_r, TRANSPOSE(CONJG(gr)))))
    
    t = 0.0_DP
    DO i = 1, n
      t = t + REAL(tmp(i,i), DP)
    END DO
    DEALLOCATE(tmp)
  END FUNCTION transmission

  ! Densidade de Estados Local (LDOS): DOS(E) = -1/pi * Im[ Tr(Gr) ]
  FUNCTION density_of_states(gr) RESULT(dos)
    COMPLEX(DP), INTENT(IN) :: gr(:,:)
    REAL(DP) :: dos
    INTEGER :: i, n
    n = SIZE(gr, 1)
    dos = 0.0_DP
    DO i = 1, n
      dos = dos - (1.0_DP/PI) * AIMAG(gr(i,i))
    END DO
  END FUNCTION density_of_states

  ! Distribuicao de Fermi-Dirac
  FUNCTION fermi_dirac(energy, mu, temperature) RESULT(f)
    REAL(DP), INTENT(IN) :: energy, mu, temperature
    REAL(DP) :: f
    REAL(DP), PARAMETER :: kB = 0.08617_DP ! Constante de Boltzmann em meV/K
    
    IF (temperature < 1.0E-6_DP) THEN
      IF (energy < mu) THEN
        f = 1.0_DP
      ELSE
        f = 0.0_DP
      END IF
    ELSE
      f = 1.0_DP / (EXP((energy - mu) / (kB * temperature)) + 1.0_DP)
    END IF
  END FUNCTION fermi_dirac

  ! Corrente de Landauer (em unidades de 2e/h)
  ! I = \int T(E) [fL(E) - fR(E)] dE
  FUNCTION landauer_current(energies, trans, mu_l, mu_r, temp) RESULT(curr)
    REAL(DP), INTENT(IN) :: energies(:), trans(:), mu_l, mu_r, temp
    REAL(DP) :: curr
    INTEGER :: i, ne
    REAL(DP) :: de, f_l, f_r
    
    ne = SIZE(energies)
    curr = 0.0_DP
    DO i = 1, ne - 1
      de = energies(i+1) - energies(i)
      f_l = fermi_dirac(energies(i), mu_l, temp)
      f_r = fermi_dirac(energies(i), mu_r, temp)
      curr = curr + trans(i) * (f_l - f_r) * de
    END DO
  END FUNCTION landauer_current

  !---------------------------------------------------------------------
  ! EXPANSÃO: Observáveis Avançados e Integração de Alta Ordem
  !---------------------------------------------------------------------

  ! TRANSMISSION_FAST: Versao otimizada que evita alocacoes densas desnecessárias
  ! e utiliza a propriedade Tr(AB) = sum(A.*B^T) para acelerar o calculo.
  FUNCTION transmission_fast(gr, gamma_l, gamma_r) RESULT(t)
    COMPLEX(DP), INTENT(IN) :: gr(:,:), gamma_l(:,:), gamma_r(:,:)
    REAL(DP) :: t
    COMPLEX(DP), ALLOCATABLE :: g_l(:,:), ga_r(:,:)
    INTEGER :: n

    n = SIZE(gr, 1)
    ALLOCATE(g_l(n,n), ga_r(n,n))
    
    ! G_L = Gr * Gamma_R
    g_l = MATMUL(gr, gamma_r)
    ! GA_R = Gr^H * Gamma_L
    ga_r = MATMUL(TRANSPOSE(CONJG(gr)), gamma_l)

    ! Tr(A*B) calculada como produto interno elemento a elemento
    t = REAL(SUM(g_l * TRANSPOSE(ga_r)), DP)
    
    DEALLOCATE(g_l, ga_r)
  END FUNCTION transmission_fast

  ! SPECTRAL_FUNCTION: A(E) = i(Gr - Gr^H)
  SUBROUTINE spectral_function(gr, spectral_mat)
    COMPLEX(DP), INTENT(IN)  :: gr(:,:)
    COMPLEX(DP), INTENT(OUT) :: spectral_mat(:,:)
    spectral_mat = (0.0_DP, 1.0_DP) * (gr - TRANSPOSE(CONJG(gr)))
  END SUBROUTINE spectral_function

  ! CURRENT_SPECTRUM: Retorna o integrando da corrente dI/dE
  FUNCTION current_spectrum(energy, transmission_val, mu_l, mu_r, temp) RESULT(die)
    REAL(DP), INTENT(IN) :: energy, transmission_val, mu_l, mu_r, temp
    REAL(DP) :: die
    die = transmission_val * (fermi_dirac(energy, mu_l, temp) - fermi_dirac(energy, mu_r, temp))
  END FUNCTION current_spectrum

  ! SIMPSON_INTEGRATION: Integração numérica de 3a ordem para corrente total
  FUNCTION simpson_current(energies, trans, mu_l, mu_r, temp) RESULT(curr)
    REAL(DP), INTENT(IN) :: energies(:), trans(:), mu_l, mu_r, temp
    REAL(DP) :: curr
    INTEGER :: i, ne
    REAL(DP) :: h, f(3)

    ne = SIZE(energies)
    curr = 0.0_DP
    h = energies(2) - energies(1)

    DO i = 1, ne - 2, 2
       f(1) = current_spectrum(energies(i),   trans(i),   mu_l, mu_r, temp)
       f(2) = current_spectrum(energies(i+1), trans(i+1), mu_l, mu_r, temp)
       f(3) = current_spectrum(energies(i+2), trans(i+2), mu_l, mu_r, temp)
       curr = curr + (h/3.0_DP) * (f(1) + 4.0_DP*f(2) + f(3))
    END DO
  END FUNCTION simpson_current

END MODULE libnegf_observables
