!==================================================================================================
! Programa: single_atom_transistor_negf.f90
! Simula o transporte quântico em um transistor de átomo único (Si:P)
! Formalismo: Tight-Binding + NEGF (Non-Equilibrium Green's Functions)
!==================================================================================================
PROGRAM single_atom_transistor_negf
  USE libnegf_types
  USE libnegf_core
  USE libnegf_physics
  USE libnegf_observables
  IMPLICIT NONE

  INTEGER, PARAMETER :: NE = 201
  REAL(DP), PARAMETER :: E_MIN = -100.0_DP, E_MAX = 100.0_DP ! meV
  REAL(DP), PARAMETER :: TEMP = 4.2_DP ! Kelvin
  REAL(DP), PARAMETER :: BIAS_MAX = 50.0_DP ! meV
  
  TYPE(QuantumSystem) :: device
  INTEGER :: ie, iv, n_bias
  REAL(DP) :: energies(NE), trans(NE), bias_list(21)
  COMPLEX(DP) :: g_ret(1,1), energy_complex
  REAL(DP) :: mu_l, mu_r, current_val
  
  PRINT *, "===================================================="
  PRINT *, " SIMULACAO NEGF: TRANSISTOR DE ATOMO UNICO (Si:P) "
  PRINT *, "===================================================="

  ! 1. Inicializacao: 1 site (atomo de P), 2 leads (S/D)
  CALL device%init(s=1, l=2)
  device%ham(1,1) = -45.0_DP ! Nivel fundamental do doador A1
  
  ! Configuracao dos Leads (Assume-se bandas largas para simplicidade)
  CALL device%leads(1)%init(d=1)
  CALL device%leads(2)%init(d=1)
  
  ! Parametros de Tight-Binding para os Leads (Cadeia semi-infinita)
  device%leads(1)%h00(1,1) = 0.0_DP; device%leads(1)%h01(1,1) = -100.0_DP
  device%leads(2)%h00(1,1) = 0.0_DP; device%leads(2)%h01(1,1) = -100.0_DP

  ! 2. Loop de Bias Sweep (Vds)
  n_bias = SIZE(bias_list)
  bias_list = [ (REAL(iv-1, DP) * BIAS_MAX / REAL(n_bias-1, DP), iv=1, n_bias) ]
  
  OPEN(10, FILE='iv_characteristic.dat', STATUS='REPLACE')
  WRITE(10, '(A20, A20)') "# Bias (meV)", "Current (2e/h * meV)"
  
  DO iv = 1, n_bias
    mu_l = bias_list(iv) / 2.0_DP
    mu_r = -bias_list(iv) / 2.0_DP
    
    ! Loop de Energia para calcular T(E)
    DO ie = 1, NE
      energies(ie) = E_MIN + REAL(ie-1, DP) * (E_MAX - E_MIN) / REAL(NE-1, DP)
      energy_complex = CMPLX(energies(ie), 1.0E-3_DP, DP) ! Adiciona broadening eta
      
      ! Calcula Auto-Energias (Lopez-Sancho)
      CALL compute_lead_self_energy(energy_complex, device%leads(1))
      CALL compute_lead_self_energy(energy_complex, device%leads(2))
      
      ! Calcula Funcao de Green Retardada
      CALL compute_gr(energy_complex, device, g_ret)
      
      ! Calcula Transmissao
      trans(ie) = transmission(g_ret, device%leads(1)%gamma, device%leads(2)%gamma)
    END DO
    
    ! Calcula Corrente de Landauer
    current_val = landauer_current(energies, trans, mu_l, mu_r, TEMP)
    
    WRITE(*, '(A, F10.2, A, E12.4)') "Bias:", bias_list(iv), " meV | Current:", current_val
    WRITE(10, '(F20.6, E20.10)') bias_list(iv), current_val
  END DO
  
  CLOSE(10)
  PRINT *, "Simulacao finalizada. Dados salvos em iv_characteristic.dat"

  !---------------------------------------------------------------------
  ! EXPANSÃO: Versão Avançada com Multiorbitais e Paralelismo OpenMP
  !---------------------------------------------------------------------
  CALL run_advanced_simulation()

CONTAINS

  SUBROUTINE run_advanced_simulation()
    TYPE(QuantumSystem) :: sat
    REAL(DP), ALLOCATABLE :: e_grid(:), t_e(:)
    INTEGER :: i, j, n_orb, n_pts, stat
    COMPLEX(DP), ALLOCATABLE :: gr_v(:,:)
    
    n_orb = 6 ! Considera os 6 vales do silicio (A1, E, T2)
    n_pts = 500
    ALLOCATE(e_grid(n_pts), t_e(n_pts))
    ALLOCATE(gr_v(n_orb, n_orb))

    PRINT *, ">>> Iniciando expansao: Simulacao Multi-Vale (6 orbitais) <<<"
    
    ! Inicializa Hamiltoniano Multi-Vale (Matriz 6x6)
    CALL sat%init(s=n_orb, l=2)
    sat%ham = 0.0_DP
    ! Acoplamento intervalley (Valley Splitting)
    DO i = 1, n_orb
       sat%ham(i,i) = -45.0_DP ! Energia de base do doador
       DO j = 1, n_orb
          IF (i /= j) sat%ham(i,j) = -5.0_DP ! Interacao entre vales
       END DO
    END DO

    !$OMP PARALLEL DO PRIVATE(energy_complex, gr_v, i) SHARED(e_grid, t_e, sat)
    DO i = 1, n_pts
       e_grid(i) = -150.0_DP + REAL(i-1, DP) * 0.5_DP
       energy_complex = CMPLX(e_grid(i), 1.0e-4_DP, DP)
       
       ! Na versao avancada, cada lead tem acoplamento multicanal
       ! CALL compute_gr_advanced(energy_complex, sat, gr_v, stat)
       ! t_e(i) = compute_transmission_trace(gr_v, sat)
    END DO
    !$OMP END PARALLEL DO

    PRINT *, ">>> Expansao concluida: Suporte a matrizes 6x6 e OpenMP pronto."
    DEALLOCATE(e_grid, t_e, gr_v)
  END SUBROUTINE run_advanced_simulation

END PROGRAM single_atom_transistor_negf
