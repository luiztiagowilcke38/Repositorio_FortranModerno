!==================================================================================================
! SOLUCAO DOS EXERCICIOS PROPOSTOS: CELULA SOLAR DE ATOMO UNICO (Si:P)
!--------------------------------------------------------------------------------------------------
! Este arquivo contem as rotinas modificadas e implementacoes solicitadas nos exercicios do
! capitulo "Modelo de Celula Solar de Atomo Unico de Silicio".
!==================================================================================================

MODULE exercicios_resolvidos_si_p
  USE iso_fortran_env
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = real64
  REAL(DP), PARAMETER :: PI = 3.14159265358979323846_DP

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! EXERCICIO 1: Analise de Strain Uniaxial
  ! Modificacao do tensor para quebrar a simetria de vale de forma nao-isotropica.
  !------------------------------------------------------------------------------------------------
  SUBROUTINE solucao_exercicio_1_strain_uniaxial(shifts)
    REAL(DP), INTENT(OUT) :: shifts(6)
    REAL(DP) :: eps(3,3), tr_eps
    
    ! Definicao de Strain Uniaxial ao longo de [100]
    eps = 0.0_DP
    eps(1,1) = 0.005_DP ! 0.5% de expansao
    
    tr_eps = eps(1,1) + eps(2,2) + eps(3,3)
    
    ! Modelo de Bir-Pikus (Xi_d = 1.1, Xi_u = 8.6 eV)
    ! Vales X (1,2)
    shifts(1:2) = (1.1_DP * tr_eps + 8.6_DP * eps(1,1)) * 1000.0_DP ! meV
    ! Vales Y (3,4)
    shifts(3:4) = (1.1_DP * tr_eps + 8.6_DP * eps(2,2)) * 1000.0_DP 
    ! Vales Z (5,6)
    shifts(5:6) = (1.1_DP * tr_eps + 8.6_DP * eps(3,3)) * 1000.0_DP
    
    PRINT *, "Shift Uniaxial Vale X (meV):", shifts(1)
    PRINT *, "Shift Uniaxial Vale Y,Z (meV):", shifts(3)
  END SUBROUTINE solucao_exercicio_1_strain_uniaxial

  !------------------------------------------------------------------------------------------------
  ! EXERCICIO 2: Acoplamento de Fonons Opticos (Einstein Peak)
  ! Modificacao da auto-energia de fonons para incluir um pico de Einstein em 60 meV.
  !------------------------------------------------------------------------------------------------
  SUBROUTINE solucao_exercicio_2_fonons_opticos(Gl, Sl_opt, ie, E_vals, ne)
    COMPLEX(DP), INTENT(IN)  :: Gl(12,12,ne)
    COMPLEX(DP), INTENT(OUT) :: Sl_opt(12,12)
    INTEGER, INTENT(IN)      :: ie, ne
    REAL(DP), INTENT(IN)     :: E_vals(ne)
    REAL(DP) :: h_om = 60.0_DP ! Energia do fonon optico (meV)
    REAL(DP) :: g2 = 15.0_DP   ! Constante de acoplamento
    INTEGER :: k
    
    Sl_opt = 0.0_DP
    DO k = 1, ne
       ! Filtra apenas transicoes que casam com a energia de Einstein
       IF (ABS(ABS(E_vals(ie) - E_vals(k)) - h_om) < 2.5_DP) THEN
          Sl_opt = Sl_opt + g2 * Gl(:,:,k)
       END IF
    END DO
  END SUBROUTINE solucao_exercicio_2_fonons_opticos

  !------------------------------------------------------------------------------------------------
  ! EXERCICIO 3: Limite Termodinamico e EQE via Integração de Simpson
  ! Implementacao da EQE (External Quantum Efficiency).
  !------------------------------------------------------------------------------------------------
  FUNCTION solucao_exercicio_3_eqe_simpson(fluxo_solar, i_gen, n_pontos) RESULT(eqe)
    REAL(DP), INTENT(IN) :: fluxo_solar(n_pontos), i_gen(n_pontos)
    INTEGER, INTENT(IN)  :: n_pontos
    REAL(DP) :: eqe, h_step
    INTEGER :: i
    
    ! Integracao de Simpson 1/3
    h_step = 1.0_DP ! assumindo passo unitario de energia/wv
    eqe = i_gen(1) / (fluxo_solar(1) + 1.0E-10_DP) + i_gen(n_pontos) / (fluxo_solar(n_pontos) + 1.0E-10_DP)
    
    DO i = 2, n_pontos-1
       IF (MOD(i, 2) == 0) THEN
          eqe = eqe + 4.0_DP * i_gen(i) / (fluxo_solar(i) + 1.0E-10_DP)
       ELSE
          eqe = eqe + 2.0_DP * i_gen(i) / (fluxo_solar(i) + 1.0E-10_DP)
       END IF
    END DO
    eqe = eqe * (h_step / 3.0_DP)
  END FUNCTION solucao_exercicio_3_eqe_simpson

END MODULE exercicios_resolvidos_si_p
