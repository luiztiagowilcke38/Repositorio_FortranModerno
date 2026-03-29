!==================================================================================================
! Modulo: libnegf_types.f90
! Biblioteca Cientifica para Transporte Quantico - Definicao de Estruturas
!==================================================================================================
MODULE libnegf_types
  USE iso_fortran_env
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: DP = real64
  REAL(DP), PARAMETER :: PI = 3.14159265358979323846_DP
  COMPLEX(DP), PARAMETER :: CI = (0.0_DP, 1.0_DP)

  TYPE :: Lead
    INTEGER :: dim
    COMPLEX(DP), ALLOCATABLE :: h00(:,:), h01(:,:)
    COMPLEX(DP), ALLOCATABLE :: sigma(:,:), gamma(:,:)
  CONTAINS
    PROCEDURE :: init => init_lead
  END TYPE Lead

  TYPE :: QuantumSystem
    INTEGER :: n_sites, n_leads
    COMPLEX(DP), ALLOCATABLE :: ham(:,:), ovlp(:,:)
    TYPE(Lead), ALLOCATABLE :: leads(:)
    REAL(DP), ALLOCATABLE :: energies(:)
  CONTAINS
    PROCEDURE :: init => init_system
  END TYPE QuantumSystem

CONTAINS

  SUBROUTINE init_lead(self, d)
    CLASS(Lead), INTENT(OUT) :: self
    INTEGER, INTENT(IN) :: d
    self%dim = d
    ALLOCATE(self%h00(d,d), self%h01(d,d), self%sigma(d,d), self%gamma(d,d))
    self%h00 = 0.0_DP; self%h01 = 0.0_DP
  END SUBROUTINE

  SUBROUTINE init_system(self, s, l)
    CLASS(QuantumSystem), INTENT(OUT) :: self
    INTEGER, INTENT(IN) :: s, l
    INTEGER :: i
    self%n_sites = s
    self%n_leads = l
    ALLOCATE(self%ham(s,s), self%ovlp(s,s), self%leads(l))
    self%ham = 0.0_DP
    self%ovlp = 0.0_DP; DO i=1,self%n_sites; self%ovlp(i,i)=1.0_DP; END DO
  END SUBROUTINE

END MODULE libnegf_types
