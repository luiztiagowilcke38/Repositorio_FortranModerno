!==================================================================================================
! Modulo: libnegf_core.f90
! Operacoes de Nucleo: Algebra Linear Complexa e Inversao de GFs
!==================================================================================================
MODULE libnegf_core
  USE libnegf_types
  IMPLICIT NONE

CONTAINS

  ! Inversao de Matriz Complexa via LAPACK (ZGETRF + ZGETRI)
  SUBROUTINE z_invert(A, Ainv)
    COMPLEX(DP), INTENT(IN) :: A(:,:)
    COMPLEX(DP), INTENT(OUT) :: Ainv(:,:)
    COMPLEX(DP), ALLOCATABLE :: work(:)
    INTEGER, ALLOCATABLE :: ipiv(:)
    INTEGER :: n, lwork, info

    n = SIZE(A, 1)
    Ainv = A
    ALLOCATE(ipiv(n))
    
    ! Fatoracao LU
    CALL zgetrf(n, n, Ainv, n, ipiv, info)
    IF (info /= 0) THEN
      PRINT *, "Erro em ZGETRF: Matriz singular ou problema de convergencia. INFO =", info
      STOP
    END IF

    ! Workspace query
    lwork = -1
    ALLOCATE(work(1))
    CALL zgetri(n, Ainv, n, ipiv, work, lwork, info)
    lwork = INT(REAL(work(1)))
    DEALLOCATE(work)
    ALLOCATE(work(lwork))

    ! Inversao real
    CALL zgetri(n, Ainv, n, ipiv, work, lwork, info)
    IF (info /= 0) THEN
      PRINT *, "Erro em ZGETRI: Falha na inversao. INFO =", info
      STOP
    END IF
  END SUBROUTINE

  ! Funcao de Green Retardada: Gr = [E*S - H - SigmaL - SigmaR]^-1
  SUBROUTINE compute_gr(energy, sys, gr)
    COMPLEX(DP), INTENT(IN) :: energy
    TYPE(QuantumSystem), INTENT(IN) :: sys
    COMPLEX(DP), INTENT(OUT) :: gr(:,:)
    COMPLEX(DP), ALLOCATABLE :: tmp(:,:)
    INTEGER :: i, n

    n = sys%n_sites
    ALLOCATE(tmp(n,n))
    tmp = energy * sys%ovlp - sys%ham
    
    ! Subtrai auto-energias (Assume-se que ja foram calculadas)
    DO i = 1, sys%n_leads
      ! Nota: Aqui assumimos que as acoplamentos estao mapeados corretamente
      tmp = tmp - sys%leads(i)%sigma
    END DO

    CALL z_invert(tmp, gr)
  END SUBROUTINE

  !---------------------------------------------------------------------
  ! EXPANSÃO: Subrotinas Avançadas e Robusta de Álgebra Linear
  !---------------------------------------------------------------------

  ! Z_INVERT_SAFE: Versao que retorna status de erro em vez de parar a execucao
  SUBROUTINE z_invert_safe(A, Ainv, stat)
    COMPLEX(DP), INTENT(IN)  :: A(:,:)
    COMPLEX(DP), INTENT(OUT) :: Ainv(:,:)
    INTEGER,     INTENT(OUT) :: stat
    COMPLEX(DP), ALLOCATABLE :: work(:)
    INTEGER,     ALLOCATABLE :: ipiv(:)
    INTEGER :: n, lwork, info

    n = SIZE(A, 1) ; stat = 0
    Ainv = A
    ALLOCATE(ipiv(n))
    
    ! Fatoracao LU via LAPACK
    CALL zgetrf(n, n, Ainv, n, ipiv, info)
    IF (info /= 0) THEN
      stat = -1; RETURN
    END IF

    ! Consulta de tamanho de buffer (LWORK)
    lwork = -1
    ALLOCATE(work(1))
    CALL zgetri(n, Ainv, n, ipiv, work, lwork, info)
    lwork = MAX(1, INT(REAL(work(1))))
    DEALLOCATE(work)
    
    ALLOCATE(work(lwork))
    ! Inversao matricial de alta performance
    CALL zgetri(n, Ainv, n, ipiv, work, lwork, info)
    IF (info /= 0) stat = -2

    DEALLOCATE(ipiv, work)
  END SUBROUTINE z_invert_safe

  ! COMPUTE_GR_ADVANCED: Versao expandida com suporte a auto-energia de muitos corpos
  SUBROUTINE compute_gr_advanced(energy, sys, sigma_int, gr, stat)
    COMPLEX(DP), INTENT(IN)          :: energy
    TYPE(QuantumSystem), INTENT(IN)  :: sys
    COMPLEX(DP), INTENT(IN), OPTIONAL :: sigma_int(:,:)
    COMPLEX(DP), INTENT(OUT)         :: gr(:,:)
    INTEGER,     INTENT(OUT)         :: stat

    COMPLEX(DP), ALLOCATABLE :: mat_inv(:,:)
    INTEGER :: i, n

    n = sys%n_sites
    ALLOCATE(mat_inv(n, n))
    
    ! Equacao de Dyson: mat_inv = E*S - H - Sigma_leads - Sigma_interaction
    mat_inv = energy * sys%ovlp - sys%ham

    IF (ALLOCATED(sys%leads)) THEN
      DO i = 1, SIZE(sys%leads)
        mat_inv = mat_inv - sys%leads(i)%sigma
      END DO
    END IF

    ! Adiciona correlações elétron-elétron ou elétron-fônon se presentes
    IF (PRESENT(sigma_int)) THEN
      mat_inv = mat_inv - sigma_int
    END IF

    ! Chama a inversao segura
    CALL z_invert_safe(mat_inv, gr, stat)
    
    DEALLOCATE(mat_inv)
  END SUBROUTINE compute_gr_advanced

END MODULE libnegf_core
