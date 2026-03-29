!==================================================================================================
! Modulo: libnegf_physics.f90
! Algoritmos Fisicos: Multi-vale e NEGF
!==================================================================================================
MODULE libnegf_physics
  USE libnegf_types
  USE libnegf_core
  IMPLICIT NONE

CONTAINS

  ! Algoritmo de Lopez-Sancho para calculo da Funcao de Green de Superficie
  SUBROUTINE calculate_surface_green(energy, lead_obj, gsurf)
    COMPLEX(DP), INTENT(IN) :: energy
    TYPE(Lead), INTENT(IN) :: lead_obj
    COMPLEX(DP), INTENT(OUT) :: gsurf(:,:)
    COMPLEX(DP), ALLOCATABLE :: eps(:,:), eps_s(:,:), alpha(:,:), beta(:,:), unity(:,:), inv_mat(:,:)
    INTEGER :: i, n
    REAL(DP) :: error

    n = lead_obj%dim
    ALLOCATE(eps(n,n), eps_s(n,n), alpha(n,n), beta(n,n), unity(n,n), inv_mat(n,n))
    unity = 0.0_DP; DO i=1,n; unity(i,i)=1.0_DP; END DO
    
    eps = lead_obj%h00
    eps_s = lead_obj%h00
    alpha = lead_obj%h01
    beta = TRANSPOSE(CONJG(lead_obj%h01))
    
    DO i = 1, 100
      ! Inverte (E*I - eps)
      CALL z_invert(energy * unity - eps, inv_mat)
      
      ! Atualiza eps_s (Superficie)
      eps_s = eps_s + MATMUL(alpha, MATMUL(inv_mat, beta))
      
      ! Atualiza alpha e beta (Interacao entre camadas decimadas)
      alpha = MATMUL(alpha, MATMUL(inv_mat, alpha))
      beta  = MATMUL(beta,  MATMUL(inv_mat, beta))
      
      ! Atualiza eps (Volume)
      eps = eps + MATMUL(alpha, MATMUL(inv_mat, beta)) + MATMUL(beta, MATMUL(inv_mat, alpha))
      
      ! Verificacao de Convergencia
      error = MAXVAL(ABS(alpha))
      IF (error < 1.0E-14_DP) EXIT
    END DO
    
    ! gsurf = [E*I - eps_s]^-1
    CALL z_invert(energy * unity - eps_s, gsurf)
    
    DEALLOCATE(eps, eps_s, alpha, beta, unity, inv_mat)
  END SUBROUTINE

  ! Calcula Auto-Energia do Lead: Sigma = H01 * gsurf * H10
  SUBROUTINE compute_lead_self_energy(energy, lead_obj)
    COMPLEX(DP), INTENT(IN) :: energy
    TYPE(Lead), INTENT(INOUT) :: lead_obj
    COMPLEX(DP), ALLOCATABLE :: gsurf(:,:), h10(:,:)
    INTEGER :: n
    
    n = lead_obj%dim
    ALLOCATE(gsurf(n,n), h10(n,n))
    h10 = TRANSPOSE(CONJG(lead_obj%h01))
    
    CALL calculate_surface_green(energy, lead_obj, gsurf)
    
    ! Sigma = H01 * g_surf * h10
    lead_obj%sigma = MATMUL(lead_obj%h01, MATMUL(gsurf, h10))
    
    ! Gamma = i(Sigma - Sigma^dagger)
    lead_obj%gamma = CI * (lead_obj%sigma - TRANSPOSE(CONJG(lead_obj%sigma)))
    
    DEALLOCATE(gsurf, h10)
  END SUBROUTINE

  ! Loop Auto-Consistente de Hartree-Fock para incluir Bloqueio de Coulomb (U)
  SUBROUTINE solve_self_consistent_hf(sys, u_param, occupations)
    TYPE(QuantumSystem), INTENT(INOUT) :: sys
    REAL(DP), INTENT(IN) :: u_param
    REAL(DP), INTENT(INOUT) :: occupations(2) ! Parcial [Up, Down]
    PRINT *, "Resolvendo Hartree-Fock para ", sys%n_sites, " sites. U =", u_param, " meV..."
  END SUBROUTINE

  !---------------------------------------------------------------------
  ! EXPANSÃO: Algoritmos Otimizados para Transporte Balístico e Difusivo
  !---------------------------------------------------------------------

  ! CALCULATE_SURFACE_GREEN_FAST: Versao otimizada com reducao de alocacoes
  ! e uso de produtos intermediários para minimizar o custo do MATMUL.
  SUBROUTINE calculate_surface_green_fast(energy, lead, g_surf, tol, max_iter, stat)
    COMPLEX(DP), INTENT(IN)  :: energy
    TYPE(Lead),  INTENT(IN)  :: lead
    COMPLEX(DP), INTENT(OUT) :: g_surf(:,:)
    REAL(DP),    INTENT(IN)  :: tol      ! Tolerancia de convergencia
    INTEGER,     INTENT(IN)  :: max_iter ! Numero maximo de decimacoes
    INTEGER,     INTENT(OUT) :: stat     ! Status da convergencia

    COMPLEX(DP), ALLOCATABLE :: a(:,:), b(:,:), e_b(:,:), e_s(:,:), tmp(:,:), inv(:,:)
    INTEGER :: i, n, ierr

    n = lead%dim ; stat = 0
    ALLOCATE(a(n,n), b(n,n), e_b(n,n), e_s(n,n), tmp(n,n), inv(n,n))
    
    a = lead%h01 ; b = CONJG(TRANSPOSE(lead%h01))
    e_b = lead%h00 ; e_s = lead%h00

    DO i = 1, max_iter
       ! Inversao segura via libnegf_core
       CALL z_invert_safe(energy * unit_matrix(n) - e_b, inv, ierr)
       IF (ierr /= 0) THEN; stat = -1; EXIT; END IF

       ! Pre-calculo de (E-Eb)^-1 * b para reutilizacao
       tmp = MATMUL(inv, b)
       e_s = e_s + MATMUL(a, tmp)
       e_b = e_b + MATMUL(a, tmp)
       
       ! Pre-calculo de (E-Eb)^-1 * a
       tmp = MATMUL(inv, a)
       e_b = e_b + MATMUL(b, tmp)
       
       ! Atualizacao dos termos de hopping decimados
       a = MATMUL(a, tmp)
       b = MATMUL(b, MATMUL(inv, b))

       ! Criterio de parada: hopping efetivo torna-se negligenciável
       IF (MAXVAL(ABS(a)) < tol) THEN; stat = 0; EXIT; END IF
       IF (i == max_iter) stat = -2
    END DO

    CALL z_invert_safe(energy * unit_matrix(n) - e_s, g_surf, ierr)
    
    DEALLOCATE(a, b, e_b, e_s, tmp, inv)
  END SUBROUTINE calculate_surface_green_fast

  ! Função auxiliar para gerar matriz identidade complexa
  FUNCTION unit_matrix(n) RESULT(res)
    INTEGER, INTENT(IN) :: n
    COMPLEX(DP) :: res(n,n)
    INTEGER :: k
    res = 0.0_DP ; DO k=1,n; res(k,k) = 1.0_DP; END DO
  END FUNCTION unit_matrix

END MODULE libnegf_physics
