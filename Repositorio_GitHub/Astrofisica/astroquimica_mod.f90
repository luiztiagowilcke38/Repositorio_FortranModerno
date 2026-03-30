MODULE astroquimica_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N_ESPECIES = 5
  REAL(real64), PARAMETER :: TEMP_ISM = 10.0_real64 ! Kelvin
  REAL(real64), PARAMETER :: DENS_TOTAL = 1.0E4_real64 ! cm^-3

CONTAINS

  SUBROUTINE evoluir_quimica_stiff(abundancias, dt_total)
    REAL(real64), DIMENSION(N_ESPECIES), INTENT(INOUT) :: abundancias
    REAL(real64), INTENT(IN) :: dt_total
    REAL(real64) :: t_atual, dt_sub, erro
    REAL(real64), DIMENSION(N_ESPECIES) :: y_velho, y_novo, taxas
    INTEGER :: iter
    
    t_atual = 0.0_real64; dt_sub = 1.0E-3_real64 * dt_total
    DO WHILE (t_atual < dt_total)
      y_velho = abundancias
      CALL resolver_passo_implicito(y_velho, y_novo, dt_sub)
      abundancias = y_novo; t_atual = t_atual + dt_sub
      ! Ajuste adaptativo simples do passo de tempo
      dt_sub = MIN(dt_sub * 1.2_real64, dt_total - t_atual)
    END DO
  END SUBROUTINE evoluir_quimica_stiff

  SUBROUTINE resolver_passo_implicito(y_v, y_n, dt)
    REAL(real64), DIMENSION(N_ESPECIES), INTENT(IN) :: y_v
    REAL(real64), DIMENSION(N_ESPECIES), INTENT(OUT) :: y_n
    REAL(real64), INTENT(IN) :: dt
    REAL(real64), DIMENSION(N_ESPECIES, N_ESPECIES) :: jacobiano, identidade
    REAL(real64), DIMENSION(N_ESPECIES) :: f_v, delta_y
    INTEGER :: i
    
    identidade = 0.0; DO i = 1, N_ESPECIES; identidade(i,i) = 1.0; END DO
    CALL calcular_taxas_producao(y_v, f_v)
    CALL calcular_jacobiano_quimico(y_v, jacobiano)
    
    ! Sistema: (I - dt*J) * delta_y = dt * f_v
    jacobiano = identidade - dt * jacobiano
    CALL resolver_linear_simples(jacobiano, dt*f_v, delta_y)
    y_n = y_v + delta_y
  END SUBROUTINE resolver_passo_implicito

  SUBROUTINE calcular_taxas_producao(y, f)
    REAL(real64), DIMENSION(N_ESPECIES), INTENT(IN) :: y
    REAL(real64), DIMENSION(N_ESPECIES), INTENT(OUT) :: f
    REAL(real64) :: k1, k2
    ! Exemplo simplificado: H + H -> H2 e H2 + CRP -> H + H
    k1 = 1.0E-17_real64; k2 = 1.0E-15_real64
    f(1) = -2.0*k1*y(1)**2 + 2.0*k2*y(2) ! H
    f(2) = k1*y(1)**2 - k2*y(2)       ! H2
    f(3:5) = 0.0 ! Outras especies
  END SUBROUTINE calcular_taxas_producao

  SUBROUTINE calcular_jacobiano_quimico(y, j)
    REAL(real64), DIMENSION(N_ESPECIES), INTENT(IN) :: y
    REAL(real64), DIMENSION(N_ESPECIES, N_ESPECIES), INTENT(OUT) :: j
    REAL(real64) :: k1, k2
    k1 = 1.0E-17_real64; k2 = 1.0E-15_real64; j = 0.0
    j(1,1) = -4.0*k1*y(1); j(1,2) = 2.0*k2
    j(2,1) = 2.0*k1*y(1); j(2,2) = -k2
  END SUBROUTINE calcular_jacobiano_quimico

  SUBROUTINE resolver_linear_simples(a, b, x)
    REAL(real64), DIMENSION(:,:), INTENT(IN) :: a
    REAL(real64), DIMENSION(:), INTENT(IN) :: b
    REAL(real64), DIMENSION(:), INTENT(OUT) :: x
    ! Eliminacao de Gauss basica para o sistema J*x = b
    INTEGER :: i, j, k, n; REAL(real64), DIMENSION(SIZE(b), SIZE(b)) :: a_tmp
    REAL(real64), DIMENSION(SIZE(b)) :: b_tmp; REAL(real64) :: fator
    n = SIZE(b); a_tmp = a; b_tmp = b
    DO k = 1, n-1
      DO i = k+1, n
        fator = a_tmp(i,k) / a_tmp(k,k)
        DO j = k+1, n; a_tmp(i,j) = a_tmp(i,j) - fator * a_tmp(k,j); END DO
        b_tmp(i) = b_tmp(i) - fator * b_tmp(k)
      END DO
    END DO
    x(n) = b_tmp(n) / a_tmp(n,n)
    DO i = n-1, 1, -1
      x(i) = (b_tmp(i) - SUM(a_tmp(i, i+1:n) * x(i+1:n))) / a_tmp(i,i)
    END DO
  END SUBROUTINE resolver_linear_simples

END MODULE astroquimica_mod
