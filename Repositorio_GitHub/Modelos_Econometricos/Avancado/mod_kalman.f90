MODULE mod_kalman
  USE iso_fortran_env
  USE mod_espaco_estados
  IMPLICIT NONE

CONTAINS

  ! Filtro de Kalman Clássico (Passe Direto)
  SUBROUTINE filtrar_kalman(modelo, y, x_filtrado, P_filtrado, log_L)
    TYPE(ModeloEspacoEstados), INTENT(IN) :: modelo
    REAL(real64), INTENT(IN) :: y(:,:) ! Serie temporal (dim_obs x T)
    REAL(real64), ALLOCATABLE, INTENT(OUT) :: x_filtrado(:,:), P_filtrado(:,:,:)
    REAL(real64), INTENT(OUT) :: log_L
    
    INTEGER :: T, i, n, l
    REAL(real64), ALLOCATABLE :: x_pred(:), P_pred(:,:), v(:), S(:,:), K(:,:)
    
    n = modelo%dim_estado
    l = modelo%dim_obs
    T = SIZE(y, 2)
    
    ALLOCATE(x_filtrado(n, T), P_filtrado(n, n, T))
    ALLOCATE(x_pred(n), P_pred(n, n), v(l), S(l, l), K(n, l))
    
    log_L = 0.0d0
    x_filtrado(:,1) = modelo%x0
    P_filtrado(:,:,1) = modelo%P0
    
    DO i = 1, T
       ! 1. Predicao
       x_pred = MATMUL(modelo%F, x_filtrado(:,MAX(1,i-1)))
       P_pred = MATMUL(modelo%F, MATMUL(P_filtrado(:,:,MAX(1,i-1)), TRANSPOSE(modelo%F))) + modelo%Q
       
       ! 2. Inovacao
       v = y(:,i) - MATMUL(modelo%H, x_pred)
       S = MATMUL(modelo%H, MATMUL(P_pred, TRANSPOSE(modelo%H))) + modelo%R
       
       ! 3. Ganho (Simplificado: assumindo inversão direta p/ exemplo, mas doutoral usa dpotrf)
       ! K = P_pred * H' * inv(S)
       ! ... (Implementação completa com LAPACK no capítulo real) ...
       
       ! 4. Atualizacao
       x_filtrado(:,i) = x_pred ! + K * v
       P_filtrado(:,:,i) = P_pred ! - K * H * P_pred
    END DO
  END SUBROUTINE filtrar_kalman

END MODULE mod_kalman
