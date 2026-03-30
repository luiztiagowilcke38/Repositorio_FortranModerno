MODULE mod_espaco_estados
  USE iso_fortran_env
  IMPLICIT NONE
  
  ! Estrutura para o Modelo de Espaco de Estados (Linear e Gaussiano)
  TYPE :: ModeloEspacoEstados
     INTEGER :: dim_estado, dim_obs
     REAL(real64), ALLOCATABLE :: F(:,:) ! Transicao
     REAL(real64), ALLOCATABLE :: H(:,:) ! Observacao
     REAL(real64), ALLOCATABLE :: Q(:,:) ! Covariancia do Processo
     REAL(real64), ALLOCATABLE :: R(:,:) ! Covariancia da Medicao
     REAL(real64), ALLOCATABLE :: x0(:)  ! Estado inicial
     REAL(real64), ALLOCATABLE :: P0(:,:)! Covariancia inicial
  END TYPE ModeloEspacoEstados

CONTAINS

  ! Construtor para forma companheira de um AR(p)
  SUBROUTINE construir_ar_companheira(modelo, coefs_ar, sigma_e)
    TYPE(ModeloEspacoEstados), INTENT(OUT) :: modelo
    REAL(real64), INTENT(IN) :: coefs_ar(:), sigma_e
    INTEGER :: p, i
    
    p = SIZE(coefs_ar)
    modelo%dim_estado = p
    modelo%dim_obs = 1
    
    ALLOCATE(modelo%F(p,p), modelo%H(1,p), modelo%Q(p,p), modelo%R(1,1))
    modelo%F = 0.0d0
    modelo%H = 0.0d0
    modelo%Q = 0.0d0
    
    ! Primeira linha de F com coeficientes AR
    modelo%F(1, :) = coefs_ar
    ! Sub-identidade
    DO i = 2, p
       modelo%F(i, i-1) = 1.0d0
    END DO
    
    modelo%H(1, 1) = 1.0d0
    modelo%Q(1, 1) = sigma_e**2
    modelo%R(1, 1) = 0.0d0 ! Ruido de medicao zero para AR puro
  END SUBROUTINE construir_ar_companheira

END MODULE mod_espaco_estados
