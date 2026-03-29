!==================================================================================================
! MODULO: MOD_GLOBAL_ANNEALING (ALGORITMO DE RECOZIMENTO SIMULADO)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Metodo stocastico para otimizacao global inspirado no resfriamento lento
! de metais. Permite "subidas" temporarias na energia para escapar de minimos locais.
!
! REFERENCIA: Numerical Recipes in Fortran 90 / 10.9. Simulated Annealing.
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE mod_global_annealing
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: simulated_annealing, sa_control_type

  INTEGER, PARAMETER :: DP = real64

  TYPE :: sa_control_type
     REAL(DP) :: t_start = 1.0_DP
     REAL(DP) :: t_min = 1.0E-6_DP
     REAL(DP) :: cooling_rate = 0.95_DP
     INTEGER  :: n_steps_per_t = 100
     REAL(DP) :: step_size = 0.1_DP
  END TYPE sa_control_type

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! ROTINA DE OTIMIZACAO GLOBAL
  !------------------------------------------------------------------------------------------------
  SUBROUTINE simulated_annealing(x_best, f_best, n, func, control)
    REAL(DP), INTENT(INOUT) :: x_best(n)
    REAL(DP), INTENT(OUT)   :: f_best
    INTEGER,  INTENT(IN)    :: n
    TYPE(sa_control_type), INTENT(IN) :: control
    
    INTERFACE
       FUNCTION func(x)
          IMPORT :: DP
          REAL(DP), INTENT(IN) :: x(:)
          REAL(DP) :: func
       END FUNCTION func
    END INTERFACE

    REAL(DP) :: x_curr(n), x_trial(n), f_curr, f_trial, t, delta_f, prob
    REAL(DP) :: rand_val
    INTEGER  :: i, step

    ! Inicializacao
    x_curr = x_best
    f_curr = func(x_curr)
    x_best = x_curr
    f_best = f_curr
    t = control%t_start
    
    ! Loop de temperatura (Esquema de Resfriamento)
    DO WHILE (t > control%t_min)
       DO step = 1, control%n_steps_per_t
          ! Proposta de novo estado (Perturbacao Aleatoria)
          DO i = 1, n
             CALL RANDOM_NUMBER(rand_val)
             x_trial(i) = x_curr(i) + (rand_val - 0.5_DP) * control%step_size
          END DO
          
          f_trial = func(x_trial)
          delta_f = f_trial - f_curr
          
          ! Criterio de Metropolis
          IF (delta_f < 0.0_DP) THEN
             x_curr = x_trial ; f_curr = f_trial
             ! Atualiza o melhor global encontrado
             IF (f_curr < f_best) THEN
                x_best = x_curr ; f_best = f_curr
             END IF
          ELSE
             CALL RANDOM_NUMBER(rand_val)
             prob = EXP(-delta_f / t)
             IF (rand_val < prob) THEN
                x_curr = x_trial ; f_curr = f_trial
             END IF
          END IF
       END DO
       
       ! Resfriamento
       t = t * control%cooling_rate
    END DO
  END SUBROUTINE simulated_annealing

END MODULE mod_global_annealing
