module mod_constantes
    implicit none
    
!==================================================================================================
! Módulo: mod_constantes
! Descrição: Define as constantes físicas fundamentais e parâmetros específicos para a simulação
!            de dispositivos semicondutores em escala nanométrica (FinFET 2nm).
!==================================================================================================

    ! Precisão numérica (Double Precision)
    integer, parameter :: dp = kind(1.0d0)
    
    ! Constantes Físicas Universais (SI)
    real(dp), parameter :: q_electron = 1.602176634d-19   ! Carga do elétron [C]
    real(dp), parameter :: h_bar     = 1.054571817d-34   ! Constante de Planck reduzida [J.s]
    real(dp), parameter :: m0        = 9.10938356d-31    ! Massa de repouso do elétron [kg]
    real(dp), parameter :: eps0      = 8.8541878128d-12  ! Permissividade do vácuo [F/m]
    real(dp), parameter :: kb        = 1.380649d-23      ! Constante de Boltzmann [J/K]
    
    ! Parâmetros do Silício (Si)
    real(dp), parameter :: eps_si   = 11.7d0 * eps0     ! Permissividade do Silício
    real(dp), parameter :: meff     = 0.19d0 * m0       ! Massa efetiva transversal (exemplo para FinFET)
    real(dp), parameter :: eg_si    = 1.12d0 * q_electron ! Bandgap do Si em Joules [J]
    
    ! Temperatura de Operação
    real(dp), parameter :: temperatura = 300.d0          ! [K]
    real(dp), parameter :: vt          = (kb * temperatura) / q_electron ! Potencial térmico [V]
    
    ! Parâmetros Geométricos do FinFET 2nm (Estimativa de Design)
    real(dp), parameter :: l_canal  = 12.0d-9           ! Comprimento do canal [m] (Gate Length)
    real(dp), parameter :: w_fin    = 2.0d-9            ! Largura do fin [m]
    real(dp), parameter :: h_fin    = 10.0d-9           ! Altura do fin [m]
    real(dp), parameter :: t_ox     = 0.8d-9            ! Espessura do óxido equivalente (EOT) [m]
    real(dp), parameter :: eps_ox   = 3.9d0 * eps0      ! Permissividade do SiO2 (ou High-k equivalente)

contains

    subroutine print_info_fisica()
        print *, "--- Parâmetros Físicos do FinFET 2nm ---"
        print *, "Comprimento do Canal: ", l_canal * 1.0d9, " nm"
        print *, "Largura do Fin:       ", w_fin * 1.0d9,   " nm"
        print *, "Altura do Fin:        ", h_fin * 1.0d9,   " nm"
        print *, "EOT:                  ", t_ox * 1.0d9,    " nm"
        print *, "Temperatura:          ", temperatura,     " K"
        print *, "----------------------------------------"
    end subroutine print_info_fisica

end module mod_constantes
