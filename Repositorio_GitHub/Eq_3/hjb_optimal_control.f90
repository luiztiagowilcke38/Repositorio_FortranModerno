module mod_hjb_optimal_control
    use iso_fortran_env
    implicit none
    
    integer, parameter :: NX = 100
    real(real64), parameter :: DX = 0.1_real64
    real(real64), parameter :: DT = 0.001_real64
    real(real64), parameter :: R_DISCOUNT = 0.05_real64 ! Taxa de Desconto
    real(real64), parameter :: SIGMA = 0.1_real64     ! Difusão (Ruído)
    
contains

    ! Hamiltoniano H(p) = max_a { -a*p - 0.5*a^2 } = 0.5*p^2 (para custo quadrático)
    function hamiltonian(p) result(h)
        real(real64), intent(in) :: p
        real(real64) :: h
        h = 0.5_real64 * p**2
    end function hamiltonian

    subroutine step_hjb(v, dv)
        real(real64), intent(in)  :: v(NX)
        real(real64), intent(out) :: dv(NX)
        real(real64) :: v_x, v_xx, h_val, utility
        integer :: i, im1, ip1
        
        do i = 1, NX
            im1 = max(i-1, 1)
            ip1 = min(i+1, NX)
            
            v_x = (v(ip1) - v(im1)) / (2.0_real64 * DX)
            v_xx = (v(ip1) - 2.0_real64*v(i) + v(im1)) / (DX**2)
            
            h_val = hamiltonian(v_x)
            utility = (real(i, real64)*DX - 5.0_real64)**2 ! Custo de estado
            
            ! dV/dt = R*V - H(grad V) - 0.5*sigma^2*Lap(V) - Utility
            dv(i) = R_DISCOUNT * v(i) - h_val - 0.5_real64 * SIGMA**2 * v_xx - utility
        end do
    end subroutine step_hjb

end module mod_hjb_optimal_control

program main_hjb
    use mod_hjb_optimal_control
    implicit none
    real(real64) :: v(NX), dv(NX)
    integer :: t, i
    
    ! Chute inicial para o valor (V=0)
    v = 0.0_real64
    
    print *, "# Simulando Equação de Hamilton-Jacobi-Bellman (Controle Ótimo Estocástico)"
    
    ! Integração no tempo (Iteração de Valor)
    do t = 1, 5000
        call step_hjb(v, dv)
        v = v - DT * dv ! Note o sinal para convergência para o estado estacionário
        
        if (mod(t, 1000) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | Max Value: ", maxval(v)
        end if
    end do
    
    print *, "# Solução de HJB Estacionária Concluída."
end program main_hjb
