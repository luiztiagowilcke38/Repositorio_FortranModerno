module mod_caputo_derivative
    use iso_fortran_env
    implicit none
    
    real(real64), parameter :: PI = 3.14159265_real64
    
contains

    ! Função Gama via aproximação de Lanczos (Doctoral Level Accuracy)
    recursive function lanczos_gamma(x) result(g)
        real(real64), intent(in) :: x
        real(real64) :: g, tmp, ser
        real(real64), dimension(7) :: coeff = [ 76.18009172947146_real64, -86.50532032941677_real64, &
                                               24.01409824083091_real64, -1.231739572450155_real64, &
                                               0.1208650973866179e-2_real64, -0.5395239384953e-5_real64, &
                                               0.0_real64 ] ! Simplified Lanczos
        integer :: j
        
        if (x < 1.0_real64) then
            g = PI / (sin(PI * x) * lanczos_gamma(1.0_real64 - x))
            return
        end if
        
        tmp = x + 5.5_real64
        tmp = (x + 0.5_real64) * log(tmp) - tmp
        ser = 1.000000000190015_real64
        do j = 1, 6
            ser = ser + coeff(j) / (x + j)
        end do
        g = exp(tmp + log(2.5066282746310005_real64 * ser / x))
    end function gamma

    ! Coeficientes L1 para a derivada de Caputo: w_j = (j+1)^(1-alpha) - j^(1-alpha)
    subroutine get_l1_weights(n, alpha, weights)
        integer, intent(in) :: n
        real(real64), intent(in) :: alpha
        real(real64), intent(out) :: weights(0:n)
        integer :: j
        
        do j = 0, n
            weights(j) = (real(j+1,real64)**(1.0_real64 - alpha)) - (real(j,real64)**(1.0_real64 - alpha))
        end do
    end subroutine get_l1_weights

end module mod_caputo_derivative

module mod_fractional_diffusion
    use mod_caputo_derivative
    use iso_fortran_env
    implicit none
    
    integer, parameter :: NX = 100
    real(real64), parameter :: DX = 1.0_real64
    real(real64), parameter :: ALPHA_ORDER = 0.7_real64 ! Ordem fracionária (Sub-difusão)
    real(real64), parameter :: K_DIFF = 1.0_real64
    
contains

    subroutine step_fractional(u_history, n_steps, dt, u_new)
        real(real64), intent(in)  :: u_history(NX, 0:n_steps)
        integer, intent(in)       :: n_steps
        real(real64), intent(in)  :: dt
        real(real64), intent(out) :: u_new(NX)
        real(real64) :: weights(0:n_steps)
        real(real64) :: memory_sum(NX), lap_u(NX)
        real(real64) :: coeff_l1
        integer :: i, j, im1, ip1
        
        call get_l1_weights(n_steps, ALPHA_ORDER, weights)
        
        ! Termo de memória: Sum_{j=0}^{n-1} w_j * (u^{n-j} - u^{n-j-1})
        memory_sum = 0.0_real64
        do j = 1, n_steps
            memory_sum = memory_sum + weights(j) * (u_history(:, n_steps-j+1) - u_history(:, n_steps-j))
        end do
        
        ! Coeficiente de normalização do esquema L1
        coeff_l1 = dt**ALPHA_ORDER * lanczos_gamma(2.0_real64 - ALPHA_ORDER)
        
        ! Resolver u_new = u_old + coeff * Lap(u_old) - memory_contribution
        do i = 1, NX
            im1 = max(i-1, 1)
            ip1 = min(i+1, NX)
            lap_u(i) = (u_history(ip1, n_steps) - 2.0_real64*u_history(i, n_steps) + u_history(im1, n_steps)) / (DX**2)
            
            ! Evolução implícita simplificada (ou explícita com compensação de memória)
            u_new(i) = u_history(i, n_steps) + coeff_l1 * K_DIFF * lap_u(i) - memory_sum(i)
        end do
    end subroutine step_fractional

end module mod_fractional_diffusion

program main_fractional
    use mod_fractional_diffusion
    implicit none
    real(real64), allocatable :: u_history(:, :)
    integer :: n_max = 500
    integer :: t, i
    real(real64) :: dt = 0.01_real64
    
    allocate(u_history(NX, 0:n_max))
    u_history = 0.0_real64
    
    ! Condição inicial: Pulso no centro
    u_history(NX/2, 0) = 10.0_real64
    
    print *, "# Simulando Difusão Fracionária (Ordem Alpha=", ALPHA_ORDER, ")"
    
    do t = 0, n_max - 1
        call step_fractional(u_history(:, 0:t), t, dt, u_history(:, t+1))
        
        if (mod(t, 100) == 0) then
            print '(A,I6,A,F12.6)', "Step: ", t, " | Peak Value: ", maxval(u_history(:, t+1))
        end if
    end do
    
    print *, "# Simulação concluída. Observe o decaimento mais lento (Long Memory)."
end program main_fractional
