module mod_sine_gordon
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 512
    real(real64), parameter :: L = 40.0_real64
    real(real64), parameter :: DX = L / N
    real(real64), parameter :: DT = 0.05_real64
    
contains

    ! d2phi/dt2 = d2phi/dx2 - sin(phi)
    ! Reformulada como sistema de 1a ordem:
    ! dphi/dt = v
    ! dv/dt = d2phi/dx2 - sin(phi)
    subroutine step_sg(phi, v, dphi, dv)
        real(real64), intent(in)  :: phi(N), v(N)
        real(real64), intent(out) :: dphi(N), dv(N)
        integer :: i, im1, ip1
        real(real64) :: phi_xx
        
        do i = 1, N
            im1 = mod(i - 2 + N, N) + 1
            ip1 = mod(i, N) + 1
            
            phi_xx = (phi(ip1) - 2.0_real64*phi(i) + phi(im1)) / (DX**2)
            
            dphi(i) = v(i)
            dv(i) = phi_xx - sin(phi(i))
        end do
    end subroutine step_sg

end module mod_sine_gordon

program main_sine_gordon
    use mod_sine_gordon
    implicit none
    real(real64) :: phi(N), v(N)
    integer :: t, i
    real(real64) :: x, energy
    
    ! Condição inicial: Colisão de Kink e Anti-kink
    ! phi(x,0) = 4 * arctan(exp((x-x0)/sqrt(1-v^2)))
    do i = 1, N
        x = (i - N/2) * DX
        phi(i) = 4.0_real64 * (atan(exp(x + 10.0_real64)) - atan(exp(x - 10.0_real64)))
        v(i) = 0.0_real64 ! Kinks inicialmente em repouso ou v baixa
    end do
    
    print *, "# Simulando Equação de Sine-Gordon (Solitons Kink-Antikink)"
    
    do t = 1, 2000
        block
            real(real64) :: dphi(N), dv(N)
            call step_sg(phi, v, dphi, dv)
            ! Integrador Semi-Implícito (Symplectic Euler) para conservar fase
            v = v + DT * dv
            phi = phi + DT * v
        end block
        
        if (mod(t, 400) == 0) then
            energy = sum(0.5_real64 * v**2 + 0.5_real64 * ((phi - phi)**2 / DX**2) + (1.0_real64 - cos(phi))) * DX
            print '(A,I6,A,F8.4)', "Step: ", t, " | Energy Proxy: ", energy
        end if
    end do
    
    print *, "# Simulação de Sine-Gordon Concluída."
end program main_sine_gordon
