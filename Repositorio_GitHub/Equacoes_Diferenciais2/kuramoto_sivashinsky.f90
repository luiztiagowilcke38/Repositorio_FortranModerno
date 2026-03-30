module mod_ks_solver
    use iso_fortran_env
    implicit none
    
    real(real64), parameter :: L = 100.0_real64 ! Comprimento do domínio
    integer, parameter :: NX = 512              ! Número de pontos
    real(real64), parameter :: DX = L / NX
    real(real64), parameter :: DT = 0.01_real64 ! Passo de tempo (estabilidade rigorosa)
    
contains

    subroutine step_ks(u, du)
        real(real64), intent(in)  :: u(NX)
        real(real64), intent(out) :: du(NX)
        integer :: i, im1, im2, ip1, ip2
        real(real64) :: u_x, u_xx, u_xxxx
        
        do i = 1, NX
            ! Índices com condições de contorno periódicas
            im1 = mod(i - 2 + NX, NX) + 1
            im2 = mod(i - 3 + NX, NX) + 1
            ip1 = mod(i, NX) + 1
            ip2 = mod(mod(i, NX) + 1, NX) + 1
            
            ! Diferenças finitas centrais de 2a ordem
            u_x    = (u(ip1) - u(im1)) / (2.0_real64 * DX)
            u_xx   = (u(ip1) - 2.0_real64*u(i) + u(im1)) / (DX**2)
            
            ! 4a derivada (Hiperdifusão) - Erro O(DX^2)
            u_xxxx = (u(ip2) - 4.0_real64*u(ip1) + 6.0_real64*u(i) - 4.0_real64*u(im1) + u(im2)) / (DX**4)
            
            ! Equação de Kuramoto-Sivashinsky: u_t = -u_xxxx - u_xx - u*u_x
            du(i) = -u_xxxx - u_xx - u(i)*u_x
        end do
    end subroutine step_ks

    subroutine rk4_ks(u)
        real(real64), intent(inout) :: u(NX)
        real(real64) :: k1(NX), k2(NX), k3(NX), k4(NX), utmp(NX)
        
        call step_ks(u, k1)
        
        utmp = u + 0.5_real64 * DT * k1
        call step_ks(utmp, k2)
        
        utmp = u + 0.5_real64 * DT * k2
        call step_ks(utmp, k3)
        
        utmp = u + DT * k3
        call step_ks(utmp, k4)
        
        u = u + (DT/6.0_real64) * (k1 + 2.0_real64*k2 + 2.0_real64*k3 + k4)
    end subroutine rk4_ks

end module mod_ks_solver

program main_ks
    use mod_ks_solver
    implicit none
    real(real64) :: u(NX)
    integer :: i, t, n_steps, plot_freq
    
    ! Condição inicial: Perturbação senoidal com ruído
    call random_seed()
    do i = 1, NX
        u(i) = sin(2.0_real64 * 3.14159265_real64 * i * DX / L) + 0.1_real64 * (rand() - 0.5_real64)
    end do
    
    n_steps = 10000
    plot_freq = 500
    
    print *, "# Iniciando Simulação Kuramoto-Sivashinsky (Caos Espaço-Temporal)"
    print *, "# L=", L, " NX=", NX, " DT=", DT
    
    do t = 1, n_steps
        call rk4_ks(u)
        
        if (mod(t, plot_freq) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | u_mean: ", sum(u)/NX
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_ks
