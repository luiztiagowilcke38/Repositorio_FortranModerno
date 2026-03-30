module mod_cahn_hilliard
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 128
    real(real64), parameter :: DX = 1.0_real64
    real(real64), parameter :: DT = 0.05_real64
    real(real64), parameter :: KAPPA = 2.0_real64 ! Energia de superfície
    real(real64), parameter :: MOBILITY = 1.0_real64
    
contains

    subroutine step_ch(phi, phi_new)
        real(real64), intent(in)  :: phi(N, N)
        real(real64), intent(out) :: phi_new(N, N)
        real(real64) :: mu(N, N), lap_mu(N, N)
        integer :: i, j, im1, ip1, jm1, jp1
        
        ! 1. Calcular o Potencial Químico mu = phi^3 - phi - kappa*Lap(phi)
        do j = 1, N
            do i = 1, N
                im1 = mod(i - 2 + N, N) + 1
                ip1 = mod(i, N) + 1
                jm1 = mod(j - 2 + N, N) + 1
                jp1 = mod(j, N) + 1
                
                mu(i, j) = phi(i, j)**3 - phi(i, j) - &
                           KAPPA * (phi(ip1, j) + phi(im1, j) + phi(i, jp1) + phi(i, jm1) - 4.0_real64*phi(i, j)) / (DX**2)
            end do
        end do
        
        ! 2. Evolução Temporal: phi_new = phi + DT * Mobility * Lap(mu)
        do j = 1, N
            do i = 1, N
                im1 = mod(i - 2 + N, N) + 1
                ip1 = mod(i, N) + 1
                jm1 = mod(j - 2 + N, N) + 1
                jp1 = mod(j, N) + 1
                
                lap_mu(i, j) = (mu(ip1, j) + mu(im1, j) + mu(i, jp1) + mu(i, jm1) - 4.0_real64*mu(i, j)) / (DX**2)
                phi_new(i, j) = phi(i, j) + DT * MOBILITY * lap_mu(i, j)
            end do
        end do
    end subroutine step_ch

end module mod_cahn_hilliard

program main_cahn_hilliard
    use mod_cahn_hilliard
    implicit none
    real(real64) :: phi(N, N), phi_next(N, N)
    integer :: t, i, j
    
    ! Condição inicial: Mistura homogênea com ruído (Spontaneous Spinodal Decomposition)
    call random_seed()
    do j = 1, N
        do i = 1, N
            call random_number(phi(i, j))
            phi(i, j) = 0.05_real64 * (phi(i, j) - 0.5_real64) ! Média zero, pequena perturbação
        end do
    end do
    
    print *, "# Iniciando Simulação Cahn-Hilliard (Separação de Fases 2D)"
    
    do t = 1, 5000
        call step_ch(phi, phi_next)
        phi = phi_next
        
        if (mod(t, 500) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | Energy Proxy: ", sum(phi**2)/N**2
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_cahn_hilliard
