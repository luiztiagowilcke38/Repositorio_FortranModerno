module mod_swift_hohenberg
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 128
    real(real64), parameter :: DX = 0.5_real64
    real(real64), parameter :: DT = 0.01_real64
    real(real64), parameter :: R_PARAM = 0.5_real64 ! Parâmetro de controle (Acima de 0 gera padrões)
    
contains

    subroutine step_sh(u, u_new)
        real(real64), intent(in)  :: u(N, N)
        real(real64), intent(out) :: u_new(N, N)
        real(real64) :: lap_u(N, N), lap_lap_u(N, N)
        integer :: i, j, im1, ip1, jm1, jp1
        
        ! 1. Laplaciano de u
        do j = 1, N
            do i = 1, N
                im1 = mod(i - 2 + N, N) + 1
                ip1 = mod(i, N) + 1
                jm1 = mod(j - 2 + N, N) + 1
                jp1 = mod(j, N) + 1
                
                lap_u(i, j) = (u(ip1, j) + u(im1, j) + u(i, jp1) + u(i, jm1) - 4.0_real64*u(i, j)) / (DX**2)
            end do
        end do
        
        ! 2. Laplaciano do Laplaciano (Bi-Laplaciano)
        do j = 1, N
            do i = 1, N
                im1 = mod(i - 2 + N, N) + 1
                ip1 = mod(i, N) + 1
                jm1 = mod(j - 2 + N, N) + 1
                jp1 = mod(j, N) + 1
                
                lap_lap_u(i, j) = (lap_u(ip1, j) + lap_u(im1, j) + lap_u(i, jp1) + lap_u(i, jm1) - 4.0_real64*lap_u(i, j)) / (DX**2)
            end do
        end do
        
        ! 3. Evolução: u_t = r*u - (1 + Lap)^2 u - u^3
        ! Operador (1 + Lap)^2 = 1 + 2*Lap + Lap^2
        do j = 1, N
            do i = 1, N
                u_new(i, j) = u(i, j) + DT * ( R_PARAM*u(i, j) - (u(i, j) + 2.0_real64*lap_u(i, j) + lap_lap_u(i, j)) - u(i, j)**3 )
            end do
        end do
    end subroutine step_sh

end module mod_swift_hohenberg

program main_swift_hohenberg
    use mod_swift_hohenberg
    implicit none
    real(real64) :: u(N, N), u_next(N, N)
    integer :: t, i, j
    
    ! Condição inicial com ruído
    call random_seed()
    do j = 1, N
        do i = 1, N
            call random_number(u(i, j))
            u(i, j) = 0.1_real64 * (u(i, j) - 0.5_real64)
        end do
    end do
    
    print *, "# Simulando Equação de Swift-Hohenberg (Formação de Padrões)"
    
    do t = 1, 10000
        call step_sh(u, u_next)
        u = u_next
        if (mod(t, 1000) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | RMS Amplitude: ", sqrt(sum(u**2)/N**2)
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_swift_hohenberg
