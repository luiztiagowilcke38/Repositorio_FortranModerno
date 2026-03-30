module mod_mhd_dynamo
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 128
    real(real64), parameter :: DX = 1.0_real64
    real(real64), parameter :: DT = 0.01_real64
    real(real64), parameter :: ETA = 0.05_real64 ! Difusividade Magnética
    
contains

    subroutine step_dynamo(bx, by, u, v, dbx, dby)
        real(real64), intent(in)  :: bx(N, N), by(N, N), u(N, N), v(N, N)
        real(real64), intent(out) :: dbx(N, N), dby(N, N)
        real(real64) :: emf_z(N, N)
        integer :: i, j, im1, ip1, jm1, jp1
        
        ! 1. Calcular EMF (Eletromotriz) Z = u*By - v*Bx
        emf_z = u * by - v * bx
        
        ! 2. Indução: dB/dt = Curl(EMF) + ETA * Lap(B)
        do j = 1, N
            do i = 1, N
                im1 = mod(i - 2 + N, N) + 1
                ip1 = mod(i, N) + 1
                jm1 = mod(j - 2 + N, N) + 1
                jp1 = mod(j, N) + 1
                
                ! dBx/dt = d(EMF_z)/dy + ETA * Lap(Bx)
                dbx(i, j) = (emf_z(i, jp1) - emf_z(i, jm1)) / (2.0_real64 * DX) + &
                            ETA * (bx(ip1, j) + bx(im1, j) + bx(i, jp1) + bx(i, jm1) - 4.0_real64*bx(i, j)) / (DX**2)
                
                ! dBy/dt = -d(EMF_z)/dx + ETA * Lap(By)
                dby(i, j) = -(emf_z(ip1, j) - emf_z(im1, j)) / (2.0_real64 * DX) + &
                            ETA * (by(ip1, j) + by(im1, j) + by(i, jp1) + by(i, jm1) - 4.0_real64*by(i, j)) / (DX**2)
            end do
        end do
    end subroutine step_dynamo

end module mod_mhd_dynamo

program main_dynamo
    use mod_mhd_dynamo
    implicit none
    real(real64) :: bx(N, N), by(N, N), u(N, N), v(N, N)
    real(real64) :: dbx(N, N), dby(N, N)
    integer :: t, i, j
    
    ! Velocidade de cisalhamento (Shear flow) que pode gerar campo (Dínamo)
    do j = 1, N
        do i = 1, N
            u(i, j) = sin(2.0_real64 * 3.14159_real64 * j / N)
            v(i, j) = cos(2.0_real64 * 3.14159_real64 * i / N)
        end do
    end do
    
    ! Semente magnética inicial
    call random_seed()
    call random_number(bx)
    bx = 0.01_real64 * (bx - 0.5_real64)
    by = 0.0_real64
    
    print *, "# Simulando Dínamo MHD (Crescimento de Campo Magnético em Escoamento)"
    
    do t = 1, 5000
        call step_dynamo(bx, by, u, v, dbx, dby)
        bx = bx + DT * dbx
        by = by + DT * dby
        
        if (mod(t, 500) == 0) then
            print '(A,I6,A,F12.6)', "Step: ", t, " | Magnetic Energy: ", sum(bx**2 + by**2) / (N**2)
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_dynamo
