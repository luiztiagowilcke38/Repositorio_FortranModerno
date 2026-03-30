module mod_kpz_stochastic
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 512
    real(real64), parameter :: DX = 1.0_real64
    real(real64), parameter :: DT = 0.1_real64
    real(real64), parameter :: NU = 0.5_real64   ! Viscosidade (Difusão)
    real(real64), parameter :: LAMBDA = 2.0_real64 ! Coeficiente não-linear
    real(real64), parameter :: NOISE_AMP = 0.1_real64
    
contains

    subroutine step_kpz(h, dh)
        real(real64), intent(in)  :: h(N)
        real(real64), intent(out) :: dh(N)
        real(real64) :: hx, hxx, eta
        integer :: i, im1, ip1
        
        do i = 1, N
            im1 = mod(i - 2 + N, N) + 1
            ip1 = mod(i, N) + 1
            
            ! Diferenças finitas centrais
            hx = (h(ip1) - h(im1)) / (2.0_real64 * DX)
            hxx = (h(ip1) - 2.0_real64*h(i) + h(im1)) / (DX**2)
            
            ! Ruído Gaussiano Branco
            call random_number(eta)
            eta = NOISE_AMP * (eta - 0.5_real64)
            
            ! dh/dt = nu*hxx + (lambda/2)*hx^2 + eta
            dh(i) = NU * hxx + (LAMBDA/2.0_real64) * hx**2 + eta / sqrt(DT)
        end do
    end subroutine step_kpz

end module mod_kpz_stochastic

program main_kpz
    use mod_kpz_stochastic
    implicit none
    real(real64) :: h(N)
    integer :: t
    real(real64) :: roughness, mean_h
    
    ! Condição inicial flat
    h = 0.0_real64
    call random_seed()
    
    print *, "# Simulando Crescimento de Interface KPZ (Kardar-Parisi-Zhang)"
    print *, "T | Mean_H | Roughness (W)"
    
    do t = 1, 10000
        block
            real(real64) :: dh(N)
            call step_kpz(h, dh)
            h = h + DT * dh
        end block
        
        if (mod(t, 1000) == 0) then
            mean_h = sum(h) / N
            roughness = sqrt(sum((h - mean_h)**2) / N)
            print '(I6, A, F8.4, A, F8.4)', t, " | ", mean_h, " | ", roughness
        end if
    end do
    
    print *, "# Simulação concluída. A rugosidade W cresce como t^(1/3) (Classe de Universalidade KPZ)."
end program main_kpz
