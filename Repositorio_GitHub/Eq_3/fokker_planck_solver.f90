module mod_fokker_planck
    use iso_fortran_env
    implicit none
    
    integer, parameter :: NX = 500
    real(real64), parameter :: L = 10.0_real64
    real(real64), parameter :: DX = L / NX
    real(real64), parameter :: DT = 0.001_real64
    real(real64), parameter :: K_DRIFT = 1.0_real64
    real(real64), parameter :: DIFFUSION = 0.5_real64
    
contains

    subroutine step_fp(p, dp)
        real(real64), intent(in)  :: p(NX)
        real(real64), intent(out) :: dp(NX)
        real(real64) :: drift(NX), drift_p(NX), lap_p(NX)
        integer :: i, im1, ip1
        real(real64) :: x
        
        do i = 1, NX
            x = (i - NX/2) * DX
            drift(i) = -K_DRIFT * x ! Drift de Ornstein-Uhlenbeck
            
            im1 = max(i-1, 1)
            ip1 = min(i+1, NX)
            
            ! Fluxo Advectivo (Upwind) + Fluxo Difusivo
            drift_p(i) = (drift(ip1)*p(ip1) - drift(im1)*p(im1)) / (2.0_real64 * DX)
            lap_p(i) = (p(ip1) - 2.0_real64*p(i) + p(im1)) / (DX**2)
            
            dp(i) = -drift_p(i) + DIFFUSION * lap_p(i)
        end do
    end subroutine step_fp

end module mod_fokker_planck

program main_fokker_planck
    use mod_fokker_planck
    implicit none
    real(real64) :: p(NX)
    integer :: t, i
    
    ! Condição inicial: Delta de Dirac aproximada (Gaussiana estreita)
    do i = 1, NX
        p(i) = exp(-((i - NX/2) * DX)**2 / 0.1_real64)
    end do
    p = p / (sum(p) * DX) ! Normalização
    
    print *, "# Simulando Equação de Fokker-Planck (Ornstein-Uhlenbeck)"
    
    do t = 1, 2000
        block
            real(real64) :: dp(NX)
            call step_fp(p, dp)
            p = p + DT * dp
        end block
        
        if (mod(t, 500) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | Variance: ", sum(((real([(i,i=1,NX)],real64))-NX/2.0_real64)**2 * DX**2 * p) * DX
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_fokker_planck
