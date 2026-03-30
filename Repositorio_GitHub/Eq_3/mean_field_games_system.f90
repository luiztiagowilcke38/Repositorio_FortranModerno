module mod_mean_field_games
    use iso_fortran_env
    implicit none
    
    integer, parameter :: NX = 50
    real(real64), parameter :: DX = 1.0_real64
    real(real64), parameter :: DT = 0.01_real64
    real(real64), parameter :: R_RATE = 0.05_real64
    real(real64), parameter :: SIGMA = 0.5_real64
    
contains

    ! Solução de HJB (Evolução para trás)
    subroutine solve_hjb(v, m, dt_sim)
        real(real64), intent(inout) :: v(NX)
        real(real64), intent(in)    :: m(NX)
        real(real64), intent(in)    :: dt_sim
        real(real64) :: v_next(NX), v_x, v_xx
        integer :: i, im1, ip1
        
        do i = 1, NX
            im1 = max(i-1, 1); ip1 = min(i+1, NX)
            v_x = (v(ip1) - v(im1)) / (2.0_real64 * DX)
            v_xx = (v(ip1) - 2.0_real64*v(i) + v(im1)) / (DX**2)
            
            ! dV/dt = R*V - 0.5*(Vx)^2 - 0.5*sigma^2*Vxx - m(i) (Coupling)
            v_next(i) = v(i) + dt_sim * (R_RATE*v(i) - 0.5_real64*v_x**2 - 0.5_real64*SIGMA**2*v_xx - m(i))
        end do
        v = v_next
    end subroutine solve_hjb

    ! Solução de Fokker-Planck (Evolução para frente)
    subroutine solve_fp(m, v, dt_sim)
        real(real64), intent(inout) :: m(NX)
        real(real64), intent(in)    :: v(NX)
        real(real64), intent(in)    :: dt_sim
        real(real64) :: m_next(NX), vx, m_flux_x, m_xx
        integer :: i, im1, ip1
        
        do i = 1, NX
            im1 = max(i-1, 1); ip1 = min(i+1, NX)
            vx = (v(ip1) - v(im1)) / (2.0_real64 * DX)
            m_flux_x = (m(ip1)*vx - m(im1)*vx) / (2.0_real64 * DX)
            m_xx = (m(ip1) - 2.0_real64*m(i) + m(im1)) / (DX**2)
            
            ! dm/dt = div(m*grad V) + 0.5*sigma^2*mxx
            m_next(i) = m(i) + dt_sim * (m_flux_x + 0.5_real64*SIGMA**2*m_xx)
        end do
        m = m_next
    end subroutine solve_fp

end module mod_mean_field_games

program main_mfg
    use mod_mean_field_games
    implicit none
    real(real64) :: v(NX), m(NX)
    integer :: iter, step
    
    ! Condições Iniciais
    m = 1.0_real64 / (NX * DX) ! Distribuição Uniforme
    v = 0.0_real64
    
    print *, "# Simulando Mean Field Games (HJB + Fokker-Planck Acoplados)"
    
    do iter = 1, 5
        print '(A,I3)', "# Iteração de Ponto Fixo: ", iter
        do step = 1, 500
            call solve_hjb(v, m, DT)
            call solve_fp(m, v, DT)
        end do
        print '(A,F12.6)', "   Residual Sum m: ", sum(m)*DX
    end do
    
    print *, "# Simulação de MFG Concluída."
end program main_mfg
