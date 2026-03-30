module mod_tov_solver
    use iso_fortran_env
    implicit none
    
    ! Constantes Astrofísicas (CGS simplificado para demonstração)
    real(real64), parameter :: G = 6.674e-8_real64
    real(real64), parameter :: C = 2.998e10_real64
    real(real64), parameter :: PI = 3.14159265_real64
    
    ! Parâmetros da Politropia (Equação de Estado P = K * rho^gamma)
    real(real64), parameter :: K_REL = 1.23e15_real64 ! K para estrelas de nêutrons (exemplo)
    real(real64), parameter :: GAMMA_REL = 2.0_real64
    
contains

    function get_rho(p) result(rho)
        real(real64), intent(in) :: p
        real(real64) :: rho
        if (p > 0.0_real64) then
            rho = (p / K_REL)**(1.0_real64 / GAMMA_REL)
        else
            rho = 0.0_real64
        end if
    end function get_rho

    subroutine tov_derivs(r, y, dy)
        real(real64), intent(in) :: r
        real(real64), intent(in) :: y(2) ! y(1) = P, y(2) = M
        real(real64), intent(out) :: dy(2)
        real(real64) :: p, m, rho
        
        p = y(1)
        m = y(2)
        rho = get_rho(p)
        
        if (p <= 0.0_real64 .or. r <= 0.0_real64) then
            dy = 0.0_real64
            return
        end if
        
        ! Equação de Continuidade da Massa: dM/dr = 4*pi*r^2*rho
        dy(2) = 4.0_real64 * PI * r**2 * rho
        
        ! Equação TOV: dP/dr
        dy(1) = -(G / r**2) * (rho + p/C**2) * (m + 4.0_real64*PI*r**3*p/(C**2)) / &
                 (1.0_real64 - 2.0_real64*G*m/(r*C**2))
    end subroutine tov_derivs

end module mod_tov_solver

program main_tov
    use mod_tov_solver
    implicit none
    real(real64) :: r, dr, y(2), k1(2), k2(2), k3(2), k4(2), ytmp(2)
    integer :: step
    
    ! Condições Iniciais no centro da estrela
    r = 10.0_real64 ! cm (pequeno offset para evitar divisão por zero)
    dr = 100.0_real64 ! cm
    y(1) = 1.0e35_real64 ! Pressão Central [Ba]
    y(2) = 4.0_real64/3.0_real64 * PI * r**3 * get_rho(y(1)) ! Massa inicial
    
    print *, "# Simulando Equações de Tolman-Oppenheimer-Volkoff (Estrelas de Nêutrons)"
    print *, "R [km] | Massa [M_sun] | Pressão [Ba]"
    
    do step = 1, 100000
        call tov_derivs(r, y, k1)
        ytmp = y + 0.5_real64 * dr * k1
        call tov_derivs(r + 0.5_real64*dr, ytmp, k2)
        ytmp = y + 0.5_real64 * dr * k2
        call tov_derivs(r + 0.5_real64*dr, ytmp, k3)
        ytmp = y + dr * k3
        call tov_derivs(r + dr, ytmp, k4)
        
        y = y + (dr/6.0_real64) * (k1 + 2.0_real64*k2 + 2.0_real64*k3 + k4)
        r = r + dr
        
        if (y(1) <= 0.0_real64) exit ! Superfície atingida
        
        if (mod(step, 2000) == 0) then
            print '(F8.2, A, F8.4, A, E12.4)', r/1.0e5_real64, " | ", y(2)/1.989e33_real64, " | ", y(1)
        end if
    end do
    
    print *, "# Superfície atingida em R = ", r/1.0e5_real64, " km"
    print *, "# Massa Final = ", y(2)/1.989e33_real64, " M_sun"
end program main_tov
