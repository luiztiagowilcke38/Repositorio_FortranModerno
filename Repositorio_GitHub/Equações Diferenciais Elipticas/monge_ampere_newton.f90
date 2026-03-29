!=======================================================================================
!  MODULO: mod_monge_ampere
!  OBJETIVO: Solucionador para a Equacao de Monge-Ampere (Totalmente Nao-Linear).
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_ma_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: ConfigMA
        integer :: nx, ny; real(real64) :: dx, dy
        real(real64), allocatable :: f(:,:), u(:,:)
    end type ConfigMA
end module mod_ma_tipos

module mod_ma_solucionador
    use mod_ma_tipos
    implicit none
contains
    subroutine resolver_ma(sis, tol, max_it)
        type(ConfigMA), intent(inout) :: sis; real(real64), intent(in) :: tol; integer, intent(in) :: max_it
        real(real64), allocatable :: res(:,:), delta(:,:); real(real64) :: uxx, uyy, uxy, h2inv, erro
        integer :: i, j, k, it
        h2inv = 1.0_8/(sis%dx**2); allocate(res(sis%nx, sis%ny), delta(sis%nx, sis%ny))
        do it = 1, max_it
            erro = 0.0
            do j = 2, sis%ny-1; do i = 2, sis%nx-1
                uxx = (sis%u(i+1,j)-2*sis%u(i,j)+sis%u(i-1,j))*h2inv
                uyy = (sis%u(i,j+1)-2*sis%u(i,j)+sis%u(i,j-1))*h2inv
                uxy = (sis%u(i+1,j+1)-sis%u(i+1,j-1)-sis%u(i-1,j+1)+sis%u(i-1,j-1))*0.25_8*h2inv
                res(i,j) = uxx*uyy-uxy**2-sis%f(i,j); erro = erro + res(i,j)**2
            end do; end do
            if (sqrt(erro/(sis%nx*sis%ny)) < tol) exit
            delta = 0.0
            do k=1,100; do j=2,sis%ny-1; do i=2,sis%nx-1
                uxx = (sis%u(i+1,j)-2*sis%u(i,j)+sis%u(i-1,j))*h2inv
                uyy = (sis%u(i,j+1)-2*sis%u(i,j)+sis%u(i,j-1))*h2inv
                delta(i,j) = ((uyy*h2inv*(delta(i+1,j)+delta(i-1,j)) + uxx*h2inv*(delta(i,j+1)+delta(i,j-1))) - res(i,j)) / (2*h2inv*(uyy+uxx))
            end do; end do; end do
            sis%u = sis%u + 0.2_8 * delta
        end do
    end subroutine resolver_ma
end module mod_ma_solucionador

program ma_final
    use mod_ma_solucionador
    implicit none
    type(ConfigMA) :: ma; integer :: i, j
    ma%nx=50; ma%ny=50; ma%dx=0.02; ma%dy=0.02; allocate(ma%u(50,50), ma%f(50,50))
    do j=1,50; do i=1,50; ma%u(i,j)=0.5*(real(i-25,8)**2+real(j-25,8)**2)*ma%dx**2; ma%f(i,j)=1.0; end do; end do
    call resolver_ma(ma, 1e-8_8, 50)
    print *, "Autor: Luiz Tiago Wilcke. U Center:", ma%u(25,25)
end program ma_final
