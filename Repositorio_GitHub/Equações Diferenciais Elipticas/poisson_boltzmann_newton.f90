!=======================================================================================
!  MODULO: mod_poisson_boltzmann
!  OBJETIVO: Solucionador Newton-Krylov para bioeletrostatica.
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_pb_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: ParametrosPB
        real(real64) :: kappa2, dx, dy
        integer :: nx, ny
    end type ParametrosPB
    type :: SistemaPB
        type(ParametrosPB) :: conf
        real(real64), allocatable :: potencial(:,:), densidade_carga(:,:)
    end type SistemaPB
end module mod_pb_tipos

module mod_pb_solucionador
    use mod_pb_tipos
    implicit none
contains
    subroutine resolver_pb(sis, tol_newton, max_it)
        type(SistemaPB), intent(inout) :: sis
        real(real64), intent(in) :: tol_newton
        integer, intent(in) :: max_it
        real(real64), allocatable :: res(:,:), delta(:,:), jac(:,:)
        real(real64) :: h2inv, erro, d_u, d_v
        integer :: i, j, k, n
        h2inv = 1.0_8/(sis%conf%dx**2); allocate(res(sis%conf%nx, sis%conf%ny), delta(sis%conf%nx, sis%conf%ny), jac(sis%conf%nx, sis%conf%ny))
        do n = 1, max_it
            erro = 0.0
            do j = 2, sis%conf%ny-1; do i = 2, sis%conf%nx-1
                res(i,j) = h2inv*(sis%potencial(i+1,j)+sis%potencial(i-1,j)+sis%potencial(i,j+1)+sis%potencial(i,j-1)-4*sis%potencial(i,j)) - &
                           sis%conf%kappa2*sinh(sis%potencial(i,j)) + sis%densidade_carga(i,j)
                jac(i,j) = -4.0_8*h2inv - sis%conf%kappa2*cosh(sis%potencial(i,j))
                erro = erro + res(i,j)**2
            end do; end do
            if (sqrt(erro/(sis%conf%nx*sis%conf%ny)) < tol_newton) exit
            delta = 0.0
            do k = 1, 100; do j = 2, sis%conf%ny-1; do i = 2, sis%conf%nx-1
                delta(i,j) = (-res(i,j)-h2inv*(delta(i+1,j)+delta(i-1,j)+delta(i,j+1)+delta(i,j-1))) / jac(i,j)
            end do; end do; end do
            sis%potencial = sis%potencial + 0.7_8 * delta
        end do
    end subroutine resolver_pb
end module mod_pb_solucionador

program pb_final
    use mod_pb_solucionador
    implicit none
    type(SistemaPB) :: pb
    pb%conf%nx=128; pb%conf%ny=128; pb%conf%dx=0.5; pb%conf%dy=0.5; pb%conf%kappa2=2.0
    allocate(pb%potencial(128,128), pb%densidade_carga(128,128))
    pb%potencial=0; pb%densidade_carga=0; pb%densidade_carga(64,64)=10.0
    call resolver_pb(pb, 1e-9_8, 20)
    print *, "Autor: Luiz Tiago Wilcke. Potencial Max:", maxval(pb%potencial)
end program pb_final
