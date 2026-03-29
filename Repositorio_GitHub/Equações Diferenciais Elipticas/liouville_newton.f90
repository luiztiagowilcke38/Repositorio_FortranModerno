!=======================================================================================
!  MODULO: mod_liouville_continuacao
!  OBJETIVO: Solucionador para a Equacao de Liouville com continuaçao de parametro.
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_liouville_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: ConfigLiouville
        real(real64) :: L_dominio, lambda_alvo
        integer :: nx, ny, n_passos
    end type ConfigLiouville
end module mod_liouville_tipos

module mod_liouville_solucionador
    use mod_liouville_tipos
    implicit none
contains
    subroutine resolver_newton(conf, u, lambda_at, tol, max_it)
        type(ConfigLiouville), intent(in) :: conf
        real(real64), intent(inout) :: u(:,:); real(real64), intent(in) :: lambda_at, tol
        integer, intent(in) :: max_it
        real(real64), allocatable :: res(:,:), delta(:,:), jac(:,:)
        real(real64) :: h2inv, erro; integer :: i, j, k, it
        h2inv = 1.0_8/(conf%L_dominio/(conf%nx-1))**2; allocate(res(conf%nx,conf%ny), delta(conf%nx,conf%ny), jac(conf%nx,conf%ny))
        do it = 1, max_it
            erro = 0.0
            do j = 2, conf%ny-1; do i = 2, conf%nx-1
                res(i,j) = h2inv*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)) + lambda_at*exp(u(i,j))
                jac(i,j) = -4.0_8*h2inv + lambda_at*exp(u(i,j))
                erro = erro + res(i,j)**2
            end do; end do
            if (sqrt(erro/(conf%nx*conf%ny)) < tol) exit
            delta = 0.0
            do k = 1, 150; do j = 2, conf%ny-1; do i = 2, conf%nx-1
                delta(i,j) = (-res(i,j)-h2inv*(delta(i+1,j)+delta(i-1,j)+delta(i,j+1)+delta(i,j-1))) / jac(i,j)
            end do; end do; end do
            u = u + 0.5_8 * delta
        end do
    end subroutine resolver_newton
    subroutine exec_continuacao(conf, u)
        type(ConfigLiouville), intent(in) :: conf; real(real64), intent(inout) :: u(:,:)
        real(real64) :: lam_p, lam_at; integer :: s
        lam_p = conf%lambda_alvo / conf%n_passos
        do s = 1, conf%n_passos
            lam_at = s * lam_p
            call resolver_newton(conf, u, lam_at, 1e-9_8, 50)
        end do
    end subroutine exec_continuacao
end module mod_liouville_solucionador

program liouville_final
    use mod_liouville_solucionador
    implicit none
    type(ConfigLiouville) :: cfg
    real(8), allocatable :: u(:,:)
    cfg%nx=64; cfg%ny=64; cfg%L_dominio=1.0; cfg%lambda_alvo=2.0; cfg%n_passos=10
    allocate(u(64,64)); u = -1.0
    call exec_continuacao(cfg, u)
    print *, "Autor: Luiz Tiago Wilcke. U Max:", maxval(u)
end program liouville_final
