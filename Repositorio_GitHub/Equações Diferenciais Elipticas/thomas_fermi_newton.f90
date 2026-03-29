!=======================================================================================
!  MODULO: mod_thomas_fermi_malha
!  OBJETIVO: Solucionador para o Potencial de Thomas-Fermi em malha logaritmica.
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_tf_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: ConfigTF
        real(real64) :: r_max, gamma
        integer :: n_pontos
        real(real64), allocatable :: r(:), dr(:)
    end type ConfigTF
end module mod_tf_tipos

module mod_tf_solucionador
    use mod_tf_tipos
    implicit none
contains
    subroutine gerar_malha(conf)
        type(ConfigTF), intent(inout) :: conf; integer :: i; real(real64) :: a, r_min=1e-5
        allocate(conf%r(conf%n_pontos), conf%dr(conf%n_pontos)); a = log(conf%r_max/r_min)/(conf%n_pontos-1)
        do i = 1, conf%n_pontos; conf%r(i) = r_min*exp(a*(i-1)); end do
        do i = 2, conf%n_pontos-1; conf%dr(i) = 0.5_8*(conf%r(i+1)-conf%r(i-1)); end do
    end subroutine gerar_malha
    subroutine resolver_tf(conf, phi, tol, max_it)
        type(ConfigTF), intent(in) :: conf; real(real64), intent(inout) :: phi(:)
        real(real64), intent(in) :: tol; integer, intent(in) :: max_it
        real(real64), allocatable :: res(:), delta(:), jac(:)
        real(real64) :: r2_i, r2_p, r2_m, erro; integer :: i, k, it
        allocate(res(conf%n_pontos), delta(conf%n_pontos), jac(conf%n_pontos))
        do it = 1, max_it
            erro = 0.0
            do i = 2, conf%n_pontos-1
                r2_i=conf%r(i)**2; r2_p=((conf%r(i+1)+conf%r(i))/2)**2; r2_m=((conf%r(i)+conf%r(i-1))/2)**2
                res(i) = (1.0_8/(r2_i*conf%dr(i))) * (r2_p*(phi(i+1)-phi(i))/(conf%r(i+1)-conf%r(i)) - &
                         r2_m*(phi(i)-phi(i-1))/(conf%r(i)-conf%r(i-1))) - conf%gamma*(phi(i)**1.5)
                jac(i) = -2.0_8/(conf%dr(i)**2) - 1.5_8*conf%gamma*sqrt(phi(i))
                erro = erro + res(i)**2
            end do
            if (sqrt(erro/conf%n_pontos) < tol) exit
            delta = 0.0
            do k=1,100; do i=2,conf%n_pontos-1
                delta(i) = (-res(i)-1.0_8/(conf%dr(i)**2)*(delta(i+1)+delta(i-1)))/jac(i)
            end do; end do
            phi = phi + 0.4_8 * delta
        end do
    end subroutine resolver_tf
end module mod_tf_solucionador

program tf_final
    use mod_tf_solucionador
    implicit none
    type(ConfigTF) :: cfg
    real(8), allocatable :: phi(:)
    cfg%n_pontos=500; cfg%r_max=20.0; cfg%gamma=1.0; call gerar_malha(cfg)
    allocate(phi(500)); phi=1.0/cfg%r; phi(1)=79.0; phi(500)=0.0
    call resolver_tf(cfg, phi, 1e-10_8, 50)
    print *, "Autor: Luiz Tiago Wilcke. Potencial r=1.0:", phi(200)
end program tf_final
