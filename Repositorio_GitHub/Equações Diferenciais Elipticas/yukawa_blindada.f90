!=======================================================================================
!  MODULO: mod_yukawa_blindada
!  OBJETIVO: Solucionador para o Potencial de Yukawa.
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_yukawa_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: ConfiguracaoYukawa
        real(real64) :: kappa, dx, dy
        integer :: nx, ny
    end type ConfiguracaoYukawa
    type :: CampoYukawa
        real(real64), allocatable :: potencial(:,:), densidade_carga(:,:)
    end type CampoYukawa
end module mod_yukawa_tipos

module mod_yukawa_solucionador
    use mod_yukawa_tipos
    implicit none
contains
    subroutine resolver_yukawa(conf, campo, tol, max_it)
        type(ConfiguracaoYukawa), intent(in) :: conf
        type(CampoYukawa), intent(inout) :: campo
        real(real64), intent(in) :: tol
        integer, intent(in) :: max_it
        real(real64) :: dx2, dy2, fator_inv, residuo, pot_ant
        integer :: i, j, it
        dx2=conf%dx**2; dy2=conf%dy**2; fator_inv=1.0_8/(2.0_8/dx2+2.0_8/dy2+conf%kappa**2)
        do it = 1, max_it
            residuo = 0.0
            !$OMP PARALLEL DO PRIVATE(i, j, pot_ant) REDUCTION(max:residuo)
            do j = 2, conf%ny - 1; do i = 2, conf%nx - 1
                pot_ant = campo%potencial(i, j)
                campo%potencial(i, j) = ((campo%potencial(i+1,j)+campo%potencial(i-1,j))/dx2 + &
                                         (campo%potencial(i,j+1)+campo%potencial(i,j-1))/dy2 + &
                                         campo%densidade_carga(i,j)) * fator_inv
                residuo = max(residuo, abs(campo%potencial(i,j) - pot_ant))
            end do; end do
            if (residuo < tol) exit
        end do
    end subroutine resolver_yukawa
end module mod_yukawa_solucionador

program yukawa_final
    use mod_yukawa_solucionador
    implicit none
    type(ConfiguracaoYukawa) :: cfg
    type(CampoYukawa) :: cam
    cfg%nx=200; cfg%ny=200; cfg%dx=0.05; cfg%dy=0.05; cfg%kappa=5.0
    allocate(cam%potencial(200,200), cam%densidade_carga(200,200))
    cam%potencial=0; cam%densidade_carga=0; cam%densidade_carga(100,100)=1000.0
    call resolver_yukawa(cfg, cam, 1e-10, 50000)
    print *, "Autor: Luiz Tiago Wilcke. Potencial Cental:", cam%potencial(100,100)
end program yukawa_final
