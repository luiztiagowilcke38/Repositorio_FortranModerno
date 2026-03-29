!=======================================================================================
!  MODULO: mod_biharmonica_13pontos
!  OBJETIVO: Solucionador de quarta ordem para a Equacao Biharmonica (Del-4) 
!            utilizando o Stencil de 13 pontos.
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_biharmonica_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: PropriedadesPlaca
        real(real64) :: espessura, modulo_E, poisson_nu, rigidez_D
    end type PropriedadesPlaca
    type :: SistemaBiharmonico
        integer :: nx, ny
        real(real64) :: dx, dy
        real(real64), allocatable :: deflexao(:,:), carga(:,:)
    end type SistemaBiharmonico
end module mod_biharmonica_tipos

module mod_biharmonica_solucionador
    use mod_biharmonica_tipos
    implicit none
contains
    subroutine resolver_sor_13p(sis, omega, toler, max_it)
        type(SistemaBiharmonico), intent(inout) :: sis
        real(real64), intent(in) :: omega, toler
        integer, intent(in) :: max_it
        real(real64) :: h4, erro_iter, val_ant, viz_1, viz_2, viz_diag
        integer :: i, j, it
        h4 = (sis%dx * sis%dy)**2
        do it = 1, max_it
            erro_iter = 0.0_8
            !$OMP PARALLEL DO PRIVATE(i, j, val_ant, viz_1, viz_2, viz_diag) REDUCTION(+:erro_iter)
            do j = 3, sis%ny - 2; do i = 3, sis%nx - 2
                val_ant = sis%deflexao(i, j)
                viz_1 = sis%deflexao(i+1,j) + sis%deflexao(i-1,j) + sis%deflexao(i,j+1) + sis%deflexao(i,j-1)
                viz_2 = sis%deflexao(i+2,j) + sis%deflexao(i-2,j) + sis%deflexao(i,j+2) + sis%deflexao(i,j-2)
                viz_diag = sis%deflexao(i+1,j+1)+sis%deflexao(i-1,j+1)+sis%deflexao(i+1,j-1)+sis%deflexao(i-1,j-1)
                sis%deflexao(i, j) = (1.0_8 - omega) * val_ant + (omega / 20.0_8) * (8.0_8*viz_1 - 2.0_8*viz_diag - viz_2 + h4*sis%carga(i,j))
                erro_iter = erro_iter + (sis%deflexao(i,j) - val_ant)**2
            end do; end do
            if (sqrt(erro_iter/(sis%nx*sis%ny)) < toler) exit
        end do
    end subroutine resolver_sor_13p
end module mod_biharmonica_solucionador

program biharmonica_final
    use mod_biharmonica_solucionador
    implicit none
    type(SistemaBiharmonico) :: placa
    placa%nx=60; placa%ny=60; placa%dx=0.01; placa%dy=0.01
    allocate(placa%deflexao(60,60), placa%carga(60,60))
    placa%deflexao=0; placa%carga=-1e5
    call resolver_sor_13p(placa, 1.4, 1e-9, 50000)
    print *, "Autor: Luiz Tiago Wilcke. Deflexao:", placa%deflexao(30,30)
end program biharmonica_final
