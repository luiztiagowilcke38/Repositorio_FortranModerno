!=======================================================================================
!  MODULO: mod_helmholtz_pml
!  OBJETIVO: Solucionador mestre para a equacao de Helmholtz complexa com Camadas de 
!            Absorcao Perfeitamente Casadas (PML).
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_helmholtz_pml_tipos
    use, intrinsic :: iso_fortran_env
    implicit none

    type :: MalhaPML
        integer :: nx, ny, npml
        real(real64) :: dx, dy, f_onda, k0
        real(real64), allocatable :: sigma_x(:), sigma_y(:)
        complex(real64), allocatable :: s_x(:), s_y(:)
    end type MalhaPML

    type :: SistemaHelmholtz
        type(MalhaPML) :: malha
        complex(real64), allocatable :: potencial(:,:)
        complex(real64), allocatable :: fonte(:,:)
    end type SistemaHelmholtz

end module mod_helmholtz_pml_tipos

module mod_helmholtz_pml_core
    use mod_helmholtz_pml_tipos
    implicit none

    interface
        module subroutine inicializar_pml(malha, nx, ny, npml, f, dx)
            type(MalhaPML), intent(inout) :: malha
            integer, intent(in) :: nx, ny, npml
            real(real64), intent(in) :: f, dx
        end subroutine inicializar_pml

        module subroutine aplicar_operador_helmholtz(sis, ent, sai)
            type(SistemaHelmholtz), intent(in) :: sis
            complex(real64), intent(in) :: ent(:,:)
            complex(real64), intent(out) :: sai(:,:)
        end subroutine aplicar_operador_helmholtz
    end interface

end module mod_helmholtz_pml_core

submodule (mod_helmholtz_pml_core) sub_helmholtz_pml_impl
    implicit none
contains
    module subroutine inicializar_pml(malha, nx, ny, npml, f, dx)
        type(MalhaPML), intent(inout) :: malha
        integer, intent(in) :: nx, ny, npml
        real(real64), intent(in) :: f, dx
        integer :: i
        real(real64) :: r_teorico, sigma_max, dist_norm
        malha%nx = nx; malha%ny = ny; malha%npml = npml
        malha%f_onda = f; malha%dx = dx; malha%dy = dx
        malha%k0 = 2.0_8 * 3.14159265358979_8 * f / 343.0_8 
        allocate(malha%sigma_x(nx), malha%sigma_y(ny))
        allocate(malha%s_x(nx), malha%s_y(ny))
        malha%sigma_x = 0.0_8; malha%sigma_y = 0.0_8
        r_teorico = 1e-8_8
        sigma_max = -(3.0_8 / 2.0_8) * log(r_teorico) / (real(npml,8)*dx)
        do i = 1, npml
            dist_norm = real(npml - i + 1, 8) / real(npml, 8)
            malha%sigma_x(i) = sigma_max * (dist_norm**2)
            malha%sigma_x(nx - i + 1) = sigma_max * (dist_norm**2)
            malha%sigma_y(i) = sigma_max * (dist_norm**2)
            malha%sigma_y(ny - i + 1) = sigma_max * (dist_norm**2)
        end do
        do i = 1, nx; malha%s_x(i) = 1.0_8 - (0.0_8, 1.0_8) * malha%sigma_x(i) / (2.0_8 * 3.14159_8 * f); end do
        do i = 1, ny; malha%s_y(i) = 1.0_8 - (0.0_8, 1.0_8) * malha%sigma_y(i) / (2.0_8 * 3.14159_8 * f); end do
    end subroutine inicializar_pml

    module subroutine aplicar_operador_helmholtz(sis, ent, sai)
        type(SistemaHelmholtz), intent(in) :: sis
        complex(real64), intent(in) :: ent(:,:)
        complex(real64), intent(out) :: sai(:,:)
        integer :: i, j
        complex(real64) :: dx2inv, dy2inv, k2, s_x, s_y, term_x, term_y
        dx2inv = 1.0_8 / (sis%malha%dx**2); dy2inv = 1.0_8 / (sis%malha%dy**2); k2 = (sis%malha%k0)**2
        !$OMP PARALLEL DO PRIVATE(i, j, s_x, s_y, term_x, term_y) SHARED(ent, sai, sis)
        do j = 2, sis%malha%ny - 1
            do i = 2, sis%malha%nx - 1
                s_x = sis%malha%s_x(i); s_y = sis%malha%s_y(j)
                term_x = ( (1.0_8/sis%malha%s_x(i+1) + 1.0_8/s_x)*(ent(i+1,j) - ent(i,j)) - &
                           (1.0_8/sis%malha%s_x(i-1) + 1.0_8/s_x)*(ent(i,j) - ent(i-1,j)) ) * 0.5_8 * dx2inv
                term_y = ( (1.0_8/sis%malha%s_y(j+1) + 1.0_8/s_y)*(ent(i,j+1) - ent(i,j)) - &
                           (1.0_8/sis%malha%s_y(j-1) + 1.0_8/s_y)*(ent(i,j) - ent(i,j-1)) ) * 0.5_8 * dy2inv
                sai(i, j) = (1.0_8/(s_x * s_y)) * (term_x + term_y) + k2 * ent(i, j)
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine aplicar_operador_helmholtz
end submodule sub_helmholtz_pml_impl

module mod_helmholtz_solucionador
    use mod_helmholtz_pml_core
    implicit none
contains
    subroutine resolver_bicgstab(sis, tolerancia, max_iteracoes)
        type(SistemaHelmholtz), intent(inout) :: sis
        real(real64), intent(in) :: tolerancia
        integer, intent(in) :: max_iteracoes
        complex(real64), allocatable :: r(:,:), r_tio(:,:), p(:,:), v(:,:), s(:,:), t(:,:)
        complex(real64) :: rho, rho_anterior, alfa, omega, beta
        real(real64) :: norma_r0, erro_relativo
        integer :: iter
        allocate(r(sis%malha%nx, sis%malha%ny), r_tio(sis%malha%nx, sis%malha%ny))
        allocate(p(sis%malha%nx, sis%malha%ny), v(sis%malha%nx, sis%malha%ny))
        allocate(s(sis%malha%nx, sis%malha%ny), t(sis%malha%nx, sis%malha%ny))
        call aplicar_operador_helmholtz(sis, sis%potencial, v)
        r = sis%fonte - v; r_tio = r
        norma_r0 = sqrt(sum(abs(r)**2))
        if (norma_r0 < 1e-15_8) return
        rho = 1.0_8; alfa = 1.0_8; omega = 1.0_8; v = 0.0_8; p = 0.0_8
        do iter = 1, max_iteracoes
            rho_anterior = rho; rho = sum(conjg(r_tio) * r)
            if (abs(rho) < 1e-25_8) exit
            if (iter == 1) then; p = r; else
                beta = (rho / rho_anterior) * (alfa / omega); p = r + beta * (p - omega * v)
            end if
            call aplicar_operador_helmholtz(sis, p, v)
            alfa = rho / sum(conjg(r_tio) * v); s = r - alfa * v
            erro_relativo = sqrt(sum(abs(s)**2)) / norma_r0
            if (erro_relativo < tolerancia) then; sis%potencial = sis%potencial + alfa * p; return; end if
            call aplicar_operador_helmholtz(sis, s, t); omega = sum(conjg(t) * s) / sum(conjg(t) * t)
            sis%potencial = sis%potencial + alfa * p + omega * s
            r = s - omega * t; erro_relativo = sqrt(sum(abs(r)**2)) / norma_r0
            if (erro_relativo < tolerancia) exit
        end do
    end subroutine resolver_bicgstab
end module mod_helmholtz_solucionador

program simulacao_pml
    use mod_helmholtz_solucionador
    implicit none
    type(SistemaHelmholtz) :: hel
    integer :: NX=256, NY=256, NPML=30
    call inicializar_pml(hel%malha, NX, NY, NPML, 500.0_8, 5.0_8/(NX-1))
    allocate(hel%potencial(NX, NY), hel%fonte(NX, NY))
    hel%potencial = (0.0_8, 0.0_8); hel%fonte = (0.0_8, 0.0_8); hel%fonte(NX/2, NY/2) = (100.0_8, 0.0_8)
    call resolver_bicgstab(hel, 1e-9_8, 2000)
    print *, "Autor: Luiz Tiago Wilcke. Potencial Central:", hel%potencial(NX/2, NY/2)
end program simulacao_pml
