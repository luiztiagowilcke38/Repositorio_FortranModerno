module mod_schrodinger
    use mod_constantes
    use mod_malha
    implicit none

!==================================================================================================
! Módulo: mod_schrodinger
! Descrição: Resolve a equação de Schrödinger 2D na seção transversal (y,z) do canal.
!            Calcula autovalores (sub-bandas) e autovetores (funções de onda).
!            Integração com LAPACK para solução do problema de autovalores.
!==================================================================================================

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character :: jobz, uplo
            integer :: n, lda, info, lwork
            real(dp) :: a(lda, *), w(*), work(*)
        end subroutine dsyev
    end interface

contains

    subroutine resolver_schrodinger_2d(m, pot_transv, energias, psi)
        type(t_malha), intent(in) :: m
        real(dp), intent(in) :: pot_transv(m%ny, m%nz) ! Potencial médio na seção
        real(dp), intent(out) :: energias(m%ny * m%nz)
        real(dp), intent(out) :: psi(m%ny * m%nz, m%ny * m%nz)
        
        real(dp), allocatable :: ham(:, :)
        real(dp), allocatable :: work(:)
        integer :: n_total, j, k, row, info, lwork
        real(dp) :: t_y, t_z

        n_total = m%ny * m%nz
        allocate(ham(n_total, n_total))
        ham = 0.0d0

        ! Termos de hopping (Energia cinética)
        t_y = h_bar**2 / (2.0d0 * meff * m%dy**2)
        t_z = h_bar**2 / (2.0d0 * meff * m%dz**2)

        ! Construção da Matriz Hamiltoniana (Tight-Binding efetivo)
        do j = 1, m%ny
            do k = 1, m%nz
                row = (k - 1) * m%ny + j
                
                ! Diagonal (Energia onsite: Potencial + Cinética 2D)
                ham(row, row) = 2.0d0 * t_y + 2.0d0 * t_z + pot_transv(j, k) * q_electron
                
                ! Acoplamentos (Vizinhos mais próximos)
                if (j > 1) ham(row, row - 1) = -t_y
                if (j < m%ny) ham(row, row + 1) = -t_y
                if (k > 1) ham(row, row - m%ny) = -t_z
                if (k < m%nz) ham(row, row + m%ny) = -t_z
            end do
        end do

        ! Chamada ao LAPACK (DSYEV)
        lwork = 3 * n_total - 1
        allocate(work(lwork))
        
        call dsyev('V', 'U', n_total, ham, n_total, energias, work, lwork, info)

        if (info /= 0) then
            print *, "Erro no LAPACK DSYEV: info = ", info
        else
            psi = ham ! DSYEV retorna os autovetores na própria matriz a
        end if

        deallocate(ham, work)
    end subroutine resolver_schrodinger_2d

    ! Calcula a densidade eletrônica a partir dos estados ocupados
    subroutine calcular_densidade(m, energias, psi, n_dens)
        type(t_malha), intent(in) :: m
        real(dp), intent(in) :: energias(:)
        real(dp), intent(in) :: psi(:, :)
        real(dp), intent(out) :: n_dens(m%nx, m%ny, m%nz)
        
        real(dp) :: ef, occ
        integer :: i, j, k, s, row
        
        ef = 0.0d0 ! Potencial químico de referência
        n_dens = 0.0d0
        
        do i = 1, m%nx
            do s = 1, size(energias)
                ! Distribuição de Fermi-Dirac simplificada para elétrons
                occ = 1.0d0 / (exp((energias(s) - ef) / (kb * temperatura)) + 1.0d0)
                
                do j = 1, m%ny
                    do k = 1, m%nz
                        row = (k - 1) * m%ny + j
                        n_dens(i, j, k) = n_dens(i, j, k) + occ * (psi(row, s)**2)
                    end do
                end do
            end do
        end do
        
        ! Normalização (Simples para demonstração)
        n_dens = n_dens / (m%dy * m%dz)
    end subroutine calcular_densidade

end module mod_schrodinger
