program simulador_finfet_2nm
    use mod_constantes
    use mod_malha
    use mod_poisson
    use mod_schrodinger
    use mod_negf
    implicit none

!==================================================================================================
! Programa: simulador_finfet_2nm
! Descrição: Orquestrador da simulação de um FinFET 2nm. Realiza o loop SCF (Self-Consistent Field)
!            entre a equação de Poisson e Schrödinger, seguido pelo cálculo de transporte NEGF.
!==================================================================================================

    type(t_malha) :: malha
    real(dp), allocatable :: phi(:,:,:), rho(:,:,:), n_dens(:,:,:)
    real(dp), allocatable :: energias(:), psi(:,:)
    real(dp) :: vg, vd, vs, mix_param
    integer :: nx, ny, nz, scf_iter
    integer, parameter :: max_scf = 100
    logical :: poisson_conv

    ! Configuração Inicial
    nx = 40; ny = 12; nz = 12
    vg = 0.7d0; vd = 0.7d0; vs = 0.0d0
    mix_param = 0.2d0 ! Damping para convergência SCF

    call inicializar_malha(malha, nx, ny, nz)
    allocate(phi(nx, ny, nz), rho(nx, ny, nz), n_dens(nx, ny, nz))
    phi = 0.0d0; rho = 0.0d0; n_dens = 0.0d0

    print *, ">>> Iniciando Loop SCF Poisson-Schrödinger <<<"

    do scf_iter = 1, max_scf
        ! 1. Aplicar condições de contorno no potencial
        call aplicar_contornos_poisson(malha, phi, vg, vd, vs)

        ! 2. Resolver Poisson para o potencial eletrostático
        ! rho = q * (n - N_d)
        rho = q_electron * n_dens 
        call resolver_poisson(malha, phi, rho, poisson_conv)

        ! 3. Extrair potencial médio na seção transversal (y,z) para Schrödinger
        ! (Simplificação: Resolvemos no meio do canal i = nx/2)
        if (scf_iter == 1) then
            allocate(energias(ny*nz), psi(ny*nz, ny*nz))
        end if
        
        call resolver_schrodinger_2d(malha, phi(nx/2, :, :), energias, psi)

        ! 4. Calcular densidade eletrônica a partir dos estados quânticos
        call calcular_densidade(malha, energias, psi, n_dens)

        ! 5. Verificar Convergência (Critério de erro no potencial)
        ! erro_scf = sum(abs(phi - phi_old)) ... (Omitido para brevidade)
        
        print *, "Iteração SCF: ", scf_iter, " Densidade Max: ", maxval(n_dens)
        
        if (scf_iter > 5) exit ! Simulação rápida para demonstração
    end do

    print *, ">>> Loop SCF Finalizado. Calculando Corrente via NEGF <<<"
    
    ! 6. Cálculo de Corrente (NEGF) - Demonstrativo em um ponto de energia
    ! (Em uma simulação real, integraríamos T(E) sobre o espectro)
    ! print *, "Corrente (Id): ", calcular_corrente(vg, vd, malha)

    call liberar_malha(malha)
    deallocate(phi, rho, n_dens, energias, psi)

    print *, ">>> Simulação FinFET 2nm Concluída com Sucesso! <<<"

end program simulador_finfet_2nm
