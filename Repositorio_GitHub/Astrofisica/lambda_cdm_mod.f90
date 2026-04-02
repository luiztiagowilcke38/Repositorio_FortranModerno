module lambda_cdm_mod
    implicit none
    private
    public :: ModeloCosmologico_t

    integer, parameter :: dp = kind(1.0d0)

    !> Classe base e modelo de dados para o universo Lambda-CDM
    type :: ModeloCosmologico_t
        real(dp) :: densidade_radiacao_atual !< Densidade de radiação atual (\Omega_R)
        real(dp) :: densidade_materia_atual  !< Densidade de matéria atual (\Omega_M)
        real(dp) :: densidade_energia_escura !< Densidade de energia escura atual (\Omega_\Lambda)
        real(dp) :: parametro_hubble         !< Parâmetro de Hubble atual em km/s/Mpc
        real(dp) :: velocidade_luz           !< Constante de velocidade da luz (km/s)
    contains
        procedure :: iniciar => inicializar_modelo_cosmologico
        procedure :: hubble => calcular_parametro_hubble
        procedure :: idade_do_universo => calcular_idade_universo
        procedure :: distancia_comovel => calcular_distancia_comovel
        procedure :: distancia_luminosidade => calcular_distancia_luminosidade
        procedure :: resolver_historico => solucionar_historico_cosmologico
    end type ModeloCosmologico_t

contains

    !> Inicializa os parâmetros da classe Cosmológica
    subroutine inicializar_modelo_cosmologico(este, rad_atual, mat_atual, lambda_atual, h_zero)
        class(ModeloCosmologico_t), intent(out) :: este
        real(dp), intent(in) :: rad_atual, mat_atual, lambda_atual, h_zero
        este%densidade_radiacao_atual = rad_atual
        este%densidade_materia_atual = mat_atual
        este%densidade_energia_escura = lambda_atual
        ! Assumimos um universo plano: \Omega_R + \Omega_M + \Omega_\Lambda = 1
        este%parametro_hubble = h_zero
        este%velocidade_luz = 299792.458_dp ! Velocidade da luz exata em km/s
    end subroutine inicializar_modelo_cosmologico

    !> Calcula o H(a) usando a equação de Friedmann local (H em km/s/Mpc)
    function calcular_parametro_hubble(este, fator_escala) result(hubble_inst)
        class(ModeloCosmologico_t), intent(in) :: este
        real(dp), intent(in) :: fator_escala
        real(dp) :: hubble_inst
        
        hubble_inst = este%parametro_hubble * sqrt(este%densidade_radiacao_atual/(fator_escala**4) + &
                                                   este%densidade_materia_atual/(fator_escala**3) +  &
                                                   este%densidade_energia_escura)
    end function calcular_parametro_hubble

    !> Calcula a idade do universo t(a) fazendo integração numérica com a Regra de Simpson 1/3
    function calcular_idade_universo(este, fator_escala) result(tempo_decorrido)
        class(ModeloCosmologico_t), intent(in) :: este
        real(dp), intent(in) :: fator_escala
        real(dp) :: tempo_decorrido
        integer, parameter :: numero_passos = 100000
        real(dp) :: a_minimo, delta_a, a_iteracao, hubble_iteracao, soma_integral
        integer :: indice
        
        a_minimo = 1.0e-6_dp
        delta_a = (fator_escala - a_minimo) / real(numero_passos, dp)
        soma_integral = 0.0_dp
        
        do indice = 0, numero_passos
            a_iteracao = a_minimo + real(indice, dp) * delta_a
            hubble_iteracao = este%hubble(a_iteracao)
            if (indice == 0 .or. indice == numero_passos) then
                soma_integral = soma_integral + 1.0_dp / (a_iteracao * hubble_iteracao)
            else if (mod(indice, 2) == 0) then
                soma_integral = soma_integral + 2.0_dp / (a_iteracao * hubble_iteracao)
            else
                soma_integral = soma_integral + 4.0_dp / (a_iteracao * hubble_iteracao)
            end if
        end do
        
        tempo_decorrido = (delta_a / 3.0_dp) * soma_integral
        ! Converter de 1/(km/s/Mpc) para Bilhões de Anos (Sendo 1/H0 ~ 9.77 Gyr / h)
        tempo_decorrido = tempo_decorrido * 977.79222168_dp
    end function calcular_idade_universo

    !> Calcula a Distância Comóvel r(a) de um observador até uma fonte no fator_escala
    function calcular_distancia_comovel(este, fator_escala) result(distancia_r)
        class(ModeloCosmologico_t), intent(in) :: este
        real(dp), intent(in) :: fator_escala
        real(dp) :: distancia_r
        integer, parameter :: numero_passos = 100000
        real(dp) :: delta_a, a_iteracao, hubble_iteracao, soma_integral
        integer :: indice
        
        delta_a = (1.0_dp - fator_escala) / real(numero_passos, dp)
        soma_integral = 0.0_dp
        
        do indice = 0, numero_passos
            a_iteracao = fator_escala + real(indice, dp) * delta_a
            hubble_iteracao = este%hubble(a_iteracao)
            if (indice == 0 .or. indice == numero_passos) then
                soma_integral = soma_integral + 1.0_dp / (a_iteracao**2 * hubble_iteracao)
            else if (mod(indice, 2) == 0) then
                soma_integral = soma_integral + 2.0_dp / (a_iteracao**2 * hubble_iteracao)
            else
                soma_integral = soma_integral + 4.0_dp / (a_iteracao**2 * hubble_iteracao)
            end if
        end do
        
        distancia_r = este%velocidade_luz * (delta_a / 3.0_dp) * soma_integral ! Resultado final em Mpc
    end function calcular_distancia_comovel

    !> Distância de Luminosidade D_L(a) baseada na distância comóvel para modelos planos
    function calcular_distancia_luminosidade(este, fator_escala) result(dist_lum)
        class(ModeloCosmologico_t), intent(in) :: este
        real(dp), intent(in) :: fator_escala
        real(dp) :: dist_lum
        dist_lum = este%distancia_comovel(fator_escala) / fator_escala
    end function calcular_distancia_luminosidade
    
    !> Resolve todo o histórico do universo e escreve resultados num arquivo formatado I/O
    subroutine solucionar_historico_cosmologico(este, nome_arquivo)
         class(ModeloCosmologico_t), intent(in) :: este
         character(len=*), intent(in) :: nome_arquivo
         integer :: ind, unidade_saida
         real(dp) :: param_a, redshift_z, tempo_vida, dist_comovel, dist_lum
         
         open(newunit=unidade_saida, file=nome_arquivo, status='replace', action='write')
         write(unidade_saida, '(A)') '# Modelo Cosmologico Lambda-CDM Dinamico: Resolucao Completa'
         write(unidade_saida, '(A,3E16.8)') '# Parametros [Rad., Mat., E.E.]: ', &
               este%densidade_radiacao_atual, este%densidade_materia_atual, &
               este%densidade_energia_escura
         write(unidade_saida, '(A,E16.8)') '# Constante de Hubble H_0 [km/s/Mpc]: ', este%parametro_hubble
         write(unidade_saida, '(A)') '# redshift_z     fator_a         Idade_Univers(Gyr)  Dist_Comovel(Mpc)    Dist_Lum(Mpc)'
         
         ! Percorre o universo do presente (z=0, a=1) até o Big Bang aproximado (z~99)
         do ind = 1000, 10, -5
            param_a = real(ind, dp) / 1000.0_dp
            redshift_z = (1.0_dp / param_a) - 1.0_dp
            tempo_vida = este%idade_do_universo(param_a)
            dist_comovel = este%distancia_comovel(param_a)
            dist_lum = este%distancia_luminosidade(param_a)
            write(unidade_saida, '(5E18.8)') redshift_z, param_a, tempo_vida, dist_comovel, dist_lum
         end do
         close(unidade_saida)
         print *, 'Histórico cósmico salvo e resolvido com sucesso.'
         print *, 'O arquivo gerado é: ', trim(nome_arquivo)
    end subroutine solucionar_historico_cosmologico

end module lambda_cdm_mod
