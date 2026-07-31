# =============================================================================
# runtests.jl -- Pruebas de regresion del proyecto
# =============================================================================
# Blindan los resultados publicados en los cuatro documentos y verifican de
# forma reproducible lo que hasta ahora solo se habia comprobado a mano:
#
#   1. Modelos       -- el modelo generico de N eslabones reproduce las dos
#                       configuraciones ya entregadas, que son implementaciones
#                       INDEPENDIENTES de la misma fisica.
#   2. Linealizacion -- los jacobianos ANALITICOS escritos a mano en
#                       linearization.jl coinciden con la diferenciacion
#                       automatica de las EOM no lineales.
#   3. Publicados    -- los espectros y las ganancias que aparecen en el
#                       informe, el resumen ejecutivo, la presentacion y el
#                       README.
#   4. Controller    -- solve_care contra MatrixEquations.arec, propiedades
#                       algebraicas de la solucion de Riccati y Cayley-Hamilton.
#   5. Metricas      -- los valores de referencia del modulo Metrics.
#
# Ejecucion:
#     julia --project=. test/runtests.jl
#
# NOTA: no se usa Pkg.test(). El repositorio es un ENTORNO de Julia, no un
# paquete: Project.toml no declara name, uuid ni version, y Pkg.test() exige
# las tres. Correr el archivo directamente da exactamente el mismo resultado
# sin obligar a reestructurar el proyecto.
# =============================================================================

using Pkg
const PROJ_ROOT = dirname(@__DIR__)
Pkg.activate(PROJ_ROOT)

include(joinpath(PROJ_ROOT, "src", "model_simple.jl"))
include(joinpath(PROJ_ROOT, "src", "model_double.jl"))
include(joinpath(PROJ_ROOT, "src", "model_nlink.jl"))
include(joinpath(PROJ_ROOT, "src", "model_triple.jl"))
include(joinpath(PROJ_ROOT, "src", "linearization.jl"))
include(joinpath(PROJ_ROOT, "src", "controller.jl"))
include(joinpath(PROJ_ROOT, "src", "metrics.jl"))
include(joinpath(PROJ_ROOT, "src", "observer.jl"))
include(joinpath(PROJ_ROOT, "src", "sweep.jl"))

using .Model, .ModelDouble, .ModelNLink, .ModelTriple
using .Linearization, .Controller, .Metrics, .Observer, .Sweep

using Test
using LinearAlgebra
using Random
using ForwardDiff
using DifferentialEquations
using MatrixEquations: arec

Random.seed!(20260731)

# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

"""
Jacobiano df/dx por diferenciacion automatica, evaluado en `x0`, para una EOM
in-place de la forma eom!(dx, x, p, t) con fuerza de control nula.
"""
function ad_jacobian(eom!, params, n; x0=zeros(n))
    f = function (x)
        dx = zeros(eltype(x), n)
        eom!(dx, x, (params=params, F=0.0), 0.0)
        return dx
    end
    return ForwardDiff.jacobian(f, x0)
end

"""
Jacobiano df/du por diferenciacion automatica, evaluado en el equilibrio.
Devuelve un vector columna n x 1, con la misma forma que la matriz B.
"""
function ad_input_jacobian(eom!, params, n; x0=zeros(n))
    f = function (u)
        dx = zeros(typeof(u), n)
        eom!(dx, x0, (params=params, F=u), 0.0)
        return dx
    end
    return reshape(ForwardDiff.derivative(f, 0.0), n, 1)
end

"""
Polinomio caracteristico de A evaluado en A, normalizado por norm(A)^n para
que la tolerancia sea comparable entre sistemas de distinta escala.
"""
function cayley_hamilton_residual(A)
    n = size(A, 1)
    phi = Matrix{ComplexF64}(I, n, n)
    for lambda in eigvals(A)
        phi = phi * (A - lambda * I)
    end
    return norm(phi) / norm(A)^n
end

"""
Comprueba que un espectro sin friccion tiene la estructura predicha:
pares +-lambda mas un cero doble.
"""
function is_symmetric_spectrum(lambda; tol=1e-8)
    re = sort(real.(lambda))
    simetrico = maximum(abs.(re .+ reverse(re))) < tol
    reales = maximum(abs.(imag.(lambda))) < tol
    ceros = count(abs.(re) .< tol) == 2
    return simetrico && reales && ceros
end

const R_STD = reshape([0.1], 1, 1)

# ---------------------------------------------------------------------------

@testset verbose = true "proyecto-pendulos" begin

    # =======================================================================
    @testset "1. Modelos" begin
        p_simple = default_params()
        p_double = default_params_double()
        p_n1 = uniform_rods(1.0, [0.3], [1.0]; b=0.1)
        p_n2 = point_masses(1.0, [0.3, 0.3], [0.5, 0.5]; b=0.0)
        p_n3 = uniform_rods(1.0, fill(0.2, 3), fill(1/3, 3); b=0.1)

        @testset "ModelNLink(N=1) reproduce Model" begin
            # Los constructores deben coincidir con los parametros del informe.
            @test p_n1.a[1] ≈ p_simple.L
            @test p_n1.Il[1] ≈ p_simple.I

            for _ in 1:20
                x = [randn(), 2randn(), pi * (2rand() - 1), 2randn()]
                F = 10randn()
                @test state_derivative(x, p_simple, F) ≈
                      state_derivative_nlink(x, p_n1, F) atol = 1e-10
            end
        end

        @testset "ModelNLink(N=2) reproduce ModelDouble" begin
            bc = link_couplings(p_n2)
            # beta_1 = (m1+m2) L1, beta_2 = m2 L2, Gam_12 = m2 L1 L2
            @test bc.beta ≈ [0.6 * 0.5, 0.3 * 0.5]
            @test bc.J ≈ [0.6 * 0.5^2, 0.3 * 0.5^2]
            @test bc.Gam[1, 2] ≈ 0.3 * 0.5 * 0.5

            for _ in 1:20
                x = [randn(), 2randn(), pi * (2rand() - 1), 2randn(),
                     pi * (2rand() - 1), 2randn()]
                F = 10randn()
                @test state_derivative_double(x, p_double, F) ≈
                      state_derivative_nlink(x, p_n2, F) atol = 1e-10
            end
        end

        @testset "M(q) simetrica definida positiva" begin
            for p in (p_n1, p_n2, p_n3)
                N = length(p.m)
                for _ in 1:50
                    Mq = mass_matrix(p, pi .* (2rand(N) .- 1))
                    @test Mq ≈ Mq'
                    @test isposdef(Symmetric(Mq))
                end
            end
        end

        @testset "Configuracion III" begin
            p3 = default_params_triple()
            @test length(p3.m) == 3
            @test p3.a ≈ p3.l ./ 2                      # barras uniformes
            @test p3.Il ≈ (1/12) .* p3.m .* p3.l .^ 2
            @test sum(p3.m) ≈ 0.6                       # misma masa que Config. II
            @test sum(p3.l) ≈ 1.0                       # misma longitud que I y II

            bc = link_couplings(p3)
            @test bc.beta ≈ [0.16667, 0.1, 0.03333] atol = 1e-5
            @test bc.J ≈ [0.05185, 0.02963, 0.00741] atol = 1e-5

            # La capa delgada debe delegar exactamente en ModelNLink.
            for _ in 1:10
                x = [randn(), randn(), pi*(2rand()-1), randn(),
                     pi*(2rand()-1), randn(), pi*(2rand()-1), randn()]
                F = 10randn()
                @test state_derivative_triple(x, p3, F) ==
                      state_derivative_nlink(x, p3, F)
            end
        end

        @testset "friccion articular disipativa" begin
            # Con u = 0 y sin friccion del carro, la energia mecanica no puede
            # crecer. Es la unica comprobacion posible del signo del par
            # relativo entre eslabones: no hay implementacion de referencia.
            p = uniform_rods(1.0, fill(0.2, 3), fill(1/3, 3);
                             b=0.0, c=fill(0.05, 3))
            bc = link_couplings(p)
            energy = function (x)
                th, om = x[3:2:end], x[4:2:end]
                qd = vcat(x[2], om)
                T = 0.5 * dot(qd, mass_matrix(p, th) * qd)
                V = p.g * sum(bc.beta[j] * cos(th[j]) for j in 1:3)
                return T + V
            end
            for _ in 1:50
                x = [randn(), randn(), pi*(2rand()-1), randn(),
                     pi*(2rand()-1), randn(), pi*(2rand()-1), randn()]
                dx = state_derivative_nlink(x, p, 0.0)
                @test dot(ForwardDiff.gradient(energy, x), dx) < 1e-8
            end
        end
    end

    # =======================================================================
    @testset "2. Linealizacion" begin
        p_simple = default_params()
        p_double = default_params_double()

        @testset "jacobiano analitico vs automatico (simple)" begin
            ss = linearize_system(p_simple)
            @test ss.A ≈ ad_jacobian(nonlinear_eom!, p_simple, 4) atol = 1e-6
            @test ss.B ≈ ad_input_jacobian(nonlinear_eom!, p_simple, 4) atol = 1e-6
        end

        @testset "jacobiano analitico vs automatico (doble)" begin
            ss = linearize_system_double(p_double)
            @test ss.A ≈ ad_jacobian(nonlinear_eom_double!, p_double, 6) atol = 1e-6
            @test ss.B ≈ ad_input_jacobian(nonlinear_eom_double!, p_double, 6) atol = 1e-6
        end

        @testset "jacobiano analitico vs automatico (triple)" begin
            p3 = default_params_triple()
            ss = linearize_system_triple(p3)
            @test size(ss.A) == (8, 8)
            @test size(ss.C) == (4, 8)
            @test ss.A ≈ ad_jacobian(nonlinear_eom_triple!, p3, 8) atol = 1e-6
            @test ss.B ≈ ad_input_jacobian(nonlinear_eom_triple!, p3, 8) atol = 1e-6
        end

        @testset "linearize_system_nlink reproduce I y II" begin
            # El linealizador generico y los dos analiticos escritos a mano son
            # implementaciones independientes: deben coincidir.
            ss1 = linearize_system(default_params())
            ss1g = linearize_system_nlink(uniform_rods(1.0, [0.3], [1.0]; b=0.1))
            @test ss1.A ≈ ss1g.A atol = 1e-10
            @test ss1.B ≈ ss1g.B atol = 1e-10

            ss2 = linearize_system_double(default_params_double())
            ss2g = linearize_system_nlink(
                point_masses(1.0, [0.3, 0.3], [0.5, 0.5]; b=0.0))
            @test ss2.A ≈ ss2g.A atol = 1e-10
            @test ss2.B ≈ ss2g.B atol = 1e-10
        end

        @testset "estructura del espectro sin friccion" begin
            # Sin friccion, buscar q(t) = v exp(lambda t) da el problema
            # generalizado G0 v = lambda^2 M0 v: los eigenvalores de A son las
            # raices cuadradas de los de M0^-1 G0 y vienen en pares +-lambda.
            # G0 tiene una fila y una columna nulas (el carro no tiene fuerza
            # recuperadora), asi que aparece un cero doble.
            p_sin_fric = SystemParams(M=1.0, m=0.3, L=0.5, g=9.81,
                                      b=0.0, I=(1/12) * 0.3 * 1.0^2)
            @test is_symmetric_spectrum(linearize_system(p_sin_fric).eigenvalues)
            @test is_symmetric_spectrum(linearize_system_double(p_double).eigenvalues)

            # Se cumple para cualquier N: es una prediccion a priori, no un
            # ajuste a posteriori.
            for N in 1:5
                p = uniform_rods(1.0, fill(0.6/N, N), fill(1.0/N, N); b=0.0)
                ss = linearize_system_nlink(p)
                @test is_symmetric_spectrum(ss.eigenvalues; tol=1e-7)
                @test dominant_mode(ss.A).n_unstable == N
            end
        end

        @testset "la friccion del carro rompe la simetria" begin
            # Con b = 0.1 el cero doble se desdobla en {0, -0.0769} y cada par
            # +-lambda queda desbalanceado en la tercera cifra.
            @test !is_symmetric_spectrum(linearize_system(p_simple).eigenvalues)
        end
    end

    # =======================================================================
    @testset "3. Valores publicados" begin
        ss1 = linearize_system(default_params())
        ss2 = linearize_system_double(default_params_double())

        @testset "espectros de lazo abierto" begin
            @test sort(real.(ss1.eigenvalues), rev=true) ≈
                  [4.2105, 0.0, -0.0769, -4.2266] atol = 1e-3
            @test sort(real.(ss2.eigenvalues), rev=true) ≈
                  [8.5726, 4.0941, 0.0, 0.0, -4.0941, -8.5726] atol = 1e-3
        end

        @testset "rangos de Kalman" begin
            @test check_controllability(ss1).rank == 4
            @test check_observability(ss1).rank == 4
            @test check_controllability(ss2).rank == 6
            @test check_observability(ss2).rank == 6
        end

        @testset "ganancias LQR" begin
            K1 = design_lqr(ss1.A, ss1.B, diagm([1.0, 0, 10, 0]), R_STD).K
            @test vec(K1) ≈ [-3.16, -4.69, -45.39, -10.93] atol = 1e-2

            K2 = design_lqr(ss2.A, ss2.B, diagm([1.0, 0, 10, 0, 10, 0]), R_STD).K
            @test vec(K2) ≈ [3.16, 5.82, -191.55, -10.99, 228.32, 36.14] atol = 1e-2
        end

        @testset "Configuracion III" begin
            ss3 = linearize_system_triple(default_params_triple())

            @test sort(real.(ss3.eigenvalues), rev=true) ≈
                  [18.0613, 9.7740, 4.3777, 0.0,
                   -0.0625, -4.3990, -9.7813, -18.0647] atol = 1e-3
            # Tres modos inestables, uno mas que el doble.
            @test dominant_mode(ss3.A).n_unstable == 3

            @test vec(ss3.B) ≈
                  [0, 0.9455, 0, -3.6, 0, 0.9818, 0, -0.3273] atol = 1e-3

            @test check_controllability(ss3).rank == 8
            @test check_observability(ss3).rank == 8

            K3 = design_lqr(ss3.A, ss3.B,
                            diagm([1.0, 0, 10, 0, 10, 0, 10, 0]), R_STD).K
            @test vec(K3) ≈ [-3.16, -6.17, -317.43, -4.37,
                             911.95, 37.73, -667.37, -61.42] rtol = 1e-2
            @test norm(K3) ≈ 1176.0 rtol = 1e-2

            polos = sort(eigvals(ss3.A - ss3.B * K3), by=real)
            @test real.(polos) ≈ [-18.14, -18.14, -10.02, -10.02,
                                  -4.57, -4.57, -0.84, -0.84] atol = 1e-2
            @test sort(abs.(imag.(polos))) ≈ [0.79, 0.79, 0.97, 0.97,
                                              1.73, 1.73, 2.32, 2.32] atol = 1e-2
        end

        @testset "Ackermann sobre el simple" begin
            polos = [-1.0, -2.0, -3.0, -4.0]
            res = design_pole_placement(ss1.A, ss1.B, polos)
            @test vec(res.K) ≈ [-1.75, -3.75, -39.01, -9.60] atol = 1e-2
            @test sort(real.(res.eigenvalues_cl)) ≈ sort(polos) atol = 1e-8
        end
    end

    # =======================================================================
    @testset "4. Controller" begin
        casos = [
            ("simple", linearize_system(default_params()), diagm([1.0, 0, 10, 0])),
            ("doble", linearize_system_double(default_params_double()),
             diagm([1.0, 0, 10, 0, 10, 0])),
        ]

        for (nombre, ss, Q) in casos
            @testset "$nombre" begin
                A, B = ss.A, ss.B
                res = design_lqr(A, B, Q, R_STD)
                P, K = res.P, res.K

                # La implementacion propia contra la de la biblioteca.
                # MatrixEquations ya estaba declarada en Project.toml y no se
                # usaba en ningun archivo del repositorio.
                X, _, _ = arec(A, B, R_STD, Q)
                @test P ≈ X atol = 1e-8

                @test P ≈ P'
                @test isposdef(Symmetric(P))

                # Residuo de la CARE: A'P + PA - P B R^-1 B' P + Q = 0
                residuo = A' * P + P * A - P * B * inv(R_STD) * B' * P + Q
                @test norm(residuo) < 1e-8

                # Ganancia optima K = R^-1 B' P
                @test K ≈ inv(R_STD) * B' * P

                # Lazo cerrado estable
                @test all(real.(eigvals(A - B * K)) .< 0)

                # Cayley-Hamilton: el polinomio caracteristico de A anula a A
                @test cayley_hamilton_residual(A) < 1e-10
            end
        end
    end

    # =======================================================================
    @testset "5. Metricas" begin
        refs = [
            ("simple", linearize_system(default_params()),
             diagm([1.0, 0, 10, 0]),
             (pbh=1.72e-2, cond=3.60e1, cstar=22.94, theta=0.962, lambda=4.2105)),
            ("doble", linearize_system_double(default_params_double()),
             diagm([1.0, 0, 10, 0, 10, 0]),
             (pbh=2.56e-3, cond=1.61e4, cstar=8.99, theta=0.234, lambda=8.5726)),
            ("triple", linearize_system_triple(default_params_triple()),
             diagm([1.0, 0, 10, 0, 10, 0, 10, 0]),
             (pbh=4.92e-4, cond=2.157e8, cstar=3.73, theta=0.152, lambda=18.0613)),
        ]

        for (nombre, ss, Q, ref) in refs
            @testset "$nombre" begin
                A, B, C = ss.A, ss.B, ss.C
                res = design_lqr(A, B, Q, R_STD)

                @test dominant_mode(A).lambda_max ≈ ref.lambda rtol = 1e-3
                @test pbh_controllability_margin_normalized(A, B) ≈ ref.pbh rtol = 0.05
                @test controllability_condition(A, B) ≈ ref.cond rtol = 0.05

                c_star = nonsaturation_level(res.K, res.P, R_STD, B, 50.0)
                @test c_star ≈ ref.cstar rtol = 0.05
                @test max_axis_angle(res.P, c_star, 3) ≈ ref.theta rtol = 0.05

                # El margen PBH solo es cero si el par no es controlable.
                @test pbh_controllability_margin(A, B) > 0
                @test pbh_observability_margin(A, C) > 0
            end
        end

        @testset "discretizacion ZOH exacta" begin
            ss = linearize_system(default_params())
            A, B = ss.A, ss.B
            for h in (0.001, 0.01, 0.1)
                d = discretize_zoh(A, B, h)
                # Ad debe ser exactamente la exponencial matricial.
                @test d.Ad ≈ exp(A * h) atol = 1e-12
                # Bd debe satisfacer la identidad de la solucion variacional:
                # A Bd + B h ... se comprueba contra la serie de la integral.
                serie = sum(A^k * h^(k + 1) / factorial(big(k + 1)) for k in 0:20)
                @test d.Bd ≈ Float64.(serie * B) atol = 1e-12
            end
        end

        @testset "criterio de Schur" begin
            ss = linearize_system(default_params())
            A, B = ss.A, ss.B
            K = design_lqr(A, B, diagm([1.0, 0, 10, 0]), R_STD).K

            # El sistema en lazo ABIERTO es inestable: nunca es de Schur.
            @test !is_schur(discretize_zoh(A, B, 0.01).Ad)

            # En lazo cerrado, con muestreo rapido es de Schur y con muestreo
            # lento deja de serlo: hay un h_max finito entre ambos.
            d_rapido = discretize_zoh(A, B, 0.01)
            d_lento = discretize_zoh(A, B, 0.5)
            @test is_schur(d_rapido.Ad - d_rapido.Bd * K)
            @test !is_schur(d_lento.Ad - d_lento.Bd * K)
        end
    end

    # =======================================================================
    @testset "6. Degradacion con N" begin
        # La tabla mas contundente del informe: el mismo algebra escala, pero
        # el CONDICIONAMIENTO no. La familia esta normalizada (longitud total
        # 1 m, masa total 0.6 kg, barras uniformes) para que la comparacion
        # entre valores de N sea limpia. OJO: esta familia NO coincide con los
        # parametros por defecto de I y II, que tienen otras masas.
        metas = [
            (N=1, lambda=4.51, pbh=1.31e-2, cond=4.8e1, normK=52.9),
            (N=2, lambda=10.59, pbh=1.79e-3, cond=6.2e4, normK=248.7),
            (N=3, lambda=18.06, pbh=4.92e-4, cond=2.2e8, normK=1176.0),
            (N=4, lambda=26.54, pbh=1.90e-4, cond=1.7e12, normK=5266.8),
            (N=5, lambda=35.64, pbh=9.03e-5, cond=3.3e16, normK=23221.9),
        ]

        for meta in metas
            N = meta.N
            @testset "N = $N" begin
                p = uniform_rods(1.0, fill(0.6/N, N), fill(1.0/N, N); b=0.1)
                ss = linearize_system_nlink(p)
                n = 2 * (N + 1)

                Q = diagm(vcat([1.0, 0.0], repeat([10.0, 0.0], N)))
                K = design_lqr(ss.A, ss.B, Q, R_STD).K

                @test dominant_mode(ss.A).lambda_max ≈ meta.lambda rtol = 0.01
                @test dominant_mode(ss.A).n_unstable == N
                @test pbh_controllability_margin_normalized(ss.A, ss.B) ≈
                      meta.pbh rtol = 0.05
                @test controllability_condition(ss.A, ss.B) ≈ meta.cond rtol = 0.15
                @test norm(K) ≈ meta.normK rtol = 0.02

                # El rango de Kalman sigue dando n aunque el condicionamiento
                # se dispare: por eso hay que reportar el margen PBH.
                @test pbh_controllability_margin(ss.A, ss.B) > 0
            end
        end

        # El margen PBH decae SUAVEMENTE mientras cond(Ctrb) estalla: esa es
        # justamente la razon de preferirlo como metrica de barrido.
        pbhs = Float64[]
        conds = Float64[]
        for N in 1:5
            p = uniform_rods(1.0, fill(0.6/N, N), fill(1.0/N, N); b=0.1)
            ss = linearize_system_nlink(p)
            push!(pbhs, pbh_controllability_margin_normalized(ss.A, ss.B))
            push!(conds, controllability_condition(ss.A, ss.B))
        end
        @test issorted(pbhs, rev=true)          # decae monotonamente
        @test issorted(conds)                   # crece monotonamente
        @test pbhs[1] / pbhs[5] < 1e3           # dos ordenes y medio
        @test conds[5] / conds[1] > 1e14        # catorce ordenes
    end

    # =======================================================================
    @testset "7. Atlas de operabilidad" begin
        lin = linearize_system_nlink
        p3 = default_params_triple()

        configs = [
            ("simple", uniform_rods(1.0, [0.3], [1.0]; b=0.1), 197.4),
            ("doble", point_masses(1.0, [0.3, 0.3], [0.5, 0.5]; b=0.0), 77.3),
            ("triple", p3, 32.0),
        ]

        @testset "periodo de muestreo maximo" begin
            for (nombre, p, meta) in configs
                ss = lin(p)
                n = size(ss.A, 1)
                K = design_lqr(ss.A, ss.B, default_weights(n), R_STD).K
                h = max_sampling_period(ss.A, ss.B, K)
                @test 1000h ≈ meta rtol = 0.05

                # A h_max el lazo esta al borde: un poco menos es de Schur y un
                # poco mas no. Si no hubiera cambio de signo, la biseccion no
                # significaria nada.
                d_ok = discretize_zoh(ss.A, ss.B, 0.95h)
                d_no = discretize_zoh(ss.A, ss.B, 1.10h)
                @test is_schur(d_ok.Ad - d_ok.Bd * K)
                @test !is_schur(d_no.Ad - d_no.Bd * K)
            end
        end

        @testset "optimo interior de la masa del carro" begin
            Ms = collect(range(0.04, 1.0, length=25))
            res = sweep_1d(lin, p3, :M, Ms)
            pbh = [r.pbh_normalized for r in res]
            i = argmax(pbh)
            # El maximo NO esta en un extremo del dominio: es un compromiso
            # entre autoridad de control y acoplamiento con los eslabones.
            @test 1 < i < length(Ms)
            @test 0.05 < Ms[i] < 0.30

            # c* y el margen PBH apuntan en direcciones OPUESTAS: no hay una
            # unica "mejor" configuracion, depende de que se optimice.
            @test issorted([r.c_star for r in res])
        end

        @testset "geometria derivada al barrer longitudes" begin
            # Al alargar una barra uniforme deben recalcularse el centro de masa
            # y la inercia. Olvidarlo da una barra que cambia de largo sin
            # cambiar de inercia, que no es ningun sistema fisico.
            @test infer_geometry(p3) == :uniform
            q = set_param(p3, (:l, 3), 1.0)
            @test q.a[3] ≈ 0.5
            @test q.Il[3] ≈ (1/12) * q.m[3] * 1.0^2

            pm = point_masses(1.0, [0.3, 0.3], [0.5, 0.5])
            @test infer_geometry(pm) == :point
            qm = set_param(pm, (:l, 2), 0.8)
            @test qm.a[2] ≈ 0.8
            @test qm.Il[2] ≈ 0.0
        end

        @testset "el eslabon superior corto es el peligroso" begin
            res = sweep_1d(lin, p3, (:l, 3), [0.02, 0.1, 1/3, 1.0])
            lambdas = [r.lambda_max for r in res]
            # Alargar el eslabon superior FACILITA el control: es el efecto de
            # la escoba larga en la palma de la mano, ahora cuantificado.
            @test issorted(lambdas, rev=true)
            @test lambdas[1] > 40.0
            @test issorted([r.pbh_normalized for r in res])
        end

        @testset "cota elipsoidal dentro de la frontera de saturacion" begin
            # c* certifica NO SATURACION mas convergencia del sistema LINEAL.
            # La comparacion honesta es contra la frontera sin restriccion de
            # riel: en el simple el riel es la restriccion activa y limita muy
            # por debajo de lo que la saturacion permitiria.
            for (nombre, p, _) in configs
                ss = lin(p)
                n = size(ss.A, 1)
                lqr = design_lqr(ss.A, ss.B, default_weights(n), R_STD)
                cs = nonsaturation_level(lqr.K, lqr.P, R_STD, ss.B, 50.0)
                cota = max_axis_angle(lqr.P, cs, 3)
                th_sat = max_recoverable_angle(p, lqr.K; umax=50.0, x_rail=1e6)
                @test cota <= th_sat
            end

            # En el simple el riel muerde antes que el motor; en el triple no.
            p1 = configs[1][2]
            ss1 = lin(p1)
            K1 = design_lqr(ss1.A, ss1.B, default_weights(4), R_STD).K
            @test max_recoverable_angle(p1, K1; umax=50.0, x_rail=1.5) <
                  max_recoverable_angle(p1, K1; umax=50.0, x_rail=1e6) - 0.05

            ss3 = lin(p3)
            K3 = design_lqr(ss3.A, ss3.B, default_weights(8), R_STD).K
            @test max_recoverable_angle(p3, K3; umax=50.0, x_rail=1.5) ≈
                  max_recoverable_angle(p3, K3; umax=50.0, x_rail=1e6) atol = 1e-3
        end

        @testset "barrido 2D" begin
            # sweep_2d debe coincidir punto a punto con viability sobre los
            # mismos parametros: es el mismo calculo, solo organizado en malla.
            l3s = [0.2, 1/3]
            m3s = [0.1, 0.2]
            malla = sweep_2d(lin, p3, (:l, 3), l3s, (:m, 3), m3s)
            @test size(malla) == (length(l3s), length(m3s))

            for i in eachindex(l3s), j in eachindex(m3s)
                @test malla[i, j].value1 == l3s[i]
                @test malla[i, j].value2 == m3s[j]
                q = set_param(set_param(p3, (:l, 3), l3s[i]), (:m, 3), m3s[j])
                @test malla[i, j].pbh_normalized ≈ viability(lin, q).pbh_normalized
            end
        end

        @testset "actuador minimo" begin
            # Se mide sobre el simple, que es el mas barato de integrar, y sin
            # restriccion de riel para aislar el efecto del actuador: es lo que
            # min_actuator_force pretende medir.
            p1 = configs[1][2]
            ss1 = lin(p1)
            K1 = design_lqr(ss1.A, ss1.B, default_weights(4), R_STD).K

            umin = min_actuator_force(p1, K1; theta0=0.2, x_rail=1e6)
            @test isfinite(umin)

            # La biseccion devuelve una frontera real: con ese limite recupera
            # y con un 20% menos ya no. Sin cambio de signo, el valor no
            # significaria nada.
            @test recovers(p1, K1, 0.2; umax=umin, x_rail=1e6)
            @test !recovers(p1, K1, 0.2; umax=0.8 * umin, x_rail=1e6)

            # Cuanto mas inclinado se arranca, mas motor hace falta.
            @test min_actuator_force(p1, K1; theta0=0.35, x_rail=1e6) > umin
        end
    end

    # =======================================================================
    @testset "8. Observador y dualidad" begin
        casos = [
            ("simple", linearize_system(default_params())),
            ("doble", linearize_system_double(default_params_double())),
            ("triple", linearize_system_triple(default_params_triple())),
        ]

        for (nombre, ss) in casos
            @testset "$nombre" begin
                A, B, C = ss.A, ss.B, ss.C
                n, p = size(A, 1), size(C, 1)
                Qo, Ro = 100.0 * Matrix(I, n, n), Matrix(I, p, p)

                K = design_lqr(A, B, default_weights(n), R_STD).K
                obs = design_observer(A, C, Qo, Ro)

                # Dualidad exacta: L es la K del sistema traspuesto.
                K_dual = design_lqr(Matrix(A'), Matrix(C'), Qo, Ro).K
                @test obs.L ≈ Matrix(K_dual') atol = 1e-12

                # El observador es estable: el error de estimacion decae.
                @test all(real.(obs.eigenvalues_obs) .< 0)
                @test size(obs.L) == (n, p)

                # PRINCIPIO DE SEPARACION: el espectro del sistema aumentado es
                # la union de los dos espectros, porque la matriz aumentada es
                # triangular por bloques.
                sep = check_separation(A, B, C, K, obs.L)
                @test sep.holds
                @test sep.max_error < 1e-8

                # La estructura de la matriz aumentada es la esperada.
                M = augmented_closed_loop(A, B, C, K, obs.L)
                @test M[1:n, 1:n] ≈ A - B * K
                @test M[1:n, n+1:2n] ≈ B * K
                @test all(M[n+1:2n, 1:n] .== 0)     # triangular por bloques
                @test M[n+1:2n, n+1:2n] ≈ A - obs.L * C
            end
        end

        @testset "subconjuntos de sensores del triple" begin
            ss = linearize_system_triple(default_params_triple())
            tabla = sensor_subset_analysis(ss.A, ss.C,
                                           ["pos", "th1", "th2", "th3"])
            @test length(tabla) == 15          # 2^4 - 1 subconjuntos no vacios

            observables = filter(f -> f.observable, tabla)
            no_obs = filter(f -> !f.observable, tabla)

            # Todo subconjunto que incluya la posicion del carro es observable;
            # ninguno que la excluya lo es. Midiendo solo angulos, la posicion
            # absoluta del carro es inobservable: la dinamica angular es
            # invariante a trasladar el carro. De ahi el rango 7/8 exacto.
            @test all(f -> 1 in f.indices, observables)
            @test all(f -> !(1 in f.indices), no_obs)
            @test all(f -> f.rank == 7, no_obs)
            @test length(observables) == 8     # los 2^3 que contienen "pos"

            # El juego minimo observable es un unico sensor: la posicion.
            minimo = observables[argmin(length(f.indices) for f in observables)]
            @test minimo.indices == [1]

            # Pero su margen es dos ordenes peor que el del juego completo:
            # observable no es lo mismo que practicable.
            completo = tabla[1]
            @test completo.indices == [1, 2, 3, 4]
            @test minimo.margin_normalized < completo.margin_normalized / 50
        end

        @testset "el observador encoge la region de atraccion" begin
            # El principio de separacion garantiza estabilidad del sistema
            # LINEAL aumentado. No dice nada sobre la region de atraccion del
            # NO LINEAL, y ahi estimar en vez de medir se paga caro: durante el
            # transitorio inicial el control actua sobre una estimacion
            # equivocada y saca a la planta del regimen donde la linealizacion
            # vale.
            p3 = default_params_triple()
            ss = linearize_system_triple(p3)
            K = design_lqr(ss.A, ss.B, default_weights(8), R_STD).K
            obs = design_observer(ss.A, ss.C, 100.0 * Matrix(I, 8, 8),
                                  Matrix(I, 4, 4))

            converge(sol) = all(abs(sol.u[end][j]) < 0.02 for j in 3:2:7)

            # Condicion inicial nominal del informe: los tres angulos a 0.1 rad.
            x0 = zeros(8)
            for j in 1:3
                x0[2j+1] = 0.10
            end

            sol_full = solve(ODEProblem(closed_loop_eom_nlink!, copy(x0),
                                        (0.0, 10.0),
                                        (params=p3, K=K, saturate=150.0)),
                             Tsit5(), saveat=0.01, verbose=false)
            @test converge(sol_full)          # con estado completo, converge

            sol_obs = solve(ODEProblem(observer_eom!, vcat(x0, zeros(8)),
                                       (0.0, 10.0),
                                       (params=p3, eom! = nonlinear_eom_triple!,
                                        A=ss.A, B=ss.B, C=ss.C, K=K, L=obs.L,
                                        saturate=150.0)),
                            Tsit5(), saveat=0.01, verbose=false)
            @test !converge(sol_obs)          # con observador, NO converge

            # Desde una desviacion mas pequena si converge, y el error de
            # estimacion decae como predice A - LC de Hurwitz.
            z0 = zeros(16)
            z0[3] = 0.10
            sol_ok = solve(ODEProblem(observer_eom!, z0, (0.0, 10.0),
                                      (params=p3, eom! = nonlinear_eom_triple!,
                                       A=ss.A, B=ss.B, C=ss.C, K=K, L=obs.L,
                                       saturate=150.0)),
                           Tsit5(), saveat=0.01, verbose=false)
            @test converge(sol_ok)
            @test norm(sol_ok.u[end][1:8] .- sol_ok.u[end][9:16]) < 1e-4
        end

        @testset "Ackermann dual solo con una salida" begin
            ss = linearize_system(default_params())
            # Con una sola salida (la posicion) el dual de Ackermann aplica y
            # debe reproducir exactamente los polos pedidos.
            C1 = reshape(ss.C[1, :], 1, 4)
            polos = [-6.0, -7.0, -8.0, -9.0]
            res = design_observer_poles(ss.A, C1, polos)
            @test sort(real.(res.eigenvalues_obs)) ≈ sort(polos) atol = 1e-6

            # Con varias salidas debe rechazar, no dar un resultado silencioso.
            @test_throws ErrorException design_observer_poles(ss.A, ss.C, polos)
        end
    end
end
