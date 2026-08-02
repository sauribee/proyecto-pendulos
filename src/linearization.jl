# =============================================================================
# linearization.jl -- Linealizacion y analisis en espacio de estados
# =============================================================================
# Calcula el Jacobiano del sistema no lineal alrededor del punto de equilibrio
# inestable (theta = 0, pendulo arriba), extrae las matrices A, B, C, D del
# espacio de estados, y realiza analisis de estabilidad via eigenvalores.
# =============================================================================

module Linearization

using LinearAlgebra
using Printf

export linearize_system, linearize_system_double, StateSpaceModel, print_analysis
export linearize_system_nlink, linearize_system_triple
export controllability_matrix, observability_matrix
export check_controllability, check_observability

"""
Modelo linealizado en espacio de estados:
    dx = A x + B u
    y  = C x + D u
"""
struct StateSpaceModel
    A::Matrix{Float64}    # Matriz de estado (n x n)
    B::Matrix{Float64}    # Matriz de entrada (n x m)
    C::Matrix{Float64}    # Matriz de salida (p x n)
    D::Matrix{Float64}    # Matriz de transmision directa (p x m)
    eigenvalues::Vector{ComplexF64}
    eigenvectors::Matrix{ComplexF64}
    state_names::Vector{String}
    output_names::Vector{String}
end

"""
    linearize_system(params) -> StateSpaceModel

Linealiza las ecuaciones de movimiento alrededor del punto de equilibrio
superior (theta = 0, omega = 0, pos = 0, vel = 0, F = 0).

El Jacobiano se calcula analiticamente a partir de las EOM.
En el punto de equilibrio: sin(theta) -> theta, cos(theta) -> 1, omega^2 -> 0.
"""
function linearize_system(params)
    M = params.M
    m = params.m
    L = params.L
    g = params.g
    b = params.b
    Ip = params.I  # momento de inercia del pendulo

    # Determinante de la matriz de masa evaluado en theta = 0 (cos(0) = 1):
    # D0 = (M + m)(I + m L^2) - (m L)^2
    D0 = (M + m) * (Ip + m * L^2) - (m * L)^2

    # ---------------------------------------------------------------
    # Matriz A (4x4): Jacobiano df/dx evaluado en el equilibrio
    # ---------------------------------------------------------------
    # Estado: [pos, vel, theta, omega]
    #
    # Las ecuaciones linealizadas son:
    #   ddx     = (-b(I+m L^2) vel - m^2 g L^2 theta) / D0 + ((I+m L^2)/D0) F
    #   ddtheta = ( b m L vel + (M+m) m g L theta) / D0    + (-m L/D0) F
    #
    # Convencion: theta desde la vertical superior. El signo POSITIVO de
    # A[4,3] = +(M+m) m g L / D0 es la firma de la inestabilidad: produce
    # un eigenvalor real positivo (el pendulo invertido cae).

    A = zeros(4, 4)
    A[1, 2] = 1.0                                       # d(pos)/dt = vel
    A[2, 2] = -b * (Ip + m * L^2) / D0                 # d(ddx)/d(vel)
    A[2, 3] = -m^2 * g * L^2 / D0                      # d(ddx)/d(theta)
    A[3, 4] = 1.0                                       # d(theta)/dt = omega
    A[4, 2] = b * m * L / D0                            # d(ddtheta)/d(vel)
    A[4, 3] = (M + m) * m * g * L / D0                 # d(ddtheta)/d(theta)

    # ---------------------------------------------------------------
    # Matriz B (4x1): Jacobiano df/du
    # ---------------------------------------------------------------
    B = zeros(4, 1)
    B[2, 1] = (Ip + m * L^2) / D0                      # d(ddx)/dF
    B[4, 1] = -m * L / D0                              # d(ddtheta)/dF

    # ---------------------------------------------------------------
    # Matriz C (2x4): salidas medibles (posicion del carro y angulo)
    # ---------------------------------------------------------------
    C = zeros(2, 4)
    C[1, 1] = 1.0    # medimos pos
    C[2, 3] = 1.0    # medimos theta

    # ---------------------------------------------------------------
    # Matriz D (2x1): transmision directa (cero para este sistema)
    # ---------------------------------------------------------------
    D_mat = zeros(2, 1)

    # ---------------------------------------------------------------
    # Analisis espectral: eigenvalores y eigenvectores de A
    # ---------------------------------------------------------------
    eig_result = eigen(A)
    lambda = eig_result.values
    V = eig_result.vectors

    state_names = ["pos (posicion)", "vel (velocidad)", "theta (angulo)", "omega (vel. angular)"]
    output_names = ["pos (posicion)", "theta (angulo)"]

    return StateSpaceModel(A, B, C, D_mat, lambda, V, state_names, output_names)
end

"""
    linearize_system_double(params) -> StateSpaceModel

Linealiza el pendulo invertido DOBLE alrededor del equilibrio superior
(theta1 = theta2 = 0, todo en reposo). El estado es de dimension 6:
    x = [pos, vel, theta1, omega1, theta2, omega2]

Las matrices A y B son la forma analitica del Jacobiano (sistema sin friccion,
masas puntuales m1 en la articulacion intermedia y m2 en el extremo). Coinciden
con la linealizacion numerica de ModelDouble.nonlinear_eom_double! y reproducen
el espectro de lazo abierto {+8.57, +4.09, 0, 0, -4.09, -8.57}: dos modos
inestables. Devuelve un StateSpaceModel, por lo que reutiliza sin cambios las
funciones de analisis (controlabilidad, observabilidad, print_analysis).
"""
function linearize_system_double(params)
    M = params.M
    m1 = params.m1
    m2 = params.m2
    L1 = params.L1
    L2 = params.L2
    g = params.g

    # ---------------------------------------------------------------
    # Matriz A (6x6): Jacobiano analitico df/dx en el equilibrio
    # ---------------------------------------------------------------
    # Estado: [pos, vel, theta1, omega1, theta2, omega2]
    A = zeros(6, 6)
    A[1, 2] = 1.0                                       # d(pos)/dt = vel
    A[2, 3] = -g * (m1 + m2) / M                        # d(ddx)/d(theta1)
    A[3, 4] = 1.0                                       # d(theta1)/dt = omega1
    A[4, 3] = g * (M + m1) * (m1 + m2) / (M * L1 * m1)  # d(ddtheta1)/d(theta1)
    A[4, 5] = -g * m2 / (L1 * m1)                       # d(ddtheta1)/d(theta2)
    A[5, 6] = 1.0                                       # d(theta2)/dt = omega2
    A[6, 3] = -g * (m1 + m2) / (L2 * m1)                # d(ddtheta2)/d(theta1)
    A[6, 5] = g * (m1 + m2) / (L2 * m1)                 # d(ddtheta2)/d(theta2)

    # ---------------------------------------------------------------
    # Matriz B (6x1): Jacobiano df/du
    # ---------------------------------------------------------------
    B = zeros(6, 1)
    B[2, 1] = 1 / M
    B[4, 1] = -1 / (M * L1)

    # ---------------------------------------------------------------
    # Matriz C (3x6): medimos posicion del carro y los dos angulos
    # ---------------------------------------------------------------
    C = zeros(3, 6)
    C[1, 1] = 1.0   # pos
    C[2, 3] = 1.0   # theta1
    C[3, 5] = 1.0   # theta2

    # Matriz D (3x1): sin transmision directa
    D_mat = zeros(3, 1)

    # Analisis espectral
    eig_result = eigen(A)
    lambda = eig_result.values
    V = eig_result.vectors

    state_names = ["pos (posicion)", "vel (velocidad)",
                   "theta1 (angulo 1)", "omega1 (vel. angular 1)",
                   "theta2 (angulo 2)", "omega2 (vel. angular 2)"]
    output_names = ["pos (posicion)", "theta1 (angulo 1)", "theta2 (angulo 2)"]

    return StateSpaceModel(A, B, C, D_mat, lambda, V, state_names, output_names)
end

"""
    linearize_system_nlink(params) -> StateSpaceModel

Linealiza el pendulo de N eslabones alrededor del equilibrio superior
(todos los theta_j = 0, todo en reposo). El estado tiene dimension n = 2(N+1),
en el orden intercalado del proyecto:

    x = [pos, vel, theta_1, omega_1, ..., theta_N, omega_N]

EL TRUCO DECISIVO. Evaluando en theta_j = 0 y qd = 0 se tiene cos -> 1,
sin(theta_j) -> theta_j, cos(theta_j - theta_k) -> 1, y los terminos
cuadraticos en las velocidades desaparecen por ser de segundo orden. Lo que
queda es un sistema lineal con MATRIZ DE MASA CONSTANTE:

    M0 qdd = G0 q - F0 qd + e1 u

con M0 = M(q) evaluada en theta = 0, G0 = g diag(0, beta_1, ..., beta_N) y F0
la matriz de friccion. Despejando qdd se obtienen los tres bloques

    Acc = M0 \\ G0        d(qdd)/dq
    Acv = -(M0 \\ F0)     d(qdd)/dqd
    Bc  = M0 \\ e1        d(qdd)/du

que se ensamblan directamente en orden intercalado. En orden por bloques la
matriz seria [0 I; Acc Acv], que es esta misma conjugada por una permutacion:
el espectro no cambia, pero se ensambla la version intercalada para ser
coherente con las otras dos configuraciones.

PREDICCION ESTRUCTURAL DEL ESPECTRO. Sin friccion, buscar soluciones de la
forma q(t) = v exp(lambda t) da el problema generalizado G0 v = lambda^2 M0 v:
los eigenvalores de A son las RAICES CUADRADAS de los de M0^-1 G0 y por tanto
vienen en pares +-lambda. Como G0 tiene una fila y una columna nulas (el carro
no tiene fuerza recuperadora), M0^-1 G0 es singular y aparece un cero doble en
un bloque de Jordan 2x2:

    spec(A) = {+-lambda_1, ..., +-lambda_N, 0, 0},   N modos inestables

La friccion del carro rompe ligeramente esa simetria y desdobla el cero.

El jacobiano es ANALITICO, no numerico. La coincidencia con la diferenciacion
automatica de las EOM no lineales se comprueba en test/runtests.jl.

`params` debe exponer los campos M, m, l, a, Il, g, b y c, es decir la
estructura de ModelNLink.SystemParamsNLink.
"""
function linearize_system_nlink(params)
    N = length(params.m)
    d = N + 1        # grados de libertad: carro + N eslabones
    n = 2 * d        # dimension del estado

    g = params.g

    # ---------------------------------------------------------------
    # Agrupaciones de parametros (informe tecnico de la 2.a entrega, seccion 2.2)
    #   beta_j = momento estatico del eslabon j y de lo que carga encima
    #   J_j    = inercia efectiva respecto a su articulacion
    #   Gam_jk = acoplamiento inercial entre los eslabones j y k
    # ---------------------------------------------------------------
    beta = zeros(N)
    Jeff = zeros(N)
    for j in 1:N
        m_above = sum(params.m[j+1:end])
        beta[j] = params.m[j] * params.a[j] + params.l[j] * m_above
        Jeff[j] = params.Il[j] + params.m[j] * params.a[j]^2 +
                  params.l[j]^2 * m_above
    end

    # ---------------------------------------------------------------
    # M0 = M(q) evaluada en theta = 0 (todos los cosenos valen 1)
    # ---------------------------------------------------------------
    M0 = zeros(d, d)
    M0[1, 1] = params.M + sum(params.m)
    for j in 1:N
        M0[1, j+1] = beta[j]
        M0[j+1, 1] = beta[j]
        M0[j+1, j+1] = Jeff[j]
        for k in 1:N
            j == k && continue
            M0[j+1, k+1] = params.l[min(j, k)] * beta[max(j, k)]
        end
    end

    # ---------------------------------------------------------------
    # G0: gravedad linealizada. El signo POSITIVO de g beta_j es la firma
    # de la inestabilidad, igual que el +(M+m) m g L / D0 de la Config. I.
    # ---------------------------------------------------------------
    G0 = zeros(d, d)
    for j in 1:N
        G0[j+1, j+1] = g * beta[j]
    end

    # ---------------------------------------------------------------
    # F0: friccion. El carro aporta b; las articulaciones, si estan
    # activas, aportan un bloque tridiagonal simetrico porque el par es
    # RELATIVO entre eslabones contiguos.
    # ---------------------------------------------------------------
    F0 = zeros(d, d)
    F0[1, 1] = params.b
    if !isempty(params.c)
        c = params.c
        for j in 1:N
            c_above = j < N ? c[j+1] : 0.0
            F0[j+1, j+1] = c[j] + c_above
            if j > 1
                F0[j+1, j] = -c[j]
            end
            if j < N
                F0[j+1, j+2] = -c_above
            end
        end
    end

    e1 = zeros(d)
    e1[1] = 1.0

    Acc = M0 \ G0
    Acv = -(M0 \ F0)
    Bc = M0 \ e1

    # ---------------------------------------------------------------
    # Ensamblaje en orden INTERCALADO
    # ---------------------------------------------------------------
    A = zeros(n, n)
    B = zeros(n, 1)
    for i in 1:d
        A[2i-1, 2i] = 1.0              # d(q_i)/dt = qd_i
        for k in 1:d
            A[2i, 2k-1] = Acc[i, k]
            A[2i, 2k] = Acv[i, k]
        end
        B[2i, 1] = Bc[i]
    end

    # ---------------------------------------------------------------
    # C: se miden la posicion del carro y todos los angulos
    # ---------------------------------------------------------------
    C = zeros(d, n)
    C[1, 1] = 1.0
    for j in 1:N
        C[j+1, 2j+1] = 1.0
    end
    D_mat = zeros(d, 1)

    eig_result = eigen(A)

    state_names = ["pos (posicion)", "vel (velocidad)"]
    output_names = ["pos (posicion)"]
    for j in 1:N
        push!(state_names, "theta$j (angulo $j)")
        push!(state_names, "omega$j (vel. angular $j)")
        push!(output_names, "theta$j (angulo $j)")
    end

    return StateSpaceModel(A, B, C, D_mat, eig_result.values, eig_result.vectors,
                           state_names, output_names)
end

"""
    linearize_system_triple(params) -> StateSpaceModel

Linealiza el pendulo invertido TRIPLE (Configuracion III) alrededor del
equilibrio superior. El estado es de dimension 8:

    x = [pos, vel, theta1, omega1, theta2, omega2, theta3, omega3]

Es `linearize_system_nlink` restringida a N = 3. Con los parametros por defecto
reproduce el espectro

    {+18.0613, +9.7740, +4.3777, 0, -0.0625, -4.3990, -9.7813, -18.0647}

es decir TRES modos inestables. El dominante implica una constante de tiempo de
caida de 1/18.06 = 55 ms: el sistema se desploma cuatro veces mas rapido que la
Configuracion I. La matriz C mide pos, theta1, theta2 y theta3 (4x8).
"""
function linearize_system_triple(params)
    length(params.m) == 3 ||
        error("linearize_system_triple espera 3 eslabones, recibio " *
              string(length(params.m)))
    return linearize_system_nlink(params)
end

"""
Matriz de controlabilidad de Kalman: C = [B  AB  A^2 B  A^3 B]
"""
function controllability_matrix(ss::StateSpaceModel)
    n = size(ss.A, 1)
    C_ctrl = hcat([ss.A^k * ss.B for k in 0:n-1]...)
    return C_ctrl
end

"""
Matriz de observabilidad de Kalman: O = [C; CA; CA^2; CA^3]
"""
function observability_matrix(ss::StateSpaceModel)
    n = size(ss.A, 1)
    O_obs = vcat([ss.C * ss.A^k for k in 0:n-1]...)
    return O_obs
end

"""
Verifica controlabilidad: rank(C) == n
"""
function check_controllability(ss::StateSpaceModel)
    C_ctrl = controllability_matrix(ss)
    n = size(ss.A, 1)
    r = rank(C_ctrl)
    return (is_controllable=r == n, rank=r, required_rank=n, matrix=C_ctrl)
end

"""
Verifica observabilidad: rank(O) == n
"""
function check_observability(ss::StateSpaceModel)
    O_obs = observability_matrix(ss)
    n = size(ss.A, 1)
    r = rank(O_obs)
    return (is_observable=r == n, rank=r, required_rank=n, matrix=O_obs)
end

"""
Imprime un resumen completo del analisis del sistema linealizado.
"""
function print_analysis(ss::StateSpaceModel)
    println("=" ^ 60)
    println("  ANALISIS DEL SISTEMA LINEALIZADO")
    println("=" ^ 60)

    println("\n  Matriz A (dinamica del sistema):")
    display(round.(ss.A, digits=4))

    println("\n\n  Matriz B (entrada de control):")
    display(round.(ss.B, digits=4))

    println("\n\n  Matriz C (salidas medibles):")
    display(round.(ss.C, digits=4))

    # Eigenvalores
    println("\n\n  Eigenvalores de A:")
    for (i, lambda) in enumerate(ss.eigenvalues)
        real_part = real(lambda)
        imag_part = imag(lambda)
        stability = real_part > 0 ? "INESTABLE" : (real_part < 0 ? "ESTABLE" : "MARGINAL")

        if abs(imag_part) < 1e-10
            @printf("  lambda_%d = %+.4f  [%s]\n", i, real_part, stability)
        else
            @printf("  lambda_%d = %+.4f %+.4fi  [%s]\n", i, real_part, imag_part, stability)
        end
    end

    any_unstable = any(real.(ss.eigenvalues) .> 1e-10)
    println("\n  Sistema: ", any_unstable ? "INESTABLE (requiere control)" : "estable")

    # Controlabilidad
    ctrl = check_controllability(ss)
    println("\n  Controlabilidad:")
    @printf("  rank(C) = %d / %d -> %s\n",
            ctrl.rank, ctrl.required_rank,
            ctrl.is_controllable ? "CONTROLABLE" : "NO CONTROLABLE")

    # Observabilidad
    obs = check_observability(ss)
    println("\n  Observabilidad:")
    @printf("  rank(O) = %d / %d -> %s\n",
            obs.rank, obs.required_rank,
            obs.is_observable ? "OBSERVABLE" : "NO OBSERVABLE")

    println("\n" * "=" ^ 60)
end

end # module Linearization
