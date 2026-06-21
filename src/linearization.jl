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
