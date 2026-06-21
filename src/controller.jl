# =============================================================================
# controller.jl -- Diseno de controladores para el pendulo invertido
# =============================================================================
# Implementa:
#   1. LQR (Linear Quadratic Regulator) via ecuacion algebraica de Riccati
#   2. Asignacion de polos via la formula de Ackermann
#   3. Funciones auxiliares para simulacion en lazo cerrado
# =============================================================================

module Controller

using LinearAlgebra
using Printf

export design_lqr, design_pole_placement, closed_loop_eom!
export solve_care, print_controller_summary

"""
    solve_care(A, B, Q, R) -> P

Resuelve la ecuacion algebraica continua de Riccati (CARE):
    A' P + P A - P B R^-1 B' P + Q = 0

Usa el metodo de eigendescomposicion del Hamiltoniano:
    H = [  A    -B R^-1 B' ]
        [ -Q      -A'      ]

Los eigenvectores asociados a eigenvalores con parte real negativa
dan la solucion P.
"""
function solve_care(A, B, Q, R)
    n = size(A, 1)

    # Construir la matriz Hamiltoniana (2n x 2n)
    R_inv = inv(R)
    S = B * R_inv * B'

    H = [A    -S;
         -Q   -A']

    # Eigendescomposicion del Hamiltoniano
    eig_result = eigen(H)
    lambda = eig_result.values
    V = eig_result.vectors

    # Seleccionar eigenvectores con eigenvalores estables (Re(lambda) < 0)
    stable_idx = findall(real.(lambda) .< 0)

    if length(stable_idx) != n
        error("No se encontraron exactamente n eigenvalores estables. " *
              "El sistema podria no ser estabilizable.")
    end

    # Particionar los eigenvectores estables:
    # V_stable = [V1; V2], donde V1, V2 son n x n
    V_stable = V[:, stable_idx]
    V1 = V_stable[1:n, :]
    V2 = V_stable[n+1:2n, :]

    # Solucion: P = V2 * V1^-1
    P = real.(V2 / V1)

    # Simetrizar (eliminar errores numericos)
    P = (P + P') / 2

    return P
end

"""
    design_lqr(A, B, Q, R) -> (K, P, eigenvalues_cl, A_cl)

Disena un controlador LQR (Linear Quadratic Regulator).

Minimiza el funcional de costo:
    J = integral_0^inf (x' Q x + u' R u) dt

Retorna:
    K  -- Ganancia de retroalimentacion de estado (u = -K x)
    P  -- Solucion de la ecuacion de Riccati
    eigenvalues_cl -- Eigenvalores del sistema en lazo cerrado (A - B K)
    A_cl -- Matriz de lazo cerrado
"""
function design_lqr(A, B, Q, R)
    # Resolver CARE
    P = solve_care(A, B, Q, R)

    # Ganancia optima: K = R^-1 B' P
    K = inv(R) * B' * P

    # Eigenvalores en lazo cerrado
    A_cl = A - B * K
    lambda_cl = eigvals(A_cl)

    return (K=K, P=P, eigenvalues_cl=lambda_cl, A_cl=A_cl)
end

"""
    design_pole_placement(A, B, desired_poles) -> (K, eigenvalues_cl, A_cl, desired)

Disena un controlador por asignacion de polos usando el metodo de Ackermann.

Para un sistema SISO de orden n con polos deseados p1, p2, ..., pn:
    K = en' * C^-1 * phi(A)

donde:
    C    = [B  AB  A^2 B  ...  A^(n-1) B]   (matriz de controlabilidad)
    en'  = [0 0 ... 0 1]
    phi(A) = prod(A - pi I)   (polinomio caracteristico deseado evaluado en A)
"""
function design_pole_placement(A, B, desired_poles)
    n = size(A, 1)

    # Verificar que B sea vector columna (SISO)
    if size(B, 2) != 1
        error("Pole placement por Ackermann solo aplica a sistemas SISO")
    end

    # Matriz de controlabilidad
    C_ctrl = hcat([A^k * B for k in 0:n-1]...)

    if rank(C_ctrl) < n
        error("El sistema no es controlable. Pole placement no es posible.")
    end

    # Polinomio caracteristico deseado evaluado en A:
    # phi(A) = (A - p1 I)(A - p2 I)...(A - pn I)
    # Se usa aritmetica compleja para admitir pares conjugados; el resultado
    # final es real cuando los polos vienen en pares conjugados.
    phi_A = Matrix{ComplexF64}(I, n, n)
    for p in desired_poles
        phi_A = phi_A * (A - p * I)
    end
    phi_A = real.(phi_A)

    # Formula de Ackermann: K = en' * C^-1 * phi(A)
    en = zeros(1, n)
    en[1, n] = 1.0
    K = en * inv(C_ctrl) * phi_A

    # Eigenvalores en lazo cerrado (verificacion)
    A_cl = A - B * K
    lambda_cl = eigvals(A_cl)

    return (K=K, eigenvalues_cl=lambda_cl, A_cl=A_cl, desired=desired_poles)
end

"""
    closed_loop_eom!(dx, x, p, t)

Ecuaciones de movimiento en lazo cerrado para simulacion con DifferentialEquations.jl.
Usa el modelo NO LINEAL con la ley de control u = -K x (retroalimentacion lineal).

Parametros esperados en p:
    p.params     -- SystemParams del modelo
    p.K          -- Matriz de ganancia (1x4)
    p.saturate   -- (opcional) limite de fuerza [N]

Las EOM coinciden con Model.nonlinear_eom! (misma convencion de signos).
"""
function closed_loop_eom!(dx, x, p, t)
    # Calcular fuerza de control: u = -K x
    K = p.K
    u = -dot(K[1, :], x)  # escalar para SISO

    # Saturacion de actuador (opcional)
    if haskey(p, :saturate)
        u = clamp(u, -p.saturate, p.saturate)
    end

    # EOM no lineales con la fuerza calculada
    sp = p.params
    M, m, L, g, b, Ip = sp.M, sp.m, sp.L, sp.g, sp.b, sp.I

    pos, vel, theta, omega = x
    sin_t = sin(theta)
    cos_t = cos(theta)

    D = (M + m) * (Ip + m * L^2) - (m * L * cos_t)^2

    rhs1 = u - b * vel + m * L * omega^2 * sin_t
    rhs2 = m * g * L * sin_t

    dx[1] = vel
    dx[2] = ((Ip + m * L^2) * rhs1 - m * L * cos_t * rhs2) / D
    dx[3] = omega
    dx[4] = ((M + m) * rhs2 - m * L * cos_t * rhs1) / D
end

"""
Imprime un resumen del controlador disenado.
"""
function print_controller_summary(result; method="LQR", labels=nothing)
    println("=" ^ 60)
    println("  DISENO DEL CONTROLADOR: $method")
    println("=" ^ 60)

    println("\n  Ganancia K:")
    K = result.K
    n = length(K)
    # Si no se dan etiquetas, generar K_1..K_n genericas (funciona para
    # cualquier dimension: 4 estados en el simple, 6 en el doble).
    if labels === nothing
        labels = ["x[$i]" for i in 1:n]
    end
    for (i, label) in enumerate(labels)
        @printf("  K_%-8s = %+.4f\n", label, K[i])
    end

    println("\n  Eigenvalores en lazo cerrado:")
    for (i, lambda) in enumerate(result.eigenvalues_cl)
        rp = real(lambda)
        ip = imag(lambda)
        if abs(ip) < 1e-10
            @printf("  lambda_%d = %+.4f  [estable]\n", i, rp)
        else
            @printf("  lambda_%d = %+.4f %+.4fi  [estable]\n", i, rp, ip)
        end
    end

    all_stable = all(real.(result.eigenvalues_cl) .< 0)
    println("\n  Sistema en lazo cerrado: ",
            all_stable ? "TODOS LOS POLOS ESTABLES" : "SISTEMA AUN INESTABLE")

    if hasfield(typeof(result), :P) && result.P !== nothing
        println("\n  Solucion de Riccati P (norma): ", round(norm(result.P), digits=4))
    end

    println("\n" * "=" ^ 60)
end

end # module Controller
