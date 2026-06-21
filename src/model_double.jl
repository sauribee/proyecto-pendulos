# =============================================================================
# model_double.jl -- Modelo fisico del pendulo invertido DOBLE sobre carro
# =============================================================================
# Configuracion II del informe: un carro con dos eslabones en serie. Las masas
# son puntuales (m1 en la articulacion intermedia, m2 en el extremo) y las
# barras se consideran de masa despreciable. Sin friccion.
#
# Variables de estado (dimension 6): x = [pos, vel, theta1, omega1, theta2, omega2]
#   pos    -> posicion horizontal del carro
#   vel    -> velocidad del carro
#   theta1 -> angulo del eslabon inferior desde la vertical superior
#   omega1 -> velocidad angular del eslabon inferior
#   theta2 -> angulo del eslabon superior desde la vertical superior
#   omega2 -> velocidad angular del eslabon superior
#
# Las ecuaciones se obtienen del formalismo de Euler-Lagrange y se escriben en
# la forma matricial estandar de la mecanica:
#     M(q) qdd + C(q, qd) qd + G(q) = F u,   con q = [pos, theta1, theta2]
# =============================================================================

module ModelDouble

using LinearAlgebra

export SystemParamsDouble, default_params_double
export nonlinear_eom_double!, closed_loop_eom_double!, state_derivative_double

"""
Parametros fisicos del sistema carro-pendulo doble.
Los valores por defecto son los del informe (Configuracion II).
"""
Base.@kwdef struct SystemParamsDouble
    M::Float64 = 1.0    # masa del carro [kg]
    m1::Float64 = 0.3   # masa en la articulacion intermedia [kg]
    m2::Float64 = 0.3   # masa en el extremo superior [kg]
    L1::Float64 = 0.5   # longitud del eslabon inferior [m]
    L2::Float64 = 0.5   # longitud del eslabon superior [m]
    g::Float64 = 9.81   # gravedad [m/s^2]
end

"""
Parametros por defecto del informe: M = 1 kg, m1 = m2 = 0.3 kg,
L1 = L2 = 0.5 m, sin friccion.
"""
default_params_double() = SystemParamsDouble()

"""
    nonlinear_eom_double!(dx, x, p, t)

Ecuaciones de movimiento no lineales del pendulo doble para
DifferentialEquations.jl. La fuerza de control se pasa como p.F (escalar)
o mediante p.controller(x, t).

Se arma la matriz de masa M(q) (simetrica y definida positiva) y el lado
derecho con los terminos centrifugos y gravitatorios; luego se resuelve el
sistema lineal M(q) qdd = rhs para obtener las aceleraciones.
"""
function nonlinear_eom_double!(dx, x, p, t)
    # Extraer estado
    pos, vel, theta1, omega1, theta2, omega2 = x

    # Parametros del sistema
    sp = p.params
    M, m1, m2, L1, L2, g = sp.M, sp.m1, sp.m2, sp.L1, sp.L2, sp.g

    # Fuerza de control
    F = haskey(p, :controller) ? p.controller(x, t) : p.F

    # Abreviaciones trigonometricas
    c1 = cos(theta1)
    s1 = sin(theta1)
    c2 = cos(theta2)
    s2 = sin(theta2)
    c12 = cos(theta1 - theta2)
    s12 = sin(theta1 - theta2)

    # Matriz de masa M(q), con q = [pos, theta1, theta2]
    Mq = [ M + m1 + m2        (m1 + m2) * L1 * c1     m2 * L2 * c2;
           (m1 + m2) * L1 * c1  (m1 + m2) * L1^2      m2 * L1 * L2 * c12;
           m2 * L2 * c2         m2 * L1 * L2 * c12     m2 * L2^2 ]

    # Lado derecho: fuerza de control + terminos centrifugos + gravedad.
    # Provienen de pasar C(q,qd) qd y G(q) al otro lado de la ecuacion.
    rhs = [ F + (m1 + m2) * L1 * s1 * omega1^2 + m2 * L2 * s2 * omega2^2,
            (m1 + m2) * g * L1 * s1 - m2 * L1 * L2 * s12 * omega2^2,
            m2 * g * L2 * s2 + m2 * L1 * L2 * s12 * omega1^2 ]

    # Resolver M(q) qdd = rhs para las aceleraciones generalizadas
    qdd = Mq \ rhs

    # Vector de derivadas del estado
    dx[1] = vel       # d(pos)/dt
    dx[2] = qdd[1]    # d(vel)/dt    = aceleracion del carro
    dx[3] = omega1    # d(theta1)/dt
    dx[4] = qdd[2]    # d(omega1)/dt = aceleracion angular eslabon 1
    dx[5] = omega2    # d(theta2)/dt
    dx[6] = qdd[3]    # d(omega2)/dt = aceleracion angular eslabon 2
end

"""
    closed_loop_eom_double!(dx, x, p, t)

Ecuaciones en lazo cerrado con la ley de control u = -K x. Reutiliza el modelo
no lineal pasando la fuerza ya calculada (asi no se duplica la fisica).

Parametros esperados en p:
    p.params    -- SystemParamsDouble del modelo
    p.K         -- Matriz de ganancia (1x6)
    p.saturate  -- (opcional) limite de fuerza [N]
"""
function closed_loop_eom_double!(dx, x, p, t)
    K = p.K
    u = -dot(K[1, :], x)  # escalar para SISO

    if haskey(p, :saturate)
        u = clamp(u, -p.saturate, p.saturate)
    end

    nonlinear_eom_double!(dx, x, (params=p.params, F=u), t)
end

"""
Version funcional (no in-place) para uso rapido.
"""
function state_derivative_double(x, params::SystemParamsDouble, F::Float64)
    dx = zeros(6)
    nonlinear_eom_double!(dx, x, (params=params, F=F), 0.0)
    return dx
end

end # module ModelDouble
