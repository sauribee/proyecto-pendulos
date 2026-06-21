# =============================================================================
# model.jl -- Modelo fisico del pendulo invertido sobre carro
# =============================================================================
# Define parametros del sistema y ecuaciones de movimiento no lineales
# derivadas del formalismo de Euler-Lagrange.
#
# Variables de estado: x = [pos, vel, theta, omega]
#   pos   -> posicion horizontal del carro
#   vel   -> velocidad del carro
#   theta -> angulo del pendulo medido desde la vertical superior (arriba = 0)
#   omega -> velocidad angular del pendulo
# =============================================================================

module Model

using LinearAlgebra

export SystemParams, default_params, nonlinear_eom!, state_derivative

"""
Parametros fisicos del sistema carro-pendulo.
"""
Base.@kwdef struct SystemParams
    M::Float64 = 1.0       # masa del carro [kg]
    m::Float64 = 0.3       # masa del pendulo [kg]
    L::Float64 = 0.5       # longitud al centro de masa [m]
    g::Float64 = 9.81      # gravedad [m/s^2]
    b::Float64 = 0.1       # coeficiente de friccion del carro [N s/m]
    I::Float64 = 0.0       # momento de inercia del pendulo respecto a su CM [kg m^2]
                            # Para masa puntual: I = 0
                            # Para barra uniforme: I = (1/12) m Lbar^2 con Lbar = 2L
end

"""
Parametros por defecto: barra uniforme de 1 m, carro de 1 kg, barra de 0.3 kg.
"""
function default_params()
    m = 0.3
    Lbar = 1.0               # longitud total de la barra
    L = Lbar / 2             # distancia al centro de masa
    Ip = (1/12) * m * Lbar^2 # momento de inercia barra uniforme
    return SystemParams(M=1.0, m=m, L=L, g=9.81, b=0.1, I=Ip)
end

"""
    nonlinear_eom!(dx, x, p, t)

Ecuaciones de movimiento no lineales en forma estandar para DifferentialEquations.jl.
Derivadas del Lagrangiano del sistema carro-pendulo.

El estado es x = [pos_carro, vel_carro, angulo, vel_angular].
La fuerza de control se pasa como p.F (escalar) o mediante p.controller(x, t).

Convencion: theta medido desde la vertical superior (theta = 0 es el equilibrio
erguido inestable). Con esta convencion el termino gravitatorio +(M+m) m g L sin(theta)
en la aceleracion angular produce un eigenvalor real positivo: el pendulo cae.
"""
function nonlinear_eom!(dx, x, p, t)
    # Extraer estado
    pos, vel, theta, omega = x

    # Parametros del sistema
    sp = p.params
    M, m, L, g, b, Ip = sp.M, sp.m, sp.L, sp.g, sp.b, sp.I

    # Fuerza de control
    F = haskey(p, :controller) ? p.controller(x, t) : p.F

    # Abreviaciones trigonometricas
    sin_t = sin(theta)
    cos_t = cos(theta)

    # Sistema 2x2 de la mecanica:
    #   (M + m) ddx + m L cos(theta) ddtheta = F - b vel + m L omega^2 sin(theta)
    #   m L cos(theta) ddx + (I + m L^2) ddtheta = m g L sin(theta)
    #
    # Matriz de masa Mm = [ M+m            m L cos(theta) ;
    #                       m L cos(theta) I + m L^2      ]
    # con determinante D = (M+m)(I+m L^2) - (m L cos(theta))^2.
    D = (M + m) * (Ip + m * L^2) - (m * L * cos_t)^2

    # Lados derechos
    rhs1 = F - b * vel + m * L * omega^2 * sin_t   # ecuacion del carro
    rhs2 = m * g * L * sin_t                        # ecuacion del pendulo

    # Resolver [ddx; ddtheta] = Mm^{-1} [rhs1; rhs2]:
    #   Mm^{-1} = (1/D) [ I+m L^2         -m L cos(theta) ;
    #                     -m L cos(theta)  M+m            ]
    x_ddot = ((Ip + m * L^2) * rhs1 - m * L * cos_t * rhs2) / D
    theta_ddot = ((M + m) * rhs2 - m * L * cos_t * rhs1) / D

    # Vector de derivadas del estado
    dx[1] = vel         # d(pos)/dt   = velocidad
    dx[2] = x_ddot      # d(vel)/dt   = aceleracion del carro
    dx[3] = omega       # d(theta)/dt = velocidad angular
    dx[4] = theta_ddot  # d(omega)/dt = aceleracion angular
end

"""
Version funcional (no in-place) para uso rapido.
"""
function state_derivative(x, params::SystemParams, F::Float64)
    dx = zeros(4)
    p = (params=params, F=F)
    nonlinear_eom!(dx, x, p, 0.0)
    return dx
end

end # module Model
