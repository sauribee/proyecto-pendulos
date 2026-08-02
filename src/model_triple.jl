# =============================================================================
# model_triple.jl -- Modelo fisico del pendulo invertido TRIPLE sobre carro
# =============================================================================
# Configuracion III del informe: un carro con tres eslabones articulados en
# serie, modelados como BARRAS UNIFORMES con inercia. Cuatro grados de libertad
# y un solo actuador: la subactuacion pasa de 2:1 (Config. I) a 3:1 (Config. II)
# a 4:1.
#
# Variables de estado (dimension 8), en orden intercalado:
#   x = [pos, vel, theta1, omega1, theta2, omega2, theta3, omega3]
#   pos    -> posicion horizontal del carro
#   vel    -> velocidad del carro
#   theta_j-> angulo ABSOLUTO del eslabon j desde la vertical superior
#   omega_j-> velocidad angular del eslabon j
#
# Este modulo es una capa delgada sobre ModelNLink: no reimplementa nada de la
# fisica, solo expone la API con los nombres del proyecto para que
# main_triple.jl se lea igual que main_double.jl. Toda la deduccion de
# Euler-Lagrange vive en model_nlink.jl.
#
# REQUISITO DE CARGA: model_nlink.jl debe incluirse ANTES que este archivo.
# =============================================================================

module ModelTriple

using ..ModelNLink

export SystemParamsTriple, default_params_triple
export nonlinear_eom_triple!, closed_loop_eom_triple!, state_derivative_triple

"""
Parametros del pendulo triple. Es exactamente `SystemParamsNLink` con N = 3:
no se define un tipo nuevo para que las tres configuraciones compartan sin
friccion los modulos genericos (Controller, Linearization, Metrics, Sweep).
"""
const SystemParamsTriple = SystemParamsNLink

"""
    default_params_triple() -> SystemParamsTriple

Parametros por defecto de la Configuracion III. Se eligen para que la
comparacion con la Configuracion II sea limpia: la MISMA masa total de pendulo
(0.6 kg) y la MISMA longitud total (1 m), redistribuidas en tres barras
uniformes en vez de dos masas puntuales.

    M          = 1.0 kg           igual que I y II
    m1=m2=m3   = 0.2 kg           total 0.6 kg, igual que la Config. II
    l1=l2=l3   = 1/3 m            total 1 m, igual que I y II
    a_j        = l_j/2 = 1/6 m    barra uniforme
    Il_j       = m_j l_j^2/12     barra uniforme
    g          = 9.81 m/s^2
    b          = 0.1 N s/m        igual que la Config. I
    c_j        = 0                sin friccion articular por defecto

Con ellos beta = (0.16667, 0.1, 0.03333) y J = (0.05185, 0.02963, 0.00741).
"""
function default_params_triple()
    return uniform_rods(1.0, fill(0.2, 3), fill(1/3, 3); g=9.81, b=0.1)
end

"""
    nonlinear_eom_triple!(dx, x, p, t)

Ecuaciones de movimiento no lineales del pendulo triple para
DifferentialEquations.jl. Delega en `ModelNLink.nonlinear_eom_nlink!`: la
matriz de masa M(q) es 4x4, simetrica y definida positiva, y las aceleraciones
se obtienen por Cholesky.

La fuerza de control se pasa como p.F (escalar) o mediante p.controller(x, t).
"""
nonlinear_eom_triple!(dx, x, p, t) = nonlinear_eom_nlink!(dx, x, p, t)

"""
    closed_loop_eom_triple!(dx, x, p, t)

Ecuaciones en lazo cerrado con u = -K x. Reutiliza la EOM no lineal, de modo
que la fisica no se duplica.

Parametros esperados en p:
    p.params    -- SystemParamsTriple del modelo
    p.K         -- Matriz de ganancia (1x8)
    p.saturate  -- (opcional) limite de fuerza [N]

Sobre la saturacion: el pipeline usa 150 N por defecto, frente a los 50 N de la
Configuracion I. No es que la condicion inicial nominal lo exija --con los tres
angulos a 0.10 rad los terminos de K alternan signo y se cancelan, y el pico
real es de 7.29 N--, sino que la frontera practica del triple se estanca en
0.252 rad a partir de unos 73 N: con 150 N el actuador deja de ser la
restriccion activa y lo que limita es la no linealidad.
"""
closed_loop_eom_triple!(dx, x, p, t) = closed_loop_eom_nlink!(dx, x, p, t)

"""
Version funcional (no in-place) para uso rapido.
"""
state_derivative_triple(x, params, F) = state_derivative_nlink(x, params, F)

end # module ModelTriple
