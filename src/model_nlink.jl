# =============================================================================
# model_nlink.jl -- Modelo generico de N eslabones articulados sobre un carro
# =============================================================================
# Generaliza las Configuraciones I y II a un numero arbitrario de eslabones en
# serie sobre el carro, con inercia arbitraria en cada uno. Las Configuraciones
# ya entregadas son casos particulares:
#
#   N = 1 con barra uniforme (a = l/2, Il = m l^2 / 12)  ->  model_simple.jl
#   N = 2 con masas puntuales (a = l, Il = 0)            ->  model_double.jl
#
# Este modulo NO reemplaza a los dos anteriores: coexiste con ellos. Asi la
# comparacion entre ambas implementaciones (independientes, escritas por
# separado) se convierte en una prueba de regresion genuina.
#
# Coordenadas generalizadas:  q = [pos, theta_1, ..., theta_N],  d = N + 1
# Variables de estado (dimension n = 2(N+1)), en orden INTERCALADO:
#   x = [pos, vel, theta_1, omega_1, theta_2, omega_2, ..., theta_N, omega_N]
#   pos      -> posicion horizontal del carro
#   vel      -> velocidad del carro
#   theta_j  -> angulo ABSOLUTO del eslabon j desde la vertical superior
#   omega_j  -> velocidad angular del eslabon j
#
# Los angulos son absolutos respecto a la vertical, no relativos al eslabon
# anterior: es la convencion de model_double.jl y hay que conservarla.
#
# Las ecuaciones salen del formalismo de Euler-Lagrange y se escriben en la
# forma matricial estandar de la mecanica:
#     M(q) qdd + C(q, qd) qd + G(q) = F u
# =============================================================================

module ModelNLink

using LinearAlgebra

export SystemParamsNLink, uniform_rods, point_masses
export link_couplings, mass_matrix
export nonlinear_eom_nlink!, closed_loop_eom_nlink!, state_derivative_nlink

"""
Parametros fisicos del sistema carro-pendulo de N eslabones.

    M    masa del carro [kg]
    m    masas de los eslabones [kg]                      (longitud N)
    l    longitudes de los eslabones [m]                  (longitud N)
    a    distancia de la articulacion inferior al centro de masa [m]
    Il   momentos de inercia respecto al centro de masa [kg m^2]
    g    gravedad [m/s^2]
    b    friccion viscosa del carro [N s/m]
    c    friccion RELATIVA en las articulaciones [N m s]  (vacio = sin friccion)

La distincion entre `l` y `a` es la que permite unificar los dos modelos
existentes: `l` es la longitud completa del eslabon (la que determina donde se
articula el siguiente) mientras que `a` es la distancia al centro de masa (la
que determina el momento estatico). Para una barra uniforme a = l/2; para una
masa puntual en el extremo a = l e Il = 0.
"""
Base.@kwdef struct SystemParamsNLink
    M::Float64
    m::Vector{Float64}
    l::Vector{Float64}
    a::Vector{Float64}
    Il::Vector{Float64}
    g::Float64 = 9.81
    b::Float64 = 0.0
    c::Vector{Float64} = Float64[]
end

"""
    uniform_rods(M, m, l; g=9.81, b=0.0, c=Float64[]) -> SystemParamsNLink

Constructor para eslabones modelados como BARRAS UNIFORMES:

    a_j = l_j / 2,      Il_j = m_j l_j^2 / 12

Es el caso de la Configuracion I y el propuesto para la Configuracion III.
"""
function uniform_rods(M, m, l; g=9.81, b=0.0, c=Float64[])
    m = collect(float.(m))
    l = collect(float.(l))
    return SystemParamsNLink(M=float(M), m=m, l=l,
                             a=l ./ 2, Il=(1/12) .* m .* l .^ 2,
                             g=float(g), b=float(b), c=collect(float.(c)))
end

"""
    point_masses(M, m, l; g=9.81, b=0.0, c=Float64[]) -> SystemParamsNLink

Constructor para eslabones modelados como MASAS PUNTUALES en el extremo, con
barras de masa despreciable:

    a_j = l_j,      Il_j = 0

Es el caso de la Configuracion II.
"""
function point_masses(M, m, l; g=9.81, b=0.0, c=Float64[])
    m = collect(float.(m))
    l = collect(float.(l))
    return SystemParamsNLink(M=float(M), m=m, l=l,
                             a=copy(l), Il=zeros(length(m)),
                             g=float(g), b=float(b), c=collect(float.(c)))
end

"""
    link_couplings(p::SystemParamsNLink) -> (beta, J, Gam)

Las tres agrupaciones de parametros que organizan todo el modelo:

    beta_j   = m_j a_j + l_j * sum_{i>j} m_i          momento estatico
    J_j      = Il_j + m_j a_j^2 + l_j^2 * sum_{i>j} m_i   inercia efectiva
    Gam_jk   = l_min(j,k) * beta_max(j,k)             acoplamiento inercial

Interpretacion fisica:

  - beta_j es el momento estatico del eslabon j y de todo lo que carga encima:
    su masa propia por su brazo, mas la longitud completa del eslabon por la
    masa de los eslabones superiores.
  - J_j es el momento de inercia del eslabon j respecto a SU articulacion,
    incluyendo la carga de los superiores (Steiner mas el resto).
  - Gam_jk mide cuanto se acelera el eslabon j cuando se acelera el k.

Observacion estructural: el mismo beta_j que acopla el carro con el eslabon j
en la energia cinetica es el que pesa la gravedad sobre ese eslabon, porque
ambos provienen del mismo momento estatico. La potencial es simplemente

    V = g * sum_j beta_j cos(theta_j)

Esa coincidencia no es casual, y es la razon estructural de que la
inestabilidad y la controlabilidad esten ligadas: el canal por el que la
gravedad desestabiliza es el mismo por el que el carro puede actuar.
"""
function link_couplings(p::SystemParamsNLink)
    N = length(p.m)
    beta = similar(p.m)
    J = similar(p.m)

    for j in 1:N
        m_above = sum(@view p.m[j+1:end])
        beta[j] = p.m[j] * p.a[j] + p.l[j] * m_above
        J[j] = p.Il[j] + p.m[j] * p.a[j]^2 + p.l[j]^2 * m_above
    end

    Gam = [p.l[min(j, k)] * beta[max(j, k)] for j in 1:N, k in 1:N]
    return (beta=beta, J=J, Gam=Gam)
end

"""
    mass_matrix(p::SystemParamsNLink, thetas) -> Matrix{Float64}

Matriz de masa M(q) de dimension (N+1) x (N+1) con q = [pos, theta_1, ...]:

    M[1,1]     = M + sum(m)
    M[1,j+1]   = beta_j cos(theta_j)
    M[j+1,j+1] = J_j
    M[j+1,k+1] = Gam_jk cos(theta_j - theta_k)      para j != k

Es simetrica y definida positiva por construccion: es la Hessiana de la energia
cinetica respecto a las velocidades generalizadas. Por eso las aceleraciones se
obtienen con una factorizacion de Cholesky en vez de una LU generica.

Evaluada en thetas = 0 da la matriz constante M_0 que necesita la
linealizacion; se expone precisamente para que linearization.jl no tenga que
duplicar la fisica.
"""
function mass_matrix(p::SystemParamsNLink, thetas)
    bc = link_couplings(p)
    return mass_matrix(p, thetas, bc)
end

function mass_matrix(p::SystemParamsNLink, thetas, bc)
    N = length(p.m)
    d = N + 1
    Mq = zeros(d, d)

    Mq[1, 1] = p.M + sum(p.m)
    for j in 1:N
        cj = cos(thetas[j])
        Mq[1, j+1] = bc.beta[j] * cj
        Mq[j+1, 1] = Mq[1, j+1]
        Mq[j+1, j+1] = bc.J[j]
        for k in 1:N
            j == k && continue
            Mq[j+1, k+1] = bc.Gam[j, k] * cos(thetas[j] - thetas[k])
        end
    end
    return Mq
end

"""
    nonlinear_eom_nlink!(dx, x, p, t)

Ecuaciones de movimiento no lineales para DifferentialEquations.jl. La fuerza
de control se pasa como p.F (escalar) o mediante p.controller(x, t).

Ecuacion del carro:

    (M + sum m_i) ddx + sum_j beta_j (cos(theta_j) ddtheta_j
                                      - sin(theta_j) omega_j^2) + b vel = u

Ecuacion del eslabon j:

    beta_j cos(theta_j) ddx + J_j ddtheta_j
      + sum_{k != j} Gam_jk [cos(theta_j - theta_k) ddtheta_k
                             + sin(theta_j - theta_k) omega_k^2]
      = g beta_j sin(theta_j) - tau_j^fric

El signo POSITIVO de g beta_j sin(theta_j) en el lado derecho es la firma de la
inestabilidad: para theta_j pequeno el termino es +g beta_j theta_j, una
realimentacion positiva. Es el mismo mecanismo que el +(M+m) m g L / D0 de la
Configuracion I.

La friccion articular, si esta activa, es un par RELATIVO entre eslabones
contiguos, que es lo fisicamente correcto:

    tau_j^fric = c_j (omega_j - omega_{j-1}) - c_{j+1} (omega_{j+1} - omega_j)

con omega_0 := 0. Produce una matriz de amortiguamiento tridiagonal simetrica.
"""
function nonlinear_eom_nlink!(dx, x, p, t)
    sp = p.params
    N = length(sp.m)
    d = N + 1

    # Fuerza de control
    F = haskey(p, :controller) ? p.controller(x, t) : p.F

    # Desempaquetar el estado intercalado
    vel = x[2]
    thetas = @view x[3:2:2d]      # theta_j en la posicion 2j+1
    omegas = @view x[4:2:2d]      # omega_j en la posicion 2j+2

    bc = link_couplings(sp)
    Mq = mass_matrix(sp, thetas, bc)

    # ---------------------------------------------------------------
    # Lado derecho: control, friccion, terminos centrifugos y gravedad
    # ---------------------------------------------------------------
    rhs = zeros(d)

    # Fila del carro
    rhs[1] = F - sp.b * vel
    for j in 1:N
        rhs[1] += bc.beta[j] * sin(thetas[j]) * omegas[j]^2
    end

    # Fila de cada eslabon
    for j in 1:N
        acc = sp.g * bc.beta[j] * sin(thetas[j])
        for k in 1:N
            j == k && continue
            acc -= bc.Gam[j, k] * sin(thetas[j] - thetas[k]) * omegas[k]^2
        end
        rhs[j+1] = acc
    end

    # Friccion articular relativa (solo si esta activada)
    if !isempty(sp.c)
        for j in 1:N
            omega_below = j > 1 ? omegas[j-1] : 0.0
            tau = sp.c[j] * (omegas[j] - omega_below)
            if j < N
                tau -= sp.c[j+1] * (omegas[j+1] - omegas[j])
            end
            rhs[j+1] -= tau
        end
    end

    # M(q) es simetrica definida positiva, asi que se resuelve por Cholesky.
    qdd = cholesky(Symmetric(Mq)) \ rhs

    # Reempaquetar en orden intercalado
    dx[1] = vel
    dx[2] = qdd[1]
    for j in 1:N
        dx[2j+1] = omegas[j]
        dx[2j+2] = qdd[j+1]
    end
    return nothing
end

"""
    closed_loop_eom_nlink!(dx, x, p, t)

Ecuaciones en lazo cerrado con la ley de control u = -K x. Reutiliza el modelo
no lineal pasando la fuerza ya calculada, para no duplicar la fisica.

Parametros esperados en p:
    p.params    -- SystemParamsNLink del modelo
    p.K         -- Matriz de ganancia (1 x n)
    p.saturate  -- (opcional) limite de fuerza [N]
"""
function closed_loop_eom_nlink!(dx, x, p, t)
    K = p.K
    u = -dot(K[1, :], x)  # escalar para SISO

    if haskey(p, :saturate)
        u = clamp(u, -p.saturate, p.saturate)
    end

    nonlinear_eom_nlink!(dx, x, (params=p.params, F=u), t)
    return nothing
end

"""
Version funcional (no in-place) para uso rapido.
"""
function state_derivative_nlink(x, params::SystemParamsNLink, F::Float64)
    dx = zeros(2 * (length(params.m) + 1))
    nonlinear_eom_nlink!(dx, x, (params=params, F=F), 0.0)
    return dx
end

end # module ModelNLink
