# =============================================================================
# sweep.jl -- Barridos del espacio de parametros y fronteras de operabilidad
# =============================================================================
# Responde tres preguntas sobre el controlador:
#
#   1. Donde funciona bien?     Existe una configuracion que maximice el margen
#                               de control? (Adelanto: si, y no esta en los
#                               extremos del dominio.)
#   2. Donde funciona apenas?   Como se degrada al alejarse de ese optimo?
#   3. Donde es imposible?      Cual es la frontera mas alla de la cual no hay
#                               controlador que sirva?
#
# Conviene separar DOS TIPOS DE IMPOSIBILIDAD, porque el informe debe
# distinguirlas:
#
#   - ESTRUCTURAL: el par (A,B) deja de ser controlable, o se acerca tanto que
#     ningun K razonable sirve. Es una propiedad DEL SISTEMA, independiente del
#     controlador. Se mide con el margen PBH de Metrics.
#   - PRACTICA: el sistema es controlable, pero el actuador satura, el riel se
#     acaba o el periodo de muestreo es muy largo. Es una propiedad DE LA
#     IMPLEMENTACION. Se mide simulando el modelo NO LINEAL y buscando la
#     frontera por biseccion.
#
# Todo es generico en la dimension del estado, igual que Controller y Metrics:
# los mismos barridos corren sobre las tres configuraciones.
#
# REQUISITO DE CARGA: model_nlink.jl, linearization.jl, controller.jl y
# metrics.jl deben incluirse ANTES que este archivo.
# =============================================================================

module Sweep

using LinearAlgebra
using Printf
using DifferentialEquations

using ..ModelNLink
using ..Linearization
using ..Controller
using ..Metrics

export viability, default_weights
export set_param, infer_geometry
export sweep_1d, sweep_2d
export max_recoverable_angle, min_actuator_force, max_sampling_period
export recovers, print_sweep_summary

# ---------------------------------------------------------------------------
# 1. Manipulacion de parametros
# ---------------------------------------------------------------------------

"""
    infer_geometry(p) -> Symbol

Deduce como estan modelados los eslabones a partir de la relacion entre `a`,
`l` e `Il`:

    :uniform  barras uniformes    a = l/2,  Il = m l^2 / 12
    :point    masas puntuales     a = l,    Il = 0
    :custom   ninguna de las dos

Sirve para que al barrer una longitud se recalculen tambien el centro de masa y
la inercia, que en una barra uniforme NO son independientes de ella. Ignorar
esto es el error silencioso mas facil de cometer en un barrido geometrico: se
obtendria una barra que cambia de largo sin cambiar de inercia.
"""
function infer_geometry(p)
    if p.a ≈ p.l ./ 2 && p.Il ≈ (1/12) .* p.m .* p.l .^ 2
        return :uniform
    elseif p.a ≈ p.l && all(abs.(p.Il) .< 1e-14)
        return :point
    else
        return :custom
    end
end

"""
    set_param(p::SystemParamsNLink, spec, value; geometry=infer_geometry(p))

Devuelve una copia de `p` con un parametro modificado. `spec` puede ser:

    :M, :g, :b            un campo escalar
    (:m, j), (:l, j)      el elemento j de un campo vectorial
    (:a, j), (:Il, j)     idem, para geometria :custom

Si la geometria es `:uniform` o `:point`, al cambiar `l[j]` o `m[j]` se
recalculan `a` e `Il` de forma coherente.
"""
function set_param(p::SystemParamsNLink, spec, value; geometry=infer_geometry(p))
    M, g, b, c = p.M, p.g, p.b, p.c
    m, l, a, Il = copy(p.m), copy(p.l), copy(p.a), copy(p.Il)

    if spec === :M
        M = value
    elseif spec === :g
        g = value
    elseif spec === :b
        b = value
    elseif spec isa Tuple
        campo, j = spec
        if campo === :m
            m[j] = value
        elseif campo === :l
            l[j] = value
        elseif campo === :a
            a[j] = value
        elseif campo === :Il
            Il[j] = value
        else
            error("campo desconocido en set_param: $campo")
        end
    else
        error("spec desconocido en set_param: $spec")
    end

    # Reconstruir la geometria derivada si corresponde
    if geometry === :uniform
        a = l ./ 2
        Il = (1/12) .* m .* l .^ 2
    elseif geometry === :point
        a = copy(l)
        Il = zeros(length(m))
    end

    return SystemParamsNLink(M=M, m=m, l=l, a=a, Il=Il, g=g, b=b, c=c)
end

"""
    default_weights(n) -> Matrix

Pesos Q del LQR usados en las tres configuraciones, extendidos a dimension n:
se penalizan la posicion del carro (peso 1) y todos los angulos (peso 10); las
velocidades no se penalizan.
"""
default_weights(n) = diagm(vcat([1.0, 0.0], repeat([10.0, 0.0], (n - 2) ÷ 2)))

# ---------------------------------------------------------------------------
# 2. Viabilidad: todas las metricas de una configuracion
# ---------------------------------------------------------------------------

"""
    viability(linearize, params; Q=nothing, R=0.1, umax=50.0, x0=nothing,
              angle_idx=3) -> NamedTuple

Evalua de una sola pasada todas las metricas de viabilidad de una
configuracion. `linearize` es una funcion `params -> StateSpaceModel`, tipicamente
`Linearization.linearize_system_nlink`.

Devuelve un NamedTuple con:

    lambda_max, tau_fall, n_unstable    escala de tiempo del problema
    rank_ctrb, cond_ctrb                criterio de Kalman y su condicionamiento
    pbh_margin, pbh_normalized          distancia a la incontrolabilidad
    pbh_obs_normalized                  idem para observabilidad
    K, norm_K, P                        controlador LQR
    J_star                              costo optimo desde x0
    c_star, theta_guaranteed            elipsoide de no saturacion
    poles_cl                            polos de lazo cerrado

Si la CARE no se puede resolver (sistema no estabilizable) devuelve `nothing` en
los campos del controlador en vez de propagar el error, para que un barrido no
se caiga entero por un punto malo.
"""
function viability(linearize, params; Q=nothing, R=0.1, umax=50.0,
                   x0=nothing, angle_idx=3)
    ss = linearize(params)
    A, B, C = ss.A, ss.B, ss.C
    n = size(A, 1)

    Qm = Q === nothing ? default_weights(n) : Q
    Rm = R isa AbstractMatrix ? R : reshape([float(R)], 1, 1)

    z0 = x0 === nothing ? (v = zeros(n); v[angle_idx] = 0.1; v) : x0

    dm = dominant_mode(A)
    ctrb = check_controllability(ss)

    base = (lambda_max=dm.lambda_max, tau_fall=dm.tau_fall,
            n_unstable=dm.n_unstable,
            rank_ctrb=ctrb.rank, cond_ctrb=controllability_condition(A, B),
            pbh_margin=pbh_controllability_margin(A, B),
            pbh_normalized=pbh_controllability_margin_normalized(A, B),
            pbh_obs_normalized=pbh_observability_margin_normalized(A, C))

    local res
    try
        res = design_lqr(A, B, Qm, Rm)
    catch
        return merge(base, (K=nothing, norm_K=NaN, P=nothing, J_star=NaN,
                            c_star=NaN, theta_guaranteed=NaN, poles_cl=nothing,
                            ss=ss))
    end

    c_star = nonsaturation_level(res.K, res.P, Rm, B, umax)

    return merge(base, (K=res.K, norm_K=norm(res.K), P=res.P,
                        J_star=optimal_cost(res.P, z0),
                        c_star=c_star,
                        theta_guaranteed=max_axis_angle(res.P, c_star, angle_idx),
                        poles_cl=res.eigenvalues_cl,
                        ss=ss))
end

# ---------------------------------------------------------------------------
# 3. Barridos
# ---------------------------------------------------------------------------

"""
    sweep_1d(linearize, base_params, spec, values; kwargs...) -> Vector

Barre un parametro y devuelve un vector de NamedTuple, uno por valor, cada uno
con el campo `value` anadido. Los `kwargs` se pasan a `viability`.
"""
function sweep_1d(linearize, base_params, spec, values; kwargs...)
    geom = infer_geometry(base_params)
    out = Vector{Any}(undef, length(values))
    Threads.@threads for i in eachindex(values)
        p = set_param(base_params, spec, values[i]; geometry=geom)
        out[i] = merge((value=values[i],), viability(linearize, p; kwargs...))
    end
    return [out[i] for i in eachindex(values)]
end

"""
    sweep_2d(linearize, base_params, spec1, vals1, spec2, vals2; kwargs...) -> Matrix

Barre dos parametros y devuelve una matriz de NamedTuple de tamano
`length(vals1) x length(vals2)`. Paraleliza con `Threads.@threads` sobre las
filas: lanzar Julia con `-t auto` reduce el tiempo por el numero de nucleos.
"""
function sweep_2d(linearize, base_params, spec1, vals1, spec2, vals2; kwargs...)
    geom = infer_geometry(base_params)
    out = Matrix{Any}(undef, length(vals1), length(vals2))
    Threads.@threads for i in eachindex(vals1)
        for j in eachindex(vals2)
            p = set_param(base_params, spec1, vals1[i]; geometry=geom)
            p = set_param(p, spec2, vals2[j]; geometry=geom)
            out[i, j] = merge((value1=vals1[i], value2=vals2[j]),
                              viability(linearize, p; kwargs...))
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# 4. Fronteras por biseccion sobre el modelo NO LINEAL
# ---------------------------------------------------------------------------

"""
    recovers(params, K, theta0; umax, T=10.0, x_rail=1.5, angle_idx=3,
             tol_final=0.02, saveat=0.01) -> Bool

Simula el modelo NO LINEAL en lazo cerrado con saturacion, partiendo del estado
con la componente `angle_idx` igual a `theta0` y el resto en cero, y decide si
el controlador recupera el equilibrio.

CRITERIO DE EXITO (hay que fijarlo y documentarlo, porque toda la frontera
depende de el):

    exito  <=>  max_j |theta_j(T)| < 0.02 rad          (converge)
           AND  max_t |x(t)| <= x_rail                 (no se sale del riel)
           AND  max_{t,j} |theta_j(t)| < pi/2          (no se cae)

El umbral de 0.02 rad es el mismo de `settling_time` en make_report_figs.jl: se
reutiliza en vez de inventar otro.

Lleva un callback de terminacion temprana que aborta en cuanto algun angulo
supera pi/2. La mayoria de los fracasos se detectan en las primeras decimas de
segundo, asi que el barrido completo baja de minutos a segundos.
"""
function recovers(params, K, theta0; umax=150.0, T=10.0, x_rail=1.5,
                  angle_idx=3, tol_final=0.02, saveat=0.01)
    n = 2 * (length(params.m) + 1)
    x0 = zeros(n)
    x0[angle_idx] = theta0

    caido = function (u, t, integrator)
        @inbounds for j in 3:2:n-1
            abs(u[j]) > pi/2 && return true
        end
        return false
    end

    prob = ODEProblem(closed_loop_eom_nlink!, x0, (0.0, T),
                      (params=params, K=K, saturate=umax))
    sol = try
        solve(prob, Tsit5(); saveat=saveat, save_everystep=false,
              callback=DiscreteCallback(caido, terminate!),
              verbose=false, maxiters=10^7)
    catch
        return false
    end

    # Aborto temprano o fallo del solver: no recupera.
    sol.retcode == ReturnCode.Success || return false

    for u in sol.u
        abs(u[1]) <= x_rail || return false
        @inbounds for j in 3:2:n-1
            abs(u[j]) < pi/2 || return false
        end
    end

    ufin = sol.u[end]
    @inbounds for j in 3:2:n-1
        abs(ufin[j]) < tol_final || return false
    end
    return true
end

"""
    max_recoverable_angle(params, K; umax=150.0, hi=1.2, iters=14, kwargs...)

Mayor angulo inicial del que el controlador todavia se recupera, por biseccion
sobre `theta0`. Responde "desde que tan inclinado puede recuperarse".

Devuelve 0.0 si ni siquiera el angulo mas pequeno de la busqueda funciona.
"""
function max_recoverable_angle(params, K; umax=150.0, hi=1.2, iters=14, kwargs...)
    lo = 0.0
    recovers(params, K, hi; umax=umax, kwargs...) && return hi
    for _ in 1:iters
        mid = (lo + hi) / 2
        recovers(params, K, mid; umax=umax, kwargs...) ? (lo = mid) : (hi = mid)
    end
    return lo
end

"""
    min_actuator_force(params, K; theta0=0.2, lo=1.0, hi=1000.0, iters=14, kwargs...)

Menor limite de actuador con el que el controlador todavia recupera desde
`theta0`, por biseccion GEOMETRICA sobre `umax` (geometrica porque el rango
util cubre tres ordenes de magnitud). Responde "que motor hay que comprar".

Devuelve `Inf` si ni con `hi` se recupera: en ese caso el actuador no es la
restriccion activa, y el limite esta en el riel o en la no linealidad.
"""
function min_actuator_force(params, K; theta0=0.2, lo=1.0, hi=1000.0,
                            iters=14, kwargs...)
    recovers(params, K, theta0; umax=hi, kwargs...) || return Inf
    for _ in 1:iters
        mid = sqrt(lo * hi)
        recovers(params, K, theta0; umax=mid, kwargs...) ? (hi = mid) : (lo = mid)
    end
    return hi
end

"""
    max_sampling_period(A, B, K; hmin=1e-4, hmax=2.0, iters=40) -> h

Mayor periodo de muestreo con el que el lazo cerrado DISCRETO sigue siendo
estable. Se discretiza con retenedor de orden cero (Metrics.discretize_zoh, que
usa el truco de Van Loan) y se exige que `Ad - Bd K` sea de SCHUR, es decir que
todos sus eigenvalores esten dentro del circulo unitario: el analogo discreto
de la condicion de Hurwitz.

Es puro exponencial matricial, sin integrar ninguna EDO. La biseccion es
geometrica porque h recorre varios ordenes de magnitud.

Valores de referencia con los parametros por defecto: simple 197 ms,
doble 77 ms, triple 32 ms.
"""
function max_sampling_period(A, B, K; hmin=1e-4, hmax=2.0, iters=40)
    es_estable = function (h)
        d = discretize_zoh(A, B, h)
        return is_schur(d.Ad - d.Bd * K)
    end
    es_estable(hmin) || return 0.0
    es_estable(hmax) && return hmax

    lo, hi = hmin, hmax
    for _ in 1:iters
        mid = sqrt(lo * hi)
        es_estable(mid) ? (lo = mid) : (hi = mid)
    end
    return lo
end

# ---------------------------------------------------------------------------
# 5. Reporte
# ---------------------------------------------------------------------------

"""
    print_sweep_summary(results; spec="parametro", label="")

Imprime una tabla de un barrido 1D con el formato de banner del proyecto.
"""
function print_sweep_summary(results; spec="parametro", label="")
    println("=" ^ 78)
    println("  BARRIDO 1D" * (isempty(label) ? "" : ": $label"))
    println("=" ^ 78)
    @printf("  %-12s %10s %12s %12s %10s %10s\n",
            string(spec), "lambda_max", "PBH norm", "cond(Ctrb)", "norm(K)", "c*")
    println("  " * "-" ^ 74)
    for r in results
        @printf("  %-12.4g %10.3f %12.3e %12.3e %10.1f %10.3f\n",
                r.value, r.lambda_max, r.pbh_normalized, r.cond_ctrb,
                r.norm_K, r.c_star)
    end
    println("=" ^ 78)
end

end # module Sweep
