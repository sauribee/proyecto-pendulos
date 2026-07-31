# =============================================================================
# observer.jl -- Estimacion de estado por observador de Luenberger
# =============================================================================
# Cierra la grieta mas visible del trabajo. El proyecto DEMUESTRA que
# rank(Obsv) = n en las tres configuraciones, y acto seguido controla con
# u = -K x, es decir con el estado COMPLETO, que en un sistema real no se mide:
# nadie tiene un sensor de velocidad angular en cada articulacion. Si el estado
# se puede reconstruir, hay que reconstruirlo.
#
# LA IDEA. Se estima el estado con una copia del modelo corregida por el error
# de prediccion de la salida:
#
#     xhat' = A xhat + B u + L (y - C xhat)
#
# y el error e = x - xhat obedece  e' = (A - L C) e.  Disenar L para que
# A - LC sea Hurwitz es EL PROBLEMA DUAL de disenar K para que A - BK lo sea:
#
#     (A, B) controlable   <=>   (A', B') observable
#
# de modo que, por la dualidad de Kalman,
#
#     L = lqr(A', C', Qo, Ro)'
#
# y NO hace falta escribir un solo algoritmo nuevo: se reutiliza design_lqr,
# con su Riccati resuelta por el Hamiltoniano. Ese es el punto elegante: la
# dualidad convierte un problema nuevo en uno ya resuelto.
#
# EL PRINCIPIO DE SEPARACION sale de una observacion de algebra lineal de dos
# lineas y esta documentado en augmented_closed_loop.
#
# REQUISITO DE CARGA: linearization.jl, controller.jl y metrics.jl deben
# incluirse ANTES que este archivo.
# =============================================================================

module Observer

using LinearAlgebra
using Printf

using ..Controller
using ..Metrics

export design_observer, design_observer_poles
export augmented_closed_loop, check_separation
export observer_eom!, sensor_subset_analysis, print_sensor_table

# ---------------------------------------------------------------------------
# 1. Diseno de la ganancia del observador
# ---------------------------------------------------------------------------

"""
    design_observer(A, C, Qo, Ro) -> (L, eigenvalues_obs, A_obs, P)

Ganancia del observador de Luenberger por DUALIDAD de Kalman:

    L = design_lqr(A', C', Qo, Ro).K'

No se reimplementa nada. El par (A', C') del sistema dual juega el papel de
(A, B) en el problema de control, asi que la misma ecuacion de Riccati que da
la ganancia optima de realimentacion da, transpuesta, la ganancia optima de
inyeccion de salida.

Interpretacion estocastica (LQG): si Qo y Ro se leen como las covarianzas del
ruido de proceso y del ruido de medicion, L es exactamente la ganancia del
FILTRO DE KALMAN en regimen permanente, y el conjunto K + L es un controlador
LQG. La sintonia entonces tiene sentido fisico: Ro grande significa sensores
ruidosos, y el observador confia menos en la medicion.

Convenciones de tamano: A es n x n, C es p x n, Qo es n x n, Ro es p x p, y la
L devuelta es n x p.
"""
function design_observer(A, C, Qo, Ro)
    res = design_lqr(Matrix(A'), Matrix(C'), Qo, Ro)
    L = Matrix(res.K')
    A_obs = A - L * C
    return (L=L, eigenvalues_obs=eigvals(A_obs), A_obs=A_obs, P=res.P)
end

"""
    design_observer_poles(A, C, desired_poles) -> (L, eigenvalues_obs, A_obs)

Asignacion de polos del observador por el dual de la formula de Ackermann,
aplicada sobre (A', C').

RESTRICCION: Ackermann solo aplica a sistemas de una entrada, asi que su dual
solo aplica a sistemas de UNA SALIDA. Para la Configuracion III, que mide
cuatro salidas, hay que usar `design_observer` (LQR dual), que no tiene esa
limitacion. Se incluye por completitud y para el paralelo con
design_pole_placement.
"""
function design_observer_poles(A, C, desired_poles)
    size(C, 1) == 1 ||
        error("design_observer_poles solo aplica con una salida (C de una " *
              "fila); para varias salidas usar design_observer")
    res = design_pole_placement(Matrix(A'), Matrix(C'), desired_poles)
    L = Matrix(res.K')
    A_obs = A - L * C
    return (L=L, eigenvalues_obs=eigvals(A_obs), A_obs=A_obs)
end

# ---------------------------------------------------------------------------
# 2. Principio de separacion
# ---------------------------------------------------------------------------

"""
    augmented_closed_loop(A, B, C, K, L) -> Matrix

Matriz del sistema aumentado en coordenadas (x, e), con e = x - xhat el error
de estimacion.

DEDUCCION. Con la planta x' = Ax + Bu, el control u = -K xhat y el observador
xhat' = A xhat + Bu + L(Cx - C xhat):

    x' = Ax - BK xhat = Ax - BK(x - e) = (A - BK) x + BK e
    e' = x' - xhat' = A(x - xhat) - LC(x - xhat) = (A - LC) e

es decir

    (x', e') = [ A - BK    BK    ] (x, e)
               [   0     A - LC  ]

EL PRINCIPIO DE SEPARACION. Esa matriz es TRIANGULAR POR BLOQUES, y el espectro
de una triangular por bloques es la union de los espectros de sus bloques
diagonales:

    spec = spec(A - BK)  union  spec(A - LC)

Por tanto el diseno del controlador y el del observador son INDEPENDIENTES: se
puede elegir K sin pensar en L y viceversa, y el lazo cerrado con estimacion es
estable si ambos lo son por separado. Un resultado de control que se demuestra
en dos lineas de algebra lineal, sin ninguna maquinaria adicional.

Notese que el bloque (1,2) NO es nulo: el error de estimacion SI perturba a la
planta. Lo que la triangularidad garantiza es que esa perturbacion no puede
desestabilizar, porque el error decae por su cuenta.
"""
function augmented_closed_loop(A, B, C, K, L)
    n = size(A, 1)
    M = zeros(2n, 2n)
    M[1:n, 1:n] = A - B * K
    M[1:n, n+1:2n] = B * K
    M[n+1:2n, n+1:2n] = A - L * C
    return M
end

"""
    check_separation(A, B, C, K, L; tol=1e-8) -> (holds, max_error, ...)

Verifica NUMERICAMENTE el principio de separacion: compara el espectro de la
matriz aumentada contra la union de spec(A - BK) y spec(A - LC), emparejando
cada eigenvalor con el mas cercano del conjunto de referencia.
"""
function check_separation(A, B, C, K, L; tol=1e-8)
    aug = eigvals(augmented_closed_loop(A, B, C, K, L))
    ref = vcat(eigvals(A - B * K), eigvals(A - L * C))

    disponibles = collect(ref)
    max_err = 0.0
    for lambda in aug
        i = argmin(abs.(disponibles .- lambda))
        max_err = max(max_err, abs(disponibles[i] - lambda))
        deleteat!(disponibles, i)
    end
    return (holds=max_err < tol, max_error=max_err,
            spectrum_augmented=aug, spectrum_reference=ref)
end

# ---------------------------------------------------------------------------
# 3. Dinamica conjunta planta no lineal + observador lineal
# ---------------------------------------------------------------------------

"""
    observer_eom!(dz, z, p, t)

Dinamica conjunta de la planta NO LINEAL y el observador LINEAL, para
DifferentialEquations.jl. El estado tiene 2n componentes:

    z[1:n]      estado real, integrado con la EOM no lineal (p.eom!)
    z[n+1:2n]   estado estimado, integrado con el observador lineal

El control es u = -K xhat, NO -K x: ese es justamente el punto. El observador
solo ve la medicion y = C x.

Parametros esperados en p:
    p.params    -- parametros fisicos del modelo
    p.eom!      -- EOM no lineal in-place del modelo (recibe (dx, x, (params, F), t))
    p.A, p.B, p.C, p.K, p.L
    p.saturate  -- (opcional) limite de fuerza [N]
    p.noise     -- (opcional) funcion t -> vector de p componentes que se suma
                   a la medicion

IMPORTANTE sobre la saturacion: el observador debe alimentarse con la fuerza
REALMENTE aplicada, es decir la ya saturada, no con la que el control pidio. Si
se le pasa la no saturada, el observador cree que la planta recibio mas fuerza
de la que recibio y el error de estimacion no converge. Es un error clasico.

NOTA sobre el ruido: `p.noise` es una funcion DETERMINISTA del tiempo, no
`randn()` evaluado dentro de la funcion. Un solver adaptativo evalua el lado
derecho varias veces por paso, con pasos que ademas rechaza y repite, asi que
llamar a `randn()` aqui no produce ruido blanco sino una perturbacion
irreproducible que rompe el control de paso. Para ruido estocastico de verdad
hay que usar la maquinaria de EDEs de DifferentialEquations.
"""
function observer_eom!(dz, z, p, t)
    n = size(p.A, 1)
    x = @view z[1:n]
    xhat = @view z[n+1:2n]

    # Control calculado con el estado ESTIMADO
    u = -dot(p.K[1, :], xhat)
    if haskey(p, :saturate)
        u = clamp(u, -p.saturate, p.saturate)
    end

    # Planta: modelo no lineal completo, con la fuerza realmente aplicada
    dx = @view dz[1:n]
    p.eom!(dx, x, (params=p.params, F=u), t)

    # Medicion
    y = p.C * x
    if haskey(p, :noise)
        y = y .+ p.noise(t)
    end

    # Observador: copia lineal del modelo corregida por el error de salida
    dz[n+1:2n] .= p.A * xhat .+ vec(p.B .* u) .+ p.L * (y .- p.C * xhat)
    return nothing
end

# ---------------------------------------------------------------------------
# 4. Que sensores hacen falta
# ---------------------------------------------------------------------------

"""
    sensor_subset_analysis(A, C_full, labels) -> Vector{NamedTuple}

Para cada uno de los 2^p - 1 subconjuntos no vacios de filas de `C_full`,
calcula el rango de la matriz de observabilidad y el margen PBH de
observabilidad normalizado. Devuelve la tabla ordenada por margen decreciente.

Responde una pregunta de ingenieria con una respuesta de algebra lineal: CUAL
ES EL JUEGO MINIMO DE SENSORES. El rango solo dice si se puede o no; el margen
dice cuanto se puede confiar en que se pueda, que es lo que importa cuando el
condicionamiento se dispara.

Para la Configuracion III, con C_full midiendo pos, theta1, theta2 y theta3,
son 15 combinaciones.

Cada entrada tiene los campos `indices`, `labels`, `rank`, `observable`,
`margin` y `margin_normalized`.
"""
function sensor_subset_analysis(A, C_full, labels)
    n = size(A, 1)
    p = size(C_full, 1)
    length(labels) == p ||
        error("se esperaban $p etiquetas, se recibieron $(length(labels))")

    filas = Vector{Any}()
    for mask in 1:(2^p - 1)
        idx = [i for i in 1:p if (mask >> (i - 1)) & 1 == 1]
        C = C_full[idx, :]
        O = vcat([C * A^k for k in 0:n-1]...)
        r = rank(O)
        push!(filas, (indices=idx, labels=labels[idx], rank=r,
                      observable=(r == n),
                      margin=pbh_observability_margin(A, C),
                      margin_normalized=pbh_observability_margin_normalized(A, C)))
    end
    return sort(filas, by=f -> -f.margin_normalized)
end

"""
    print_sensor_table(tabla, n; label="")

Imprime la tabla de subconjuntos de sensores con el formato del proyecto y
senala cual es el juego minimo que conserva la observabilidad.
"""
function print_sensor_table(tabla, n; label="")
    println("=" ^ 78)
    println("  SUBCONJUNTOS DE SENSORES" * (isempty(label) ? "" : ": $label"))
    println("=" ^ 78)
    @printf("  %-34s %6s %8s %14s\n", "sensores", "rank", "observ", "margen PBH")
    println("  " * "-" ^ 74)
    for f in tabla
        @printf("  %-34s %4d/%d %8s %14.4e\n",
                join(f.labels, "+"), f.rank, n,
                f.observable ? "SI" : "NO", f.margin_normalized)
    end

    observables = filter(f -> f.observable, tabla)
    if isempty(observables)
        println("\n  Ningun subconjunto conserva la observabilidad.")
    else
        minimos = filter(f -> length(f.indices) ==
                              minimum(length(g.indices) for g in observables),
                         observables)
        mejor = minimos[1]
        @printf("\n  Juego MINIMO observable: %s (%d sensor(es), margen %.4e)\n",
                join(mejor.labels, "+"), length(mejor.indices),
                mejor.margin_normalized)
        @printf("  Mejor margen en total:   %s (%.4e)\n",
                join(observables[1].labels, "+"), observables[1].margin_normalized)
    end
    println("=" ^ 78)
end

end # module Observer
