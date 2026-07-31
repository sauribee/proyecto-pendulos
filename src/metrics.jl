# =============================================================================
# metrics.jl -- Metricas cuantitativas de controlabilidad y viabilidad
# =============================================================================
# El criterio de Kalman rank(C) = n es binario y se evalua con una tolerancia
# numerica implicita. La matriz de controlabilidad C = [B AB ... A^(n-1)B] es
# una base de Krylov, notoriamente mal condicionada: su numero de condicion
# crece unos cuatro ordenes de magnitud por cada eslabon que se agrega, y para
# cinco eslabones alcanza 1/eps_maquina. En ese punto el rango numerico deja de
# significar algo, aunque el sistema siga siendo controlable en aritmetica
# exacta.
#
# Este modulo sustituye (o acompana) ese criterio binario por medidas continuas:
#   1. Margen de Popov-Belevitch-Hautus: distancia a la incontrolabilidad.
#   2. Numero de condicion de la matriz de Kalman: mide cuando dejar de confiar
#      en el rango.
#   3. Costo optimo J* = x0' P x0 y norma de la ganancia.
#   4. Elipsoide de no saturacion: garantia dura sobre el sistema lineal.
#   5. Discretizacion exacta con retenedor de orden cero y criterio de Schur.
#
# Todas las rutinas son genericas en la dimension del estado, igual que
# Controller y las funciones de analisis de Linearization: el mismo codigo
# opera sobre el sistema 4x4 del pendulo simple, el 6x6 del doble y el 8x8
# del triple.
# =============================================================================

module Metrics

using LinearAlgebra
using Printf

export dominant_mode
export pbh_controllability_margin, pbh_controllability_margin_normalized
export pbh_observability_margin, pbh_observability_margin_normalized
export pbh_controllability_profile
export controllability_condition, observability_condition
export optimal_cost, nonsaturation_level, max_axis_angle
export discretize_zoh, is_schur
export print_metrics_summary

# ---------------------------------------------------------------------------
# 1. Escala de tiempo del problema
# ---------------------------------------------------------------------------

"""
    dominant_mode(A) -> (lambda_max, tau_fall, n_unstable)

Modo inestable dominante del sistema en lazo abierto:

    lambda_max = max_i Re(lambda_i(A)),    tau_fall = 1 / lambda_max

`tau_fall` es la constante de tiempo con la que el sistema se desploma: fija la
escala temporal del problema, porque el controlador debe actuar en una fraccion
de ella. Si no hay modos inestables, `tau_fall` es `Inf`.

`n_unstable` cuenta los eigenvalores con parte real estrictamente positiva
(tolerancia 1e-10 para no contar el modo marginal del carro, que es cero).
"""
function dominant_mode(A)
    lambda = eigvals(A)
    lambda_max = maximum(real.(lambda))
    tau_fall = lambda_max > 1e-10 ? 1.0 / lambda_max : Inf
    n_unstable = count(real.(lambda) .> 1e-10)
    return (lambda_max=lambda_max, tau_fall=tau_fall, n_unstable=n_unstable)
end

# ---------------------------------------------------------------------------
# 2. Margen de Popov-Belevitch-Hautus
# ---------------------------------------------------------------------------

"""
    pbh_controllability_margin(A, B) -> mu

Version cuantitativa del test de Popov-Belevitch-Hautus.

El test clasico dice que el par (A, B) es controlable si y solo si

    rank[ A - lambda I | B ] = n   para todo lambda en C

y basta comprobarlo en los eigenvalores de A: fuera del espectro, A - lambda I
ya tiene rango n por si sola. Sustituyendo el rango por el menor valor singular
se obtiene una medida continua:

    mu = min_{lambda en spec(A)} sigma_min([ A - lambda I | B ])

Por el teorema de Eckart-Young, sigma_min de una matriz es exactamente su
distancia en norma 2 al conjunto de matrices de rango deficiente. Es decir: mu
mide cuanto habria que perturbar el sistema para volverlo incontrolable. Un mu
pequeno significa que el sistema es controlable "por poco".

Nota: evaluar el minimo solo sobre spec(A) da una cota SUPERIOR de la verdadera
distancia a la incontrolabilidad de Eising,

    mu_E = min_{s en C} sigma_min([ A - s I | B ])

cuyo minimo puede caer fuera del espectro. Para comparar configuraciones la
version sobre el espectro es suficiente y cuesta solo n descomposiciones SVD.
"""
function pbh_controllability_margin(A, B)
    n = size(A, 1)
    Ac = Matrix{ComplexF64}(A)
    Bc = Matrix{ComplexF64}(B)
    Ic = Matrix{ComplexF64}(I, n, n)

    mu = Inf
    for lambda in eigvals(A)
        s_min = minimum(svdvals(hcat(Ac - lambda * Ic, Bc)))
        mu = min(mu, s_min)
    end
    return mu
end

"""
    pbh_controllability_margin_normalized(A, B) -> mu_hat

Margen PBH dividido por `opnorm([A B])`:

    mu_hat = mu_PBH / || [A  B] ||_2

La normalizacion vuelve la metrica adimensional y comparable entre sistemas de
distinta escala (distinta dimension, distintas unidades de masa y longitud).
Sin ella, un sistema con matrices de entradas grandes parece tener mas margen
solo por su escala.
"""
function pbh_controllability_margin_normalized(A, B)
    return pbh_controllability_margin(A, B) / opnorm(hcat(Matrix(A), Matrix(B)))
end

"""
    pbh_controllability_profile(A, B) -> Vector{NamedTuple}

Margen PBH desglosado por eigenvalor, ordenado de menor a mayor margen. El
primer elemento identifica el modo que esta mas cerca de perderse: es el cuello
de botella de la controlabilidad y suele ser el modo mas rapido.

Cada entrada tiene los campos `lambda` y `sigma_min`.
"""
function pbh_controllability_profile(A, B)
    n = size(A, 1)
    Ac = Matrix{ComplexF64}(A)
    Bc = Matrix{ComplexF64}(B)
    Ic = Matrix{ComplexF64}(I, n, n)

    profile = [(lambda=lambda,
                sigma_min=minimum(svdvals(hcat(Ac - lambda * Ic, Bc))))
               for lambda in eigvals(A)]
    return sort(profile, by=p -> p.sigma_min)
end

"""
    pbh_observability_margin(A, C) -> mu

Analogo dual del margen de controlabilidad. El par (A, C) es observable si y
solo si

    rank[ A - lambda I ; C ] = n   para todo lambda en C

(la matriz se apila VERTICALMENTE, no horizontalmente: es la transpuesta
estructural del caso de controlabilidad). La version cuantitativa es

    mu = min_{lambda en spec(A)} sigma_min([ A - lambda I ; C ])

Mide cuanto habria que perturbar el sistema para perder observabilidad, y sirve
para decidir que juego de sensores es suficiente.
"""
function pbh_observability_margin(A, C)
    n = size(A, 1)
    Ac = Matrix{ComplexF64}(A)
    Cc = Matrix{ComplexF64}(C)
    Ic = Matrix{ComplexF64}(I, n, n)

    mu = Inf
    for lambda in eigvals(A)
        s_min = minimum(svdvals(vcat(Ac - lambda * Ic, Cc)))
        mu = min(mu, s_min)
    end
    return mu
end

"""
    pbh_observability_margin_normalized(A, C) -> mu_hat

Margen PBH de observabilidad dividido por `opnorm([A; C])`.
"""
function pbh_observability_margin_normalized(A, C)
    return pbh_observability_margin(A, C) / opnorm(vcat(Matrix(A), Matrix(C)))
end

# ---------------------------------------------------------------------------
# 3. Condicionamiento de las matrices de Kalman
# ---------------------------------------------------------------------------

"""
    controllability_condition(A, B) -> cond

Numero de condicion (sigma_max / sigma_min) de la matriz de controlabilidad de
Kalman C = [B  AB  ...  A^(n-1)B].

Esta cantidad es el argumento experimental de por que el criterio de rango,
siendo un teorema correcto, es un mal algoritmo. Para la familia de pendulos de
N eslabones con longitud total 1 m y masa total 0.6 kg:

    N = 1  ->  4.8e1
    N = 2  ->  6.2e4
    N = 3  ->  2.2e8
    N = 4  ->  1.7e12
    N = 5  ->  3.0e16

Con N = 5 el condicionamiento alcanza 1/eps_maquina: el rango numerico de C es
ruido, y sin embargo el sistema sigue siendo controlable en aritmetica exacta.
A partir de N = 3 conviene reportar el margen PBH en lugar del rango.
"""
function controllability_condition(A, B)
    n = size(A, 1)
    C_ctrl = hcat([A^k * B for k in 0:n-1]...)
    s = svdvals(C_ctrl)
    return s[1] / s[end]
end

"""
    observability_condition(A, C) -> cond

Numero de condicion de la matriz de observabilidad O = [C; CA; ...; CA^(n-1)].
Sufre la misma degradacion que la de controlabilidad, por el mismo motivo: es
una base de Krylov.
"""
function observability_condition(A, C)
    n = size(A, 1)
    O_obs = vcat([C * A^k for k in 0:n-1]...)
    s = svdvals(O_obs)
    return s[1] / s[end]
end

# ---------------------------------------------------------------------------
# 4. Costo optimo y esfuerzo de control
# ---------------------------------------------------------------------------

"""
    optimal_cost(P, x0) -> J

Costo optimo alcanzable desde el estado inicial x0 bajo la ley LQR:

    J*(x0) = x0' P x0

donde P es la solucion de la CARE que devuelve Controller.solve_care. Es una
medida directa y sin ambiguedad de "que tan dificil es" estabilizar desde x0:
no depende de ninguna tolerancia ni de ningun criterio de exito arbitrario.

Nota: el gramiano de controlabilidad de horizonte infinito NO sirve aqui como
alternativa, porque A es inestable y la integral que lo define diverge. Si se
quisiera un gramiano habria que usar horizonte finito o el de lazo cerrado, con
A - BK, que si es Hurwitz.
"""
optimal_cost(P, x0) = dot(x0, P * x0)

"""
    nonsaturation_level(K, P, R, B, umax) -> c_star

Mayor conjunto de nivel de la funcion de Lyapunov V(x) = x' P x que cabe entero
dentro de la franja donde el actuador no satura, {x : |K x| <= umax}.

El maximo de una forma lineal sobre un elipsoide es su norma dual (ejercicio de
Cauchy-Schwarz):

    max_{x' P x <= c} |K x| = sqrt( c * K P^-1 K' )

Exigiendo que ese maximo no pase de umax se despeja

    c* = umax^2 / (K P^-1 K')

Dentro del elipsoide E = {x : x' P x <= c*} el control nunca satura y, para el
sistema lineal, Vdot = -x'(Q + K'RK)x < 0: el elipsoide es invariante y toda
trayectoria que empieza dentro converge al origen.

Como la ganancia optima es K = R^-1 B' P, se tiene

    K P^-1 K' = R^-1 B' P P^-1 P B R^-1 = R^-1 B' P B R^-1

que para R escalar es B' P B / R^2. Se usa esta forma para no invertir P.

ADVERTENCIA sobre el alcance: c* garantiza no saturacion y convergencia DEL
SISTEMA LINEAL. Para el sistema no lineal hace falta ademas acotar el residuo
de la linealizacion, lo que da un c menor. Presentar c* como si fuera la region
de atraccion no lineal seria un error.
"""
function nonsaturation_level(K, P, R, B, umax)
    r = R isa AbstractMatrix ? R[1, 1] : R
    BPB = (B' * P * B)
    q = (BPB isa AbstractMatrix ? BPB[1, 1] : BPB) / r^2
    return umax^2 / q
end

"""
    max_axis_angle(P, c, idx) -> value

Maximo valor que puede tomar la componente `idx` del estado, sobre su propio
eje (con el resto del estado en cero), sin salirse del elipsoide
{x : x' P x <= c}:

    x = v * e_idx  =>  v^2 P[idx,idx] <= c  =>  |v| <= sqrt(c / P[idx,idx])

Particularizado a idx = 3 (el angulo del primer eslabon en el orden intercalado
del proyecto) responde: desde que inclinacion inicial se garantiza que el
control no satura.
"""
max_axis_angle(P, c, idx) = sqrt(c / P[idx, idx])

# ---------------------------------------------------------------------------
# 5. Discretizacion y criterio de Schur
# ---------------------------------------------------------------------------

"""
    discretize_zoh(A, B, h) -> (Ad, Bd)

Discretizacion EXACTA del sistema continuo con retenedor de orden cero (ZOH),
es decir suponiendo que la entrada se mantiene constante en cada intervalo de
muestreo de duracion h. La solucion exacta de la EDO lineal da

    Ad = exp(A h),    Bd = (integral_0^h exp(A tau) dtau) B

Ambas se obtienen de UNA sola exponencial matricial mediante el truco de
Van Loan: se construye la matriz aumentada de tamano (n+m) x (n+m)

    Maug = [ A  B ]
           [ 0  0 ]

y su exponencial tiene la estructura por bloques

    exp(Maug * h) = [ Ad  Bd ]
                    [ 0   I  ]

de modo que Ad y Bd se leen directamente de los bloques superiores. Esto evita
integrar numericamente y evita invertir A, que aqui seria singular (el carro
aporta un eigenvalor cero).
"""
function discretize_zoh(A, B, h)
    n = size(A, 1)
    m = size(B, 2)

    Maug = zeros(n + m, n + m)
    Maug[1:n, 1:n] = A
    Maug[1:n, n+1:n+m] = B

    Phi = exp(Maug * h)
    Ad = Phi[1:n, 1:n]
    Bd = Phi[1:n, n+1:n+m]
    return (Ad=Ad, Bd=Bd)
end

"""
    is_schur(M; tol=1.0) -> Bool

Criterio de estabilidad en TIEMPO DISCRETO: una matriz es de Schur si todos sus
eigenvalores estan dentro del circulo unitario, |lambda_i| < 1. Es el analogo
discreto de la condicion de Hurwitz (Re lambda_i < 0) del caso continuo, y sale
de que la solucion discreta es x_k = M^k x_0, que decae si y solo si el radio
espectral es menor que uno.
"""
is_schur(M; tol=1.0) = maximum(abs.(eigvals(M))) < tol

# ---------------------------------------------------------------------------
# 6. Reporte
# ---------------------------------------------------------------------------

"""
    print_metrics_summary(A, B, C, K, P, R; umax=50.0, x0=nothing, angle_idx=3, label="")

Imprime todas las metricas del sistema con el formato de banner del proyecto.
`x0` por defecto es el estado con la componente `angle_idx` igual a 0.1 rad,
que es la condicion inicial estandar de los pipelines.
"""
function print_metrics_summary(A, B, C, K, P, R;
                               umax=50.0, x0=nothing, angle_idx=3, label="")
    n = size(A, 1)
    if x0 === nothing
        x0 = zeros(n)
        x0[angle_idx] = 0.1
    end

    println("=" ^ 60)
    println("  METRICAS DE VIABILIDAD" * (isempty(label) ? "" : ": $label"))
    println("=" ^ 60)

    dm = dominant_mode(A)
    println("\n  Escala de tiempo:")
    @printf("  lambda_max          = %+.4f s^-1\n", dm.lambda_max)
    @printf("  tau_caida = 1/l_max = %.1f ms\n", 1000 * dm.tau_fall)
    @printf("  modos inestables    = %d\n", dm.n_unstable)

    println("\n  Controlabilidad:")
    mu_c = pbh_controllability_margin(A, B)
    @printf("  margen PBH          = %.4e\n", mu_c)
    @printf("  margen PBH normaliz = %.4e\n",
            pbh_controllability_margin_normalized(A, B))
    @printf("  cond(Ctrb)          = %.4e", controllability_condition(A, B))
    println(controllability_condition(A, B) > 1e12 ?
            "   [rango numerico NO confiable]" : "")
    worst = pbh_controllability_profile(A, B)[1]
    @printf("  modo mas debil      = %+.4f %+.4fi\n",
            real(worst.lambda), imag(worst.lambda))

    println("\n  Observabilidad:")
    @printf("  margen PBH          = %.4e\n", pbh_observability_margin(A, C))
    @printf("  margen PBH normaliz = %.4e\n",
            pbh_observability_margin_normalized(A, C))
    @printf("  cond(Obsv)          = %.4e\n", observability_condition(A, C))

    println("\n  Esfuerzo de control:")
    @printf("  norm(K)             = %.4f\n", norm(K))
    @printf("  J*(x0)              = %.4f\n", optimal_cost(P, x0))

    println("\n  Elipsoide de no saturacion (umax = $(umax) N):")
    c_star = nonsaturation_level(K, P, R, B, umax)
    theta_g = max_axis_angle(P, c_star, angle_idx)
    @printf("  c*                  = %.4f\n", c_star)
    @printf("  theta_%d garantizado = %.4f rad = %.1f grados\n",
            angle_idx, theta_g, rad2deg(theta_g))

    println("\n" * "=" ^ 60)
end

end # module Metrics
