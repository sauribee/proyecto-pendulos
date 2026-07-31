### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 00000001-0000-4000-8000-000000000000
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))

    using LinearAlgebra
    using Printf
    using DifferentialEquations
    using CairoMakie
    using PlutoUI
end

# ╔═╡ 00000002-0000-4000-8000-000000000000
md"""
# Pendulo Invertido Triple -- Exploracion Interactiva

**Proyecto de Algebra Lineal Aplicada -- Configuracion III**

Un carro con **tres** eslabones articulados en serie, modelados como barras
uniformes con inercia. Cuatro grados de libertad y **un solo actuador**: la
subactuacion pasa de ``2\!:\!1`` en la Configuracion I a ``4\!:\!1`` aqui. El
estado tiene dimension 8:

``\mathbf{x} = [x,\ \dot{x},\ \theta_1,\ \omega_1,\ \theta_2,\ \omega_2,\ \theta_3,\ \omega_3]^T``

Este notebook recorre todo lo implementado en la segunda entrega: el modelo
generico de ``N`` eslabones, las metricas de viabilidad (margen PBH,
condicionamiento, elipsoide de no saturacion, periodo de muestreo), el
observador de Luenberger por dualidad y los barridos del atlas de operabilidad.

Mueve los sliders y todas las celdas dependientes se recalculan automaticamente.

---
"""

# ╔═╡ 00000003-0000-4000-8000-000000000000
begin
    include(joinpath(@__DIR__, "..", "src", "model_nlink.jl"))
    include(joinpath(@__DIR__, "..", "src", "model_triple.jl"))
    include(joinpath(@__DIR__, "..", "src", "linearization.jl"))
    include(joinpath(@__DIR__, "..", "src", "controller.jl"))
    include(joinpath(@__DIR__, "..", "src", "metrics.jl"))
    include(joinpath(@__DIR__, "..", "src", "observer.jl"))
    include(joinpath(@__DIR__, "..", "src", "sweep.jl"))

    using .ModelNLink
    using .ModelTriple
    using .Linearization
    using .Controller
    using .Metrics
    using .Observer
    using .Sweep

    md"""
    Modulos cargados: `ModelNLink`, `ModelTriple`, `Linearization`,
    `Controller`, `Metrics`, `Observer`, `Sweep`.
    """
end

# ╔═╡ 00000004-0000-4000-8000-000000000000
md"""
## 1. Parametros del sistema

Los tres eslabones son **barras uniformes**: ``a_j = \ell_j/2`` e
``I_j = \tfrac{1}{12}m_j\ell_j^2``. Los valores por defecto reparten la misma
masa total (0.6 kg) y la misma longitud total (1 m) de la Configuracion II en
tres barras, para que la comparacion sea limpia.
"""

# ╔═╡ 00000005-0000-4000-8000-000000000000
md"**Masa del carro M [kg]:**"

# ╔═╡ 00000006-0000-4000-8000-000000000000
@bind M_val PlutoUI.Slider(0.05:0.05:5.0, default=1.0, show_value=true)

# ╔═╡ 00000007-0000-4000-8000-000000000000
md"**Masa del eslabon 1 (inferior) m1 [kg]:**"

# ╔═╡ 00000008-0000-4000-8000-000000000000
@bind m1_val PlutoUI.Slider(0.02:0.02:1.0, default=0.2, show_value=true)

# ╔═╡ 00000009-0000-4000-8000-000000000000
md"**Masa del eslabon 2 (intermedio) m2 [kg]:**"

# ╔═╡ 00000010-0000-4000-8000-000000000000
@bind m2_val PlutoUI.Slider(0.02:0.02:1.0, default=0.2, show_value=true)

# ╔═╡ 00000011-0000-4000-8000-000000000000
md"**Masa del eslabon 3 (superior) m3 [kg]:**"

# ╔═╡ 00000012-0000-4000-8000-000000000000
@bind m3_val PlutoUI.Slider(0.02:0.02:1.0, default=0.2, show_value=true)

# ╔═╡ 00000013-0000-4000-8000-000000000000
md"**Longitud del eslabon 1 l1 [m]:**"

# ╔═╡ 00000014-0000-4000-8000-000000000000
@bind l1_val PlutoUI.Slider(0.05:0.01:1.5, default=1/3, show_value=true)

# ╔═╡ 00000015-0000-4000-8000-000000000000
md"**Longitud del eslabon 2 l2 [m]:**"

# ╔═╡ 00000016-0000-4000-8000-000000000000
@bind l2_val PlutoUI.Slider(0.05:0.01:1.5, default=1/3, show_value=true)

# ╔═╡ 00000017-0000-4000-8000-000000000000
md"**Longitud del eslabon 3 l3 [m]** -- el eslabon superior *corto* es el peligroso:"

# ╔═╡ 00000018-0000-4000-8000-000000000000
@bind l3_val PlutoUI.Slider(0.05:0.01:1.5, default=1/3, show_value=true)

# ╔═╡ 00000019-0000-4000-8000-000000000000
md"**Gravedad g [m/s^2]:**"

# ╔═╡ 00000020-0000-4000-8000-000000000000
@bind g_val PlutoUI.Slider(1.0:0.1:20.0, default=9.81, show_value=true)

# ╔═╡ 00000021-0000-4000-8000-000000000000
md"**Friccion viscosa del carro b [N s/m]:**"

# ╔═╡ 00000022-0000-4000-8000-000000000000
@bind b_val PlutoUI.Slider(0.0:0.02:1.0, default=0.1, show_value=true)

# ╔═╡ 00000023-0000-4000-8000-000000000000
begin
    params = uniform_rods(M_val, [m1_val, m2_val, m3_val],
                          [l1_val, l2_val, l3_val]; g=g_val, b=b_val)
    bc = link_couplings(params)

    md"""
    ### Parametros activos

    | | Eslabon 1 | Eslabon 2 | Eslabon 3 |
    |---|---|---|---|
    | masa ``m_j`` [kg] | $(m1_val) | $(m2_val) | $(m3_val) |
    | longitud ``\ell_j`` [m] | $(round(l1_val, digits=4)) | $(round(l2_val, digits=4)) | $(round(l3_val, digits=4)) |
    | brazo ``a_j = \ell_j/2`` [m] | $(round(params.a[1], digits=4)) | $(round(params.a[2], digits=4)) | $(round(params.a[3], digits=4)) |
    | inercia ``I_j`` [kg m^2] | $(round(params.Il[1], digits=6)) | $(round(params.Il[2], digits=6)) | $(round(params.Il[3], digits=6)) |

    Carro: ``M`` = $(M_val) kg, ``b`` = $(b_val) N s/m, ``g`` = $(g_val) m/s^2.

    Masa total del pendulo: **$(round(sum(params.m), digits=3)) kg** ·
    Longitud total: **$(round(sum(params.l), digits=3)) m**
    """
end

# ╔═╡ 00000024-0000-4000-8000-000000000000
md"""
## 2. Las tres agrupaciones que organizan el modelo

Al expandir las energias aparecen tres agrupaciones de parametros. Con ellas la
matriz de masa y el lado derecho se escriben de una vez para cualquier ``N``:

``\beta_j = m_j a_j + \ell_j\!\!\sum_{i>j}m_i``  (momento estatico)

``J_j = I_j + m_j a_j^2 + \ell_j^2\!\!\sum_{i>j}m_i``  (inercia efectiva)

``\Gamma_{jk} = \ell_{\min(j,k)}\,\beta_{\max(j,k)}``  (acoplamiento inercial)

El mismo ``\beta_j`` que acopla el carro con el eslabon ``j`` en la energia
**cinetica** es el que pesa la gravedad sobre ese eslabon en la **potencial**,
``V = g\sum_j \beta_j\cos\theta_j``. Por eso inestabilidad y controlabilidad
estan ligadas: el canal por el que la gravedad desestabiliza es el mismo por el
que el carro puede actuar.
"""

# ╔═╡ 00000025-0000-4000-8000-000000000000
begin
    Mq0 = mass_matrix(params, zeros(3))
    Markdown.parse("""
    | ``j`` | ``\\beta_j`` | ``J_j`` |
    |---|---|---|
    | 1 | $(round(bc.beta[1], digits=5)) | $(round(bc.J[1], digits=5)) |
    | 2 | $(round(bc.beta[2], digits=5)) | $(round(bc.J[2], digits=5)) |
    | 3 | $(round(bc.beta[3], digits=5)) | $(round(bc.J[3], digits=5)) |

    Matriz de masa evaluada en el equilibrio, ``\\mathbf{M}_0 = \\mathbf{M}(0)``
    (simetrica definida positiva; por eso las aceleraciones se obtienen por
    Cholesky y no por una LU generica):
    """)
end

# ╔═╡ 00000026-0000-4000-8000-000000000000
round.(Mq0, digits=5)

# ╔═╡ 00000027-0000-4000-8000-000000000000
md"""
## 3. Linealizacion -- Espacio de estados

Evaluando en ``\theta_j = 0`` y ``\dot q = 0`` queda un sistema lineal con
**matriz de masa constante**, ``\mathbf{M}_0\ddot q = \mathbf{G}_0 q -
\mathbf{F}_0\dot q + \mathbf{e}_1 u``, de donde salen los tres bloques
``A_{cc} = \mathbf{M}_0^{-1}\mathbf{G}_0``,
``A_{cv} = -\mathbf{M}_0^{-1}\mathbf{F}_0`` y
``B_c = \mathbf{M}_0^{-1}\mathbf{e}_1``. El jacobiano es **analitico**, no
numerico.
"""

# ╔═╡ 00000028-0000-4000-8000-000000000000
begin
    ss = linearize_system_triple(params)
    md"Sistema linealizado: A es $(size(ss.A, 1))x$(size(ss.A, 2)), B es $(size(ss.B, 1))x$(size(ss.B, 2)), C es $(size(ss.C, 1))x$(size(ss.C, 2))"
end

# ╔═╡ 00000029-0000-4000-8000-000000000000
md"### Matriz A (8x8):"

# ╔═╡ 00000030-0000-4000-8000-000000000000
round.(ss.A, digits=3)

# ╔═╡ 00000031-0000-4000-8000-000000000000
md"### Matriz B (8x1):"

# ╔═╡ 00000032-0000-4000-8000-000000000000
round.(ss.B, digits=4)

# ╔═╡ 00000033-0000-4000-8000-000000000000
md"### Matriz C (4x8) -- medimos la posicion del carro y los tres angulos:"

# ╔═╡ 00000034-0000-4000-8000-000000000000
Int.(ss.C)

# ╔═╡ 00000035-0000-4000-8000-000000000000
md"""
## 4. Estabilidad -- Eigenvalores de A

**Prediccion estructural.** Sin friccion, buscar ``q(t) = v e^{\lambda t}`` da
el problema generalizado ``\mathbf{G}_0 v = \lambda^2\mathbf{M}_0 v``: los
eigenvalores de ``A`` son las *raices cuadradas* de los de
``\mathbf{M}_0^{-1}\mathbf{G}_0`` y vienen en pares ``\pm\lambda``. Como
``\mathbf{G}_0`` tiene una fila y una columna nulas (el carro no tiene fuerza
recuperadora), aparece un cero doble:

``\sigma(A) = \{\pm\lambda_1, \pm\lambda_2, \pm\lambda_3,\ 0,\ 0\}``, con
**tres modos inestables**.

Sube ``b`` por encima de cero y veras como la friccion del carro rompe esa
simetria y desdobla el cero.
"""

# ╔═╡ 00000036-0000-4000-8000-000000000000
begin
    dm = dominant_mode(ss.A)
    eig_sorted = sort(ss.eigenvalues, by=real, rev=true)
    eig_rows = map(enumerate(eig_sorted)) do (i, lam)
        rp = round(real(lam), digits=4)
        ip = round(imag(lam), digits=4)
        estab = rp > 1e-9 ? "**INESTABLE**" : (rp < -1e-9 ? "estable" : "marginal")
        val = abs(ip) < 1e-9 ? "$rp" : "$rp + $(ip)i"
        "| ``\\lambda_{$i}`` | $val | $estab |"
    end
    Markdown.parse("""
    | Eigenvalor | Valor | Estabilidad |
    |---|---|---|
    $(join(eig_rows, "\n"))

    **Modos inestables:** $(dm.n_unstable) ·
    **``\\lambda_{\\max}``** = $(round(dm.lambda_max, digits=4)) s⁻¹ ·
    **``\\tau_{caida} = 1/\\lambda_{\\max}``** = $(round(1000 * dm.tau_fall, digits=1)) ms

    El controlador debe actuar en una fraccion de ``\\tau_{caida}``: esa es la
    escala de tiempo del problema.
    """)
end

# ╔═╡ 00000037-0000-4000-8000-000000000000
md"""
## 5. Controlabilidad y observabilidad (Kalman)

El criterio clasico: el sistema es completamente controlable y observable si
``\text{rank}\,\mathcal{C} = \text{rank}\,\mathcal{O} = 8``.
"""

# ╔═╡ 00000038-0000-4000-8000-000000000000
begin
    ctrb = check_controllability(ss)
    obsv = check_observability(ss)
    Markdown.parse("""
    | Propiedad | Rango | Requerido | Resultado |
    |---|---|---|---|
    | Controlabilidad | $(ctrb.rank) | $(ctrb.required_rank) | $(ctrb.is_controllable ? "**CONTROLABLE**" : "NO CONTROLABLE") |
    | Observabilidad | $(obsv.rank) | $(obsv.required_rank) | $(obsv.is_observable ? "**OBSERVABLE**" : "NO OBSERVABLE") |
    """)
end

# ╔═╡ 00000039-0000-4000-8000-000000000000
md"""
## 6. Cuando el rango deja de servir -- el margen PBH

Aqui esta el resultado central de la segunda entrega. El criterio de Kalman es
un teorema correcto pero un **mal algoritmo**: ``\mathcal{C}`` es una base de
Krylov, notoriamente mal condicionada. Con tres eslabones
``\text{cond}\,\mathcal{C} \approx 2\times10^8``, y con cinco alcanza
``1/\varepsilon_{maq}``: el rango numerico deja de significar algo.

La version **cuantitativa** del test de Popov-Belevitch-Hautus sustituye el
rango por el menor valor singular:

``\mu_{PBH} = \min_{\lambda\in\sigma(A)}\ \sigma_{\min}[\,A - \lambda I \mid B\,]``

Por el teorema de Eckart-Young eso es exactamente la distancia en norma 2 al
conjunto de matrices de rango deficiente: mide **cuanto habria que perturbar el
sistema para volverlo incontrolable**. Normalizado,
``\hat\mu = \mu_{PBH}/\lVert[A\ B]\rVert_2``.
"""

# ╔═╡ 00000040-0000-4000-8000-000000000000
begin
    mu_c = pbh_controllability_margin_normalized(ss.A, ss.B)
    mu_o = pbh_observability_margin_normalized(ss.A, ss.C)
    cond_c = controllability_condition(ss.A, ss.B)
    cond_o = observability_condition(ss.A, ss.C)
    peor = pbh_controllability_profile(ss.A, ss.B)[1]
    confiable = cond_c < 1e12

    Markdown.parse("""
    | Metrica | Valor |
    |---|---|
    | Margen PBH normalizado (control) | **$(@sprintf("%.4e", mu_c))** |
    | Margen PBH normalizado (observacion) | $(@sprintf("%.4e", mu_o)) |
    | ``\\text{cond}\\,\\mathcal{C}`` | $(@sprintf("%.3e", cond_c)) $(confiable ? "" : "— **rango numerico NO confiable**") |
    | ``\\text{cond}\\,\\mathcal{O}`` | $(@sprintf("%.3e", cond_o)) |
    | Modo mas debil (cuello de botella) | ``\\lambda`` = $(round(real(peor.lambda), digits=3)) $(abs(imag(peor.lambda)) < 1e-9 ? "" : "+ $(round(imag(peor.lambda), digits=3))i") |

    El minimo del perfil PBH se alcanza en el modo inestable mas rapido: la
    direccion mas cercana a perderse es exactamente la que hay que controlar.
    """)
end

# ╔═╡ 00000041-0000-4000-8000-000000000000
md"""
## 7. Diseno del controlador LQR

Minimiza ``J = \int_0^\infty (\mathbf{x}^T Q\mathbf{x} + u^T R u)\,dt``
resolviendo la ecuacion algebraica de Riccati por eigendescomposicion de la
matriz **hamiltoniana**. La ganancia optima es ``K = R^{-1}B^T P``.
"""

# ╔═╡ 00000042-0000-4000-8000-000000000000
md"**Peso Q en la posicion del carro:**"

# ╔═╡ 00000043-0000-4000-8000-000000000000
@bind q_x PlutoUI.Slider(0.1:0.1:50.0, default=1.0, show_value=true)

# ╔═╡ 00000044-0000-4000-8000-000000000000
md"**Peso Q en los tres angulos (igual para los tres):**"

# ╔═╡ 00000045-0000-4000-8000-000000000000
@bind q_th PlutoUI.Slider(1.0:1.0:200.0, default=10.0, show_value=true)

# ╔═╡ 00000046-0000-4000-8000-000000000000
md"**Peso R del esfuerzo de control** (escala logaritmica; el nominal es 0.1):"

# ╔═╡ 00000047-0000-4000-8000-000000000000
@bind log_r PlutoUI.Slider(-3.0:0.1:1.0, default=-1.0, show_value=true)

# ╔═╡ 00000048-0000-4000-8000-000000000000
begin
    R_val = 10.0^log_r
    Qm = diagm([q_x, 0.0, q_th, 0.0, q_th, 0.0, q_th, 0.0])
    Rm = reshape([R_val], 1, 1)
    lqr = design_lqr(ss.A, ss.B, Qm, Rm)

    etiquetas = ["pos", "vel", "th1", "w1", "th2", "w2", "th3", "w3"]
    k_rows = ["| ``K_{$(etiquetas[i])}`` | $(round(lqr.K[i], digits=3)) |" for i in 1:8]
    all_stable = all(real.(lqr.eigenvalues_cl) .< 0)

    Markdown.parse("""
    ``R`` = $(round(R_val, digits=5)) · ``\\lVert K\\rVert_2`` = **$(round(norm(lqr.K), digits=2))**

    | Componente | Ganancia |
    |---|---|
    $(join(k_rows, "\n"))

    **Lazo cerrado:** $(all_stable ? "todos los polos estables" : "AUN INESTABLE")
    """)
end

# ╔═╡ 00000049-0000-4000-8000-000000000000
begin
    polos_sorted = sort(lqr.eigenvalues_cl, by=real)
    polo_rows = map(enumerate(polos_sorted)) do (i, lam)
        rp = round(real(lam), digits=3)
        ip = round(imag(lam), digits=3)
        "| ``\\mu_{$i}`` | $(abs(ip) < 1e-9 ? "$rp" : "$rp + $(ip)i") |"
    end
    Markdown.parse("""
    ### Polos de lazo cerrado ``\\sigma(A - BK)``

    | Polo | Valor |
    |---|---|
    $(join(polo_rows, "\n"))

    Con ``R`` pequeno el LQR tiende a **reflejar** los polos inestables al
    semiplano izquierdo: compara el polo mas rapido de aqui con
    ``+\\lambda_{\\max} =`` $(round(dm.lambda_max, digits=2)) del lazo abierto.
    Es la propiedad del *regulador barato*. Baja el slider de ``R`` y observa
    que ``\\lVert K\\rVert`` apenas se mueve: la inestabilidad de la planta ya
    fija la escala de la ganancia.
    """)
end

# ╔═╡ 00000050-0000-4000-8000-000000000000
md"""
## 8. El elipsoide de no saturacion -- una garantia dura

Con ``V(\mathbf{x}) = \mathbf{x}^T P\mathbf{x}`` la funcion de Lyapunov del LQR,
el mayor conjunto de nivel contenido en la franja de no saturacion
``\{|K\mathbf{x}| \le u_{\max}\}`` sale de que **el maximo de una forma lineal
sobre un elipsoide es su norma dual** (Cauchy-Schwarz):

``c^\star = \dfrac{u_{\max}^2}{K P^{-1} K^T}``

Dentro de ``\mathcal{E}_{c^\star}`` el control **nunca satura** y, para el
sistema lineal, ``\dot V < 0``: el elipsoide es invariante. Sobre el eje del
primer angulo la garantia es ``|\theta_0| \le \sqrt{c^\star/P_{33}}``.

> **Honestidad sobre el alcance.** Esto certifica no saturacion mas convergencia
> del sistema **lineal**. Para el no lineal haria falta acotar el residuo de la
> linealizacion. Presentar ``c^\star`` como la region de atraccion no lineal
> seria un error.
"""

# ╔═╡ 00000051-0000-4000-8000-000000000000
md"**Limite del actuador u_max [N]:**"

# ╔═╡ 00000052-0000-4000-8000-000000000000
@bind umax_val PlutoUI.Slider(5.0:5.0:300.0, default=150.0, show_value=true)

# ╔═╡ 00000053-0000-4000-8000-000000000000
begin
    c_star = nonsaturation_level(lqr.K, lqr.P, Rm, ss.B, umax_val)
    th_guar = max_axis_angle(lqr.P, c_star, 3)
    x0_ref = zeros(8); x0_ref[3] = 0.1
    J_star = optimal_cost(lqr.P, x0_ref)

    Markdown.parse("""
    | Cantidad | Valor |
    |---|---|
    | ``P_{33}`` | $(round(lqr.P[3,3], digits=4)) |
    | ``c^\\star`` | $(round(c_star, digits=4)) |
    | ``\\theta_0`` garantizado sin saturar | **$(round(th_guar, digits=4)) rad = $(round(rad2deg(th_guar), digits=2))°** |
    | ``J^\\star`` desde ``\\theta_1 = 0.1`` rad | $(round(J_star, digits=4)) |

    Como ``c^\\star \\propto u_{\\max}^2``, el angulo garantizado crece
    **linealmente** con el actuador. Pero la frontera real se estanca: mas alla
    de cierto punto lo que limita no es el motor sino la no linealidad.
    """)
end

# ╔═╡ 00000054-0000-4000-8000-000000000000
md"""
## 9. Control digital -- periodo de muestreo maximo

Todo controlador real muestrea. Con retenedor de orden cero, la solucion exacta
da ``A_d = e^{Ah}`` y ``B_d = \left(\int_0^h e^{A\tau}d\tau\right)B``, y **ambas
salen de una sola exponencial matricial** por el truco de Van Loan: se exponencia
la matriz aumentada ``\begin{pmatrix}A & B\\ 0 & 0\end{pmatrix}`` y los bloques
superiores son ``A_d`` y ``B_d``. Esto evita integrar y evita invertir ``A``,
que aqui es singular.

La estabilidad en tiempo discreto pide que ``A_d - B_d K`` sea de **Schur**
(todos los eigenvalores dentro del circulo unitario): el analogo discreto de
Hurwitz.
"""

# ╔═╡ 00000055-0000-4000-8000-000000000000
begin
    h_max = max_sampling_period(ss.A, ss.B, lqr.K)
    producto = h_max * dm.lambda_max
    Markdown.parse("""
    | Cantidad | Valor |
    |---|---|
    | ``h_{\\max}`` | **$(round(1000 * h_max, digits=1)) ms** (``\\ge`` $(round(1/h_max, digits=1)) Hz) |
    | ``1/\\lambda_{\\max}`` | $(round(1000 * dm.tau_fall, digits=1)) ms |
    | ``h_{\\max}\\,\\lambda_{\\max}`` | $(round(producto, digits=3)) |

    El margen relativo respecto al tiempo de caida se estrecha al crecer la
    complejidad (0.83 en el simple, 0.66 en el doble, 0.58 en el triple): no
    basta con escalar la frecuencia proporcionalmente al modo dominante.
    """)
end

# ╔═╡ 00000056-0000-4000-8000-000000000000
md"**Periodo de muestreo h a probar [ms]** -- mueve hasta cruzar ``h_{\max}``:"

# ╔═╡ 00000057-0000-4000-8000-000000000000
@bind h_test_ms PlutoUI.Slider(1.0:1.0:300.0, default=20.0, show_value=true)

# ╔═╡ 00000058-0000-4000-8000-000000000000
let
    h = h_test_ms / 1000
    d = discretize_zoh(ss.A, ss.B, h)
    lam_d = eigvals(d.Ad - d.Bd * lqr.K)
    radio = maximum(abs.(lam_d))
    es_schur = is_schur(d.Ad - d.Bd * lqr.K)

    fig = Figure(size=(640, 480))
    ax = Axis(fig[1, 1],
        title=@sprintf("Polos discretos con h = %.0f ms  (radio espectral %.4f)",
                       h_test_ms, radio),
        xlabel="Re(z)", ylabel="Im(z)", aspect=DataAspect())

    tt = range(0, 2pi, length=200)
    lines!(ax, cos.(tt), sin.(tt), color=:gray50, linewidth=1.5, linestyle=:dash)
    hlines!(ax, [0.0], color=:gray85, linewidth=0.5)
    vlines!(ax, [0.0], color=:gray85, linewidth=0.5)
    scatter!(ax, real.(lam_d), imag.(lam_d), markersize=13,
        color=es_schur ? RGBf(0.17, 0.45, 0.71) : RGBf(0.84, 0.15, 0.16))
    text!(ax, 0.0, -1.35,
        text=es_schur ? "DENTRO del circulo: estable" : "FUERA: el lazo digital diverge",
        align=(:center, :center), fontsize=14,
        color=es_schur ? :gray25 : RGBf(0.7, 0.1, 0.1))
    limits!(ax, -1.5, 1.5, -1.5, 1.5)
    fig
end

# ╔═╡ 00000059-0000-4000-8000-000000000000
md"""
## 10. Simulacion no lineal en lazo cerrado

Se integra el modelo **no lineal completo** con ``u = -K\mathbf{x}`` y
saturacion, partiendo de una inclinacion en los tres eslabones, y se compara con
la caida libre. El triple sin control da varias vueltas completas en tres
segundos.
"""

# ╔═╡ 00000060-0000-4000-8000-000000000000
md"**Angulo inicial de los tres eslabones [rad]:**"

# ╔═╡ 00000061-0000-4000-8000-000000000000
@bind th0 PlutoUI.Slider(0.01:0.01:0.35, default=0.1, show_value=true)

# ╔═╡ 00000062-0000-4000-8000-000000000000
begin
    x0 = [0.0, 0.0, th0, 0.0, th0, 0.0, th0, 0.0]

    sol_free = solve(ODEProblem(nonlinear_eom_triple!, copy(x0), (0.0, 3.0),
                     (params=params, F=0.0)), Tsit5(), saveat=0.005)
    sol_lqr = solve(ODEProblem(closed_loop_eom_triple!, copy(x0), (0.0, 10.0),
                    (params=params, K=lqr.K, saturate=umax_val)),
                    Tsit5(), saveat=0.005)

    th_lqr = [[u[2j+1] for u in sol_lqr.u] for j in 1:3]
    x_lqr = [u[1] for u in sol_lqr.u]
    F_lqr = [clamp(-dot(lqr.K[1, :], u), -umax_val, umax_val) for u in sol_lqr.u]
    pico = maximum(abs.(F_lqr))

    function t_asent(t, y; tol=0.02)
        idx = findlast(>(tol), abs.(y))
        return idx === nothing ? t[1] : t[min(idx + 1, length(t))]
    end
    ts = [t_asent(sol_lqr.t, th_lqr[j]) for j in 1:3]
    converge = all(abs(sol_lqr.u[end][2j+1]) < 0.02 for j in 1:3)

    Markdown.parse("""
    | Cantidad | Valor |
    |---|---|
    | ``t_s(\\theta_1)`` | $(round(ts[1], digits=2)) s |
    | ``t_s(\\theta_2)`` | $(round(ts[2], digits=2)) s |
    | ``t_s(\\theta_3)`` | $(round(ts[3], digits=2)) s |
    | Fuerza pico | $(round(pico, digits=2)) N de $(umax_val) N disponibles ($(round(100*pico/umax_val, digits=1)) %) |
    | Excursion del carro | $(round(maximum(abs.(x_lqr)), digits=3)) m |
    | Se satura el actuador | $(pico >= umax_val - 1e-6 ? "**SI**" : "no") |

    **$(converge ? "El controlador recupera el equilibrio." : "NO converge: fuera de la region de atraccion.")**
    """)
end

# ╔═╡ 00000063-0000-4000-8000-000000000000
begin
    COL_LQR = RGBf(0.17, 0.45, 0.71)
    COL_FREE = RGBf(0.84, 0.15, 0.16)
    COL_2 = RGBf(0.85, 0.50, 0.11)
    COL_3 = RGBf(0.18, 0.49, 0.20)
    COL_OBS = RGBf(0.55, 0.24, 0.62)
    md"Paleta configurada"
end

# ╔═╡ 00000064-0000-4000-8000-000000000000
let
    fig = Figure(size=(920, 660))

    ax1 = Axis(fig[1, 1], title="Angulos con LQR",
        xlabel="t [s]", ylabel="theta [rad]")
    for (j, col) in enumerate((COL_LQR, COL_2, COL_3))
        lines!(ax1, sol_lqr.t, th_lqr[j], color=col, linewidth=2, label="theta$j")
    end
    hlines!(ax1, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="Sin control: el triple se desploma",
        xlabel="t [s]", ylabel="theta [rad]")
    for (j, col) in enumerate((COL_LQR, COL_2, COL_3))
        lines!(ax2, sol_free.t, [u[2j+1] for u in sol_free.u],
            color=col, linewidth=2, label="theta$j")
    end
    hlines!(ax2, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)
    axislegend(ax2, position=:lt)

    ax3 = Axis(fig[2, 1], title="Posicion del carro",
        xlabel="t [s]", ylabel="x [m]")
    lines!(ax3, sol_lqr.t, x_lqr, color=COL_LQR, linewidth=2)
    hlines!(ax3, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)

    ax4 = Axis(fig[2, 2],
        title=@sprintf("Fuerza de control (pico %.2f N)", pico),
        xlabel="t [s]", ylabel="F [N]")
    lines!(ax4, sol_lqr.t, F_lqr, color=COL_3, linewidth=2)
    hlines!(ax4, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)
    if pico > 0.5 * umax_val
        hlines!(ax4, [umax_val, -umax_val], color=:orange,
            linestyle=:dot, linewidth=1.5)
    end

    Label(fig[0, :], "Configuracion III: sin control vs LQR",
        fontsize=15, font=:bold)
    fig
end

# ╔═╡ 00000065-0000-4000-8000-000000000000
let
    fig = Figure(size=(660, 450))
    ax = Axis(fig[1, 1], title="Eigenvalores: lazo abierto vs cerrado",
        xlabel="Re(lambda) [1/s]", ylabel="Im(lambda) [1/s]")
    vspan!(ax, 0.0, 40.0, color=(COL_FREE, 0.07))
    vlines!(ax, [0.0], color=:gray50, linewidth=1.2, linestyle=:dash)
    hlines!(ax, [0.0], color=:gray85, linewidth=0.5)
    scatter!(ax, real.(ss.eigenvalues), imag.(ss.eigenvalues),
        color=COL_FREE, markersize=14, marker=:xcross, label="Lazo abierto (A)")
    scatter!(ax, real.(lqr.eigenvalues_cl), imag.(lqr.eigenvalues_cl),
        color=COL_LQR, markersize=13, label="Lazo cerrado (A - BK)")
    axislegend(ax, position=:lt)
    fig
end

# ╔═╡ 00000066-0000-4000-8000-000000000000
md"""
## 11. El observador de Luenberger -- controlar sin medir todo

Hasta aqui el control ha sido ``u = -K\mathbf{x}``, con el estado **completo**.
En un sistema real eso no se mide: nadie pone un sensor de velocidad angular en
cada articulacion. Como ``\text{rank}\,\mathcal{O} = 8``, el estado se puede
reconstruir:

``\dot{\hat{\mathbf{x}}} = A\hat{\mathbf{x}} + Bu + L(\mathbf{y} - C\hat{\mathbf{x}})``,
con el error ``\mathbf{e} = \mathbf{x} - \hat{\mathbf{x}}`` obedeciendo
``\dot{\mathbf{e}} = (A - LC)\mathbf{e}``.

**La dualidad de Kalman** convierte el problema nuevo en uno ya resuelto:
``L = \text{lqr}(A^T, C^T, Q_o, R_o)^T``. No hace falta escribir un solo
algoritmo nuevo.
"""

# ╔═╡ 00000067-0000-4000-8000-000000000000
md"**Confianza en los sensores** -- ``Q_o = 10^{\alpha} I``; mas alto = observador mas rapido:"

# ╔═╡ 00000068-0000-4000-8000-000000000000
@bind log_qo PlutoUI.Slider(0.0:0.5:5.0, default=2.0, show_value=true)

# ╔═╡ 00000069-0000-4000-8000-000000000000
begin
    Qo = (10.0^log_qo) * Matrix(I, 8, 8)
    Ro = Matrix(I, 4, 4)
    obs_res = design_observer(ss.A, ss.C, Qo, Ro)

    sep = check_separation(ss.A, ss.B, ss.C, lqr.K, obs_res.L)
    obs_stable = all(real.(obs_res.eigenvalues_obs) .< 0)
    obs_pol = sort(obs_res.eigenvalues_obs, by=real)
    pol_rows = map(enumerate(obs_pol)) do (i, lam)
        rp = round(real(lam), digits=3)
        ip = round(imag(lam), digits=3)
        "| ``\\nu_{$i}`` | $(abs(ip) < 1e-9 ? "$rp" : "$rp + $(ip)i") |"
    end

    Markdown.parse("""
    ### Polos del observador ``\\sigma(A - LC)``

    | Polo | Valor |
    |---|---|
    $(join(pol_rows, "\n"))

    Observador $(obs_stable ? "**estable**: el error de estimacion decae" : "INESTABLE") ·
    ``\\lVert L\\rVert`` = $(round(norm(obs_res.L), digits=2))
    """)
end

# ╔═╡ 00000070-0000-4000-8000-000000000000
md"""
### El principio de separacion

En coordenadas ``(\mathbf{x}, \mathbf{e})`` la matriz de lazo cerrado es
**triangular por bloques**:

``\begin{pmatrix}A - BK & BK\\ 0 & A - LC\end{pmatrix}
\ \Longrightarrow\ \sigma = \sigma(A-BK)\ \cup\ \sigma(A-LC)``

y por tanto controlador y observador se disenan **por separado**. Es un
resultado de control que se demuestra en dos lineas de algebra lineal: el
espectro de una triangular por bloques es la union de los espectros diagonales.

Notese que el bloque ``(1,2)`` **no** es nulo: el error de estimacion si perturba
a la planta. Lo que la triangularidad garantiza es que esa perturbacion no puede
desestabilizar, porque el error decae por su cuenta.
"""

# ╔═╡ 00000071-0000-4000-8000-000000000000
begin
    Maug = augmented_closed_loop(ss.A, ss.B, ss.C, lqr.K, obs_res.L)
    bloque_nulo = all(Maug[9:16, 1:8] .== 0)
    Markdown.parse("""
    | Verificacion | Resultado |
    |---|---|
    | Bloque inferior izquierdo nulo (triangular) | $(bloque_nulo ? "SI" : "no") |
    | Error maximo al emparejar espectros | $(@sprintf("%.3e", sep.max_error)) |
    | Principio de separacion | $(sep.holds ? "**SE CUMPLE**" : "NO se cumple") |
    """)
end

# ╔═╡ 00000072-0000-4000-8000-000000000000
md"""
### Que sensores hacen falta de verdad

De los ``2^4 - 1 = 15`` subconjuntos no vacios de
``\{x, \theta_1, \theta_2, \theta_3\}``, ¿cuales conservan la observabilidad?
La respuesta es una particion limpia, y el **margen PBH de observabilidad**
distingue lo que el rango no puede: *observable* no es lo mismo que
*practicable*.
"""

# ╔═╡ 00000073-0000-4000-8000-000000000000
begin
    tabla_sensores = sensor_subset_analysis(ss.A, ss.C, ["x", "th1", "th2", "th3"])
    filas_s = map(tabla_sensores) do f
        "| $(join(f.labels, " + ")) | $(f.rank)/8 | $(f.observable ? "SI" : "**NO**") | $(@sprintf("%.3e", f.margin_normalized)) |"
    end
    n_obs = count(f -> f.observable, tabla_sensores)
    Markdown.parse("""
    | Sensores | rank ``\\mathcal{O}`` | Observable | Margen PBH normalizado |
    |---|---|---|---|
    $(join(filas_s, "\n"))

    **$n_obs de $(length(tabla_sensores)) subconjuntos son observables.** Todos
    los que incluyen la posicion del carro lo son; ninguno que la excluya. La
    razon es geometrica: midiendo solo angulos, la dinamica angular es invariante
    a trasladar el carro, de modo que la posicion absoluta es inobservable.
    Falta exactamente una direccion.
    """)
end

# ╔═╡ 00000074-0000-4000-8000-000000000000
md"""
### Simulacion con observador

La planta arranca inclinada y el estimador arranca en **cero**: el maximo
desconocimiento posible sobre el estado inicial. El control es ``u = -K\hat{\mathbf{x}}``,
no ``-K\mathbf{x}``: ese es justamente el punto.
"""

# ╔═╡ 00000075-0000-4000-8000-000000000000
md"**Angulo inicial para la prueba del observador [rad]** (solo ``\theta_1``):"

# ╔═╡ 00000076-0000-4000-8000-000000000000
@bind th0_obs PlutoUI.Slider(0.01:0.01:0.30, default=0.1, show_value=true)

# ╔═╡ 00000077-0000-4000-8000-000000000000
begin
    z0 = zeros(16)
    z0[3] = th0_obs
    p_obs = (params=params, eom! = nonlinear_eom_triple!,
             A=ss.A, B=ss.B, C=ss.C, K=lqr.K, L=obs_res.L, saturate=umax_val)
    sol_obs = solve(ODEProblem(observer_eom!, z0, (0.0, 10.0), p_obs),
                    Tsit5(), saveat=0.005)

    err_obs = [norm(u[1:8] .- u[9:16]) for u in sol_obs.u]
    u_obs = [clamp(-dot(lqr.K[1, :], u[9:16]), -umax_val, umax_val)
             for u in sol_obs.u]
    i_1pct = findfirst(e -> e < 0.01 * max(err_obs[1], eps()), err_obs)
    conv_obs = all(abs(sol_obs.u[end][2j+1]) < 0.02 for j in 1:3)

    Markdown.parse("""
    | Cantidad | Valor |
    |---|---|
    | ``\\lVert e(0)\\rVert`` | $(round(err_obs[1], digits=4)) |
    | ``\\lVert e(T)\\rVert`` | $(@sprintf("%.3e", err_obs[end])) |
    | Cae al 1 % del inicial | $(i_1pct === nothing ? "no alcanza" : "t = $(round(sol_obs.t[i_1pct], digits=2)) s") |
    | Fuerza pico con observador | $(round(maximum(abs.(u_obs)), digits=2)) N |
    | Converge la planta | $(conv_obs ? "**SI**" : "**NO**") |

    Compara la fuerza pico con la del estado completo
    ($(round(pico, digits=2)) N): durante el primer transitorio el controlador
    actua sobre una estimacion equivocada, y eso se paga.
    """)
end

# ╔═╡ 00000078-0000-4000-8000-000000000000
let
    fig = Figure(size=(920, 640))

    ax1 = Axis(fig[1, 1], title="theta1 -- variable MEDIDA",
        xlabel="t [s]", ylabel="theta1 [rad]")
    lines!(ax1, sol_obs.t, [u[3] for u in sol_obs.u], color=:black,
        linewidth=2.5, label="real")
    lines!(ax1, sol_obs.t, [u[11] for u in sol_obs.u], color=COL_2,
        linewidth=2, linestyle=:dash, label="estimado")
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="omega1 -- variable NO medida, reconstruida",
        xlabel="t [s]", ylabel="omega1 [rad/s]")
    lines!(ax2, sol_obs.t, [u[4] for u in sol_obs.u], color=:black,
        linewidth=2.5, label="real")
    lines!(ax2, sol_obs.t, [u[12] for u in sol_obs.u], color=COL_2,
        linewidth=2, linestyle=:dash, label="estimado")
    axislegend(ax2, position=:rb)

    ax3 = Axis(fig[2, 1], title="Error de estimacion ||x - xhat||",
        xlabel="t [s]", ylabel="||e||", yscale=log10)
    lines!(ax3, sol_obs.t, max.(err_obs, 1e-14), color=COL_OBS, linewidth=2)

    ax4 = Axis(fig[2, 2], title="Fuerza: u = -K xhat vs u = -K x",
        xlabel="t [s]", ylabel="F [N]")
    lines!(ax4, sol_obs.t, u_obs, color=COL_OBS, linewidth=2,
        label="con observador")
    lines!(ax4, sol_lqr.t, F_lqr, color=COL_LQR, linewidth=2,
        linestyle=:dash, label="estado completo")
    axislegend(ax4, position=:rt)

    Label(fig[0, :], "Observador de Luenberger: el estimador arranca en cero",
        fontsize=15, font=:bold)
    fig
end

# ╔═╡ 00000079-0000-4000-8000-000000000000
md"""
> **El observador no sale gratis.** El principio de separacion garantiza
> estabilidad del sistema **lineal** aumentado. No dice nada sobre la region de
> atraccion del **no lineal**, y ahi estimar en vez de medir se paga: durante el
> transitorio inicial el control actua sobre una estimacion equivocada y saca a
> la planta del regimen donde la linealizacion vale. Sube el slider de
> ``\theta_1`` y compara desde donde converge con observador frente a con estado
> completo.

---
"""

# ╔═╡ 00000080-0000-4000-8000-000000000000
md"""
## 12. Barridos: el atlas de operabilidad en vivo

El atlas responde tres preguntas: ¿donde funciona bien el controlador?, ¿donde
apenas?, ¿donde es imposible? Y separa **dos tipos de imposibilidad**:

- **Estructural:** el par ``(A,B)`` deja de ser controlable, o se acerca tanto
  que ningun ``K`` razonable sirve. Es una propiedad *del sistema*; se mide con
  el margen PBH.
- **Practica:** el sistema es controlable, pero el actuador satura, el riel se
  acaba o el muestreo es muy lento. Es una propiedad *de la implementacion*; se
  mide integrando el modelo no lineal y biseccionando.

Elige que parametro barrer:
"""

# ╔═╡ 00000081-0000-4000-8000-000000000000
@bind spec_sel PlutoUI.Select([
    "M" => "Masa del carro M",
    "l3" => "Longitud del eslabon superior l3",
    "m3" => "Masa del eslabon superior m3",
    "l1" => "Longitud del eslabon inferior l1",
], default="M")

# ╔═╡ 00000082-0000-4000-8000-000000000000
begin
    spec_map = Dict("M" => (:M, 10 .^ range(log10(0.03), log10(20.0), length=40)),
                    "l3" => ((:l, 3), range(0.02, 2.0, length=40)),
                    "m3" => ((:m, 3), 10 .^ range(log10(0.02), log10(3.2), length=40)),
                    "l1" => ((:l, 1), range(0.05, 2.0, length=40)))
    spec_sw, vals_sw = spec_map[spec_sel]
    res_sw = sweep_1d(linearize_system_nlink, params, spec_sw, collect(vals_sw))

    pbh_sw = [r.pbh_normalized for r in res_sw]
    normk_sw = [r.norm_K for r in res_sw]
    lam_sw = [r.lambda_max for r in res_sw]
    cstar_sw = [r.c_star for r in res_sw]
    i_best = argmax(pbh_sw)
    interior = 1 < i_best < length(vals_sw)

    Markdown.parse("""
    Barrido de **$(spec_sel)** en $(length(vals_sw)) puntos.

    | | Valor |
    |---|---|
    | Mejor margen PBH en | $(spec_sel) = $(round(vals_sw[i_best], digits=4)) |
    | Margen ahi | $(@sprintf("%.4e", pbh_sw[i_best])) |
    | ¿Optimo interior? | $(interior ? "**SI** — el optimo no esta en un extremo del dominio" : "no, esta en un extremo") |
    | Rango de ``\\lVert K\\rVert`` | $(round(minimum(normk_sw), digits=1)) a $(round(maximum(normk_sw), digits=1)) |
    """)
end

# ╔═╡ 00000083-0000-4000-8000-000000000000
let
    escala_log = spec_sel in ("M", "m3")
    fig = Figure(size=(920, 620))

    ax1 = Axis(fig[1, 1], title="Margen PBH normalizado",
        xlabel=spec_sel, ylabel="PBH norm", yscale=log10,
        xscale=escala_log ? log10 : identity)
    lines!(ax1, collect(vals_sw), pbh_sw, color=COL_LQR, linewidth=2.5)
    scatter!(ax1, [vals_sw[i_best]], [pbh_sw[i_best]], color=COL_FREE,
        markersize=14, marker=:star5)

    ax2 = Axis(fig[1, 2], title="Esfuerzo de control ||K||",
        xlabel=spec_sel, ylabel="norm(K)", yscale=log10,
        xscale=escala_log ? log10 : identity)
    lines!(ax2, collect(vals_sw), normk_sw, color=COL_2, linewidth=2.5)

    ax3 = Axis(fig[2, 1], title="Modo inestable dominante",
        xlabel=spec_sel, ylabel="lambda_max [1/s]",
        xscale=escala_log ? log10 : identity)
    lines!(ax3, collect(vals_sw), lam_sw, color=COL_FREE, linewidth=2.5)

    ax4 = Axis(fig[2, 2], title="Elipsoide de no saturacion c*",
        xlabel=spec_sel, ylabel="c*",
        xscale=escala_log ? log10 : identity)
    lines!(ax4, collect(vals_sw), cstar_sw, color=COL_3, linewidth=2.5)

    Label(fig[0, :], "Barrido de $(spec_sel): la estrella marca el mejor margen PBH",
        fontsize=15, font=:bold)
    fig
end

# ╔═╡ 00000084-0000-4000-8000-000000000000
md"""
> Al barrer la masa del carro veras el hallazgo mas interesante del atlas: el
> margen PBH **no es monotono**, tiene un maximo interior. Un carro muy pesado se
> mueve poco ante la misma fuerza y pierde autoridad de control; uno muy liviano
> se mueve mucho pero degrada el acoplamiento con los eslabones. Y fijate en que
> ``c^\star`` y ``\hat\mu`` apuntan en **direcciones opuestas**: no existe *la*
> mejor configuracion sin decir antes que se optimiza.
"""

# ╔═╡ 00000085-0000-4000-8000-000000000000
md"""
### La frontera practica: hasta que angulo se recupera

Bisecciona sobre el modelo **no lineal** con saturacion. El criterio de exito
es: todos los angulos por debajo de 0.02 rad al final, el carro dentro del riel
y ningun angulo pasando de ``\pi/2`` en el camino.
"""

# ╔═╡ 00000086-0000-4000-8000-000000000000
md"**Longitud del riel [m]:**"

# ╔═╡ 00000087-0000-4000-8000-000000000000
@bind x_rail PlutoUI.Slider(0.5:0.1:5.0, default=1.5, show_value=true)

# ╔═╡ 00000088-0000-4000-8000-000000000000
begin
    th_max_riel = max_recoverable_angle(params, lqr.K; umax=umax_val, x_rail=x_rail)
    th_max_libre = max_recoverable_angle(params, lqr.K; umax=umax_val, x_rail=1e6)
    manda_riel = abs(th_max_libre - th_max_riel) > 0.01
    cota_lin = max_axis_angle(lqr.P, c_star, 3)

    Markdown.parse("""
    | Frontera | Valor |
    |---|---|
    | Cota elipsoidal ``\\sqrt{c^\\star/P_{33}}`` (garantia lineal) | $(round(cota_lin, digits=4)) rad |
    | Frontera real con riel de $(x_rail) m | **$(round(th_max_riel, digits=4)) rad** |
    | Frontera real sin restriccion de riel | $(round(th_max_libre, digits=4)) rad |
    | Restriccion activa | **$(manda_riel ? "EL RIEL" : "LA SATURACION")** |

    La cota elipsoidal hay que compararla contra la columna **sin riel**: solo
    certifica no saturacion mas convergencia lineal, no dice nada sobre la
    carrera del carro. Con ``u_{\\max}`` grande la cota crece linealmente y
    acaba prometiendo mas de lo que el sistema no lineal cumple.
    """)
end

# ╔═╡ 00000089-0000-4000-8000-000000000000
md"""
## 13. Animacion del pendulo triple

Visualizacion mecanica del carro y los tres eslabones. Como los angulos son
**absolutos** respecto a la vertical, la posicion del extremo del eslabon ``j``
es la suma acumulada ``x_j = x + \sum_{i\le j}\ell_i\sin\theta_i``.
"""

# ╔═╡ 00000090-0000-4000-8000-000000000000
begin
    function link_points(state, p)
        N = length(p.l)
        pts = Vector{Tuple{Float64,Float64}}(undef, N + 1)
        pts[1] = (state[1], 0.0)
        for j in 1:N
            th = state[2j+1]
            pts[j+1] = (pts[j][1] + p.l[j] * sin(th),
                        pts[j][2] + p.l[j] * cos(th))
        end
        return pts
    end

    function draw_triple_frame!(ax, state, p; time_label=nothing)
        pts = link_points(state, p)
        N = length(p.l)
        reach = sum(p.l)

        empty!(ax)
        lines!(ax, [-4.0, 4.0], [0.0, 0.0], color=:gray70,
            linewidth=1.5, linestyle=:dash)
        poly!(ax, Rect(state[1] - 0.15, -0.075, 0.3, 0.15),
            color=(COL_LQR, 0.85), strokecolor=:gray20, strokewidth=1.5)
        for dx in (-0.1, 0.1)
            scatter!(ax, [state[1] + dx], [-0.09], markersize=9, color=:gray45)
        end
        tonos = [:gray20, :gray40, :gray55]
        for j in 1:N
            lines!(ax, [pts[j][1], pts[j+1][1]], [pts[j][2], pts[j+1][2]],
                color=tonos[min(j, 3)], linewidth=3.5)
        end
        for j in 1:N
            es_extremo = j == N
            scatter!(ax, [pts[j+1][1]], [pts[j+1][2]],
                markersize=es_extremo ? 18 : 12,
                color=es_extremo ? :red : :orange,
                strokecolor=:black, strokewidth=1)
        end
        if !isnothing(time_label)
            text!(ax, 3.9, reach + 0.35, text=time_label, fontsize=12,
                color=:gray40, align=(:right, :top))
        end
    end
    md"Funciones de dibujo listas: `link_points`, `draw_triple_frame!`"
end

# ╔═╡ 00000091-0000-4000-8000-000000000000
md"**Tiempo t [s] -- trayectoria controlada (LQR):**"

# ╔═╡ 00000092-0000-4000-8000-000000000000
@bind t_frame PlutoUI.Slider(range(0.0, 10.0, length=501), default=0.0, show_value=true)

# ╔═╡ 00000093-0000-4000-8000-000000000000
let
    reach = sum(params.l)
    lim_x = max(1.2, maximum(abs.(x_lqr)) + reach * 0.6)
    fig = Figure(size=(760, 580))

    ax = Axis(fig[1, 1], title="Pendulo triple -- trayectoria LQR",
        xlabel="x [m]", ylabel="y [m]", aspect=DataAspect(),
        limits=(-lim_x, lim_x, -0.35, reach + 0.45))
    draw_triple_frame!(ax, sol_lqr(t_frame), params,
        time_label=@sprintf("t = %.2f s", t_frame))

    ax2 = Axis(fig[2, 1], xlabel="t [s]", ylabel="angulo [rad]", height=130,
        title="theta1 (azul), theta2 (naranja), theta3 (verde)")
    for (j, col) in enumerate((COL_LQR, COL_2, COL_3))
        lines!(ax2, sol_lqr.t, th_lqr[j], color=col, linewidth=1.5)
    end
    hlines!(ax2, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)
    vlines!(ax2, [t_frame], color=:red, linewidth=1.8)

    rowgap!(fig.layout, 12)
    fig
end

# ╔═╡ 00000094-0000-4000-8000-000000000000
md"""
## 14. Exportar animacion (GIF)

Marca la casilla para generar `figures/12_triple_exploracion.gif` (los primeros
6 s de la trayectoria controlada). Desmarcala para no regenerar el archivo en
cada cambio de parametro.
"""

# ╔═╡ 00000095-0000-4000-8000-000000000000
@bind do_export PlutoUI.CheckBox(default=false)

# ╔═╡ 00000096-0000-4000-8000-000000000000
begin
    if do_export
        fig_dir = joinpath(@__DIR__, "..", "figures")
        mkpath(fig_dir)
        gif_path = joinpath(fig_dir, "12_triple_exploracion.gif")

        reach_e = sum(params.l)
        lim_e = max(1.2, maximum(abs.(x_lqr)) + reach_e * 0.6)
        fig_e = Figure(size=(660, 540))
        ax_e = Axis(fig_e[1, 1], title="Pendulo triple -- LQR",
            xlabel="x [m]", ylabel="y [m]", aspect=DataAspect(),
            limits=(-lim_e, lim_e, -0.35, reach_e + 0.45))

        ts_e = range(0.0, 6.0, step=1/30)
        record(fig_e, gif_path, eachindex(ts_e); framerate=30) do i
            draw_triple_frame!(ax_e, sol_lqr(ts_e[i]), params,
                time_label=@sprintf("t = %.2f s", ts_e[i]))
        end
        md"GIF exportado en: `$gif_path`"
    else
        md"Marca la casilla de arriba para exportar el GIF."
    end
end

# ╔═╡ 00000097-0000-4000-8000-000000000000
md"""
---

## Resumen: que mirar en este notebook

1. **Sube ``m_3`` o acorta ``\ell_3``** y observa como se hunde el margen PBH:
   el eslabon superior corto y pesado es el peligroso. Es el efecto de la escoba
   larga en la palma de la mano, cuantificado.
2. **Baja ``M`` hasta ~0.13 kg** en el barrido: el margen PBH tiene un maximo
   interior que no esta en ningun extremo del dominio.
3. **Mueve ``R`` cuatro ordenes de magnitud**: ``\lVert K\rVert`` apenas cambia
   un factor 2. La planta es tan inestable que cualquier ``K`` que estabilice
   necesita ganancia grande.
4. **Sube ``h``** por encima de ``h_{\max}`` y mira los polos discretos salir
   del circulo unitario.
5. **Sube ``\theta_1``** en la prueba del observador hasta que deje de
   converger, y compara con la region de atraccion del estado completo.

El mismo algebra lineal escala de ``\mathbb{R}^4`` a ``\mathbb{R}^8`` sin cambios
conceptuales. Lo que **no** escala es el condicionamiento: ahi esta la leccion.
"""

# ╔═╡ Cell order:
# ╠═00000001-0000-4000-8000-000000000000
# ╟─00000002-0000-4000-8000-000000000000
# ╠═00000003-0000-4000-8000-000000000000
# ╟─00000004-0000-4000-8000-000000000000
# ╟─00000005-0000-4000-8000-000000000000
# ╠═00000006-0000-4000-8000-000000000000
# ╟─00000007-0000-4000-8000-000000000000
# ╠═00000008-0000-4000-8000-000000000000
# ╟─00000009-0000-4000-8000-000000000000
# ╠═00000010-0000-4000-8000-000000000000
# ╟─00000011-0000-4000-8000-000000000000
# ╠═00000012-0000-4000-8000-000000000000
# ╟─00000013-0000-4000-8000-000000000000
# ╠═00000014-0000-4000-8000-000000000000
# ╟─00000015-0000-4000-8000-000000000000
# ╠═00000016-0000-4000-8000-000000000000
# ╟─00000017-0000-4000-8000-000000000000
# ╠═00000018-0000-4000-8000-000000000000
# ╟─00000019-0000-4000-8000-000000000000
# ╠═00000020-0000-4000-8000-000000000000
# ╟─00000021-0000-4000-8000-000000000000
# ╠═00000022-0000-4000-8000-000000000000
# ╠═00000023-0000-4000-8000-000000000000
# ╟─00000024-0000-4000-8000-000000000000
# ╠═00000025-0000-4000-8000-000000000000
# ╠═00000026-0000-4000-8000-000000000000
# ╟─00000027-0000-4000-8000-000000000000
# ╠═00000028-0000-4000-8000-000000000000
# ╟─00000029-0000-4000-8000-000000000000
# ╠═00000030-0000-4000-8000-000000000000
# ╟─00000031-0000-4000-8000-000000000000
# ╠═00000032-0000-4000-8000-000000000000
# ╟─00000033-0000-4000-8000-000000000000
# ╠═00000034-0000-4000-8000-000000000000
# ╟─00000035-0000-4000-8000-000000000000
# ╠═00000036-0000-4000-8000-000000000000
# ╟─00000037-0000-4000-8000-000000000000
# ╠═00000038-0000-4000-8000-000000000000
# ╟─00000039-0000-4000-8000-000000000000
# ╠═00000040-0000-4000-8000-000000000000
# ╟─00000041-0000-4000-8000-000000000000
# ╟─00000042-0000-4000-8000-000000000000
# ╠═00000043-0000-4000-8000-000000000000
# ╟─00000044-0000-4000-8000-000000000000
# ╠═00000045-0000-4000-8000-000000000000
# ╟─00000046-0000-4000-8000-000000000000
# ╠═00000047-0000-4000-8000-000000000000
# ╠═00000048-0000-4000-8000-000000000000
# ╠═00000049-0000-4000-8000-000000000000
# ╟─00000050-0000-4000-8000-000000000000
# ╟─00000051-0000-4000-8000-000000000000
# ╠═00000052-0000-4000-8000-000000000000
# ╠═00000053-0000-4000-8000-000000000000
# ╟─00000054-0000-4000-8000-000000000000
# ╠═00000055-0000-4000-8000-000000000000
# ╟─00000056-0000-4000-8000-000000000000
# ╠═00000057-0000-4000-8000-000000000000
# ╠═00000058-0000-4000-8000-000000000000
# ╟─00000059-0000-4000-8000-000000000000
# ╟─00000060-0000-4000-8000-000000000000
# ╠═00000061-0000-4000-8000-000000000000
# ╠═00000062-0000-4000-8000-000000000000
# ╠═00000063-0000-4000-8000-000000000000
# ╠═00000064-0000-4000-8000-000000000000
# ╠═00000065-0000-4000-8000-000000000000
# ╟─00000066-0000-4000-8000-000000000000
# ╟─00000067-0000-4000-8000-000000000000
# ╠═00000068-0000-4000-8000-000000000000
# ╠═00000069-0000-4000-8000-000000000000
# ╟─00000070-0000-4000-8000-000000000000
# ╠═00000071-0000-4000-8000-000000000000
# ╟─00000072-0000-4000-8000-000000000000
# ╠═00000073-0000-4000-8000-000000000000
# ╟─00000074-0000-4000-8000-000000000000
# ╟─00000075-0000-4000-8000-000000000000
# ╠═00000076-0000-4000-8000-000000000000
# ╠═00000077-0000-4000-8000-000000000000
# ╠═00000078-0000-4000-8000-000000000000
# ╟─00000079-0000-4000-8000-000000000000
# ╟─00000080-0000-4000-8000-000000000000
# ╠═00000081-0000-4000-8000-000000000000
# ╠═00000082-0000-4000-8000-000000000000
# ╠═00000083-0000-4000-8000-000000000000
# ╟─00000084-0000-4000-8000-000000000000
# ╟─00000085-0000-4000-8000-000000000000
# ╟─00000086-0000-4000-8000-000000000000
# ╠═00000087-0000-4000-8000-000000000000
# ╠═00000088-0000-4000-8000-000000000000
# ╟─00000089-0000-4000-8000-000000000000
# ╠═00000090-0000-4000-8000-000000000000
# ╟─00000091-0000-4000-8000-000000000000
# ╠═00000092-0000-4000-8000-000000000000
# ╠═00000093-0000-4000-8000-000000000000
# ╟─00000094-0000-4000-8000-000000000000
# ╠═00000095-0000-4000-8000-000000000000
# ╠═00000096-0000-4000-8000-000000000000
# ╟─00000097-0000-4000-8000-000000000000
