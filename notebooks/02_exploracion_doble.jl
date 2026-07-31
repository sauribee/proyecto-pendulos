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
# Pendulo Invertido Doble -- Exploracion Interactiva

**Proyecto de Algebra Lineal Aplicada -- Configuracion II**

Un carro con dos eslabones en serie (masas puntuales ``m_1`` en la articulacion
intermedia y ``m_2`` en el extremo). El estado tiene dimension 6:

``\mathbf{x} = [x,\ \dot{x},\ \theta_1,\ \dot{\theta}_1,\ \theta_2,\ \dot{\theta}_2]^T``

Mueve los sliders y todas las celdas dependientes se recalculan automaticamente.

---
"""

# ╔═╡ 00000003-0000-4000-8000-000000000000
begin
    include(joinpath(@__DIR__, "..", "src", "model_double.jl"))
    include(joinpath(@__DIR__, "..", "src", "model_nlink.jl"))
    include(joinpath(@__DIR__, "..", "src", "linearization.jl"))
    include(joinpath(@__DIR__, "..", "src", "controller.jl"))
    include(joinpath(@__DIR__, "..", "src", "metrics.jl"))
    include(joinpath(@__DIR__, "..", "src", "observer.jl"))
    include(joinpath(@__DIR__, "..", "src", "sweep.jl"))

    using .ModelDouble
    using .ModelNLink
    using .Linearization
    using .Controller
    using .Metrics
    using .Observer
    using .Sweep

    md"""
    Modulos cargados: `ModelDouble`, `ModelNLink`, `Linearization`,
    `Controller`, `Metrics`, `Observer`, `Sweep`.
    """
end

# ╔═╡ 00000004-0000-4000-8000-000000000000
md"""
## 1. Parametros del sistema

El pendulo doble se modela sin friccion (igual que en el informe). Ajusta las
masas, las longitudes y la gravedad.
"""

# ╔═╡ 00000005-0000-4000-8000-000000000000
md"**Masa del carro M [kg]:**"

# ╔═╡ 00000006-0000-4000-8000-000000000000
@bind M_val PlutoUI.Slider(0.2:0.1:5.0, default=1.0, show_value=true)

# ╔═╡ 00000007-0000-4000-8000-000000000000
md"**Masa de la articulacion intermedia m1 [kg]:**"

# ╔═╡ 00000008-0000-4000-8000-000000000000
@bind m1_val PlutoUI.Slider(0.05:0.05:2.0, default=0.3, show_value=true)

# ╔═╡ 00000009-0000-4000-8000-000000000000
md"**Masa del extremo superior m2 [kg]:**"

# ╔═╡ 00000010-0000-4000-8000-000000000000
@bind m2_val PlutoUI.Slider(0.05:0.05:2.0, default=0.3, show_value=true)

# ╔═╡ 00000011-0000-4000-8000-000000000000
md"**Longitud del eslabon inferior L1 [m]:**"

# ╔═╡ 00000012-0000-4000-8000-000000000000
@bind L1_val PlutoUI.Slider(0.2:0.1:2.0, default=0.5, show_value=true)

# ╔═╡ 00000013-0000-4000-8000-000000000000
md"**Longitud del eslabon superior L2 [m]:**"

# ╔═╡ 00000014-0000-4000-8000-000000000000
@bind L2_val PlutoUI.Slider(0.2:0.1:2.0, default=0.5, show_value=true)

# ╔═╡ 00000015-0000-4000-8000-000000000000
md"**Gravedad g [m/s^2]:**"

# ╔═╡ 00000016-0000-4000-8000-000000000000
@bind g_val PlutoUI.Slider(1.0:0.1:20.0, default=9.81, show_value=true)

# ╔═╡ 00000017-0000-4000-8000-000000000000
begin
    params = SystemParamsDouble(M=M_val, m1=m1_val, m2=m2_val,
                                L1=L1_val, L2=L2_val, g=g_val)
    md"""
    ### Parametros activos
    - Carro: M = $(M_val) kg
    - Eslabon inferior: m1 = $(m1_val) kg, L1 = $(L1_val) m
    - Eslabon superior: m2 = $(m2_val) kg, L2 = $(L2_val) m
    - Gravedad: g = $(g_val) m/s^2
    """
end

# ╔═╡ 00000018-0000-4000-8000-000000000000
md"""
## 2. Linealizacion -- Espacio de estados

Linealizando las ecuaciones de Euler-Lagrange alrededor del equilibrio superior
``\theta_1 = \theta_2 = 0`` se obtiene ``\dot{\mathbf{x}} = A\mathbf{x} + B u``,
con salida ``\mathbf{y} = C\mathbf{x}``.
"""

# ╔═╡ 00000019-0000-4000-8000-000000000000
begin
    ss = linearize_system_double(params)
    md"Sistema linealizado (estado de dimension 6)"
end

# ╔═╡ 00000020-0000-4000-8000-000000000000
md"### Matriz A (6x6):"

# ╔═╡ 00000021-0000-4000-8000-000000000000
round.(ss.A, digits=3)

# ╔═╡ 00000022-0000-4000-8000-000000000000
md"### Matriz B (6x1):"

# ╔═╡ 00000023-0000-4000-8000-000000000000
round.(ss.B, digits=3)

# ╔═╡ 00000024-0000-4000-8000-000000000000
md"### Matriz C (3x6) -- medimos posicion del carro y ambos angulos:"

# ╔═╡ 00000025-0000-4000-8000-000000000000
ss.C

# ╔═╡ 00000026-0000-4000-8000-000000000000
md"""
## 3. Estabilidad -- Eigenvalores de A

El pendulo doble tiene **dos** modos inestables (dos eigenvalores con parte real
positiva): por eso es mas exigente de controlar que el simple.
"""

# ╔═╡ 00000027-0000-4000-8000-000000000000
begin
    eig_rows = map(enumerate(ss.eigenvalues)) do (i, lam)
        rp = round(real(lam), digits=4)
        ip = round(imag(lam), digits=4)
        estab = rp > 1e-9 ? "INESTABLE" : (rp < -1e-9 ? "estable" : "marginal")
        val = abs(ip) < 1e-9 ? "$rp" : "$rp + $(ip)i"
        "| lambda_$i | $val | $estab |"
    end
    n_unstable = count(>(1e-9), real.(ss.eigenvalues))
    Markdown.parse("""
    | Eigenvalor | Valor | Estabilidad |
    |---|---|---|
    $(join(eig_rows, "\n"))

    **Modos inestables:** $n_unstable
    """)
end

# ╔═╡ 00000028-0000-4000-8000-000000000000
md"""
## 4. Controlabilidad y observabilidad (Kalman)

Si ``\text{rank}(\mathcal{C}) = \text{rank}(\mathcal{O}) = 6``, el sistema es
completamente controlable y observable, y se puede disenar el control.
"""

# ╔═╡ 00000029-0000-4000-8000-000000000000
begin
    ctrl = check_controllability(ss)
    obs = check_observability(ss)
    Markdown.parse("""
    | Propiedad | Rango | Requerido | Resultado |
    |---|---|---|---|
    | Controlabilidad | $(ctrl.rank) | $(ctrl.required_rank) | $(ctrl.is_controllable ? "CONTROLABLE" : "NO") |
    | Observabilidad | $(obs.rank) | $(obs.required_rank) | $(obs.is_observable ? "OBSERVABLE" : "NO") |
    """)
end

# ╔═╡ 00000030-0000-4000-8000-000000000000
md"""
## 5. Diseno del controlador LQR

Minimiza ``J = \int_0^\infty (\mathbf{x}^T Q \mathbf{x} + u^T R u)\, dt``.
La matriz ``Q`` penaliza posicion y ambos angulos; ``R`` penaliza el esfuerzo
de control. Ajusta los pesos y observa como cambian ``K`` y los polos.
"""

# ╔═╡ 00000031-0000-4000-8000-000000000000
md"**Peso Q en posicion del carro:**"

# ╔═╡ 00000032-0000-4000-8000-000000000000
@bind q_x PlutoUI.Slider(0.1:0.1:50.0, default=1.0, show_value=true)

# ╔═╡ 00000033-0000-4000-8000-000000000000
md"**Peso Q en angulo del eslabon 1 (theta1):**"

# ╔═╡ 00000034-0000-4000-8000-000000000000
@bind q_th1 PlutoUI.Slider(0.1:0.5:100.0, default=10.0, show_value=true)

# ╔═╡ 00000035-0000-4000-8000-000000000000
md"**Peso Q en angulo del eslabon 2 (theta2):**"

# ╔═╡ 00000036-0000-4000-8000-000000000000
@bind q_th2 PlutoUI.Slider(0.1:0.5:100.0, default=10.0, show_value=true)

# ╔═╡ 00000037-0000-4000-8000-000000000000
md"**Peso R en esfuerzo de control:**"

# ╔═╡ 00000038-0000-4000-8000-000000000000
@bind r_val PlutoUI.Slider(0.01:0.01:5.0, default=0.1, show_value=true)

# ╔═╡ 00000039-0000-4000-8000-000000000000
begin
    Q_mat = diagm([q_x, 0.0, q_th1, 0.0, q_th2, 0.0])
    R_mat = reshape([r_val], 1, 1)
    lqr = design_lqr(ss.A, ss.B, Q_mat, R_mat)

    labels6 = ["pos", "vel", "th1", "w1", "th2", "w2"]
    k_rows = ["| K_$(labels6[i]) | $(round(lqr.K[i], digits=4)) |" for i in 1:6]
    cl_rows = map(enumerate(lqr.eigenvalues_cl)) do (i, lam)
        rp = round(real(lam), digits=4)
        ip = round(imag(lam), digits=4)
        val = abs(ip) < 1e-9 ? "$rp" : "$rp + $(ip)i"
        "| lambda_$i | $val |"
    end
    all_stable = all(real.(lqr.eigenvalues_cl) .< 0)
    Markdown.parse("""
    ### Ganancia K  (ley de control u = -K x)
    | Componente | Valor |
    |---|---|
    $(join(k_rows, "\n"))

    ### Polos en lazo cerrado (A - B K)
    | Eigenvalor | Valor |
    |---|---|
    $(join(cl_rows, "\n"))

    **Lazo cerrado:** $(all_stable ? "todos los polos estables" : "aun inestable")
    """)
end

# ╔═╡ 00000040-0000-4000-8000-000000000000
md"""
## 6. Simulacion no lineal

Se integra el modelo **no lineal** con la ley ``u = -K\mathbf{x}``, partiendo de
una pequena inclinacion en ambos eslabones, y se compara con la respuesta sin
control (caida libre).
"""

# ╔═╡ 00000041-0000-4000-8000-000000000000
md"**Angulo inicial del eslabon 1, theta1_0 [rad]:**"

# ╔═╡ 00000042-0000-4000-8000-000000000000
@bind th1_0 PlutoUI.Slider(0.01:0.01:0.4, default=0.1, show_value=true)

# ╔═╡ 00000043-0000-4000-8000-000000000000
md"**Angulo inicial del eslabon 2, theta2_0 [rad]:**"

# ╔═╡ 00000044-0000-4000-8000-000000000000
@bind th2_0 PlutoUI.Slider(0.01:0.01:0.4, default=0.1, show_value=true)

# ╔═╡ 00000045-0000-4000-8000-000000000000
begin
    x0 = [0.0, 0.0, th1_0, 0.0, th2_0, 0.0]

    prob_free = ODEProblem(nonlinear_eom_double!, x0, (0.0, 3.0),
        (params=params, F=0.0))
    sol_free = solve(prob_free, Tsit5(), saveat=0.02)

    prob_lqr = ODEProblem(closed_loop_eom_double!, x0, (0.0, 10.0),
        (params=params, K=lqr.K, saturate=100.0))
    sol_lqr = solve(prob_lqr, Tsit5(), saveat=0.02)

    md"""
    Simulacion lista: libre (3 s) y con LQR (10 s).
    - theta1 final con LQR = $(round(sol_lqr.u[end][3], digits=5)) rad
    - theta2 final con LQR = $(round(sol_lqr.u[end][5], digits=5)) rad
    """
end

# ╔═╡ 00000046-0000-4000-8000-000000000000
md"""
## 7. Graficas
"""

# ╔═╡ 00000047-0000-4000-8000-000000000000
begin
    COL_LQR = RGBf(0.17, 0.45, 0.71)
    COL_FREE = RGBf(0.84, 0.15, 0.16)
    COL_FORCE = RGBf(0.12, 0.47, 0.71)
    md"Colores configurados"
end

# ╔═╡ 00000048-0000-4000-8000-000000000000
begin
    t_free = sol_free.t
    t_lqr = sol_lqr.t
    th1_free = [u[3] for u in sol_free.u]
    th2_free = [u[5] for u in sol_free.u]
    x_free = [u[1] for u in sol_free.u]
    th1_lqr = [u[3] for u in sol_lqr.u]
    th2_lqr = [u[5] for u in sol_lqr.u]
    x_lqr = [u[1] for u in sol_lqr.u]
    F_lqr = [clamp(-dot(lqr.K[1, :], u), -100.0, 100.0) for u in sol_lqr.u]
    md"Trayectorias extraidas ($(length(t_lqr)) puntos con LQR)"
end

# ╔═╡ 00000049-0000-4000-8000-000000000000
let
    fig = Figure(size=(880, 640))

    ax1 = Axis(fig[1, 1], title="Angulo eslabon 1 (theta1)",
        xlabel="t [s]", ylabel="theta1 [rad]")
    lines!(ax1, t_lqr, th1_lqr, color=COL_LQR, linewidth=2, label="Con LQR")
    lines!(ax1, t_free, th1_free, color=COL_FREE, linewidth=2, linestyle=:dash, label="Sin control")
    hlines!(ax1, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="Angulo eslabon 2 (theta2)",
        xlabel="t [s]", ylabel="theta2 [rad]")
    lines!(ax2, t_lqr, th2_lqr, color=COL_LQR, linewidth=2)
    lines!(ax2, t_free, th2_free, color=COL_FREE, linewidth=2, linestyle=:dash)
    hlines!(ax2, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)

    ax3 = Axis(fig[2, 1], title="Posicion del carro",
        xlabel="t [s]", ylabel="x [m]")
    lines!(ax3, t_lqr, x_lqr, color=COL_LQR, linewidth=2)
    lines!(ax3, t_free, x_free, color=COL_FREE, linewidth=2, linestyle=:dash)
    hlines!(ax3, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)

    ax4 = Axis(fig[2, 2], title="Fuerza de control u = -K x",
        xlabel="t [s]", ylabel="F [N]")
    lines!(ax4, t_lqr, F_lqr, color=COL_FORCE, linewidth=2)
    hlines!(ax4, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)

    Label(fig[0, :], "Respuesta temporal: sin control vs LQR",
        fontsize=15, font=:bold)
    fig
end

# ╔═╡ 00000050-0000-4000-8000-000000000000
let
    fig = Figure(size=(640, 430))
    ax = Axis(fig[1, 1],
        title="Eigenvalores: lazo abierto vs cerrado",
        xlabel="Re(lambda)", ylabel="Im(lambda)")
    vlines!(ax, [0.0], color=:gray50, linewidth=1.2, linestyle=:dash)
    hlines!(ax, [0.0], color=:gray80, linewidth=0.5)
    scatter!(ax, real.(ss.eigenvalues), imag.(ss.eigenvalues),
        color=COL_FREE, markersize=14, label="Lazo abierto (A)")
    scatter!(ax, real.(lqr.eigenvalues_cl), imag.(lqr.eigenvalues_cl),
        color=COL_LQR, markersize=14, marker=:utriangle, label="Lazo cerrado (A - BK)")
    axislegend(ax, position=:lt)
    fig
end

# ╔═╡ 00000051-0000-4000-8000-000000000000
md"""
## 8. Animacion del pendulo doble

Visualizacion mecanica del carro y los dos eslabones. Usa el slider de tiempo
para recorrer la trayectoria controlada cuadro a cuadro.
"""

# ╔═╡ 00000052-0000-4000-8000-000000000000
begin
    function draw_double_frame!(ax, state, p; time_label=nothing)
        x_cart = state[1]
        t1 = state[3]
        t2 = state[5]
        jx = x_cart + p.L1 * sin(t1)
        jy = p.L1 * cos(t1)
        ex = jx + p.L2 * sin(t2)
        ey = jy + p.L2 * cos(t2)

        empty!(ax)
        lines!(ax, [-3.0, 3.0], [0.0, 0.0], color=:gray70, linewidth=1.5, linestyle=:dash)
        poly!(ax, Rect(x_cart - 0.2, -0.1, 0.4, 0.2),
            color=(COL_LQR, 0.85), strokecolor=:gray20, strokewidth=1.5)
        for dx in (-0.12, 0.12)
            scatter!(ax, [x_cart + dx], [-0.12], markersize=10, color=:gray45)
        end
        lines!(ax, [x_cart, jx], [0.0, jy], color=:gray25, linewidth=3.5)
        lines!(ax, [jx, ex], [jy, ey], color=:gray45, linewidth=3.5)
        scatter!(ax, [jx], [jy], markersize=15, color=:orange,
            strokecolor=:black, strokewidth=1)
        scatter!(ax, [ex], [ey], markersize=19, color=:red,
            strokecolor=:black, strokewidth=1)
        scatter!(ax, [x_cart], [0.0], markersize=6, color=:gray30)
        if !isnothing(time_label)
            text!(ax, 2.9, p.L1 + p.L2 + 0.3, text=time_label,
                fontsize=12, color=:gray40, align=(:right, :top))
        end
    end
    md"Funcion de dibujo draw_double_frame! lista"
end

# ╔═╡ 00000053-0000-4000-8000-000000000000
md"**Tiempo t [s] -- trayectoria controlada (LQR):**"

# ╔═╡ 00000054-0000-4000-8000-000000000000
@bind t_frame PlutoUI.Slider(range(0.0, 10.0, length=501), default=0.0, show_value=true)

# ╔═╡ 00000055-0000-4000-8000-000000000000
let
    reach = params.L1 + params.L2
    fig = Figure(size=(720, 540))

    ax = Axis(fig[1, 1], title="Pendulo doble -- trayectoria LQR",
        xlabel="x [m]", ylabel="y [m]", aspect=DataAspect(),
        limits=(-2.5, 2.5, -0.4, reach + 0.4))
    state_now = sol_lqr(t_frame)
    draw_double_frame!(ax, state_now, params,
        time_label=@sprintf("t = %.2f s", t_frame))

    ax2 = Axis(fig[2, 1], xlabel="t [s]", ylabel="angulo [rad]", height=120,
        title="theta1 (azul) y theta2 (naranja)")
    lines!(ax2, t_lqr, th1_lqr, color=COL_LQR, linewidth=1.5)
    lines!(ax2, t_lqr, th2_lqr, color=:darkorange, linewidth=1.5)
    hlines!(ax2, [0.0], color=:gray60, linewidth=0.5, linestyle=:dot)
    vlines!(ax2, [t_frame], color=:red, linewidth=1.8)

    rowgap!(fig.layout, 12)
    fig
end

# ╔═╡ 00000056-0000-4000-8000-000000000000
md"""
## 9. Exportar animacion (GIF)

Marca la casilla para generar `figures/08_doble_exploracion.gif` (los primeros
6 s de la trayectoria controlada). Desmarcala para no regenerar el archivo en
cada cambio de parametro.
"""

# ╔═╡ 00000057-0000-4000-8000-000000000000
@bind do_export PlutoUI.CheckBox(default=false)

# ╔═╡ 00000058-0000-4000-8000-000000000000
begin
    if do_export
        fig_dir = joinpath(@__DIR__, "..", "figures")
        mkpath(fig_dir)
        gif_path = joinpath(fig_dir, "08_doble_exploracion.gif")

        reach_e = params.L1 + params.L2
        fig_e = Figure(size=(640, 520))
        ax_e = Axis(fig_e[1, 1], title="Pendulo doble -- LQR",
            xlabel="x [m]", ylabel="y [m]", aspect=DataAspect(),
            limits=(-2.5, 2.5, -0.4, reach_e + 0.4))

        ts = range(0.0, 6.0, step=1/30)
        record(fig_e, gif_path, eachindex(ts); framerate=30) do i
            draw_double_frame!(ax_e, sol_lqr(ts[i]), params,
                time_label=@sprintf("t = %.2f s", ts[i]))
        end
        md"GIF exportado en: `$gif_path`"
    else
        md"Marca la casilla de arriba para exportar el GIF."
    end
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000001
md"""
---

# Parte II -- Lo que agrega la segunda entrega

Los modulos `Metrics`, `Observer` y `Sweep` son **genericos en la dimension del
estado**: el mismo codigo corre sobre el simple de ``\mathbb{R}^4``, este de
``\mathbb{R}^6`` y el triple de ``\mathbb{R}^8``. Todo lo que sigue reacciona a
los sliders de arriba.

## 10. Metricas de viabilidad -- mas alla del criterio de rango

``\text{rank}\,\mathcal{C} = n`` es **binario** y se evalua con una tolerancia
implicita. La version cuantitativa del test de Popov-Belevitch-Hautus sustituye
el rango por el menor valor singular:

``\mu_{PBH} = \min_{\lambda\in\sigma(A)}\ \sigma_{\min}[\,A-\lambda I \mid B\,]``

que por Eckart-Young es la distancia en norma 2 a perder la controlabilidad:
mide **cuanto habria que perturbar el sistema** para volverlo incontrolable.
"""

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000002
begin
    met_dm = dominant_mode(ss.A)
    met_pbh = pbh_controllability_margin_normalized(ss.A, ss.B)
    met_pbh_o = pbh_observability_margin_normalized(ss.A, ss.C)
    met_cond = controllability_condition(ss.A, ss.B)
    met_worst = pbh_controllability_profile(ss.A, ss.B)[1]

    Markdown.parse("""
    | Metrica | Valor |
    |---|---|
    | ``\\lambda_{\\max}`` | $(round(met_dm.lambda_max, digits=4)) s⁻¹ |
    | ``\\tau_{caida} = 1/\\lambda_{\\max}`` | $(met_dm.n_unstable == 0 ? "—" : "$(round(1000*met_dm.tau_fall, digits=1)) ms") |
    | Modos inestables | $(met_dm.n_unstable) |
    | Margen PBH normalizado (control) | **$(@sprintf("%.4e", met_pbh))** |
    | Margen PBH normalizado (observacion) | $(@sprintf("%.4e", met_pbh_o)) |
    | ``\\text{cond}\\,\\mathcal{C}`` | $(@sprintf("%.3e", met_cond)) |
    | Modo mas debil | ``\\lambda`` = $(round(real(met_worst.lambda), digits=3)) |

    El margen del doble esta un orden por debajo del simple
    (``1.7\\times10^{-2}``) y uno por encima del triple
    (``4.9\\times10^{-4}``): cae **suavemente**, mientras que
    ``\\text{cond}\\,\\mathcal{C}`` crece cuatro ordenes por eslabon.
    """)
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000003
md"""
## 11. El elipsoide de no saturacion

Con ``V(\mathbf{x}) = \mathbf{x}^T P\mathbf{x}`` la funcion de Lyapunov del LQR,
el mayor conjunto de nivel contenido en ``\{|K\mathbf{x}|\le u_{\max}\}`` sale de
Cauchy-Schwarz: ``c^\star = u_{\max}^2/(KP^{-1}K^T)``. Dentro de el, el control
**nunca satura** y el sistema lineal converge. Sobre el eje de ``\theta_1`` la
garantia es ``|\theta_0| \le \sqrt{c^\star/P_{33}}``.
"""

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000004
md"**Limite del actuador u_max [N]:**"

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000005
@bind met_umax PlutoUI.Slider(5.0:5.0:250.0, default=100.0, show_value=true)

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000006
begin
    met_cstar = nonsaturation_level(lqr.K, lqr.P, R_mat, ss.B, met_umax)
    met_theta = max_axis_angle(lqr.P, met_cstar, 3)
    met_x0 = zeros(6); met_x0[3] = 0.1
    met_J = optimal_cost(lqr.P, met_x0)

    Markdown.parse("""
    | Cantidad | Valor |
    |---|---|
    | ``P_{33}`` | $(round(lqr.P[3,3], digits=4)) |
    | ``c^\\star`` | $(round(met_cstar, digits=4)) |
    | ``\\theta_0`` garantizado sin saturar | **$(round(met_theta, digits=4)) rad = $(round(rad2deg(met_theta), digits=2))°** |
    | ``J^\\star`` desde ``\\theta_1 = 0.1`` rad | $(round(met_J, digits=4)) |

    > **Honestidad sobre el alcance.** ``c^\\star`` certifica no saturacion mas
    > convergencia del sistema **lineal**. No es la region de atraccion no
    > lineal: para eso haria falta acotar el residuo de la linealizacion.
    """)
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000007
md"""
## 12. Control digital -- periodo de muestreo maximo

Con retenedor de orden cero, ``A_d = e^{Ah}`` y
``B_d = \left(\int_0^h e^{A\tau}d\tau\right)B`` salen de **una sola**
exponencial matricial por el truco de Van Loan: se exponencia
``\begin{pmatrix}A & B\\ 0 & 0\end{pmatrix}`` y los bloques superiores son
``A_d`` y ``B_d``. La estabilidad discreta pide que ``A_d - B_dK`` sea de
**Schur**.
"""

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000008
begin
    met_hmax = max_sampling_period(ss.A, ss.B, lqr.K)
    Markdown.parse("""
    | Cantidad | Valor |
    |---|---|
    | ``h_{\\max}`` | **$(round(1000*met_hmax, digits=1)) ms** (``\\ge`` $(round(1/met_hmax, digits=1)) Hz) |
    | ``1/\\lambda_{\\max}`` | $(met_dm.n_unstable == 0 ? "—" : "$(round(1000*met_dm.tau_fall, digits=1)) ms") |
    | ``h_{\\max}\\,\\lambda_{\\max}`` | $(round(met_hmax*met_dm.lambda_max, digits=3)) |

    El margen relativo se estrecha al crecer la complejidad: 0.83 en el simple,
    0.66 aqui, 0.58 en el triple.
    """)
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000009
md"**Periodo de muestreo h a probar [ms]** -- cruza ``h_{\max}`` y mira los polos salir del circulo:"

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000010
@bind met_h_ms PlutoUI.Slider(1.0:1.0:200.0, default=40.0, show_value=true)

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000011
let
    h = met_h_ms / 1000
    d = discretize_zoh(ss.A, ss.B, h)
    lam_d = eigvals(d.Ad - d.Bd * lqr.K)
    radio = maximum(abs.(lam_d))
    ok = is_schur(d.Ad - d.Bd * lqr.K)

    fig = Figure(size=(620, 470))
    ax = Axis(fig[1, 1],
        title=@sprintf("Polos discretos con h = %.0f ms  (radio %.4f)", met_h_ms, radio),
        xlabel="Re(z)", ylabel="Im(z)", aspect=DataAspect())
    tt = range(0, 2pi, length=200)
    lines!(ax, cos.(tt), sin.(tt), color=:gray50, linewidth=1.5, linestyle=:dash)
    hlines!(ax, [0.0], color=:gray85, linewidth=0.5)
    vlines!(ax, [0.0], color=:gray85, linewidth=0.5)
    scatter!(ax, real.(lam_d), imag.(lam_d), markersize=14,
        color=ok ? COL_LQR : COL_FREE)
    text!(ax, 0.0, -1.35,
        text=ok ? "DENTRO del circulo: el lazo digital es estable" : "FUERA: el lazo digital diverge",
        align=(:center, :center), fontsize=13,
        color=ok ? :gray25 : RGBf(0.7, 0.1, 0.1))
    limits!(ax, -1.5, 1.5, -1.5, 1.5)
    fig
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000012
md"""
## 13. El observador de Luenberger

Hasta aqui el control ha sido ``u = -K\mathbf{x}``, con el estado **completo**:
nadie mide la velocidad angular de cada eslabon. Como
``\text{rank}\,\mathcal{O} = 6``, el estado se puede reconstruir:

``\dot{\hat{\mathbf{x}}} = A\hat{\mathbf{x}} + Bu + L(\mathbf{y} - C\hat{\mathbf{x}})``,
con ``\dot{\mathbf{e}} = (A - LC)\mathbf{e}``

Por la **dualidad de Kalman**, ``L = \text{lqr}(A^T, C^T, Q_o, R_o)^T``: la misma
ecuacion de Riccati que dio ``K`` da, transpuesta, la ganancia del observador.
"""

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000013
md"**Confianza en los sensores** -- ``Q_o = 10^{\alpha}I``; mas alto = observador mas rapido:"

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000014
@bind obs_log_qo PlutoUI.Slider(0.0:0.5:5.0, default=2.0, show_value=true)

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000015
begin
    obs_Qo = (10.0^obs_log_qo) * Matrix(I, 6, 6)
    obs_Ro = Matrix(I, size(ss.C, 1), size(ss.C, 1))
    obs_res = design_observer(ss.A, ss.C, obs_Qo, obs_Ro)
    obs_sep = check_separation(ss.A, ss.B, ss.C, lqr.K, obs_res.L)
    obs_M = augmented_closed_loop(ss.A, ss.B, ss.C, lqr.K, obs_res.L)

    obs_rows = map(enumerate(sort(obs_res.eigenvalues_obs, by=real))) do (i, lam)
        rp = round(real(lam), digits=3); ip = round(imag(lam), digits=3)
        "| ``\\nu_{$i}`` | $(abs(ip) < 1e-9 ? "$rp" : "$rp + $(ip)i") |"
    end

    Markdown.parse("""
    ### Polos del observador ``\\sigma(A - LC)``

    | Polo | Valor |
    |---|---|
    $(join(obs_rows, "\n"))

    ``\\lVert L\\rVert`` = $(round(norm(obs_res.L), digits=3)) ·
    Observador $(all(real.(obs_res.eigenvalues_obs) .< 0) ? "**estable**" : "INESTABLE")

    ### Principio de separacion

    En coordenadas ``(\\mathbf{x}, \\mathbf{e})`` la matriz aumentada es
    **triangular por bloques**, y el espectro de una triangular por bloques es la
    union de los espectros diagonales:
    ``\\sigma = \\sigma(A-BK)\\cup\\sigma(A-LC)``. Controlador y observador se
    disenan por separado.

    | Verificacion | Resultado |
    |---|---|
    | Bloque inferior izquierdo nulo | $(all(obs_M[7:12, 1:6] .== 0) ? "SI" : "no") |
    | Error maximo al emparejar espectros | $(@sprintf("%.3e", obs_sep.max_error)) |
    | Principio de separacion | $(obs_sep.holds ? "**SE CUMPLE**" : "NO se cumple") |
    """)
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000016
md"""
### Que sensores hacen falta

De los ``2^3 - 1 = 7`` subconjuntos no vacios de
``\{x, \theta_1, \theta_2\}``, ¿cuales conservan la observabilidad? El margen
PBH de observabilidad distingue lo que el rango no puede: *observable* no es lo
mismo que *practicable*.
"""

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000017
begin
    obs_tabla = sensor_subset_analysis(ss.A, ss.C, ["x", "th1", "th2"])
    obs_filas = map(obs_tabla) do f
        "| $(join(f.labels, " + ")) | $(f.rank)/6 | $(f.observable ? "SI" : "**NO**") | $(@sprintf("%.3e", f.margin_normalized)) |"
    end
    Markdown.parse("""
    | Sensores | rank ``\\mathcal{O}`` | Observable | Margen PBH normalizado |
    |---|---|---|---|
    $(join(obs_filas, "\n"))

    $(count(f -> f.observable, obs_tabla)) de $(length(obs_tabla)) subconjuntos
    son observables. Midiendo solo angulos, la posicion absoluta del carro es
    inobservable: la dinamica angular es invariante a trasladarlo.
    """)
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000018
md"""
### Simulacion con observador

La planta arranca inclinada y el estimador arranca en **cero**: el maximo
desconocimiento posible sobre el estado inicial. El control es
``u = -K\hat{\mathbf{x}}``, no ``-K\mathbf{x}``.
"""

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000019
md"**Angulo inicial de theta1 para la prueba del observador [rad]:**"

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000020
@bind obs_th0 PlutoUI.Slider(0.01:0.01:0.35, default=0.1, show_value=true)

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000021
begin
    obs_z0 = zeros(12)
    obs_z0[3] = obs_th0
    obs_p = (params=params, eom! = nonlinear_eom_double!, A=ss.A, B=ss.B, C=ss.C,
        K=lqr.K, L=obs_res.L, saturate=met_umax)
    obs_sol = solve(ODEProblem(observer_eom!, obs_z0, (0.0, 10.0), obs_p),
        Tsit5(), saveat=0.01)

    obs_err = [norm(u[1:6] .- u[7:12]) for u in obs_sol.u]
    obs_u = [clamp(-dot(lqr.K[1, :], u[7:12]), -met_umax, met_umax) for u in obs_sol.u]
    obs_i1 = findfirst(e -> e < 0.01 * max(obs_err[1], eps()), obs_err)
    obs_conv = all(abs(obs_sol.u[end][j]) < 0.02 for j in (3, 5))

    Markdown.parse("""
    | Cantidad | Valor |
    |---|---|
    | ``\\lVert e(0)\\rVert`` | $(round(obs_err[1], digits=4)) |
    | ``\\lVert e(T)\\rVert`` | $(@sprintf("%.3e", obs_err[end])) |
    | Cae al 1 % del inicial | $(obs_i1 === nothing ? "no alcanza" : "t = $(round(obs_sol.t[obs_i1], digits=2)) s") |
    | Fuerza pico con observador | $(round(maximum(abs.(obs_u)), digits=2)) N |
    | Converge la planta | $(obs_conv ? "**SI**" : "**NO**") |

    > **El observador no sale gratis.** El principio de separacion garantiza
    > estabilidad del sistema **lineal**. Sobre el modelo no lineal la region de
    > atraccion se encoge: durante el transitorio el control actua sobre una
    > estimacion equivocada. Sube el slider hasta que deje de converger.
    """)
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000022
let
    fig = Figure(size=(900, 620))

    ax1 = Axis(fig[1, 1], title="theta1 -- variable MEDIDA",
        xlabel="t [s]", ylabel="theta1 [rad]")
    lines!(ax1, obs_sol.t, [u[3] for u in obs_sol.u], color=:black,
        linewidth=2.5, label="real")
    lines!(ax1, obs_sol.t, [u[9] for u in obs_sol.u], color=:darkorange,
        linewidth=2, linestyle=:dash, label="estimado")
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="omega1 -- variable NO medida, reconstruida",
        xlabel="t [s]", ylabel="omega1 [rad/s]")
    lines!(ax2, obs_sol.t, [u[4] for u in obs_sol.u], color=:black,
        linewidth=2.5, label="real")
    lines!(ax2, obs_sol.t, [u[10] for u in obs_sol.u], color=:darkorange,
        linewidth=2, linestyle=:dash, label="estimado")
    axislegend(ax2, position=:rb)

    ax3 = Axis(fig[2, 1], title="Error de estimacion ||x - xhat||",
        xlabel="t [s]", ylabel="||e||", yscale=log10)
    lines!(ax3, obs_sol.t, max.(obs_err, 1e-14), color=RGBf(0.55, 0.24, 0.62),
        linewidth=2)

    ax4 = Axis(fig[2, 2], title="Fuerza: u = -K xhat vs u = -K x",
        xlabel="t [s]", ylabel="F [N]")
    lines!(ax4, obs_sol.t, obs_u, color=RGBf(0.55, 0.24, 0.62), linewidth=2,
        label="con observador")
    lines!(ax4, t_lqr, F_lqr, color=COL_LQR, linewidth=2, linestyle=:dash,
        label="estado completo")
    axislegend(ax4, position=:rt)

    Label(fig[0, :], "Observador de Luenberger: el estimador arranca en cero",
        fontsize=15, font=:bold)
    fig
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000023
md"""
## 14. Barrido de parametros

El mismo sistema visto como instancia del modelo generico de ``N`` eslabones
(``N = 2`` con masas puntuales). Elige que parametro barrer.
"""

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000024
@bind sw_sel PlutoUI.Select([
    "M" => "Masa del carro M",
    "m2" => "Masa del extremo m2",
    "L2" => "Longitud del eslabon superior L2",
    "L1" => "Longitud del eslabon inferior L1",
], default="M")

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000025
begin
    sw_base = point_masses(M_val, [m1_val, m2_val], [L1_val, L2_val]; g=g_val, b=0.0)
    sw_map = Dict(
        "M" => (:M, collect(10 .^ range(log10(0.05), log10(20.0), length=40))),
        "m2" => ((:m, 2), collect(10 .^ range(log10(0.02), log10(3.0), length=40))),
        "L2" => ((:l, 2), collect(range(0.05, 2.0, length=40))),
        "L1" => ((:l, 1), collect(range(0.05, 2.0, length=40))),
    )
    sw_spec, sw_vals = sw_map[sw_sel]
    sw_res = sweep_1d(linearize_system_nlink, sw_base, sw_spec, sw_vals)
    sw_pbh = [r.pbh_normalized for r in sw_res]
    sw_normk = [r.norm_K for r in sw_res]
    sw_lam = [r.lambda_max for r in sw_res]
    sw_best = argmax(sw_pbh)

    Markdown.parse("""
    Barrido de **$(sw_sel)** en $(length(sw_vals)) puntos.
    Mejor margen PBH en $(sw_sel) = **$(round(sw_vals[sw_best], digits=4))**
    ($(@sprintf("%.4e", sw_pbh[sw_best]))).
    $(1 < sw_best < length(sw_vals) ? "El optimo es **interior**: no esta en ningun extremo del dominio." : "El optimo cae en un extremo del dominio.")
    """)
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000026
let
    logx = sw_sel in ("M", "m2")
    fig = Figure(size=(900, 320))

    ax1 = Axis(fig[1, 1], title="Margen PBH normalizado", xlabel=sw_sel,
        ylabel="PBH norm", yscale=log10, xscale=logx ? log10 : identity)
    lines!(ax1, sw_vals, sw_pbh, color=COL_LQR, linewidth=2.5)
    scatter!(ax1, [sw_vals[sw_best]], [sw_pbh[sw_best]], color=COL_FREE,
        markersize=14, marker=:star5)

    ax2 = Axis(fig[1, 2], title="Esfuerzo de control ||K||", xlabel=sw_sel,
        ylabel="norm(K)", yscale=log10, xscale=logx ? log10 : identity)
    lines!(ax2, sw_vals, sw_normk, color=:darkorange, linewidth=2.5)

    ax3 = Axis(fig[1, 3], title="Modo dominante", xlabel=sw_sel,
        ylabel="lambda_max [1/s]", xscale=logx ? log10 : identity)
    lines!(ax3, sw_vals, sw_lam, color=COL_FREE, linewidth=2.5)

    fig
end

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000027
md"""
### La frontera practica

Biseccion sobre el modelo **no lineal** con saturacion. Criterio de exito:
``\max_j|\theta_j(T)| < 0.02`` rad, el carro dentro del riel y sin que ningun
angulo pase de ``\pi/2``. La cota elipsoidal debe compararse contra la columna
**sin riel**: solo certifica saturacion mas convergencia lineal.
"""

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000028
md"**Longitud del riel [m]:**"

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000029
@bind sw_rail PlutoUI.Slider(0.5:0.1:5.0, default=1.5, show_value=true)

# ╔═╡ f6a7b8c9-5555-4000-d000-000000000030
begin
    sw_th_riel = max_recoverable_angle(sw_base, lqr.K; umax=met_umax, x_rail=sw_rail)
    sw_th_libre = max_recoverable_angle(sw_base, lqr.K; umax=met_umax, x_rail=1e6)
    sw_manda = abs(sw_th_libre - sw_th_riel) > 0.01

    Markdown.parse("""
    | Frontera | Valor |
    |---|---|
    | Cota elipsoidal (garantia lineal) | $(round(met_theta, digits=4)) rad |
    | Frontera real con riel de $(sw_rail) m | **$(round(sw_th_riel, digits=4)) rad** |
    | Frontera real sin restriccion de riel | $(round(sw_th_libre, digits=4)) rad |
    | Restriccion activa | **$(sw_manda ? "EL RIEL" : "LA SATURACION")** |

    En el doble la restriccion activa suele ser la **saturacion**: el riel nunca
    llega a morder, al reves que en el simple. Cual de los limites manda cambia
    con la configuracion, y eso obliga a matizar la idea de "imposibilidad
    practica" como un bloque homogeneo.
    """)
end

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
# ╠═00000017-0000-4000-8000-000000000000
# ╟─00000018-0000-4000-8000-000000000000
# ╠═00000019-0000-4000-8000-000000000000
# ╟─00000020-0000-4000-8000-000000000000
# ╠═00000021-0000-4000-8000-000000000000
# ╟─00000022-0000-4000-8000-000000000000
# ╠═00000023-0000-4000-8000-000000000000
# ╟─00000024-0000-4000-8000-000000000000
# ╠═00000025-0000-4000-8000-000000000000
# ╟─00000026-0000-4000-8000-000000000000
# ╠═00000027-0000-4000-8000-000000000000
# ╟─00000028-0000-4000-8000-000000000000
# ╠═00000029-0000-4000-8000-000000000000
# ╟─00000030-0000-4000-8000-000000000000
# ╟─00000031-0000-4000-8000-000000000000
# ╠═00000032-0000-4000-8000-000000000000
# ╟─00000033-0000-4000-8000-000000000000
# ╠═00000034-0000-4000-8000-000000000000
# ╟─00000035-0000-4000-8000-000000000000
# ╠═00000036-0000-4000-8000-000000000000
# ╟─00000037-0000-4000-8000-000000000000
# ╠═00000038-0000-4000-8000-000000000000
# ╠═00000039-0000-4000-8000-000000000000
# ╟─00000040-0000-4000-8000-000000000000
# ╟─00000041-0000-4000-8000-000000000000
# ╠═00000042-0000-4000-8000-000000000000
# ╟─00000043-0000-4000-8000-000000000000
# ╠═00000044-0000-4000-8000-000000000000
# ╠═00000045-0000-4000-8000-000000000000
# ╟─00000046-0000-4000-8000-000000000000
# ╠═00000047-0000-4000-8000-000000000000
# ╠═00000048-0000-4000-8000-000000000000
# ╠═00000049-0000-4000-8000-000000000000
# ╠═00000050-0000-4000-8000-000000000000
# ╟─00000051-0000-4000-8000-000000000000
# ╠═00000052-0000-4000-8000-000000000000
# ╟─00000053-0000-4000-8000-000000000000
# ╠═00000054-0000-4000-8000-000000000000
# ╠═00000055-0000-4000-8000-000000000000
# ╟─00000056-0000-4000-8000-000000000000
# ╠═00000057-0000-4000-8000-000000000000
# ╠═00000058-0000-4000-8000-000000000000
# ╟─f6a7b8c9-5555-4000-d000-000000000001
# ╠═f6a7b8c9-5555-4000-d000-000000000002
# ╟─f6a7b8c9-5555-4000-d000-000000000003
# ╟─f6a7b8c9-5555-4000-d000-000000000004
# ╠═f6a7b8c9-5555-4000-d000-000000000005
# ╠═f6a7b8c9-5555-4000-d000-000000000006
# ╟─f6a7b8c9-5555-4000-d000-000000000007
# ╠═f6a7b8c9-5555-4000-d000-000000000008
# ╟─f6a7b8c9-5555-4000-d000-000000000009
# ╠═f6a7b8c9-5555-4000-d000-000000000010
# ╠═f6a7b8c9-5555-4000-d000-000000000011
# ╟─f6a7b8c9-5555-4000-d000-000000000012
# ╟─f6a7b8c9-5555-4000-d000-000000000013
# ╠═f6a7b8c9-5555-4000-d000-000000000014
# ╠═f6a7b8c9-5555-4000-d000-000000000015
# ╟─f6a7b8c9-5555-4000-d000-000000000016
# ╠═f6a7b8c9-5555-4000-d000-000000000017
# ╟─f6a7b8c9-5555-4000-d000-000000000018
# ╟─f6a7b8c9-5555-4000-d000-000000000019
# ╠═f6a7b8c9-5555-4000-d000-000000000020
# ╠═f6a7b8c9-5555-4000-d000-000000000021
# ╠═f6a7b8c9-5555-4000-d000-000000000022
# ╟─f6a7b8c9-5555-4000-d000-000000000023
# ╠═f6a7b8c9-5555-4000-d000-000000000024
# ╠═f6a7b8c9-5555-4000-d000-000000000025
# ╠═f6a7b8c9-5555-4000-d000-000000000026
# ╟─f6a7b8c9-5555-4000-d000-000000000027
# ╟─f6a7b8c9-5555-4000-d000-000000000028
# ╠═f6a7b8c9-5555-4000-d000-000000000029
# ╠═f6a7b8c9-5555-4000-d000-000000000030
