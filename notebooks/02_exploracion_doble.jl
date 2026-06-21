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
    include(joinpath(@__DIR__, "..", "src", "linearization.jl"))
    include(joinpath(@__DIR__, "..", "src", "controller.jl"))

    using .ModelDouble
    using .Linearization
    using .Controller

    md"Modulos cargados: ModelDouble, Linearization, Controller"
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
