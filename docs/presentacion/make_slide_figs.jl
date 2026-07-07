# =============================================================================
# make_slide_figs.jl -- Genera las figuras de la presentacion (Beamer)
# =============================================================================
# Produce, en docs/presentacion/figs/, graficas con la paleta de la
# presentacion (azul marino, azul petroleo, verde gris, salvia, limon).
# Reutiliza los modulos del proyecto; no duplica fisica ni algoritmos.
#
# Ejecutar desde la raiz del proyecto:
#   julia --project=. docs/presentacion/make_slide_figs.jl
# =============================================================================

using Pkg
const PROJ_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(PROJ_ROOT)

include(joinpath(PROJ_ROOT, "src", "model_simple.jl"))
include(joinpath(PROJ_ROOT, "src", "model_double.jl"))
include(joinpath(PROJ_ROOT, "src", "linearization.jl"))
include(joinpath(PROJ_ROOT, "src", "controller.jl"))

using .Model
using .ModelDouble
using .Linearization
using .Controller
using LinearAlgebra
using Printf
using DifferentialEquations
using CairoMakie

const FIGS = joinpath(@__DIR__, "figs")
isdir(FIGS) || mkdir(FIGS)

# Paleta de la presentacion
const MARINO    = RGBf(27/255, 64/255, 121/255)    # 1B4079
const PETROLEO  = RGBf(77/255, 124/255, 138/255)   # 4D7C8A
const GRISVERDE = RGBf(127/255, 156/255, 150/255)  # 7F9C96
const SALVIA    = RGBf(143/255, 173/255, 136/255)  # 8FAD88
const LIMON     = RGBf(199/255, 219/255, 148/255)  # C7DB94

set_theme!(fontsize=15)

"Primer instante tras el cual |y| permanece por debajo de tol (banda absoluta)."
function settling_time(t, y; tol=0.02)
    idx = findlast(>(tol), abs.(y))
    return idx === nothing ? t[1] : t[min(idx + 1, length(t))]
end

# =============================================================================
# Disenos (mismos pesos y polos del informe)
# =============================================================================

ps = default_params()
ss = linearize_system(ps)
Q = diagm([1.0, 0.0, 10.0, 0.0]); R = reshape([0.1], 1, 1)
lqr = design_lqr(ss.A, ss.B, Q, R)
ack = design_pole_placement(ss.A, ss.B, [-1.0, -2.0, -3.0, -4.0])

pd = default_params_double()
ssd = linearize_system_double(pd)
Qd = diagm([1.0, 0.0, 10.0, 0.0, 10.0, 0.0]); Rd = reshape([0.1], 1, 1)
lqrd = design_lqr(ssd.A, ssd.B, Qd, Rd)

# =============================================================================
# Figura 1: mapa de polos, lazo abierto frente a lazo cerrado
# =============================================================================

let
    # Esta figura se muestra a media columna en la diapositiva: fuentes grandes
    fig = Figure(size=(1000, 420), fontsize=21)

    ax1 = Axis(fig[1, 1], title="Péndulo simple", xlabel="Re(λ)", ylabel="Im(λ)")
    vspan!(ax1, 0.0, 6.5, color=(LIMON, 0.35))
    vlines!(ax1, [0.0], color=:gray50, linewidth=1)
    hlines!(ax1, [0.0], color=:gray80, linewidth=0.6)
    scatter!(ax1, real.(ss.eigenvalues), imag.(ss.eigenvalues),
             marker=:xcross, markersize=22, color=GRISVERDE, label="lazo abierto")
    scatter!(ax1, real.(lqr.eigenvalues_cl), imag.(lqr.eigenvalues_cl),
             marker=:circle, markersize=18, color=MARINO, label="LQR")
    scatter!(ax1, real.(ack.eigenvalues_cl), imag.(ack.eigenvalues_cl),
             marker=:diamond, markersize=18, color=PETROLEO, label="Ackermann")
    xlims!(ax1, -6.5, 6.5); ylims!(ax1, -2.9, 2.9)
    axislegend(ax1, position=:lt, framevisible=true,
               backgroundcolor=(:white, 0.85), labelsize=17)
    text!(ax1, 3.3, -2.45, text="Re(λ) > 0", align=(:center, :center),
          color=:gray40, fontsize=17)

    ax2 = Axis(fig[1, 2], title="Péndulo doble", xlabel="Re(λ)", ylabel="Im(λ)")
    vspan!(ax2, 0.0, 10.5, color=(LIMON, 0.35))
    vlines!(ax2, [0.0], color=:gray50, linewidth=1)
    hlines!(ax2, [0.0], color=:gray80, linewidth=0.6)
    scatter!(ax2, real.(ssd.eigenvalues), imag.(ssd.eigenvalues),
             marker=:xcross, markersize=22, color=GRISVERDE, label="lazo abierto")
    scatter!(ax2, real.(lqrd.eigenvalues_cl), imag.(lqrd.eigenvalues_cl),
             marker=:circle, markersize=18, color=MARINO, label="LQR")
    xlims!(ax2, -10.5, 10.5); ylims!(ax2, -2.9, 2.9)
    axislegend(ax2, position=:lt, framevisible=true,
               backgroundcolor=(:white, 0.85), labelsize=17)

    save(joinpath(FIGS, "slides_polos.png"), fig, px_per_unit=2)
end

# =============================================================================
# Figura 2: respuesta del pendulo simple (libre, LQR, Ackermann)
# =============================================================================

x0s = [0.0, 0.0, 0.15, 0.0]
sat = 50.0

sol_free = solve(ODEProblem(nonlinear_eom!, x0s, (0.0, 1.5), (params=ps, F=0.0)),
                 Tsit5(), saveat=0.005)
sol_lqr = solve(ODEProblem(closed_loop_eom!, copy(x0s), (0.0, 6.0),
                (params=ps, K=lqr.K, saturate=sat)), Tsit5(), saveat=0.005)
sol_ack = solve(ODEProblem(closed_loop_eom!, copy(x0s), (0.0, 6.0),
                (params=ps, K=ack.K, saturate=sat)), Tsit5(), saveat=0.005)

th_lqr = [u[3] for u in sol_lqr.u]; x_lqr = [u[1] for u in sol_lqr.u]
th_ack = [u[3] for u in sol_ack.u]; x_ack = [u[1] for u in sol_ack.u]
th_free = [u[3] for u in sol_free.u]
u_lqr = [clamp(-dot(lqr.K[1, :], u), -sat, sat) for u in sol_lqr.u]
u_ack = [clamp(-dot(ack.K[1, :], u), -sat, sat) for u in sol_ack.u]

let
    fig = Figure(size=(1040, 320))
    ax1 = Axis(fig[1, 1], title="Ángulo del péndulo", xlabel="t [s]", ylabel="θ [rad]")
    lines!(ax1, sol_free.t, th_free, color=GRISVERDE, linewidth=2.5,
           linestyle=:dash, label="sin control")
    lines!(ax1, sol_lqr.t, th_lqr, color=MARINO, linewidth=2.5, label="LQR")
    lines!(ax1, sol_ack.t, th_ack, color=PETROLEO, linewidth=2.5, label="Ackermann")
    hlines!(ax1, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    # Limitar el eje: la respuesta libre diverge y aplastaria a las controladas
    ylims!(ax1, -0.17, 0.55)
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="Posición del carro", xlabel="t [s]", ylabel="x [m]")
    lines!(ax2, sol_lqr.t, x_lqr, color=MARINO, linewidth=2.5, label="LQR")
    lines!(ax2, sol_ack.t, x_ack, color=PETROLEO, linewidth=2.5, label="Ackermann")
    hlines!(ax2, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax2, position=:rt)

    ax3 = Axis(fig[1, 3], title="Fuerza de control", xlabel="t [s]", ylabel="u [N]")
    lines!(ax3, sol_lqr.t, u_lqr, color=MARINO, linewidth=2.5, label="LQR")
    lines!(ax3, sol_ack.t, u_ack, color=PETROLEO, linewidth=2.5, label="Ackermann")
    hlines!(ax3, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax3, position=:rt)

    save(joinpath(FIGS, "slides_simple_respuesta.png"), fig, px_per_unit=2)
end

println("== SIMPLE ==")
@printf("LQR:       ts=%.2f s, max|u|=%.1f N, max|x|=%.3f m\n",
        settling_time(sol_lqr.t, th_lqr), maximum(abs.(u_lqr)), maximum(abs.(x_lqr)))
@printf("Ackermann: ts=%.2f s, max|u|=%.1f N\n",
        settling_time(sol_ack.t, th_ack), maximum(abs.(u_ack)))

# =============================================================================
# Figura 3: distintos casos, LQR frente a condiciones iniciales crecientes
# =============================================================================

theta0s = [0.10, 0.20, 0.30, 0.40]
colores = [MARINO, PETROLEO, GRISVERDE, SALVIA]

let
    fig = Figure(size=(1040, 320))
    ax1 = Axis(fig[1, 1], title="Ángulo del péndulo", xlabel="t [s]", ylabel="θ [rad]")
    ax2 = Axis(fig[1, 2], title="Fuerza de control", xlabel="t [s]", ylabel="u [N]")
    ax3 = Axis(fig[1, 3], title="Plano de fase", xlabel="θ [rad]", ylabel="ω [rad/s]")

    println("== CASOS theta0 (LQR, saturacion 50 N) ==")
    for (th0, col) in zip(theta0s, colores)
        x0 = [0.0, 0.0, th0, 0.0]
        sol = solve(ODEProblem(closed_loop_eom!, x0, (0.0, 6.0),
                    (params=ps, K=lqr.K, saturate=sat)), Tsit5(), saveat=0.005)
        th = [u[3] for u in sol.u]
        om = [u[4] for u in sol.u]
        uc = [clamp(-dot(lqr.K[1, :], u), -sat, sat) for u in sol.u]
        lab = @sprintf("θ₀ = %.2f rad", th0)
        lines!(ax1, sol.t, th, color=col, linewidth=2.5, label=lab)
        lines!(ax2, sol.t, uc, color=col, linewidth=2.5)
        lines!(ax3, th, om, color=col, linewidth=2.0)
        @printf("theta0=%.2f rad (%4.1f deg): ts=%.2f s, max|u|=%.1f N\n",
                th0, rad2deg(th0), settling_time(sol.t, th), maximum(abs.(uc)))
    end
    hlines!(ax1, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    hlines!(ax2, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    scatter!(ax3, [0.0], [0.0], marker=:star5, markersize=16, color=LIMON,
             strokecolor=:black, strokewidth=0.8)
    axislegend(ax1, position=:rt)

    save(joinpath(FIGS, "slides_simple_casos.png"), fig, px_per_unit=2)
end

# =============================================================================
# Figura 4: respuesta del pendulo doble (libre frente a LQR)
# =============================================================================

x0d = [0.0, 0.0, 0.10, 0.0, 0.10, 0.0]
satd = 100.0

solf_d = solve(ODEProblem(nonlinear_eom_double!, x0d, (0.0, 1.5), (params=pd, F=0.0)),
               Tsit5(), saveat=0.005)
soll_d = solve(ODEProblem(closed_loop_eom_double!, copy(x0d), (0.0, 6.0),
               (params=pd, K=lqrd.K, saturate=satd)), Tsit5(), saveat=0.005)

th1_l = [u[3] for u in soll_d.u]; th2_l = [u[5] for u in soll_d.u]
x_l = [u[1] for u in soll_d.u]
th1_f = [u[3] for u in solf_d.u]; th2_f = [u[5] for u in solf_d.u]
u_l = [clamp(-dot(lqrd.K[1, :], u), -satd, satd) for u in soll_d.u]

let
    fig = Figure(size=(1040, 520))
    ax1 = Axis(fig[1, 1], title="Ángulo eslabón 1", xlabel="t [s]", ylabel="θ₁ [rad]")
    lines!(ax1, solf_d.t, th1_f, color=GRISVERDE, linewidth=2.5,
           linestyle=:dash, label="sin control")
    lines!(ax1, soll_d.t, th1_l, color=MARINO, linewidth=2.5, label="LQR")
    hlines!(ax1, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="Ángulo eslabón 2", xlabel="t [s]", ylabel="θ₂ [rad]")
    lines!(ax2, solf_d.t, th2_f, color=GRISVERDE, linewidth=2.5,
           linestyle=:dash, label="sin control")
    lines!(ax2, soll_d.t, th2_l, color=MARINO, linewidth=2.5, label="LQR")
    hlines!(ax2, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax2, position=:rt)

    ax3 = Axis(fig[2, 1], title="Posición del carro", xlabel="t [s]", ylabel="x [m]")
    lines!(ax3, soll_d.t, x_l, color=MARINO, linewidth=2.5)
    hlines!(ax3, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)

    ax4 = Axis(fig[2, 2], title="Fuerza de control", xlabel="t [s]", ylabel="u [N]")
    lines!(ax4, soll_d.t, u_l, color=PETROLEO, linewidth=2.5)
    hlines!(ax4, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)

    save(joinpath(FIGS, "slides_doble_respuesta.png"), fig, px_per_unit=2)
end

println("== DOBLE ==")
@printf("LQR: ts(θ1)=%.2f s, ts(θ2)=%.2f s, max|u|=%.1f N, max|x|=%.3f m\n",
        settling_time(soll_d.t, th1_l), settling_time(soll_d.t, th2_l),
        maximum(abs.(u_l)), maximum(abs.(x_l)))

println("== Figuras guardadas en docs/presentacion/figs/ ==")
foreach(println, readdir(FIGS))
