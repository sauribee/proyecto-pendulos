# =============================================================================
# make_report_figs.jl -- Genera las figuras de respuesta temporal del informe
# =============================================================================
# Produce, en resumen_tecnico/figs/, las graficas que se incluyen en el .tex,
# y reporta metricas (tiempo de asentamiento, esfuerzo de control pico) que se
# citan en la discusion de resultados. Usa CairoMakie (salida estatica).
#
# Ejecutar desde la carpeta del proyecto:
#   julia --project=. docs/resumen_tecnico/make_report_figs.jl
# =============================================================================

using Pkg
# Este script vive en docs/resumen_tecnico/; la raiz del proyecto esta dos
# niveles arriba (docs/resumen_tecnico -> docs -> raiz).
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

# Colores
const C_LQR = RGBf(0.122, 0.349, 0.553)   # azul
const C_FREE = RGBf(0.706, 0.165, 0.165)  # rojo
const C_ACK = RGBf(0.851, 0.498, 0.114)   # naranja
const C_U = RGBf(0.180, 0.490, 0.196)     # verde

# Tiempo de asentamiento con un umbral ABSOLUTO sobre el angulo (en rad).
# Nota: tol = 0.02 rad (~1.1 grados) es una banda fija, no el 2% del valor
# inicial. Se reporta asi en el informe (umbral |theta| < 0.02 rad).
"Primer instante tras el cual |y| permanece por debajo de tol (banda absoluta)."
function settling_time(t, y; tol=0.02)
    idx = findlast(>(tol), abs.(y))
    return idx === nothing ? t[1] : t[min(idx + 1, length(t))]
end

# =============================================================================
# CONFIGURACION I: PENDULO SIMPLE
# =============================================================================

ps = default_params()
ss = linearize_system(ps)
Q = diagm([1.0, 0.0, 10.0, 0.0]); R = reshape([0.1], 1, 1)
lqr = design_lqr(ss.A, ss.B, Q, R)
ack = design_pole_placement(ss.A, ss.B, [-1.0, -2.0, -3.0, -4.0])

x0s = [0.0, 0.0, 0.15, 0.0]
sat = 50.0

sol_free = solve(ODEProblem(nonlinear_eom!, x0s, (0.0, 1.5), (params=ps, F=0.0)),
                 Tsit5(), saveat=0.005)
sol_lqr = solve(ODEProblem(closed_loop_eom!, copy(x0s), (0.0, 6.0),
                (params=ps, K=lqr.K, saturate=sat)), Tsit5(), saveat=0.005)
sol_ack = solve(ODEProblem(closed_loop_eom!, copy(x0s), (0.0, 6.0),
                (params=ps, K=ack.K, saturate=sat)), Tsit5(), saveat=0.005)

th_lqr = [u[3] for u in sol_lqr.u]; x_lqr = [u[1] for u in sol_lqr.u]
th_free = [u[3] for u in sol_free.u]
u_lqr = [clamp(-dot(lqr.K[1, :], u), -sat, sat) for u in sol_lqr.u]
th_ack = [u[3] for u in sol_ack.u]
u_ack = [clamp(-dot(ack.K[1, :], u), -sat, sat) for u in sol_ack.u]

# Figura 1: respuesta temporal en lazo cerrado (simple, LQR)
let
    fig = Figure(size=(960, 290), fontsize=13)
    ax1 = Axis(fig[1, 1], title="Angulo del pendulo", xlabel="t [s]", ylabel="theta [rad]")
    lines!(ax1, sol_free.t, th_free, color=C_FREE, linewidth=2, linestyle=:dash, label="sin control")
    lines!(ax1, sol_lqr.t, th_lqr, color=C_LQR, linewidth=2, label="con LQR")
    hlines!(ax1, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="Posicion del carro", xlabel="t [s]", ylabel="x [m]")
    lines!(ax2, sol_lqr.t, x_lqr, color=C_LQR, linewidth=2)
    hlines!(ax2, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)

    ax3 = Axis(fig[1, 3], title="Fuerza de control", xlabel="t [s]", ylabel="u [N]")
    lines!(ax3, sol_lqr.t, u_lqr, color=C_U, linewidth=2)
    hlines!(ax3, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    save(joinpath(FIGS, "simple_respuesta.png"), fig, px_per_unit=2)
end

# Figura 2: LQR vs Ackermann (simple)
let
    fig = Figure(size=(760, 300), fontsize=13)
    ax1 = Axis(fig[1, 1], title="Angulo: LQR vs Ackermann", xlabel="t [s]", ylabel="theta [rad]")
    lines!(ax1, sol_lqr.t, th_lqr, color=C_LQR, linewidth=2, label="LQR")
    lines!(ax1, sol_ack.t, th_ack, color=C_ACK, linewidth=2, label="Ackermann")
    hlines!(ax1, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="Fuerza de control", xlabel="t [s]", ylabel="u [N]")
    lines!(ax2, sol_lqr.t, u_lqr, color=C_LQR, linewidth=2, label="LQR")
    lines!(ax2, sol_ack.t, u_ack, color=C_ACK, linewidth=2, label="Ackermann")
    hlines!(ax2, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax2, position=:rt)
    save(joinpath(FIGS, "simple_lqr_vs_acker.png"), fig, px_per_unit=2)
end

println("== METRICAS SIMPLE (umbral |theta|<0.02 rad) ==")
@printf("LQR:       ts=%.2f s, max|u|=%.1f N, max|x|=%.3f m\n",
        settling_time(sol_lqr.t, th_lqr), maximum(abs.(u_lqr)), maximum(abs.(x_lqr)))
@printf("Ackermann: ts=%.2f s, max|u|=%.1f N\n",
        settling_time(sol_ack.t, th_ack), maximum(abs.(u_ack)))

# =============================================================================
# CONFIGURACION II: PENDULO DOBLE
# =============================================================================

pd = default_params_double()
ssd = linearize_system_double(pd)
Qd = diagm([1.0, 0.0, 10.0, 0.0, 10.0, 0.0]); Rd = reshape([0.1], 1, 1)
lqrd = design_lqr(ssd.A, ssd.B, Qd, Rd)

x0d = [0.0, 0.0, 0.10, 0.0, 0.10, 0.0]
satd = 100.0

solf_d = solve(ODEProblem(nonlinear_eom_double!, x0d, (0.0, 1.5), (params=pd, F=0.0)),
               Tsit5(), saveat=0.005)
soll_d = solve(ODEProblem(closed_loop_eom_double!, copy(x0d), (0.0, 6.0),
               (params=pd, K=lqrd.K, saturate=satd)), Tsit5(), saveat=0.005)

th1_l = [u[3] for u in soll_d.u]; th2_l = [u[5] for u in soll_d.u]; x_l = [u[1] for u in soll_d.u]
th1_f = [u[3] for u in solf_d.u]; th2_f = [u[5] for u in solf_d.u]
u_l = [clamp(-dot(lqrd.K[1, :], u), -satd, satd) for u in soll_d.u]

# Figura 3: respuesta temporal en lazo cerrado (doble, LQR)
let
    fig = Figure(size=(960, 560), fontsize=13)
    ax1 = Axis(fig[1, 1], title="Angulo eslabon 1", xlabel="t [s]", ylabel="theta_1 [rad]")
    lines!(ax1, solf_d.t, th1_f, color=C_FREE, linewidth=2, linestyle=:dash, label="sin control")
    lines!(ax1, soll_d.t, th1_l, color=C_LQR, linewidth=2, label="con LQR")
    hlines!(ax1, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax1, position=:rt)

    ax2 = Axis(fig[1, 2], title="Angulo eslabon 2", xlabel="t [s]", ylabel="theta_2 [rad]")
    lines!(ax2, solf_d.t, th2_f, color=C_FREE, linewidth=2, linestyle=:dash, label="sin control")
    lines!(ax2, soll_d.t, th2_l, color=C_LQR, linewidth=2, label="con LQR")
    hlines!(ax2, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    axislegend(ax2, position=:rt)

    ax3 = Axis(fig[2, 1], title="Posicion del carro", xlabel="t [s]", ylabel="x [m]")
    lines!(ax3, soll_d.t, x_l, color=C_LQR, linewidth=2)
    hlines!(ax3, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)

    ax4 = Axis(fig[2, 2], title="Fuerza de control", xlabel="t [s]", ylabel="u [N]")
    lines!(ax4, soll_d.t, u_l, color=C_U, linewidth=2)
    hlines!(ax4, [0.0], color=:gray70, linewidth=0.6, linestyle=:dot)
    save(joinpath(FIGS, "doble_respuesta.png"), fig, px_per_unit=2)
end

println("== METRICAS DOBLE (umbral |theta|<0.02 rad) ==")
@printf("LQR: ts(th1)=%.2f s, ts(th2)=%.2f s, max|u|=%.1f N, max|x|=%.3f m\n",
        settling_time(soll_d.t, th1_l), settling_time(soll_d.t, th2_l),
        maximum(abs.(u_l)), maximum(abs.(x_l)))

println("== Figuras guardadas en resumen_tecnico/figs/ ==")
foreach(println, readdir(FIGS))
