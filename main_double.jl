# =============================================================================
# main_double.jl -- Pipeline del pendulo invertido DOBLE (Configuracion II)
# =============================================================================
# Ejecuta el flujo completo para el pendulo doble:
#   1. Definir parametros del sistema
#   2. Simular respuesta libre (sin control)
#   3. Linealizar y analizar (eigenvalores, controlabilidad, observabilidad)
#   4. Disenar controlador LQR (via Riccati)
#   5. Simular respuesta controlada (lazo cerrado)
#   6. Generar graficas comparativas
#   7. Animar el pendulo doble
#
# Reutiliza el modulo generico Controller (LQR/Riccati) y las funciones de
# analisis de Linearization, que no dependen de la dimension del estado.
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)

include("src/model_double.jl")
include("src/linearization.jl")
include("src/controller.jl")
include("src/animation_simple.jl")
include("src/animation_double.jl")

using .ModelDouble
using .Linearization
using .Controller
using .Animation
using .AnimationDouble

using DifferentialEquations
using LinearAlgebra
using Printf
using GLMakie

const FIG_DIR = joinpath(@__DIR__, "figures")
isdir(FIG_DIR) || mkdir(FIG_DIR)

# ===========================================================================
# PASO 1: Parametros del sistema
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PENDULO DOBLE -- PASO 1: PARAMETROS DEL SISTEMA")
println("=" ^ 60)

params = default_params_double()

@printf("  M  (masa carro)            = %.2f kg\n", params.M)
@printf("  m1 (masa articulacion)     = %.2f kg\n", params.m1)
@printf("  m2 (masa extremo)          = %.2f kg\n", params.m2)
@printf("  L1 (eslabon inferior)      = %.2f m\n", params.L1)
@printf("  L2 (eslabon superior)      = %.2f m\n", params.L2)
@printf("  g  (gravedad)              = %.2f m/s^2\n", params.g)

# ===========================================================================
# PASO 2: Simulacion libre (sin control)
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 2: SIMULACION LIBRE (SIN CONTROL)")
println("=" ^ 60)

# Estado inicial: ambos eslabones ligeramente desviados de la vertical
x0_free = [0.0,    # posicion del carro
           0.0,    # velocidad del carro
           0.10,   # theta1 (aprox 5.7 grados)
           0.0,    # omega1
           0.10,   # theta2
           0.0]    # omega2

tspan_free = (0.0, 3.0)
p_free = (params=params, F=0.0)

prob_free = ODEProblem(nonlinear_eom_double!, x0_free, tspan_free, p_free)
sol_free = solve(prob_free, Tsit5(), saveat=0.01)

println("  Condicion inicial: theta1 = theta2 = $(x0_free[3]) rad")
println("  Solver: Tsit5 (Runge-Kutta explicito de orden 5)")
println("  Simulacion libre completada")

# ===========================================================================
# PASO 3: Linealizacion y analisis
# ===========================================================================

println("\n")
ss = linearize_system_double(params)
print_analysis(ss)

# ===========================================================================
# PASO 4: Diseno del controlador LQR
# ===========================================================================

# Pesos del informe: se penalizan posicion y ambos angulos; el control es barato.
#   Q = diag(1, 0, 10, 0, 10, 0),  R = 0.1
Q = diagm([1.0, 0.0, 10.0, 0.0, 10.0, 0.0])
R = reshape([0.1], 1, 1)

lqr_result = design_lqr(ss.A, ss.B, Q, R)
print_controller_summary(lqr_result, method="LQR (pendulo doble)",
                         labels=["pos", "vel", "th1", "w1", "th2", "w2"])

# Nota: la asignacion de polos por Ackermann tambien aplica (el par (A,B) es
# controlable), pero elegir a mano 6 polos no es evidente; por eso, igual que en
# el informe, para el doble se usa el LQR, que selecciona los polos de forma
# optima y sistematica.

# ===========================================================================
# PASO 5: Simulacion controlada (lazo cerrado)
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 5: SIMULACION CONTROLADA (LAZO CERRADO)")
println("=" ^ 60)

tspan_ctrl = (0.0, 10.0)
saturation = 100.0  # limite del actuador [N]

p_lqr = (params=params, K=lqr_result.K, saturate=saturation)
prob_lqr = ODEProblem(closed_loop_eom_double!, copy(x0_free), tspan_ctrl, p_lqr)
sol_lqr = solve(prob_lqr, Tsit5(), saveat=0.01)

println("  Control: u = -K x (retroalimentacion de estado completo)")
println("  Saturacion del actuador: +/- $(saturation) N")
println("  Simulacion controlada completada")

# ===========================================================================
# PASO 6: Graficas comparativas
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 6: GENERANDO GRAFICAS")
println("=" ^ 60)

fig = Figure(size=(1200, 800))

# --- Angulo theta1 ---
ax1 = Axis(fig[1, 1], title="Angulo eslabon inferior (theta1)",
           xlabel="Tiempo [s]", ylabel="theta1 [rad]")
lines!(ax1, sol_free.t, [u[3] for u in sol_free.u], color=:red, linewidth=2, label="Sin control")
lines!(ax1, sol_lqr.t, [u[3] for u in sol_lqr.u], color=:blue, linewidth=2, label="Con LQR")
hlines!(ax1, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax1, position=:rt)

# --- Angulo theta2 ---
ax2 = Axis(fig[1, 2], title="Angulo eslabon superior (theta2)",
           xlabel="Tiempo [s]", ylabel="theta2 [rad]")
lines!(ax2, sol_free.t, [u[5] for u in sol_free.u], color=:red, linewidth=2, label="Sin control")
lines!(ax2, sol_lqr.t, [u[5] for u in sol_lqr.u], color=:blue, linewidth=2, label="Con LQR")
hlines!(ax2, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax2, position=:rt)

# --- Posicion del carro ---
ax3 = Axis(fig[2, 1], title="Posicion del carro",
           xlabel="Tiempo [s]", ylabel="x [m]")
lines!(ax3, sol_free.t, [u[1] for u in sol_free.u], color=:red, linewidth=2, label="Sin control")
lines!(ax3, sol_lqr.t, [u[1] for u in sol_lqr.u], color=:blue, linewidth=2, label="Con LQR")
hlines!(ax3, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax3, position=:rt)

# --- Fuerza de control ---
ax4 = Axis(fig[2, 2], title="Fuerza de control",
           xlabel="Tiempo [s]", ylabel="F [N]")
F_lqr = [clamp(-dot(lqr_result.K[1, :], sol_lqr(t)), -saturation, saturation) for t in sol_lqr.t]
lines!(ax4, sol_lqr.t, F_lqr, color=:green, linewidth=2, label="u = -K x")
hlines!(ax4, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax4, position=:rt)

path_06 = joinpath(FIG_DIR, "06_doble_comparativa_lqr.png")
save(path_06, fig, px_per_unit=2)
println("  Guardada: $path_06")

# ===========================================================================
# PASO 7: Animacion del pendulo doble
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 7: ANIMACION DEL PENDULO DOBLE")
println("=" ^ 60)

anim_data = animate_pendulum_double(sol_lqr, params,
                                    title="Pendulo invertido doble -- Control LQR",
                                    fps=30)

path_07 = joinpath(FIG_DIR, "07_doble_animacion_lqr.mp4")
save_animation(anim_data, path_07, fps=30)

println("\nPendulo doble completo. Revisa la carpeta figures/")
