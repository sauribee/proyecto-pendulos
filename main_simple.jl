# =============================================================================
# main.jl -- Script principal del proyecto: pendulo invertido
# =============================================================================
# Ejecuta el pipeline completo:
#   1. Definir parametros del sistema
#   2. Simular respuesta libre (sin control)
#   3. Linealizar y analizar (eigenvalores, controlabilidad, observabilidad)
#   4. Disenar controlador LQR (via Riccati)
#   5. Disenar controlador por asignacion de polos (Ackermann)
#   6. Simular respuesta controlada (lazo cerrado)
#   7. Generar graficas comparativas
#   8. Animar el pendulo
# =============================================================================

# Activar el proyecto local
using Pkg
Pkg.activate(@__DIR__)

# Cargar modulos del proyecto
include("src/model.jl")
include("src/linearization.jl")
include("src/controller.jl")
include("src/animation.jl")

using .Model
using .Linearization
using .Controller
using .Animation

using DifferentialEquations
using LinearAlgebra
using Printf
using GLMakie

# Carpeta de salida para figuras
const FIG_DIR = joinpath(@__DIR__, "figures")
isdir(FIG_DIR) || mkdir(FIG_DIR)

# ===========================================================================
# PASO 1: Parametros del sistema
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 1: PARAMETROS DEL SISTEMA")
println("=" ^ 60)

params = default_params()

@printf("  M (masa carro)      = %.2f kg\n", params.M)
@printf("  m (masa pendulo)    = %.2f kg\n", params.m)
@printf("  L (dist. al CM)     = %.2f m\n", params.L)
@printf("  g (gravedad)        = %.2f m/s^2\n", params.g)
@printf("  b (friccion)        = %.2f N s/m\n", params.b)
@printf("  I (inercia pendulo) = %.4f kg m^2\n", params.I)

# ===========================================================================
# PASO 2: Simulacion libre (sin control)
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 2: SIMULACION LIBRE (SIN CONTROL)")
println("=" ^ 60)

# Estado inicial: pendulo ligeramente desviado de la vertical
x0_free = [0.0,    # posicion del carro
           0.0,    # velocidad del carro
           0.15,   # angulo (aprox 8.6 grados, desviacion pequena)
           0.0]    # velocidad angular

tspan_free = (0.0, 5.0)
p_free = (params=params, F=0.0)

prob_free = ODEProblem(nonlinear_eom!, x0_free, tspan_free, p_free)
sol_free = solve(prob_free, Tsit5(), saveat=0.01)

println("  Condicion inicial: theta0 = $(x0_free[3]) rad ($(round(rad2deg(x0_free[3]), digits=1)) grados)")
println("  Tiempo de simulacion: $(tspan_free[2]) s")
println("  Solver: Tsit5 (Runge-Kutta explicito de orden 5)")
println("  Simulacion libre completada")

# ===========================================================================
# PASO 3: Linealizacion y analisis
# ===========================================================================

println("\n")
ss = linearize_system(params)
print_analysis(ss)

# ===========================================================================
# PASO 4: Diseno del controlador LQR
# ===========================================================================

# Matrices de peso para el LQR:
#   Q penaliza las desviaciones del estado
#   R penaliza el esfuerzo de control
#
# Interpretacion:
#   Q[1,1] = 1.0  -> penaliza desviacion de posicion (moderada)
#   Q[2,2] = 0.0  -> no penaliza velocidad directamente
#   Q[3,3] = 10.0 -> penaliza fuertemente la desviacion angular
#   Q[4,4] = 0.0  -> no penaliza velocidad angular directamente
#   R = 0.1        -> esfuerzo de control relativamente barato

Q = diagm([1.0, 0.0, 10.0, 0.0])
R = reshape([0.1], 1, 1)  # debe ser matriz para la formulacion general

lqr_result = design_lqr(ss.A, ss.B, Q, R)
print_controller_summary(lqr_result, method="LQR", labels=["pos", "vel", "theta", "omega"])

# ===========================================================================
# PASO 5: Diseno por asignacion de polos (Ackermann)
# ===========================================================================

# Polos deseados en lazo cerrado (todos en el semiplano izquierdo)
desired_poles = [-1.0, -2.0, -3.0, -4.0]

acker_result = design_pole_placement(ss.A, ss.B, desired_poles)
print_controller_summary(acker_result, method="Ackermann (polos -1,-2,-3,-4)",
                         labels=["pos", "vel", "theta", "omega"])

# ===========================================================================
# PASO 6: Simulacion controlada (lazo cerrado)
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 6: SIMULACION CONTROLADA (LAZO CERRADO)")
println("=" ^ 60)

tspan_ctrl = (0.0, 10.0)
saturation = 50.0  # limite del actuador [N]

# Lazo cerrado con LQR
p_lqr = (params=params, K=lqr_result.K, saturate=saturation)
prob_lqr = ODEProblem(closed_loop_eom!, copy(x0_free), tspan_ctrl, p_lqr)
sol_lqr = solve(prob_lqr, Tsit5(), saveat=0.01)

# Lazo cerrado con Ackermann
p_acker = (params=params, K=acker_result.K, saturate=saturation)
prob_acker = ODEProblem(closed_loop_eom!, copy(x0_free), tspan_ctrl, p_acker)
sol_acker = solve(prob_acker, Tsit5(), saveat=0.01)

println("  Control: u = -K x (retroalimentacion de estado completo)")
println("  Saturacion del actuador: +/- $(saturation) N")
println("  Simulacion controlada completada (LQR y Ackermann)")

# ===========================================================================
# PASO 7: Graficas comparativas
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 7: GENERANDO GRAFICAS")
println("=" ^ 60)

# --- Figura 01: comparativa LQR (4 paneles: libre vs LQR) ---
fig_lqr = Figure(size=(1200, 800))

ax1 = Axis(fig_lqr[1, 1], title="Angulo del pendulo",
           xlabel="Tiempo [s]", ylabel="theta [rad]")
lines!(ax1, sol_free.t, [u[3] for u in sol_free.u], color=:red, linewidth=2, label="Sin control")
lines!(ax1, sol_lqr.t, [u[3] for u in sol_lqr.u], color=:blue, linewidth=2, label="Con LQR")
hlines!(ax1, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax1, position=:rt)

ax2 = Axis(fig_lqr[1, 2], title="Posicion del carro",
           xlabel="Tiempo [s]", ylabel="x [m]")
lines!(ax2, sol_free.t, [u[1] for u in sol_free.u], color=:red, linewidth=2, label="Sin control")
lines!(ax2, sol_lqr.t, [u[1] for u in sol_lqr.u], color=:blue, linewidth=2, label="Con LQR")
hlines!(ax2, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax2, position=:rt)

ax3 = Axis(fig_lqr[2, 1], title="Fuerza de control",
           xlabel="Tiempo [s]", ylabel="F [N]")
F_lqr = [clamp(-dot(lqr_result.K[1, :], sol_lqr(t)), -saturation, saturation) for t in sol_lqr.t]
lines!(ax3, sol_lqr.t, F_lqr, color=:green, linewidth=2, label="u = -K x")
hlines!(ax3, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax3, position=:rt)

ax4 = Axis(fig_lqr[2, 2], title="Velocidades (lazo cerrado LQR)",
           xlabel="Tiempo [s]", ylabel="Velocidad")
lines!(ax4, sol_lqr.t, [u[2] for u in sol_lqr.u], color=:purple, linewidth=2, label="vel [m/s]")
lines!(ax4, sol_lqr.t, [u[4] for u in sol_lqr.u], color=:orange, linewidth=2, label="omega [rad/s]")
hlines!(ax4, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax4, position=:rt)

path_01 = joinpath(FIG_DIR, "01_comparativa_lqr.png")
save(path_01, fig_lqr, px_per_unit=2)
println("  Guardada: $path_01")

# --- Figura 02: comparativa LQR vs Ackermann ---
fig_cmp = Figure(size=(1200, 500))

axA = Axis(fig_cmp[1, 1], title="Angulo del pendulo (LQR vs Ackermann)",
           xlabel="Tiempo [s]", ylabel="theta [rad]")
lines!(axA, sol_lqr.t, [u[3] for u in sol_lqr.u], color=:blue, linewidth=2, label="LQR")
lines!(axA, sol_acker.t, [u[3] for u in sol_acker.u], color=:darkorange, linewidth=2, label="Ackermann")
hlines!(axA, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(axA, position=:rt)

axB = Axis(fig_cmp[1, 2], title="Posicion del carro (LQR vs Ackermann)",
           xlabel="Tiempo [s]", ylabel="x [m]")
lines!(axB, sol_lqr.t, [u[1] for u in sol_lqr.u], color=:blue, linewidth=2, label="LQR")
lines!(axB, sol_acker.t, [u[1] for u in sol_acker.u], color=:darkorange, linewidth=2, label="Ackermann")
hlines!(axB, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(axB, position=:rt)

path_02 = joinpath(FIG_DIR, "02_comparativa_ackermann.png")
save(path_02, fig_cmp, px_per_unit=2)
println("  Guardada: $path_02")

# ===========================================================================
# PASO 8: Animacion del pendulo
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 8: ANIMACION DEL PENDULO")
println("=" ^ 60)

anim_data = animate_pendulum(sol_lqr, params,
                             title="Pendulo invertido -- Control LQR",
                             fps=30)

path_03 = joinpath(FIG_DIR, "03_animacion_lqr.mp4")
save_animation(anim_data, path_03, fps=30)

println("\nProyecto completo. Revisa la carpeta figures/")
