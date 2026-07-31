# =============================================================================
# main_triple.jl -- Pipeline del pendulo invertido TRIPLE (Configuracion III)
# =============================================================================
# Ejecuta el flujo completo para el pendulo triple:
#   1. Definir parametros del sistema
#   2. Simular respuesta libre (sin control)
#   3. Linealizar y analizar (eigenvalores, controlabilidad, observabilidad)
#   4. Metricas de viabilidad (margen PBH, condicionamiento, elipsoide)
#   5. Disenar controlador LQR (via Riccati)
#   6. Simular respuesta controlada (lazo cerrado)
#   7. Generar graficas comparativas
#   8. Animar el pendulo triple
#
# Reutiliza sin cambios los modulos genericos Controller (LQR/Riccati), las
# funciones de analisis de Linearization y las metricas de Metrics: ninguno
# depende de la dimension del estado, asi que el mismo codigo que analizo el
# sistema de R^4 y el de R^6 analiza el de R^8.
#
# La fisica no se duplica: ModelTriple es una capa delgada sobre ModelNLink.
# =============================================================================

using Pkg
Pkg.activate(@__DIR__)

include("src/model_nlink.jl")
include("src/model_triple.jl")
include("src/linearization.jl")
include("src/controller.jl")
include("src/metrics.jl")
include("src/animation_simple.jl")
include("src/animation_triple.jl")

using .ModelNLink
using .ModelTriple
using .Linearization
using .Controller
using .Metrics
using .Animation
using .AnimationTriple

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
println("  PENDULO TRIPLE -- PASO 1: PARAMETROS DEL SISTEMA")
println("=" ^ 60)

params = default_params_triple()
bc = link_couplings(params)

@printf("  M  (masa carro)            = %.2f kg\n", params.M)
for j in 1:3
    @printf("  m%d (masa eslabon %d)        = %.2f kg\n", j, j, params.m[j])
end
for j in 1:3
    @printf("  l%d (longitud eslabon %d)    = %.4f m\n", j, j, params.l[j])
end
@printf("  g  (gravedad)              = %.2f m/s^2\n", params.g)
@printf("  b  (friccion del carro)    = %.2f N s/m\n", params.b)
println()
@printf("  masa total del pendulo     = %.2f kg   (igual que la Config. II)\n",
        sum(params.m))
@printf("  longitud total             = %.2f m    (igual que I y II)\n",
        sum(params.l))
println()
println("  beta (momentos estaticos)  = ", round.(bc.beta, digits=5))
println("  J    (inercias efectivas)  = ", round.(bc.J, digits=5))

# ===========================================================================
# PASO 2: Simulacion libre (sin control)
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 2: SIMULACION LIBRE (SIN CONTROL)")
println("=" ^ 60)

# Estado inicial: los tres eslabones ligeramente desviados de la vertical
x0_free = [0.0,    # posicion del carro
           0.0,    # velocidad del carro
           0.10,   # theta1 (aprox 5.7 grados)
           0.0,    # omega1
           0.10,   # theta2
           0.0,    # omega2
           0.10,   # theta3
           0.0]    # omega3

tspan_free = (0.0, 3.0)
p_free = (params=params, F=0.0)

prob_free = ODEProblem(nonlinear_eom_triple!, x0_free, tspan_free, p_free)
sol_free = solve(prob_free, Tsit5(), saveat=0.002)

println("  Condicion inicial: theta1 = theta2 = theta3 = $(x0_free[3]) rad")
println("  Solver: Tsit5 (Runge-Kutta explicito de orden 5)")
println("  Paso de guardado: 0.002 s (el triple es mas rigido que I y II)")
println("  Simulacion libre completada")

# ===========================================================================
# PASO 3: Linealizacion y analisis
# ===========================================================================

println("\n")
ss = linearize_system_triple(params)
print_analysis(ss)

# ===========================================================================
# PASO 4: Diseno del controlador LQR
# ===========================================================================

# Pesos del informe, extendidos: se penalizan la posicion y los tres angulos.
#   Q = diag(1, 0, 10, 0, 10, 0, 10, 0),  R = 0.1
Q = diagm([1.0, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0, 0.0])
R = reshape([0.1], 1, 1)

lqr_result = design_lqr(ss.A, ss.B, Q, R)
print_controller_summary(lqr_result, method="LQR (pendulo triple)",
                         labels=["pos", "vel", "th1", "w1", "th2", "w2",
                                 "th3", "w3"])

# Con R pequeno el LQR tiende a REFLEJAR los polos inestables al semiplano
# izquierdo: el polo mas rapido de lazo cerrado (-18.14) es practicamente el
# espejo del modo inestable mas rapido de lazo abierto (+18.06). Es una
# propiedad conocida del regulador barato.

# ===========================================================================
# PASO 5: Metricas de viabilidad
# ===========================================================================

println()
print_metrics_summary(ss.A, ss.B, ss.C, lqr_result.K, lqr_result.P, R;
                      umax=50.0, label="Configuracion III (triple)")

println("  Nota: el angulo garantizado se calcula con umax = 50 N, el limite")
println("  de la Configuracion I, para que la comparacion entre las tres")
println("  configuraciones sea directa. El pipeline simula con 150 N.")

# ===========================================================================
# PASO 6: Simulacion controlada (lazo cerrado)
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 6: SIMULACION CONTROLADA (LAZO CERRADO)")
println("=" ^ 60)

tspan_ctrl = (0.0, 10.0)
saturation = 150.0  # limite del actuador [N]

# Los 50 N de la Configuracion I son insuficientes aqui: con norm(K) del orden
# de 1176, una desviacion de theta2 = 0.05 rad ya pide unos 46 N solo por ese
# termino.

p_lqr = (params=params, K=lqr_result.K, saturate=saturation)
prob_lqr = ODEProblem(closed_loop_eom_triple!, copy(x0_free), tspan_ctrl, p_lqr)
sol_lqr = solve(prob_lqr, Tsit5(), saveat=0.002)

println("  Control: u = -K x (retroalimentacion de estado completo)")
println("  Saturacion del actuador: +/- $(saturation) N")
println("  Simulacion controlada completada")

# Metricas de la respuesta
function settling_time(t, y; tol=0.02)
    idx = findlast(abs.(y) .> tol)
    return idx === nothing ? 0.0 : t[idx]
end

th_lqr = [[u[2j+1] for u in sol_lqr.u] for j in 1:3]
x_lqr = [u[1] for u in sol_lqr.u]
F_lqr = [clamp(-dot(lqr_result.K[1, :], u), -saturation, saturation)
         for u in sol_lqr.u]

println()
for j in 1:3
    @printf("  t_asentamiento theta%d = %.2f s\n", j,
            settling_time(sol_lqr.t, th_lqr[j]))
end
@printf("  fuerza pico            = %.2f N\n", maximum(abs.(F_lqr)))
@printf("  excursion del carro    = %.3f m\n", maximum(abs.(x_lqr)))
@printf("  se satura el actuador  = %s\n",
        maximum(abs.(F_lqr)) >= saturation - 1e-6 ? "SI" : "no")

# ===========================================================================
# PASO 7: Graficas comparativas
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 7: GENERANDO GRAFICAS")
println("=" ^ 60)

fig = Figure(size=(1400, 900))

for j in 1:3
    fila = j <= 2 ? 1 : 2
    col = j <= 2 ? j : 1
    ax = Axis(fig[fila, col],
              title="Angulo eslabon $j (theta$j) -- la respuesta libre sale del marco",
              xlabel="Tiempo [s]", ylabel="theta$j [rad]")
    lines!(ax, sol_free.t, [u[2j+1] for u in sol_free.u],
           color=:red, linewidth=2, label="Sin control")
    lines!(ax, sol_lqr.t, th_lqr[j], color=:blue, linewidth=2, label="Con LQR")
    hlines!(ax, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
    axislegend(ax, position=:rt)

    # La escala se fija a la respuesta CONTROLADA. El pendulo triple libre da
    # varias vueltas completas en 3 s (theta3 llega a 24 rad), asi que con la
    # escala automatica la curva del LQR quedaria aplastada contra el cero y la
    # grafica no mostraria nada.
    ylim = max(1.8 * maximum(abs.(th_lqr[j])), 0.25)
    ylims!(ax, -ylim, ylim)
end

# --- Posicion del carro ---
ax_pos = Axis(fig[2, 2], title="Posicion del carro",
              xlabel="Tiempo [s]", ylabel="x [m]")
lines!(ax_pos, sol_free.t, [u[1] for u in sol_free.u],
       color=:red, linewidth=2, label="Sin control")
lines!(ax_pos, sol_lqr.t, x_lqr, color=:blue, linewidth=2, label="Con LQR")
hlines!(ax_pos, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)
axislegend(ax_pos, position=:rt)

# --- Fuerza de control ---
F_peak = maximum(abs.(F_lqr))
ax_F = Axis(fig[3, 1:2],
            title=@sprintf("Fuerza de control -- pico %.2f N de %.0f N disponibles (%.1f%%)",
                           F_peak, saturation, 100 * F_peak / saturation),
            xlabel="Tiempo [s]", ylabel="F [N]")
lines!(ax_F, sol_lqr.t, F_lqr, color=:green, linewidth=2, label="u = -K x")
hlines!(ax_F, [0.0], color=:gray, linestyle=:dash, linewidth=0.5)

# Las lineas de saturacion solo se dibujan si la senal se acerca a ellas. Con
# la condicion inicial nominal el pico es de 7.29 N contra un limite de 150 N:
# forzar el eje hasta +/-150 dejaria la senal de control como una linea plana.
if F_peak > 0.5 * saturation
    hlines!(ax_F, [saturation, -saturation], color=:orange, linestyle=:dot,
            linewidth=1.5, label="saturacion")
end
axislegend(ax_F, position=:rt)

path_09 = joinpath(FIG_DIR, "09_triple_comparativa_lqr.png")
save(path_09, fig, px_per_unit=2)
println("  Guardada: $path_09")

# ===========================================================================
# PASO 8: Animacion del pendulo triple
# ===========================================================================

println("\n" * "=" ^ 60)
println("  PASO 8: ANIMACION DEL PENDULO TRIPLE")
println("=" ^ 60)

anim_data = animate_pendulum_triple(sol_lqr, params,
                                    title="Pendulo invertido triple -- Control LQR",
                                    fps=30)

path_10 = joinpath(FIG_DIR, "10_triple_animacion_lqr.mp4")
save_animation(anim_data, path_10, fps=30)

println("\nPendulo triple completo. Revisa la carpeta figures/")
