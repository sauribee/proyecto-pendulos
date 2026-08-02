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
#   9. Observador de Luenberger: controlar sin medir el estado completo
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
include("src/observer.jl")
include("src/animation_simple.jl")
include("src/animation_triple.jl")

using .ModelNLink
using .ModelTriple
using .Linearization
using .Controller
using .Metrics
using .Observer
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

# Por que 150 N y no los 50 N de la Configuracion I. NO es porque la condicion
# inicial nominal lo exija: con los tres angulos a 0.10 rad los terminos de K
# alternan signo (-317, +912, -667) y se cancelan en buena medida, de modo que
# el pico real es de 7.29 N. (El argumento de que una desviacion de theta2 =
# 0.05 rad pediria unos 46 N solo por ese termino supone que una UNICA
# componente se desvia, y resulto ser una mala guia: ver la seccion "Lo que la
# implementacion corrigio" del informe tecnico de la segunda entrega.)
#
# El motivo real es que la frontera practica del triple se estanca en 0.252 rad
# a partir de unos 73 N (mapa B del atlas): con 150 N el actuador deja de ser la
# restriccion activa en todo el rango barrido, y lo que limita es la no
# linealidad. El limite es un parametro del barrido, no una constante escondida.

p_lqr = (params=params, K=lqr_result.K, saturate=saturation)
prob_lqr = ODEProblem(closed_loop_eom_triple!, copy(x0_free), tspan_ctrl, p_lqr)
sol_lqr = solve(prob_lqr, Tsit5(), saveat=0.002)

println("  Control: u = -K x (retroalimentacion de estado completo)")
println("  Saturacion del actuador: +/- $(saturation) N")
println("  Simulacion controlada completada")

# Metricas de la respuesta
# Primer instante TRAS EL CUAL |y| ya no vuelve a salirse de la banda. Es la
# misma definicion que usan los tres scripts make_*_figs.jl: devuelve t[idx+1],
# no t[idx], porque t[idx] es el ultimo instante en que la senal TODAVIA esta
# fuera de la banda. Tenerlas distintas hacia que el pipeline y el informe
# difirieran en una muestra sobre el mismo numero.
function settling_time(t, y; tol=0.02)
    idx = findlast(>(tol), abs.(y))
    return idx === nothing ? t[1] : t[min(idx + 1, length(t))]
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

# ===========================================================================
# PASO 9: Observador de Luenberger
# ===========================================================================
# Hasta aqui el control ha sido u = -K x, con el estado COMPLETO. En un sistema
# real eso no se mide: nadie pone un sensor de velocidad angular en cada
# articulacion. Como rank(Obsv) = 8, el estado se puede reconstruir, y por la
# dualidad de Kalman la ganancia del observador sale de la MISMA rutina que
# diseno K, sin escribir un algoritmo nuevo.

println("\n" * "=" ^ 60)
println("  PASO 9: OBSERVADOR DE LUENBERGER")
println("=" ^ 60)

# Que sensores hacen falta de verdad
tabla_sensores = sensor_subset_analysis(ss.A, ss.C, ["pos", "th1", "th2", "th3"])
println()
print_sensor_table(tabla_sensores, 8; label="Configuracion III")

# Diseno por dualidad. Qo grande frente a Ro significa confiar poco en el
# modelo y mucho en los sensores: el observador se hace rapido.
Qo = 100.0 * Matrix(I, 8, 8)
Ro = Matrix(I, 4, 4)
obs = design_observer(ss.A, ss.C, Qo, Ro)

println("\n  Polos del observador (A - L C):")
for lambda in sort(obs.eigenvalues_obs, by=real)
    if abs(imag(lambda)) < 1e-10
        @printf("  %+.4f\n", real(lambda))
    else
        @printf("  %+.4f %+.4fi\n", real(lambda), imag(lambda))
    end
end

# Principio de separacion, comprobado numericamente
sep = check_separation(ss.A, ss.B, ss.C, lqr_result.K, obs.L)
@printf("\n  Principio de separacion: error maximo = %.3e -> %s\n",
        sep.max_error, sep.holds ? "SE CUMPLE" : "NO SE CUMPLE")
println("  El espectro del sistema aumentado es la union de spec(A-BK) y")
println("  spec(A-LC), porque la matriz aumentada es triangular por bloques.")

# Simulacion: la planta arranca inclinada y el estimador arranca en CERO, es
# decir con el maximo desconocimiento posible sobre el estado inicial.
x0_obs = zeros(8)
x0_obs[3] = 0.10
z0 = vcat(x0_obs, zeros(8))

p_obs = (params=params, eom! = nonlinear_eom_triple!,
         A=ss.A, B=ss.B, C=ss.C, K=lqr_result.K, L=obs.L, saturate=saturation)
sol_obs = solve(ODEProblem(observer_eom!, z0, tspan_ctrl, p_obs),
                Tsit5(), saveat=0.002)

err_obs = [norm(u[1:8] .- u[9:16]) for u in sol_obs.u]
@printf("\n  Error inicial de estimacion: ||e(0)|| = %.4f\n", err_obs[1])
@printf("  Error final:                 ||e(T)|| = %.3e\n", err_obs[end])
i_1pct = findfirst(e -> e < 0.01 * err_obs[1], err_obs)
if i_1pct !== nothing
    @printf("  Cae al 1%% del inicial en t = %.2f s\n", sol_obs.t[i_1pct])
end

fig_obs = Figure(size=(1300, 800))

# --- theta1: real contra estimado (variable MEDIDA) ---
ax_o1 = Axis(fig_obs[1, 1], title="theta1 -- variable medida",
             xlabel="Tiempo [s]", ylabel="theta1 [rad]")
lines!(ax_o1, sol_obs.t, [u[3] for u in sol_obs.u], color=:black,
       linewidth=2.5, label="real")
lines!(ax_o1, sol_obs.t, [u[11] for u in sol_obs.u], color=:orange,
       linewidth=2, linestyle=:dash, label="estimado")
axislegend(ax_o1, position=:rt)

# --- omega1: real contra estimado (variable NO medida) ---
ax_o2 = Axis(fig_obs[1, 2], title="omega1 -- variable NO medida (reconstruida)",
             xlabel="Tiempo [s]", ylabel="omega1 [rad/s]")
lines!(ax_o2, sol_obs.t, [u[4] for u in sol_obs.u], color=:black,
       linewidth=2.5, label="real")
lines!(ax_o2, sol_obs.t, [u[12] for u in sol_obs.u], color=:orange,
       linewidth=2, linestyle=:dash, label="estimado")
axislegend(ax_o2, position=:rt)

# --- Norma del error de estimacion, en escala log ---
ax_o3 = Axis(fig_obs[2, 1], title="Error de estimacion ||e(t)|| = ||x - xhat||",
             xlabel="Tiempo [s]", ylabel="||e|| ", yscale=log10)
lines!(ax_o3, sol_obs.t, max.(err_obs, 1e-12), color=:purple, linewidth=2)

# --- Fuerza de control: con observador contra estado completo ---
ax_o4 = Axis(fig_obs[2, 2],
             title="Fuerza de control: u = -K xhat contra u = -K x",
             xlabel="Tiempo [s]", ylabel="F [N]")
u_obs = [clamp(-dot(lqr_result.K[1, :], u[9:16]), -saturation, saturation)
         for u in sol_obs.u]
lines!(ax_o4, sol_obs.t, u_obs, color=:green, linewidth=2,
       label="con observador")
lines!(ax_o4, sol_lqr.t, F_lqr, color=:blue, linewidth=2, linestyle=:dash,
       label="con estado completo")
axislegend(ax_o4, position=:rt)

path_11 = joinpath(FIG_DIR, "11_triple_observador.png")
save(path_11, fig_obs, px_per_unit=2)
println("\n  Guardada: $path_11")

println("\nPendulo triple completo. Revisa la carpeta figures/")
