# =============================================================================
# make_report_figs.jl -- Figuras del informe tecnico de la SEGUNDA ENTREGA
# =============================================================================
# Produce, en Entrega_2/resumen_tecnico/figs/, las graficas propias de este
# informe, y reporta por consola las metricas que se citan en la discusion.
# Usa CairoMakie (salida estatica).
#
# Los cuatro mapas del atlas de operabilidad NO se regeneran aqui: los produce
# docs/Entrega_2/atlas/make_atlas_figs.jl y el .tex los toma de esa carpeta
# mediante \graphicspath. Duplicar ese calculo no aportaria nada.
#
# Ejecutar desde la carpeta del proyecto:
#   julia --project=. -t auto docs/Entrega_2/resumen_tecnico/make_report_figs.jl
# =============================================================================

using Pkg
# Este script vive en docs/Entrega_2/resumen_tecnico/; la raiz del proyecto
# esta tres niveles arriba (resumen_tecnico -> Entrega_2 -> docs -> raiz).
const PROJ_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
Pkg.activate(PROJ_ROOT)

include(joinpath(PROJ_ROOT, "src", "model_nlink.jl"))
include(joinpath(PROJ_ROOT, "src", "model_triple.jl"))
include(joinpath(PROJ_ROOT, "src", "linearization.jl"))
include(joinpath(PROJ_ROOT, "src", "controller.jl"))
include(joinpath(PROJ_ROOT, "src", "metrics.jl"))
include(joinpath(PROJ_ROOT, "src", "observer.jl"))
include(joinpath(PROJ_ROOT, "src", "sweep.jl"))

using .ModelNLink, .ModelTriple, .Linearization, .Controller
using .Metrics, .Observer, .Sweep
using LinearAlgebra
using Printf
using DifferentialEquations
using CairoMakie

const FIGS = joinpath(@__DIR__, "figs")
isdir(FIGS) || mkpath(FIGS)

# Paleta del informe (la misma de Entrega_1/resumen_tecnico/make_report_figs.jl)
const C_LQR  = RGBf(0.122, 0.349, 0.553)   # azul
const C_FREE = RGBf(0.706, 0.165, 0.165)   # rojo
const C_ACK  = RGBf(0.851, 0.498, 0.114)   # naranja
const C_U    = RGBf(0.180, 0.490, 0.196)   # verde
const C_GRIS = RGBf(0.45, 0.45, 0.45)

"Primer instante tras el cual |y| permanece por debajo de tol (banda absoluta)."
function settling_time(t, y; tol=0.02)
    idx = findlast(>(tol), abs.(y))
    return idx === nothing ? t[1] : t[min(idx + 1, length(t))]
end

const R_STD = reshape([0.1], 1, 1)

# Las tres configuraciones, todas como instancias del modelo generico.
const CONFIGS = [
    ("Simple", uniform_rods(1.0, [0.3], [1.0]; b=0.1)),
    ("Doble",  point_masses(1.0, [0.3, 0.3], [0.5, 0.5]; b=0.0)),
    ("Triple", default_params_triple()),
]

println("=" ^ 74)
println("  FIGURAS DEL INFORME TECNICO -- SEGUNDA ENTREGA")
println("=" ^ 74)

# =============================================================================
# FIGURA 1: respuesta temporal de la Configuracion III en lazo cerrado
# =============================================================================

p3 = default_params_triple()
ss3 = linearize_system_triple(p3)
lqr3 = design_lqr(ss3.A, ss3.B, default_weights(8), R_STD)
K3 = lqr3.K
const SAT = 150.0

x0 = [0.0, 0.0, 0.10, 0.0, 0.10, 0.0, 0.10, 0.0]

sol_free = solve(ODEProblem(nonlinear_eom_triple!, copy(x0), (0.0, 3.0),
                            (params=p3, F=0.0)), Tsit5(), saveat=0.002)
sol_lqr = solve(ODEProblem(closed_loop_eom_triple!, copy(x0), (0.0, 10.0),
                           (params=p3, K=K3, saturate=SAT)), Tsit5(), saveat=0.002)

th_lqr = [[u[2j+1] for u in sol_lqr.u] for j in 1:3]
x_lqr = [u[1] for u in sol_lqr.u]
u_lqr = [clamp(-dot(K3[1, :], u), -SAT, SAT) for u in sol_lqr.u]

fig1 = Figure(size=(1000, 720))

ax = Axis(fig1[1, 1], title="Angulos de los tres eslabones (LQR)",
          xlabel="Tiempo [s]", ylabel="theta [rad]")
for (j, col) in enumerate((C_LQR, C_ACK, C_U))
    lines!(ax, sol_lqr.t, th_lqr[j], color=col, linewidth=2, label="theta$j")
end
hlines!(ax, [0.0], color=C_GRIS, linestyle=:dash, linewidth=0.6)
axislegend(ax, position=:rt)

ax2 = Axis(fig1[1, 2], title="Sin control: el triple libre se desploma",
           xlabel="Tiempo [s]", ylabel="theta [rad]")
for (j, col) in enumerate((C_LQR, C_ACK, C_U))
    lines!(ax2, sol_free.t, [u[2j+1] for u in sol_free.u], color=col,
           linewidth=2, label="theta$j")
end
hlines!(ax2, [0.0], color=C_GRIS, linestyle=:dash, linewidth=0.6)
axislegend(ax2, position=:lt)

ax3 = Axis(fig1[2, 1], title="Posicion del carro",
           xlabel="Tiempo [s]", ylabel="x [m]")
lines!(ax3, sol_lqr.t, x_lqr, color=C_LQR, linewidth=2)
hlines!(ax3, [0.0], color=C_GRIS, linestyle=:dash, linewidth=0.6)

pico = maximum(abs.(u_lqr))
ax4 = Axis(fig1[2, 2],
           title=@sprintf("Fuerza de control (pico %.2f N de %.0f N)", pico, SAT),
           xlabel="Tiempo [s]", ylabel="F [N]")
lines!(ax4, sol_lqr.t, u_lqr, color=C_U, linewidth=2)
hlines!(ax4, [0.0], color=C_GRIS, linestyle=:dash, linewidth=0.6)

save(joinpath(FIGS, "triple_respuesta.png"), fig1, px_per_unit=2)
println("\n  Guardada: triple_respuesta.png")
for j in 1:3
    @printf("    t_s(theta%d) = %.2f s\n", j, settling_time(sol_lqr.t, th_lqr[j]))
end
@printf("    fuerza pico = %.2f N   excursion del carro = %.3f m\n",
        pico, maximum(abs.(x_lqr)))

# =============================================================================
# FIGURA 2: degradacion con N -- el rango sobrevive, el condicionamiento no
# =============================================================================

Ns = 1:5
pbhs = Float64[]
conds = Float64[]
lams = Float64[]
normKs = Float64[]
for N in Ns
    p = uniform_rods(1.0, fill(0.6/N, N), fill(1.0/N, N); b=0.1)
    v = viability(linearize_system_nlink, p)
    push!(pbhs, v.pbh_normalized)
    push!(conds, v.cond_ctrb)
    push!(lams, v.lambda_max)
    push!(normKs, v.norm_K)
end

fig2 = Figure(size=(1000, 400))

axa = Axis(fig2[1, 1], title="El criterio de rango deja de ser confiable",
           xlabel="Numero de eslabones N", ylabel="cond(C) de Kalman",
           yscale=log10, xticks=1:5)
lines!(axa, collect(Ns), conds, color=C_FREE, linewidth=2.5)
scatter!(axa, collect(Ns), conds, color=C_FREE, markersize=12)
hlines!(axa, [1 / eps(Float64)], color=C_GRIS, linestyle=:dash, linewidth=1.5)
text!(axa, 1.1, 1.3 / eps(Float64), text="1 / eps de maquina", fontsize=11,
      color=C_GRIS)

axb = Axis(fig2[1, 2], title="El margen PBH si es utilizable",
           xlabel="Numero de eslabones N",
           ylabel="margen PBH normalizado", yscale=log10, xticks=1:5)
lines!(axb, collect(Ns), pbhs, color=C_LQR, linewidth=2.5)
scatter!(axb, collect(Ns), pbhs, color=C_LQR, markersize=12)

save(joinpath(FIGS, "degradacion_N.png"), fig2, px_per_unit=2)
println("\n  Guardada: degradacion_N.png")
@printf("    cond(C): %.1e (N=1) -> %.1e (N=5), factor %.1e\n",
        conds[1], conds[5], conds[5]/conds[1])
@printf("    PBH:     %.2e (N=1) -> %.2e (N=5), factor %.0f\n",
        pbhs[1], pbhs[5], pbhs[1]/pbhs[5])

# =============================================================================
# FIGURA 3: espectros de las tres configuraciones en el plano complejo
# =============================================================================

fig3 = Figure(size=(1000, 400))
axc = Axis(fig3[1, 1], title="Lazo abierto: N modos inestables",
           xlabel="Re(lambda) [1/s]", ylabel="Im(lambda) [1/s]")
axd = Axis(fig3[1, 2], title="Lazo cerrado con LQR: todos a la izquierda",
           xlabel="Re(lambda) [1/s]", ylabel="Im(lambda) [1/s]")

marcas = (:circle, :rect, :utriangle)
colores = (C_LQR, C_ACK, C_FREE)
for (i, (nombre, p)) in enumerate(CONFIGS)
    ss = linearize_system_nlink(p)
    n = size(ss.A, 1)
    K = design_lqr(ss.A, ss.B, default_weights(n), R_STD).K
    ol = eigvals(ss.A)
    cl = eigvals(ss.A - ss.B * K)
    scatter!(axc, real.(ol), imag.(ol), marker=marcas[i], markersize=13,
             color=colores[i], label=nombre)
    scatter!(axd, real.(cl), imag.(cl), marker=marcas[i], markersize=13,
             color=colores[i], label=nombre)
end
for ax in (axc, axd)
    vlines!(ax, [0.0], color=C_GRIS, linestyle=:dash, linewidth=1)
    axislegend(ax, position=:rt)
end

save(joinpath(FIGS, "espectros_tres.png"), fig3, px_per_unit=2)
println("\n  Guardada: espectros_tres.png")

# =============================================================================
# FIGURA 4: observador -- reconstruir lo que no se mide
# =============================================================================

obs = design_observer(ss3.A, ss3.C, 100.0 * Matrix(I, 8, 8), Matrix(I, 4, 4))

# ATENCION a la condicion inicial. Aqui se desvia SOLO theta1, no los tres
# angulos como en la figura de respuesta. El motivo es un resultado por derecho
# propio: el lazo con observador tiene una region de atraccion bastante menor
# que el de estado completo, y desde los tres angulos a 0.10 rad DIVERGE. El
# principio de separacion garantiza estabilidad del sistema LINEAL; no dice
# nada sobre la region de atraccion del no lineal. Ver la comparacion que este
# script imprime mas abajo.
x0_obs = zeros(8)
x0_obs[3] = 0.10

const P_OBS = (params=p3, eom! = nonlinear_eom_triple!,
               A=ss3.A, B=ss3.B, C=ss3.C, K=K3, L=obs.L, saturate=SAT)
sol_obs = solve(ODEProblem(observer_eom!, vcat(x0_obs, zeros(8)), (0.0, 10.0),
                           P_OBS), Tsit5(), saveat=0.002)
err = [norm(u[1:8] .- u[9:16]) for u in sol_obs.u]

fig4 = Figure(size=(1000, 400))

axe = Axis(fig4[1, 1], title="omega1: variable NO medida, reconstruida",
           xlabel="Tiempo [s]", ylabel="omega1 [rad/s]")
lines!(axe, sol_obs.t, [u[4] for u in sol_obs.u], color=:black, linewidth=2.5,
       label="real")
lines!(axe, sol_obs.t, [u[12] for u in sol_obs.u], color=C_ACK, linewidth=2,
       linestyle=:dash, label="estimado")
axislegend(axe, position=:rb)

axf = Axis(fig4[1, 2], title="El error de estimacion decae exponencialmente",
           xlabel="Tiempo [s]", ylabel="||x - xhat||", yscale=log10)
lines!(axf, sol_obs.t, max.(err, 1e-12), color=C_LQR, linewidth=2.5)

save(joinpath(FIGS, "observador.png"), fig4, px_per_unit=2)
println("\n  Guardada: observador.png")
i1 = findfirst(e -> e < 0.01 * err[1], err)
@printf("    ||e(0)|| = %.4f  ->  ||e(T)|| = %.2e\n", err[1], err[end])
@printf("    cae al 1%% en t = %.2f s\n", i1 === nothing ? NaN : sol_obs.t[i1])

# ---------------------------------------------------------------------------
# El precio de estimar: cuanto se encoge la region de atraccion
# ---------------------------------------------------------------------------
"Recupera el lazo CON observador, arrancando el estimador en cero?"
function recovers_with_observer(theta0; n_ang=1, T=10.0)
    z0 = zeros(16)
    for j in 1:n_ang
        z0[2j+1] = theta0
    end
    s = try
        solve(ODEProblem(observer_eom!, z0, (0.0, T), P_OBS), Tsit5();
              saveat=0.01, save_everystep=false, verbose=false, maxiters=10^7)
    catch
        return false
    end
    s.retcode == ReturnCode.Success || return false
    for u in s.u
        all(abs(u[j]) < pi/2 for j in 3:2:7) || return false
        abs(u[1]) <= 1.5 || return false
    end
    return all(abs(s.u[end][j]) < 0.02 for j in 3:2:7)
end

function bisect_theta(pred; lo=0.0, hi=0.5, iters=16)
    pred(hi) && return hi
    for _ in 1:iters
        mid = (lo + hi) / 2
        pred(mid) ? (lo = mid) : (hi = mid)
    end
    return lo
end

println("\n  Region de atraccion: estado completo contra observador")
@printf("  %-22s %16s %16s %8s\n", "angulos desviados", "estado completo",
        "con observador", "ratio")
for n_ang in 1:3
    th_obs = bisect_theta(t -> recovers_with_observer(t; n_ang=n_ang))
    th_full = bisect_theta(function (t)
        x = zeros(8)
        for j in 1:n_ang
            x[2j+1] = t
        end
        s = try
            solve(ODEProblem(closed_loop_eom_nlink!, x, (0.0, 10.0),
                             (params=p3, K=K3, saturate=SAT)), Tsit5();
                  saveat=0.01, save_everystep=false, verbose=false,
                  maxiters=10^7)
        catch
            return false
        end
        s.retcode == ReturnCode.Success || return false
        for u in s.u
            all(abs(u[j]) < pi/2 for j in 3:2:7) || return false
            abs(u[1]) <= 1.5 || return false
        end
        return all(abs(s.u[end][j]) < 0.02 for j in 3:2:7)
    end)
    @printf("  %-22d %16.4f %16.4f %7.0f%%\n", n_ang, th_full, th_obs,
            100 * th_obs / th_full)
end
println("  El principio de separacion garantiza estabilidad del sistema LINEAL;")
println("  no dice nada sobre la region de atraccion del no lineal, que se")
println("  encoge notablemente al estimar en vez de medir.")

# ---------------------------------------------------------------------------
# POR QUE se encoge: la comparacion tiene que ser PAREADA
# ---------------------------------------------------------------------------
# OJO: comparar el pico de fuerza del lazo con observador contra el del lazo con
# estado completo SOLO tiene sentido desde la MISMA condicion inicial. Mezclar
# la CI de tres angulos (donde los terminos de K se cancelan y el pico baja a
# 7.29 N) con la de un angulo lleva a concluir que el observador dispara la
# fuerza, que es justo lo contrario de lo que pasa: el estimador arranca en
# cero, luego u(0) = 0, y el pico resulta MENOR. Lo que crece es la EXCURSION.

"Pico de fuerza y excursion angular maxima, con estado completo."
function full_metrics(theta0; n_ang=1, T=10.0)
    x0 = zeros(8)
    for j in 1:n_ang
        x0[2j+1] = theta0
    end
    s = solve(ODEProblem(closed_loop_eom_nlink!, x0, (0.0, T),
                         (params=p3, K=K3, saturate=SAT)), Tsit5(), saveat=0.002)
    u = [clamp(-dot(K3[1, :], v), -SAT, SAT) for v in s.u]
    return (pico=maximum(abs.(u)),
            exc=maximum(maximum(abs(v[j]) for j in 3:2:7) for v in s.u))
end

"Idem con observador arrancando en cero."
function obs_metrics(theta0; n_ang=1, T=10.0)
    z0 = zeros(16)
    for j in 1:n_ang
        z0[2j+1] = theta0
    end
    s = solve(ODEProblem(observer_eom!, z0, (0.0, T), P_OBS), Tsit5(), saveat=0.002)
    u = [clamp(-dot(K3[1, :], v[9:16]), -SAT, SAT) for v in s.u]
    return (pico=maximum(abs.(u)),
            exc=maximum(maximum(abs(v[j]) for j in 3:2:7) for v in s.u))
end

println("\n  Comparacion PAREADA (misma condicion inicial en ambas ramas):")
@printf("  %-24s %9s %9s %9s %9s %8s\n", "condicion inicial",
        "pico[N]", "pico_o[N]", "exc[rad]", "exc_o[rad]", "factor")
for (n_ang, th) in [(1, 0.10), (1, 0.1452), (2, 0.0450), (3, 0.0743)]
    f = full_metrics(th; n_ang=n_ang)
    o = obs_metrics(th; n_ang=n_ang)
    @printf("  %-24s %9.2f %9.2f %9.4f %9.4f %8.2f\n",
            "$(n_ang) angulo(s) a $(th)", f.pico, o.pico, f.exc, o.exc,
            o.exc / f.exc)
end
@printf("\n  u(0) con estado completo = -K[3]*0.10 = %.2f N ; con observador = 0\n",
        -K3[1, 3] * 0.10)
i_1 = findfirst(e -> e < 0.01 * err[1], err)
@printf("  el error de estimacion cae al 1%% en t = %.2f s\n",
        i_1 === nothing ? NaN : sol_obs.t[i_1])
println("  El observador NO dispara la fuerza: la baja. Lo que dispara es la")
println("  excursion angular, y es eso lo que saca a la planta del regimen")
println("  donde la linealizacion vale.")

# =============================================================================
# TABLA COMPARATIVA que se cita en el informe
# =============================================================================

println("\n" * "=" ^ 74)
println("  TABLA COMPARATIVA DE LAS TRES CONFIGURACIONES")
println("=" ^ 74)
@printf("  %-8s %3s %6s %9s %9s %11s %11s %9s %8s %8s\n",
        "config", "n", "inest", "lam_max", "tau[ms]", "cond(C)", "PBH",
        "norm(K)", "grados", "h[ms]")
for (nombre, p) in CONFIGS
    ss = linearize_system_nlink(p)
    n = size(ss.A, 1)
    lqr = design_lqr(ss.A, ss.B, default_weights(n), R_STD)
    dm = dominant_mode(ss.A)
    cs = nonsaturation_level(lqr.K, lqr.P, R_STD, ss.B, 50.0)
    h = 1000 * max_sampling_period(ss.A, ss.B, lqr.K)
    @printf("  %-8s %3d %6d %9.2f %9.1f %11.2e %11.2e %9.1f %8.1f %8.1f\n",
            nombre, n, dm.n_unstable, dm.lambda_max, 1000dm.tau_fall,
            controllability_condition(ss.A, ss.B),
            pbh_controllability_margin_normalized(ss.A, ss.B),
            norm(lqr.K), rad2deg(max_axis_angle(lqr.P, cs, 3)), h)
end

println("\n  Fronteras practicas (biseccion sobre el modelo no lineal):")
@printf("  %-8s %10s %12s %12s %12s\n",
        "config", "cota c*", "frontera", "sin riel", "restriccion")
for (nombre, p) in CONFIGS
    ss = linearize_system_nlink(p)
    n = size(ss.A, 1)
    lqr = design_lqr(ss.A, ss.B, default_weights(n), R_STD)
    cs = nonsaturation_level(lqr.K, lqr.P, R_STD, ss.B, 50.0)
    cota = max_axis_angle(lqr.P, cs, 3)
    fr = max_recoverable_angle(p, lqr.K; umax=50.0, x_rail=1.5)
    sr = max_recoverable_angle(p, lqr.K; umax=50.0, x_rail=1e6)
    @printf("  %-8s %10.4f %12.4f %12.4f %12s\n", nombre, cota, fr, sr,
            abs(sr - fr) > 0.01 ? "RIEL" : "SATURACION")
end

println("\n  Subconjuntos de sensores de la Configuracion III:")
tabla = sensor_subset_analysis(ss3.A, ss3.C, ["pos", "th1", "th2", "th3"])
for f in tabla
    @printf("    %-24s %d/8  %s  %.4e\n", join(f.labels, "+"), f.rank,
            f.observable ? "SI" : "NO", f.margin_normalized)
end

println("\nFiguras del informe de la segunda entrega guardadas en Entrega_2/resumen_tecnico/figs/")
