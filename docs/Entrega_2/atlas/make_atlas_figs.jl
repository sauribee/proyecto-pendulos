# =============================================================================
# make_atlas_figs.jl -- Atlas de operabilidad de la Configuracion III
# =============================================================================
# Genera los cuatro mapas del atlas y las tablas 1D que van al informe.
#
#   Mapa A  (l3/l1, m3/m1)   -> margen PBH normalizado, escala log
#   Mapa B  (umax, theta0)   -> exito/fracaso, con la frontera numerica y la
#                               cota elipsoidal superpuesta   [EL PRINCIPAL]
#   Mapa C  (M, masa total)  -> margen PBH y norm(K), dos paneles
#   Mapa D  (R, theta0)      -> exito/fracaso, frontera theta_max(R)
#
# Los mapas A y C son algebra lineal pura (rapidos). Los mapas B y D exigen
# integrar el modelo NO LINEAL con saturacion en cada punto de la malla, asi
# que se paralelizan con hilos y usan terminacion temprana. Los resultados se
# cachean en CSV para no recomputar al regenerar las figuras.
#
# Ejecucion (los hilos importan mucho aqui):
#     julia --project=. -t auto docs/Entrega_2/atlas/make_atlas_figs.jl
#
# Para forzar el recomputo, borrar docs/Entrega_2/atlas/data/ o exportar
# ATLAS_FORCE=1.
# =============================================================================

using Pkg
const PROJ_ROOT = dirname(dirname(dirname(@__DIR__)))
Pkg.activate(PROJ_ROOT)

include(joinpath(PROJ_ROOT, "src", "model_nlink.jl"))
include(joinpath(PROJ_ROOT, "src", "model_triple.jl"))
include(joinpath(PROJ_ROOT, "src", "linearization.jl"))
include(joinpath(PROJ_ROOT, "src", "controller.jl"))
include(joinpath(PROJ_ROOT, "src", "metrics.jl"))
include(joinpath(PROJ_ROOT, "src", "sweep.jl"))

using .ModelNLink, .ModelTriple, .Linearization, .Controller, .Metrics, .Sweep
using LinearAlgebra
using Printf
using CairoMakie

# ---------------------------------------------------------------------------
# Paleta del proyecto (la misma de docs/Entrega_1/presentacion/make_slide_figs.jl)
# ---------------------------------------------------------------------------
const MARINO    = RGBf(27/255, 64/255, 121/255)    # 1B4079
const PETROLEO  = RGBf(77/255, 124/255, 138/255)   # 4D7C8A
const GRISVERDE = RGBf(127/255, 156/255, 150/255)  # 7F9C96
const SALVIA    = RGBf(143/255, 173/255, 136/255)  # 8FAD88
const LIMON     = RGBf(199/255, 219/255, 148/255)  # C7DB94

const PALETA = cgrad([MARINO, PETROLEO, GRISVERDE, SALVIA, LIMON])

const FIG_DIR = joinpath(@__DIR__, "figs")
const DATA_DIR = joinpath(@__DIR__, "data")
isdir(FIG_DIR) || mkpath(FIG_DIR)
isdir(DATA_DIR) || mkpath(DATA_DIR)

const FORZAR = get(ENV, "ATLAS_FORCE", "0") == "1"

# ---------------------------------------------------------------------------
# Cache en CSV. Se evita una dependencia nueva escribiendo y leyendo a mano.
# ---------------------------------------------------------------------------

"""
Guarda una matriz con sus ejes en un CSV legible: la primera fila lleva los
valores del eje y (columnas) y la primera columna los del eje x (filas).
"""
function save_grid(path, xs, ys, Z)
    open(path, "w") do io
        println(io, "x\\y," * join(string.(ys), ","))
        for (i, x) in enumerate(xs)
            println(io, string(x) * "," * join(string.(Z[i, :]), ","))
        end
    end
end

function load_grid(path)
    lineas = readlines(path)
    ys = parse.(Float64, split(lineas[1], ",")[2:end])
    xs = Float64[]
    filas = Vector{Vector{Float64}}()
    for l in lineas[2:end]
        campos = split(l, ",")
        push!(xs, parse(Float64, campos[1]))
        push!(filas, parse.(Float64, campos[2:end]))
    end
    return xs, ys, reduce(vcat, [reshape(f, 1, :) for f in filas])
end

"""
Calcula `f()` o lo recupera del cache si el CSV ya existe.
"""
function cached(nombre, f)
    path = joinpath(DATA_DIR, nombre * ".csv")
    if isfile(path) && !FORZAR
        println("  cache: $nombre.csv")
        return load_grid(path)
    end
    print("  calculando $nombre ... ")
    t0 = time()
    xs, ys, Z = f()
    save_grid(path, xs, ys, Z)
    @printf("%.1f s\n", time() - t0)
    return xs, ys, Z
end

# ---------------------------------------------------------------------------
# Configuracion base
# ---------------------------------------------------------------------------
const LIN = linearize_system_nlink
const R_STD = reshape([0.1], 1, 1)
const P_BASE = default_params_triple()

@printf("Hilos: %d\n\n", Threads.nthreads())

# ===========================================================================
# TABLAS 1D PARA EL INFORME
# ===========================================================================

println("=" ^ 78)
println("  TABLAS 1D")
println("=" ^ 78)

println("\n--- Escalamiento con N (longitud total 1 m, masa total 0.6 kg) ---")
@printf("  %2s %4s %11s %9s %12s %12s %11s\n",
        "N", "n", "lambda_max", "tau[ms]", "PBH norm", "cond(Ctrb)", "norm(K)")
for N in 1:5
    p = uniform_rods(1.0, fill(0.6/N, N), fill(1.0/N, N); b=0.1)
    v = viability(LIN, p)
    @printf("  %2d %4d %11.2f %9.0f %12.3e %12.3e %11.1f\n",
            N, 2(N+1), v.lambda_max, 1000v.tau_fall, v.pbh_normalized,
            v.cond_ctrb, v.norm_K)
end

println("\n--- Masa del carro (triple) ---")
res_M = sweep_1d(LIN, P_BASE, :M, [0.05, 0.10, 0.25, 1.00, 5.00, 20.00])
print_sweep_summary(res_M; spec="M [kg]", label="masa del carro")

println("\n--- Longitud del eslabon superior (triple) ---")
res_l3 = sweep_1d(LIN, P_BASE, (:l, 3), [0.02, 0.10, 1/3, 1.00, 2.00])
print_sweep_summary(res_l3; spec="l3 [m]", label="longitud del eslabon 3")

println("\n--- Masa del eslabon superior (triple) ---")
res_m3 = sweep_1d(LIN, P_BASE, (:m, 3), [0.02, 0.20, 0.80, 3.20])
print_sweep_summary(res_m3; spec="m3 [kg]", label="masa del eslabon 3")

# ===========================================================================
# MAPA A: (l3/l1, m3/m1) -> margen PBH normalizado
# ===========================================================================

println("\n" * "=" ^ 78)
println("  MAPA A: region estructuralmente sana")
println("=" ^ 78)

ratios_l = collect(range(0.1, 3.0, length=60))
ratios_m = collect(range(0.1, 3.0, length=60))

xs_a, ys_a, Z_a = cached("atlas_a", function ()
    Z = zeros(length(ratios_l), length(ratios_m))
    Threads.@threads for i in eachindex(ratios_l)
        for j in eachindex(ratios_m)
            p = set_param(P_BASE, (:l, 3), ratios_l[i] * P_BASE.l[1])
            p = set_param(p, (:m, 3), ratios_m[j] * P_BASE.m[1])
            ss = LIN(p)
            Z[i, j] = pbh_controllability_margin_normalized(ss.A, ss.B)
        end
    end
    return ratios_l, ratios_m, Z
end)

fig_a = Figure(size=(800, 640))
ax_a = Axis(fig_a[1, 1],
            title="Mapa A -- margen PBH normalizado (escala log)",
            xlabel="l3 / l1", ylabel="m3 / m1")
hm_a = heatmap!(ax_a, xs_a, ys_a, log10.(Z_a), colormap=PALETA)
contour!(ax_a, xs_a, ys_a, log10.(Z_a), levels=8, color=(:black, 0.35),
         linewidth=0.8)
scatter!(ax_a, [1.0], [1.0], marker=:star5, markersize=22, color=:white,
         strokecolor=:black, strokewidth=1.5)
text!(ax_a, 1.05, 1.05, text="nominal", color=:black, fontsize=13)
Colorbar(fig_a[1, 2], hm_a, label="log10 del margen PBH normalizado")
save(joinpath(FIG_DIR, "atlas_a.png"), fig_a, px_per_unit=2)
println("  Guardada: atlas_a.png")

@printf("  rango del margen: %.3e (peor) a %.3e (mejor)\n",
        minimum(Z_a), maximum(Z_a))
i_best = argmax(Z_a)
@printf("  mejor punto: l3/l1 = %.2f, m3/m1 = %.2f\n",
        xs_a[i_best[1]], ys_a[i_best[2]])

# ===========================================================================
# MAPA B: (umax, theta0) -> exito/fracaso    [EL PRINCIPAL]
# ===========================================================================

println("\n" * "=" ^ 78)
println("  MAPA B: frontera de operabilidad practica")
println("=" ^ 78)

ss3 = LIN(P_BASE)
lqr3 = design_lqr(ss3.A, ss3.B, default_weights(8), R_STD)
K3, P3 = lqr3.K, lqr3.P

umaxs = collect(range(5.0, 200.0, length=50))
thetas = collect(range(0.01, 0.45, length=50))

xs_b, ys_b, Z_b = cached("atlas_b", function ()
    Z = zeros(length(umaxs), length(thetas))
    Threads.@threads for i in eachindex(umaxs)
        for j in eachindex(thetas)
            Z[i, j] = recovers(P_BASE, K3, thetas[j]; umax=umaxs[i]) ? 1.0 : 0.0
        end
    end
    return umaxs, thetas, Z
end)

# Cota elipsoidal: c* = umax^2 R^2 / (B'PB) y theta0 = sqrt(c*/P33), de modo
# que la cota es LINEAL en umax. Es una recta por el origen.
cota_theta = [max_axis_angle(P3, nonsaturation_level(K3, P3, R_STD, ss3.B, u), 3)
              for u in xs_b]

# Frontera numerica: mayor theta0 con exito para cada umax
frontera = [(idx = findlast(Z_b[i, :] .> 0.5);
             idx === nothing ? NaN : ys_b[idx]) for i in eachindex(xs_b)]

# Punto donde la cota deja de ser valida: la cota crece LINEALMENTE con umax
# (porque c* es proporcional a umax^2), pero la frontera real se estanca en
# cuanto la saturacion deja de ser la restriccion activa. A partir de ahi la
# cota promete mas de lo que el sistema no lineal cumple.
i_cruce = findfirst(i -> cota_theta[i] > frontera[i], eachindex(xs_b))
u_cruce = i_cruce === nothing ? NaN : xs_b[i_cruce]

fig_b = Figure(size=(950, 680))
ax_b = Axis(fig_b[1, 1],
            title="Mapa B -- frontera de operabilidad del pendulo triple",
            xlabel="Limite del actuador u_max [N]",
            ylabel="Angulo inicial theta1 [rad]")
heatmap!(ax_b, xs_b, ys_b, Z_b, colormap=[RGBf(0.93, 0.93, 0.93), LIMON])
lines!(ax_b, xs_b, frontera, color=MARINO, linewidth=3,
       label="frontera numerica (modelo no lineal)")

if i_cruce === nothing
    lines!(ax_b, xs_b, cota_theta, color=PETROLEO, linewidth=2.5,
           linestyle=:dash, label="cota elipsoidal  theta = sqrt(c*/P33)")
else
    # Tramo valido: la cota queda por dentro de la frontera
    lines!(ax_b, xs_b[1:i_cruce-1], cota_theta[1:i_cruce-1], color=PETROLEO,
           linewidth=2.5, linestyle=:dash,
           label="cota elipsoidal (valida: garantiza no saturar)")
    # Tramo invalido: la cota supera la frontera
    lines!(ax_b, xs_b[i_cruce-1:end], cota_theta[i_cruce-1:end],
           color=(GRISVERDE, 0.75), linewidth=2, linestyle=:dot,
           label="cota elipsoidal (ya no acota: manda la no linealidad)")
    vlines!(ax_b, [u_cruce], color=:gray40, linestyle=:dashdot, linewidth=1.5)
    text!(ax_b, u_cruce + 3, 0.40,
          text=@sprintf("u_max = %.0f N\nla cota deja de ser valida", u_cruce),
          fontsize=12, color=:gray20, align=(:left, :top))
end

ylims!(ax_b, 0.0, 0.45)
axislegend(ax_b, position=:rb)
save(joinpath(FIG_DIR, "atlas_b.png"), fig_b, px_per_unit=2)
println("  Guardada: atlas_b.png")

i_sat = findfirst(i -> i > 1 && frontera[i] <= frontera[i-1] + 1e-6,
                  eachindex(frontera))
if i_sat !== nothing
    @printf("  la frontera se estanca a partir de u_max = %.0f N (theta = %.4f):\n",
            xs_b[i_sat], frontera[i_sat])
    println("  mas actuador ya no compra mas angulo, porque la restriccion")
    println("  activa deja de ser el motor y pasa a ser la no linealidad")
end
if i_cruce !== nothing
    @printf("  la cota elipsoidal es valida hasta u_max = %.0f N y a partir de\n",
            u_cruce)
    println("  ahi supera la frontera real: c* solo certifica NO SATURACION mas")
    println("  convergencia del sistema LINEAL, no convergencia del no lineal")
else
    println("  la cota elipsoidal queda por dentro en todo el rango barrido")
end

# ===========================================================================
# MAPA C: (M, masa total del pendulo) -> margen PBH y norm(K)
# ===========================================================================

println("\n" * "=" ^ 78)
println("  MAPA C: el optimo interior de la masa del carro")
println("=" ^ 78)

Ms = collect(10 .^ range(log10(0.03), log10(20.0), length=60))
mtots = collect(range(0.1, 2.0, length=60))

xs_c, ys_c, Z_c = cached("atlas_c_pbh", function ()
    Z = zeros(length(Ms), length(mtots))
    Threads.@threads for i in eachindex(Ms)
        for j in eachindex(mtots)
            p = SystemParamsNLink(M=Ms[i], m=fill(mtots[j]/3, 3),
                                  l=copy(P_BASE.l), a=copy(P_BASE.a),
                                  Il=(1/12) .* fill(mtots[j]/3, 3) .* P_BASE.l .^ 2,
                                  g=P_BASE.g, b=P_BASE.b, c=Float64[])
            ss = LIN(p)
            Z[i, j] = pbh_controllability_margin_normalized(ss.A, ss.B)
        end
    end
    return Ms, mtots, Z
end)

_, _, Z_c_K = cached("atlas_c_normk", function ()
    Z = zeros(length(Ms), length(mtots))
    Threads.@threads for i in eachindex(Ms)
        for j in eachindex(mtots)
            p = SystemParamsNLink(M=Ms[i], m=fill(mtots[j]/3, 3),
                                  l=copy(P_BASE.l), a=copy(P_BASE.a),
                                  Il=(1/12) .* fill(mtots[j]/3, 3) .* P_BASE.l .^ 2,
                                  g=P_BASE.g, b=P_BASE.b, c=Float64[])
            ss = LIN(p)
            Z[i, j] = norm(design_lqr(ss.A, ss.B, default_weights(8), R_STD).K)
        end
    end
    return Ms, mtots, Z
end)

fig_c = Figure(size=(1250, 560))

ax_c1 = Axis(fig_c[1, 1], title="Mapa C1 -- margen PBH normalizado",
             xlabel="Masa del carro M [kg]", ylabel="Masa total del pendulo [kg]",
             xscale=log10)
hm_c1 = heatmap!(ax_c1, xs_c, ys_c, log10.(Z_c), colormap=PALETA)
contour!(ax_c1, xs_c, ys_c, log10.(Z_c), levels=10, color=(:black, 0.3),
         linewidth=0.8)
# Cresta: para cada masa total, la M que maximiza el margen
cresta = [xs_c[argmax(Z_c[:, j])] for j in eachindex(ys_c)]
lines!(ax_c1, cresta, ys_c, color=:white, linewidth=3)
lines!(ax_c1, cresta, ys_c, color=:black, linewidth=1.5, linestyle=:dash)
Colorbar(fig_c[1, 2], hm_c1, label="log10 PBH normalizado")

ax_c2 = Axis(fig_c[1, 3], title="Mapa C2 -- esfuerzo de control norm(K)",
             xlabel="Masa del carro M [kg]", ylabel="Masa total del pendulo [kg]",
             xscale=log10)
# Paleta invertida para que en AMBOS paneles el color claro signifique "mejor":
# margen PBH alto en C1, esfuerzo de control bajo en C2.
hm_c2 = heatmap!(ax_c2, xs_c, ys_c, log10.(Z_c_K), colormap=reverse(PALETA))
contour!(ax_c2, xs_c, ys_c, log10.(Z_c_K), levels=10, color=(:black, 0.3),
         linewidth=0.8)
# La misma cresta de C1 superpuesta: el optimo del margen PBH NO coincide con
# el del esfuerzo de control, que es monotono y no tiene optimo interior.
lines!(ax_c2, cresta, ys_c, color=:white, linewidth=3)
lines!(ax_c2, cresta, ys_c, color=:black, linewidth=1.5, linestyle=:dash)
Colorbar(fig_c[1, 4], hm_c2, label="log10 norm(K)")

save(joinpath(FIG_DIR, "atlas_c.png"), fig_c, px_per_unit=2)
println("  Guardada: atlas_c.png")

j_nom = argmin(abs.(ys_c .- 0.6))
@printf("  con masa total 0.6 kg, el margen PBH se maximiza en M = %.3f kg\n",
        cresta[j_nom])
println("  el optimo es INTERIOR: un carro muy pesado se mueve poco ante la")
println("  misma fuerza y uno muy liviano pierde acoplamiento con los eslabones")

# ===========================================================================
# MAPA D: (R, theta0) -> exito/fracaso
# ===========================================================================

println("\n" * "=" ^ 78)
println("  MAPA D: compromiso agresividad-robustez")
println("=" ^ 78)

Rs = collect(10 .^ range(-3, 1, length=40))
thetas_d = collect(range(0.01, 0.32, length=45))

# Se barre R para DOS limites de actuador. La razon es el hallazgo del mapa B:
# por encima de unos 53 N la saturacion deja de ser la restriccion activa, y
# entonces el peso del control no puede influir en el angulo recuperable. El
# compromiso agresividad-robustez solo EXISTE cuando el actuador es el cuello
# de botella; con un panel solo, a 150 N, la grafica saldria plana y el
# resultado se leeria como "R no importa", que seria una conclusion falsa.
const UMAX_D = [25.0, 150.0]

resultados_d = map(UMAX_D) do umax
    xs, ys, Z = cached(@sprintf("atlas_d_%d", Int(umax)), function ()
        Zl = zeros(length(Rs), length(thetas_d))
        Threads.@threads for i in eachindex(Rs)
            Ri = reshape([Rs[i]], 1, 1)
            Ki = design_lqr(ss3.A, ss3.B, default_weights(8), Ri).K
            for j in eachindex(thetas_d)
                Zl[i, j] = recovers(P_BASE, Ki, thetas_d[j]; umax=umax) ? 1.0 : 0.0
            end
        end
        return Rs, thetas_d, Zl
    end)
    frontera = [(idx = findlast(Z[i, :] .> 0.5);
                 idx === nothing ? NaN : ys[idx]) for i in eachindex(xs)]
    return (umax=umax, xs=xs, ys=ys, Z=Z, frontera=frontera)
end

fig_d = Figure(size=(1250, 560))
for (k, r) in enumerate(resultados_d)
    rango = extrema(filter(!isnan, r.frontera))
    ax = Axis(fig_d[1, k],
              title=@sprintf("Mapa D%d -- u_max = %d N   (theta_max de %.3f a %.3f)",
                             k, Int(r.umax), rango[1], rango[2]),
              xlabel="Peso del control R",
              ylabel="Angulo inicial theta1 [rad]", xscale=log10)
    heatmap!(ax, r.xs, r.ys, r.Z, colormap=[RGBf(0.93, 0.93, 0.93), SALVIA])
    lines!(ax, r.xs, r.frontera, color=MARINO, linewidth=3, label="theta_max(R)")
    vlines!(ax, [0.1], color=PETROLEO, linestyle=:dash, linewidth=2,
            label="R nominal = 0.1")
    ylims!(ax, 0.0, 0.32)
    axislegend(ax, position=:lb)
end
save(joinpath(FIG_DIR, "atlas_d.png"), fig_d, px_per_unit=2)
println("  Guardada: atlas_d.png")

function reportar_d(r)
    f = replace(r.frontera, NaN => -Inf)
    ib = argmax(f)
    rango = extrema(filter(!isnan, r.frontera))
    th_nom = r.frontera[argmin(abs.(r.xs .- 0.1))]
    @printf("  u_max = %3d N: theta_max entre %.4f y %.4f (variacion %.1f%%)\n",
            Int(r.umax), rango[1], rango[2],
            100 * (rango[2] - rango[1]) / rango[2])
    @printf("                 optimo en R = %.3f (%.4f); con R nominal 0.1: %.4f",
            r.xs[ib], r.frontera[ib], th_nom)
    @printf("  (a %.1f%% del optimo)\n", 100 * (1 - th_nom / r.frontera[ib]))
end

for r in resultados_d
    reportar_d(r)
end

# El esfuerzo de control si depende fuertemente de R, aunque el angulo no.
normas_K = [norm(design_lqr(ss3.A, ss3.B, default_weights(8),
                            reshape([R], 1, 1)).K) for R in (1e-3, 0.1, 10.0)]
@printf("\n  norm(K) para R = 1e-3, 0.1, 10: %.0f, %.0f, %.0f (solo un factor %.1fx)\n",
        normas_K[1], normas_K[2], normas_K[3], normas_K[1] / normas_K[3])
println("  Lectura: cuatro ordenes de magnitud en R mueven norm(K) apenas un")
println("  factor 2 y el angulo recuperable un 11-13%. El motivo es que la")
println("  planta es MUY inestable: cualquier K que estabilice tiene que")
println("  cancelar tres modos con parte real hasta +18.06, y eso ya fija la")
println("  escala de la ganancia. El peso del control no puede abaratar lo que")
println("  la inestabilidad exige, asi que ajustar R aqui compra muy poco.")

# ===========================================================================
# CIERRE: fronteras de imposibilidad, separadas por tipo
# ===========================================================================

println("\n" * "=" ^ 78)
println("  FRONTERAS DE IMPOSIBILIDAD")
println("=" ^ 78)

configs = [("simple", uniform_rods(1.0, [0.3], [1.0]; b=0.1)),
           ("doble", point_masses(1.0, [0.3, 0.3], [0.5, 0.5]; b=0.0)),
           ("triple", P_BASE)]

println("\n  Imposibilidad ESTRUCTURAL (propiedad del sistema):")
for (nombre, p) in configs
    v = viability(LIN, p)
    @printf("    %-8s PBH normalizado = %.3e   cond(Ctrb) = %.2e   rank = %d\n",
            nombre, v.pbh_normalized, v.cond_ctrb, v.rank_ctrb)
end
println("    Ninguna es estructuralmente incontrolable: el rango es pleno en")
println("    las tres. Pero con cond(Ctrb) ~ 1e8 en el triple, el rango ya no")
println("    es un criterio confiable, y el margen PBH cae dos ordenes de I a III.")

println("\n  Imposibilidad PRACTICA (propiedad de la implementacion):")
@printf("    %-8s %10s %12s %12s %11s\n",
        "config", "cota c*", "frontera", "sin riel", "restriccion")
for (nombre, p) in configs
    ss = LIN(p)
    n = size(ss.A, 1)
    lqr = design_lqr(ss.A, ss.B, default_weights(n), R_STD)
    cs = nonsaturation_level(lqr.K, lqr.P, R_STD, ss.B, 50.0)
    cota = max_axis_angle(lqr.P, cs, 3)
    th_nom = max_recoverable_angle(p, lqr.K; umax=50.0, x_rail=1.5)
    th_sat = max_recoverable_angle(p, lqr.K; umax=50.0, x_rail=1e6)
    restriccion = abs(th_sat - th_nom) > 0.01 ? "RIEL" : "SATURACION"
    @printf("    %-8s %10.4f %12.4f %12.4f %11s\n",
            nombre, cota, th_nom, th_sat, restriccion)
end
println("    La cota elipsoidal solo acota SATURACION mas convergencia lineal:")
println("    hay que compararla contra la columna sin riel. En el simple la")
println("    restriccion activa es el RIEL, no el motor ni la no linealidad.")

println("\n  Periodo de muestreo maximo:")
for (nombre, p) in configs
    ss = LIN(p)
    n = size(ss.A, 1)
    K = design_lqr(ss.A, ss.B, default_weights(n), R_STD).K
    h = max_sampling_period(ss.A, ss.B, K)
    lam = dominant_mode(ss.A).lambda_max
    @printf("    %-8s h_max = %6.1f ms   1/lambda_max = %6.1f ms   " ,
            nombre, 1000h, 1000/lam)
    @printf("h_max*lambda_max = %.2f\n", h * lam)
end
println("    El margen relativo se estrecha al crecer la complejidad.")

println("\nAtlas completo. Figuras en docs/Entrega_2/atlas/figs/")
