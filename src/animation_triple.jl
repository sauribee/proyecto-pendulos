# =============================================================================
# animation_triple.jl -- Animacion del pendulo invertido TRIPLE con Makie
# =============================================================================
# Dibuja el carro y los tres eslabones en serie usando GLMakie. Devuelve la
# misma estructura que Animation.animate_pendulum (fig, ax, frames, update), de
# modo que se puede guardar con Animation.save_animation.
#
# Generaliza link_positions de animation_double.jl a un numero arbitrario de
# eslabones: la cinematica de una cadena articulada es la misma suma acumulada
# para cualquier N, asi que no hay razon para fijarla en tres.
#
# Ajustes de dibujo respecto al doble: con tres eslabones de 1/3 m las
# articulaciones quedan mucho mas juntas que con dos de 1/2 m, asi que cada
# eslabon se dibuja de un color distinto y las masas intermedias se marcan mas
# pequenas que el extremo, para que la cadena se lea sin ambiguedad.
# =============================================================================

module AnimationTriple

using Printf
import GLMakie

export animate_pendulum_triple, link_positions_nlink

"""
    link_positions_nlink(x_cart, thetas, l) -> Vector{Tuple{Float64,Float64}}

Coordenadas de todos los puntos de la cadena: el pivote sobre el carro y el
extremo superior de cada eslabon. Devuelve N+1 puntos.

Como los angulos son ABSOLUTOS respecto a la vertical superior (no relativos al
eslabon anterior), la posicion del extremo del eslabon j es simplemente la suma
acumulada

    x_j = x_carro + sum_{i<=j} l_i sin(theta_i)
    y_j =           sum_{i<=j} l_i cos(theta_i)
"""
function link_positions_nlink(x_cart, thetas, l)
    N = length(thetas)
    pts = Vector{Tuple{Float64,Float64}}(undef, N + 1)
    pts[1] = (x_cart, 0.0)
    for j in 1:N
        pts[j+1] = (pts[j][1] + l[j] * sin(thetas[j]),
                    pts[j][2] + l[j] * cos(thetas[j]))
    end
    return pts
end

"""
    animate_pendulum_triple(sol, params; kwargs...) -> NamedTuple

Anima la solucion de una simulacion del pendulo triple (o de cualquier numero
de eslabones: N se deduce de `params.l`).

Argumentos:
    sol     -- Solucion de DifferentialEquations.jl (estado de dimension 2(N+1))
    params  -- SystemParamsTriple / SystemParamsNLink del modelo

Keyword arguments:
    fps      -- Cuadros por segundo (default: 30)
    title    -- Titulo de la animacion
    trail    -- Mostrar estela del extremo (default: true)
"""
function animate_pendulum_triple(sol, params;
                                 fps=30,
                                 title="Pendulo invertido triple",
                                 trail=true,
                                 fig=nothing,
                                 ax=nothing)
    l = params.l
    N = length(l)

    # Muestreo temporal uniforme
    t_span = sol.t
    dt = 1.0 / fps
    t_anim = range(t_span[1], t_span[end], step=dt)

    xs = [sol(t)[1] for t in t_anim]
    thetas = [[sol(t)[2j+1] for j in 1:N] for t in t_anim]

    # Rango de visualizacion
    reach = sum(l)
    x_range = maximum(abs.(xs)) + reach + 0.5
    y_range = reach + 0.5

    if isnothing(fig)
        fig = GLMakie.Figure(size=(800, 600))
    end
    if isnothing(ax)
        ax = GLMakie.Axis(fig[1, 1],
                          title=title,
                          xlabel="x [m]",
                          ylabel="y [m]",
                          aspect=GLMakie.DataAspect(),
                          limits=(-x_range, x_range, -0.5, y_range))
    end

    # Estado inicial para los observables
    pts0 = link_positions_nlink(xs[1], thetas[1], l)

    cart_x = GLMakie.Observable(xs[1])
    rod_xs = [GLMakie.Observable([pts0[j][1], pts0[j+1][1]]) for j in 1:N]
    rod_ys = [GLMakie.Observable([pts0[j][2], pts0[j+1][2]]) for j in 1:N]
    joint_pts = [GLMakie.Observable(GLMakie.Point2f(pts0[j+1][1], pts0[j+1][2]))
                 for j in 1:N]
    time_text = GLMakie.Observable("t = 0.00 s")

    # Suelo
    GLMakie.hlines!(ax, [0.0], color=:gray, linewidth=1, linestyle=:dash)

    # Carro
    cart_w, cart_h = 0.3, 0.15
    cart_rect = GLMakie.@lift GLMakie.Rect($cart_x - cart_w/2, -cart_h/2, cart_w, cart_h)
    GLMakie.poly!(ax, cart_rect, color=(:steelblue, 0.8), strokecolor=:black, strokewidth=1)

    # Ruedas
    for offset in [-0.1, 0.1]
        wheel_center = GLMakie.@lift GLMakie.Point2f($cart_x + offset, -cart_h/2)
        GLMakie.scatter!(ax, wheel_center, markersize=12, color=:gray30)
    end

    # Eslabones: un tono distinto por eslabon, del mas oscuro (abajo) al mas
    # claro (arriba), para distinguirlos cuando quedan casi alineados.
    rod_colors = [GLMakie.RGBf(0.15 + 0.20 * (j - 1) / max(N - 1, 1),
                               0.15 + 0.20 * (j - 1) / max(N - 1, 1),
                               0.15 + 0.20 * (j - 1) / max(N - 1, 1)) for j in 1:N]
    for j in 1:N
        GLMakie.lines!(ax, rod_xs[j], rod_ys[j], color=rod_colors[j], linewidth=3)
    end

    # Masas: las articulaciones intermedias mas pequenas que el extremo
    for j in 1:N
        es_extremo = j == N
        GLMakie.scatter!(ax, joint_pts[j],
                         markersize=es_extremo ? 20 : 14,
                         color=es_extremo ? :red : :orange,
                         strokecolor=:black, strokewidth=1)
    end

    # Estela del extremo
    if trail
        trail_xs = GLMakie.Observable(Float64[])
        trail_ys = GLMakie.Observable(Float64[])
        GLMakie.lines!(ax, trail_xs, trail_ys, color=(:red, 0.2), linewidth=1)
    end

    # Texto de tiempo
    GLMakie.text!(ax, time_text, position=(x_range - 0.3, y_range - 0.2),
                  align=(:right, :top), fontsize=14)

    return (fig=fig, ax=ax,
            frames=length(t_anim),
            update=function(i)
                pts = link_positions_nlink(xs[i], thetas[i], l)
                cart_x[] = xs[i]
                for j in 1:N
                    rod_xs[j][] = [pts[j][1], pts[j+1][1]]
                    rod_ys[j][] = [pts[j][2], pts[j+1][2]]
                    joint_pts[j][] = GLMakie.Point2f(pts[j+1][1], pts[j+1][2])
                end
                time_text[] = @sprintf("t = %.2f s", t_anim[i])
                if trail
                    push!(trail_xs.val, pts[N+1][1])
                    push!(trail_ys.val, pts[N+1][2])
                    GLMakie.notify(trail_xs)
                    GLMakie.notify(trail_ys)
                end
            end)
end

end # module AnimationTriple
