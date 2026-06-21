# =============================================================================
# animation_double.jl -- Animacion del pendulo invertido DOBLE con Makie
# =============================================================================
# Dibuja el carro y los dos eslabones en serie usando GLMakie. Devuelve la
# misma estructura que Animation.animate_pendulum (fig, ax, frames, update),
# de modo que se puede guardar con Animation.save_animation.
# =============================================================================

module AnimationDouble

using Printf
import GLMakie

export animate_pendulum_double

"""
    link_positions(x_cart, theta1, theta2, L1, L2)

Devuelve las coordenadas del pivote (sobre el carro), la articulacion
intermedia y el extremo superior del pendulo doble.
"""
function link_positions(x_cart, theta1, theta2, L1, L2)
    pivot = (x_cart, 0.0)
    joint = (x_cart + L1 * sin(theta1), L1 * cos(theta1))
    tip = (joint[1] + L2 * sin(theta2), joint[2] + L2 * cos(theta2))
    return pivot, joint, tip
end

"""
    animate_pendulum_double(sol, params; kwargs...) -> NamedTuple

Anima la solucion de una simulacion del pendulo doble.

Argumentos:
    sol     -- Solucion de DifferentialEquations.jl (estado de dimension 6)
    params  -- SystemParamsDouble del modelo

Keyword arguments:
    fps      -- Cuadros por segundo (default: 30)
    title    -- Titulo de la animacion
    trail    -- Mostrar estela del extremo (default: true)
"""
function animate_pendulum_double(sol, params;
                                 fps=30,
                                 title="Pendulo invertido doble",
                                 trail=true,
                                 fig=nothing,
                                 ax=nothing)
    L1 = params.L1
    L2 = params.L2

    # Muestreo temporal uniforme
    t_span = sol.t
    dt = 1.0 / fps
    t_anim = range(t_span[1], t_span[end], step=dt)

    xs = [sol(t)[1] for t in t_anim]
    th1 = [sol(t)[3] for t in t_anim]
    th2 = [sol(t)[5] for t in t_anim]

    # Rango de visualizacion
    reach = L1 + L2
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
    p0, j0, e0 = link_positions(xs[1], th1[1], th2[1], L1, L2)

    cart_x = GLMakie.Observable(xs[1])
    joint_pt = GLMakie.Observable(GLMakie.Point2f(j0[1], j0[2]))
    tip_pt = GLMakie.Observable(GLMakie.Point2f(e0[1], e0[2]))
    rod1_xs = GLMakie.Observable([p0[1], j0[1]])
    rod1_ys = GLMakie.Observable([p0[2], j0[2]])
    rod2_xs = GLMakie.Observable([j0[1], e0[1]])
    rod2_ys = GLMakie.Observable([j0[2], e0[2]])
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

    # Eslabon inferior y superior
    GLMakie.lines!(ax, rod1_xs, rod1_ys, color=:gray20, linewidth=3)
    GLMakie.lines!(ax, rod2_xs, rod2_ys, color=:gray40, linewidth=3)

    # Masas: articulacion intermedia (m1) y extremo (m2)
    GLMakie.scatter!(ax, joint_pt, markersize=16, color=:orange,
                     strokecolor=:black, strokewidth=1)
    GLMakie.scatter!(ax, tip_pt, markersize=20, color=:red,
                     strokecolor=:black, strokewidth=1)

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
                p, j, e = link_positions(xs[i], th1[i], th2[i], L1, L2)
                cart_x[] = xs[i]
                rod1_xs[] = [p[1], j[1]]
                rod1_ys[] = [p[2], j[2]]
                rod2_xs[] = [j[1], e[1]]
                rod2_ys[] = [j[2], e[2]]
                joint_pt[] = GLMakie.Point2f(j[1], j[2])
                tip_pt[] = GLMakie.Point2f(e[1], e[2])
                time_text[] = @sprintf("t = %.2f s", t_anim[i])
                if trail
                    push!(trail_xs.val, e[1])
                    push!(trail_ys.val, e[2])
                    GLMakie.notify(trail_xs)
                    GLMakie.notify(trail_ys)
                end
            end)
end

end # module AnimationDouble
