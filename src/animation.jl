# =============================================================================
# animation.jl -- Animacion del pendulo invertido con Makie
# =============================================================================
# Genera animaciones del carro-pendulo usando GLMakie.
# Soporta:
#   - Animacion de la respuesta libre (caida sin control)
#   - Animacion de la respuesta controlada (estabilizacion LQR)
#   - Exportacion a MP4/GIF
# =============================================================================

module Animation

using Printf
import GLMakie

export animate_pendulum, save_animation

"""
    pendulum_coords(x_cart, theta, L_draw)

Calcula las coordenadas del extremo del pendulo a partir
de la posicion del carro y el angulo.
"""
function pendulum_coords(x_cart, theta, L_draw)
    x_tip = x_cart + L_draw * sin(theta)
    y_tip = L_draw * cos(theta)  # y positivo hacia arriba
    return (x_tip, y_tip)
end

"""
    animate_pendulum(sol, params; kwargs...) -> NamedTuple

Anima la solucion de una simulacion del pendulo.

Argumentos:
    sol     -- Solucion de DifferentialEquations.jl
    params  -- SystemParams del modelo

Keyword arguments:
    L_draw   -- Longitud visual del pendulo (default: 2*params.L)
    fps      -- Cuadros por segundo (default: 30)
    title    -- Titulo de la animacion
    trail    -- Mostrar estela del pendulo (default: true)
"""
function animate_pendulum(sol, params;
                          L_draw=nothing,
                          fps=30,
                          title="Pendulo invertido",
                          trail=true,
                          fig=nothing,
                          ax=nothing)
    if isnothing(L_draw)
        L_draw = 2 * params.L
    end

    # Extraer trayectoria
    t_span = sol.t
    dt = 1.0 / fps
    t_anim = range(t_span[1], t_span[end], step=dt)

    xs = [sol(t)[1] for t in t_anim]
    thetas = [sol(t)[3] for t in t_anim]

    # Rango de visualizacion
    x_range = maximum(abs.(xs)) + L_draw + 0.5
    y_range = L_draw + 0.5

    # Crear figura si no se proporciona
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

    # Observables para animacion reactiva
    cart_x = GLMakie.Observable(xs[1])
    bob_x = GLMakie.Observable(xs[1] + L_draw * sin(thetas[1]))
    bob_y = GLMakie.Observable(L_draw * cos(thetas[1]))
    time_text = GLMakie.Observable("t = 0.00 s")

    # Dibujar suelo
    GLMakie.hlines!(ax, [0.0], color=:gray, linewidth=1, linestyle=:dash)

    # Carro (rectangulo)
    cart_w, cart_h = 0.3, 0.15
    cart_rect = GLMakie.@lift GLMakie.Rect($cart_x - cart_w/2, -cart_h/2, cart_w, cart_h)
    GLMakie.poly!(ax, cart_rect, color=(:steelblue, 0.8), strokecolor=:black, strokewidth=1)

    # Ruedas
    for offset in [-0.1, 0.1]
        wheel_center = GLMakie.@lift GLMakie.Point2f($cart_x + offset, -cart_h/2)
        GLMakie.scatter!(ax, wheel_center, markersize=12, color=:gray30)
    end

    # Barra del pendulo
    rod_xs = GLMakie.@lift [$cart_x, $bob_x]
    rod_ys = GLMakie.@lift [0.0, $bob_y]
    GLMakie.lines!(ax, rod_xs, rod_ys, color=:gray20, linewidth=3)

    # Masa del pendulo
    bob_point = GLMakie.@lift GLMakie.Point2f($bob_x, $bob_y)
    GLMakie.scatter!(ax, bob_point, markersize=20, color=:red, strokecolor=:black, strokewidth=1)

    # Estela del pendulo
    if trail
        trail_xs = GLMakie.Observable(Float64[])
        trail_ys = GLMakie.Observable(Float64[])
        GLMakie.lines!(ax, trail_xs, trail_ys, color=(:red, 0.2), linewidth=1)
    end

    # Texto de tiempo
    GLMakie.text!(ax, time_text, position=(x_range - 0.3, y_range - 0.2),
                  align=(:right, :top), fontsize=14)

    # Retornar datos para la animacion
    return (fig=fig, ax=ax,
            frames=length(t_anim),
            update=function(i)
                cart_x[] = xs[i]
                bx, by = pendulum_coords(xs[i], thetas[i], L_draw)
                bob_x[] = bx
                bob_y[] = by
                time_text[] = @sprintf("t = %.2f s", t_anim[i])
                if trail
                    push!(trail_xs.val, bx)
                    push!(trail_ys.val, by)
                    GLMakie.notify(trail_xs)
                    GLMakie.notify(trail_ys)
                end
            end)
end

"""
    save_animation(anim_data, filename; fps=30)

Guarda la animacion como MP4 o GIF.
"""
function save_animation(anim_data, filename; fps=30)
    GLMakie.record(anim_data.fig, filename, 1:anim_data.frames; framerate=fps) do i
        anim_data.update(i)
    end

    println("Animacion guardada en: $filename")
end

end # module Animation
