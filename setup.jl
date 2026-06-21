#!/usr/bin/env julia
# =============================================================================
# setup.jl -- Instala todas las dependencias del proyecto
# =============================================================================
# Ejecutar: julia setup.jl
# =============================================================================

println("=" ^ 60)
println("  SETUP: Pendulo Invertido -- Proyecto de Algebra Lineal")
println("=" ^ 60)

using Pkg

# Activar el proyecto local
println("\n  Activando proyecto...")
Pkg.activate(@__DIR__)

# Instalar dependencias
println("\n  Instalando paquetes...")
deps = [
    "DifferentialEquations",   # Solvers de EDOs
    "CairoMakie",               # Graficas inline (notebooks Pluto)
    "GLMakie",                  # Visualizacion y animacion interactiva
    "ControlSystems",           # Analisis de sistemas de control (apoyo y verificacion)
    "MatrixEquations",          # Ecuaciones matriciales (Riccati, Lyapunov)
    "Symbolics",                # Calculo simbolico (derivaciones)
    "PlutoUI",                  # Widgets interactivos para Pluto
]

for dep in deps
    println("  -> Instalando $dep...")
    try
        Pkg.add(dep)
        println("     $dep instalado")
    catch e
        println("     Error con $dep: $(e.msg)")
    end
end

# Precompilar todo
println("\n  Precompilando paquetes (esto puede tomar unos minutos)...")
Pkg.precompile()

println("\n" * "=" ^ 60)
println("  Setup completado")
println("=" ^ 60)
println("""

Estructura del proyecto:
  Proyecto/
    main_simple.jl       Pipeline del pendulo simple (Configuracion I)
    main_double.jl       Pipeline del pendulo doble (Configuracion II)
    setup.jl             Este archivo
    Project.toml         Dependencias
    src/
      model_simple.jl       Pendulo simple: parametros y EOM
      model_double.jl       Pendulo doble: parametros y EOM
      linearization.jl      Linealizacion (simple y doble), eigenvalores, Kalman
      controller.jl         LQR, asignacion de polos (Ackermann), Riccati
      animation_simple.jl   Animacion del pendulo simple (GLMakie)
      animation_double.jl   Animacion del pendulo doble (GLMakie)
    notebooks/           Pluto notebooks
    figures/             Graficas y animaciones generadas

Proximo paso:
  julia main_simple.jl     (o julia main_double.jl)

O para uso interactivo:
  julia --project=.
  julia> include("main_simple.jl")
""")
