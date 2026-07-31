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
    "CairoMakie",               # Graficas estaticas (informes y notebooks Pluto)
    "GLMakie",                  # Visualizacion y animacion interactiva
    "ForwardDiff",              # Diferenciacion automatica (pruebas de jacobianos)
    "MatrixEquations",          # Ecuaciones matriciales (Riccati, Lyapunov)
    "ControlSystems",           # Analisis de sistemas de control (apoyo y verificacion)
    "Symbolics",                # Calculo simbolico (derivaciones)
    "Pluto",                    # Servidor de notebooks interactivos
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
    main_triple.jl       Pipeline del pendulo triple (Configuracion III)
    setup.jl             Este archivo
    Project.toml         Dependencias
    src/
      model_simple.jl       Pendulo simple: parametros y EOM
      model_double.jl       Pendulo doble: parametros y EOM
      model_nlink.jl        Modelo GENERICO de N eslabones (contiene a I, II y III)
      model_triple.jl       Pendulo triple: capa delgada sobre ModelNLink
      linearization.jl      Linealizacion (simple, doble, generica en N), Kalman
      controller.jl         LQR, asignacion de polos (Ackermann), Riccati
      metrics.jl            Margen PBH, condicionamiento, elipsoide, ZOH, Schur
      sweep.jl              Barridos y fronteras de operabilidad
      observer.jl           Observador de Luenberger por dualidad, separacion
      animation_simple.jl   Animacion del pendulo simple (GLMakie)
      animation_double.jl   Animacion del pendulo doble (GLMakie)
      animation_triple.jl   Animacion de N eslabones (GLMakie)
    test/
      runtests.jl        Pruebas de regresion (618, en ocho conjuntos)
    notebooks/           Pluto notebooks
    figures/             Graficas y animaciones generadas
    docs/
      guia_maestra.md    Documento maestro (cubre las dos entregas)
      Entrega_1/         Informe tecnico, resumen ejecutivo y presentacion
      Entrega_2/         Informe tecnico, resumen ejecutivo y atlas de operabilidad

Proximo paso:
  julia main_simple.jl     (o main_double.jl, o main_triple.jl)

Para correr las pruebas:
  julia --project=. -t auto test/runtests.jl

O para uso interactivo:
  julia --project=.
  julia> include("main_simple.jl")
""")
