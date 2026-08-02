#!/usr/bin/env julia
# =============================================================================
# setup.jl -- Instancia el entorno del proyecto
# =============================================================================
# Ejecutar: julia setup.jl
#
# Reproduce el entorno EXACTO fijado en Manifest.toml y lo precompila. No
# actualiza versiones: si lo que se quiere es actualizar, hay que hacerlo a
# mano con Pkg.update() y volver a validar con test/runtests.jl, porque los
# valores publicados en los informes se calcularon con estas versiones.
# =============================================================================

println("=" ^ 60)
println("  SETUP: Pendulo Invertido -- Proyecto de Algebra Lineal")
println("=" ^ 60)

using Pkg

# Activar el proyecto local
println("\n  Activando proyecto...")
Pkg.activate(@__DIR__)

# ---------------------------------------------------------------------------
# Instalar las dependencias EXACTAS que fija Manifest.toml.
#
# Se usa Pkg.instantiate() y NO Pkg.add() en bucle, y la diferencia importa:
# Pkg.add vuelve a resolver el grafo de dependencias y puede reescribir
# Manifest.toml con versiones mas nuevas, que es justo lo contrario de lo que un
# Manifest versionado pretende garantizar. Pkg.instantiate lee el Manifest y
# descarga esas versiones y no otras, de modo que el entorno queda identico al
# que produjo los resultados de los informes.
#
# Las dependencias directas se declaran en Project.toml: esa es la unica fuente
# de verdad. Listarlas tambien aqui solo crearia una copia que se desincroniza.
# ---------------------------------------------------------------------------
println("\n  Instalando las dependencias fijadas en Manifest.toml...")
println("  (la primera vez descarga bastantes paquetes y puede tardar)")
Pkg.instantiate()

println("\n  Dependencias directas del entorno:")
Pkg.status()

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
