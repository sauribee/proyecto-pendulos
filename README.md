# Péndulo invertido sobre un carro — Álgebra Lineal Aplicada

Análisis y control del péndulo invertido sobre un carro mediante herramientas de
álgebra lineal: representación en espacio de estados, análisis espectral,
controlabilidad y observabilidad de Kalman, asignación de polos (Ackermann) y
regulación lineal cuadrática (LQR). El proyecto estudia **dos configuraciones**
—un péndulo simple y un péndulo doble— que comparten la misma maquinaria genérica
de análisis y control.

Universidad Nacional de Colombia, Sede Medellín — Facultad de Ciencias.
Autores: Mateo Bedoya Rojas, Camilo Alejandro Patiño Osorio, Santiago Uribe Echavarría.

## Estructura del repositorio

```
proyecto-pendulos/
├── Project.toml                  Entorno Julia (dependencias)
├── Manifest.toml                 Versiones exactas resueltas (entorno reproducible)
├── setup.jl                      Instala e instancia las dependencias
├── main_simple.jl                Pipeline del pendulo SIMPLE (Configuracion I)
├── main_double.jl                Pipeline del pendulo DOBLE (Configuracion II)
├── README.md
├── src/
│   ├── model_simple.jl           Pendulo simple: parametros, EOM y lazo cerrado
│   ├── model_double.jl           Pendulo doble: parametros, EOM y lazo cerrado
│   ├── linearization.jl          Linealizacion (simple y doble), espectro, Kalman
│   ├── controller.jl             LQR (Riccati via Hamiltoniano) y Ackermann (genericos)
│   ├── animation_simple.jl       Animacion del pendulo simple (GLMakie)
│   └── animation_double.jl       Animacion del pendulo doble (GLMakie)
├── notebooks/                    Exploradores interactivos (Pluto)
│   ├── 01_exploracion_simple.jl
│   └── 02_exploracion_doble.jl                  
└── docs/
    ├── resumen_ejecutivo/      
    │   ├── resumen_ejecutivo.tex
    │   └── resumen_ejecutivo.pdf
    ├── resumen_tecnico/        
    │   ├── resumen_tecnico.tex
    │   ├── resumen_tecnico.pdf
    │   ├── make_report_figs.jl     Genera las figuras de respuesta del informe
    │   └── figs/                   Figuras que el PDF necesita para compilar
    └── presentacion/           
        ├── presentacion.tex        Diapositivas (Beamer, 16:9, 20 minutos)
        ├── presentacion.pdf
        ├── make_slide_figs.jl      Genera las figuras de las diapositivas
        └── figs/                   Figuras que el PDF necesita para compilar
```

Las dos configuraciones comparten el módulo `Controller` (LQR, Ackermann,
Riccati) y las funciones de análisis de `Linearization`, que son **genéricas en
la dimensión del estado**: el mismo código analiza el sistema de $\mathbb{R}^4$ y
el de $\mathbb{R}^6$ sin cambios.


## Requisitos

- [Julia](https://julialang.org/) (probado con la versión del `Manifest.toml`).
- Las dependencias se declaran en `Project.toml` y se fijan en `Manifest.toml`.

Paquetes principales: `DifferentialEquations` (solver `Tsit5`),
`CairoMakie` / `GLMakie` (gráficas y animación), `ControlSystems` y
`MatrixEquations` (verificación), `Symbolics` (derivaciones) y
`Pluto` / `PlutoUI` (notebooks interactivos).

## Uso rápido

```bash
# Instalar dependencias (una sola vez)
julia setup.jl

# Pendulo simple (Configuracion I, estado en R^4)
julia main_simple.jl

# Pendulo doble (Configuracion II, estado en R^6)
julia main_double.jl
```

O de forma interactiva, activando el entorno del proyecto:

```bash
julia --project=.
julia> include("main_simple.jl")   # o include("main_double.jl")
```

Cada pipeline ejecuta el flujo completo: define parámetros, simula la respuesta
libre (sin control), linealiza y analiza (eigenvalores, controlabilidad,
observabilidad), diseña el o los controladores, simula el lazo cerrado, genera
las gráficas comparativas en `figures/` y produce la animación.

## Exploración interactiva con Pluto (paso a paso)

Los notebooks de `notebooks/` permiten correr todo el análisis de forma
interactiva: se mueve un slider (una masa, una longitud, un peso del LQR) y
Pluto recalcula automáticamente las matrices, los eigenvalores, la ganancia $K$,
las gráficas y la animación.

### Paso 0. Requisito (una sola vez)

Instalar las dependencias del proyecto. Desde la carpeta del proyecto:

```bash
julia setup.jl
```

(Alternativa equivalente: `julia --project=. -e "using Pkg; Pkg.instantiate()"`.)

### Paso 1. Abrir Julia en la carpeta del proyecto

```bash
cd ruta/al/proyecto-pendulos
julia --project=.
```

El flag `--project=.` es importante: activa el entorno del proyecto (el de
`Project.toml`), de modo que el notebook encuentre todos los paquetes.

### Paso 2. Lanzar el servidor de Pluto

Dentro de la sesión de Julia (en el prompt `julia>`):

```julia
import Pluto
Pluto.run()
```

Esto abre Pluto en el navegador (normalmente en `http://localhost:1234`). La
primera vez puede tardar mientras precompila. Deja esa terminal abierta: es el
servidor; si la cierras, se cierra Pluto.

### Paso 3. Abrir un notebook

En la pantalla de inicio de Pluto, en el campo **"Open a notebook"**, escribe o
pega la ruta del notebook y pulsa **Open**:

- `notebooks/01_exploracion_simple.jl` — péndulo simple (estado en $\mathbb{R}^4$)
- `notebooks/02_exploracion_doble.jl` — péndulo doble (estado en $\mathbb{R}^6$)

La primera celda activa el proyecto y la segunda carga los módulos de `src/`. La
primera apertura precompila CairoMakie y DifferentialEquations (puede tardar
varios minutos); las siguientes son rápidas.

### Paso 4. Interactuar con los sliders

Al cambiar cualquier slider, todas las celdas dependientes se recalculan solas.

| Notebook | Sliders de parámetros | Sliders de control | Condición inicial |
|---|---|---|---|
| `01_exploracion_simple.jl` | `M`, `m`, `Lbar`, `g`, `b` | `Q11`, `Q33`, `R` | `theta0` |
| `02_exploracion_doble.jl` | `M`, `m1`, `m2`, `L1`, `L2`, `g` | `Q` (pos, `theta1`, `theta2`), `R` | `theta1_0`, `theta2_0` |

Verás reaccionar en vivo: el eigenvalor inestable, la ganancia $K$, los polos de
lazo cerrado, las gráficas de respuesta y la animación.

### Paso 5. Exportar la animación (opcional)

Al final de cada notebook hay una casilla (checkbox). Al marcarla se genera el
archivo de animación en `figures/`:

- simple: `04_comparacion_libre_vs_lqr.gif` (y `05_..._.mp4` si hay `ffmpeg`)
- doble: `08_doble_exploracion.gif`

Desmárcala para no regenerar el archivo en cada cambio de slider.

### Paso 6. Cerrar

Guarda el notebook (Pluto guarda solo el `.jl`), cierra la pestaña del navegador
y detén el servidor con `Ctrl-C` en la terminal de Julia.

## Resultados esperados

- **Simple:** espectro de lazo abierto $\{+4.21,\ 0,\ -0.077,\ -4.23\}$ (1 modo
  inestable); LQR $K=(-3.16,\,-4.69,\,-45.39,\,-10.93)$; Ackermann con polos
  $\{-1,-2,-3,-4\}$ da $K=(-1.75,\,-3.75,\,-39.01,\,-9.60)$.
- **Doble:** espectro de lazo abierto
  $\{+8.57,\ +4.09,\ 0,\ 0,\ -4.09,\ -8.57\}$ (2 modos inestables); rangos de
  controlabilidad y observabilidad $6/6$;
  LQR $K=(3.16,\,5.82,\,-191.55,\,-10.99,\,228.32,\,36.14)$.

Ambas configuraciones reproducen los valores del informe técnico
(`docs/resumen_tecnico/`).

## Regenerar las figuras del informe

Las figuras de respuesta temporal que aparecen en el informe técnico se generan
con un script aparte (usa CairoMakie, salida estática):

```bash
julia --project=. docs/resumen_tecnico/make_report_figs.jl
```

Esto reescribe `docs/resumen_tecnico/figs/` y reporta las métricas (tiempo de
asentamiento, esfuerzo de control pico) que se citan en la discusión.

Las figuras de las diapositivas se regeneran de forma análoga:

```bash
julia --project=. docs/presentacion/make_slide_figs.jl
```

## Nota sobre el entorno

`Manifest.toml` fija las versiones exactas para reproducibilidad. Tras clonar el
repositorio conviene ejecutar una vez:

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
```
