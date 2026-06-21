# Péndulo invertido sobre un carro — Álgebra Lineal Aplicada

Análisis y control del péndulo invertido sobre un carro mediante herramientas
de álgebra lineal: representación en espacio de estados, análisis espectral,
controlabilidad y observabilidad de Kalman, asignación de polos (Ackermann) y
regulación lineal cuadrática (LQR).

Universidad Nacional de Colombia, Sede Medellín — Facultad de Ciencias.
Autores: Mateo Bedoya Rojas, Santiago Uribe Echavarría, Camilo Alejandro Patiño Osorio.

## Estructura del proyecto

```
Proyecto/
├── Project.toml        Entorno Julia (dependencias)
├── Manifest.toml       Versiones exactas resueltas
├── main_simple.jl      Pipeline del pendulo SIMPLE (Configuracion I)
├── main_double.jl      Pipeline del pendulo DOBLE (Configuracion II)
├── setup.jl            Instala las dependencias
├── src/
│   ├── model.jl              Pendulo simple: parametros y EOM no lineales
│   ├── model_double.jl       Pendulo doble: parametros y EOM no lineales
│   ├── linearization.jl      Linealizacion (simple y doble), eigenvalores, Kalman
│   ├── controller.jl         LQR (Riccati via Hamiltoniano) y Ackermann (genericos)
│   ├── animation.jl          Animacion del pendulo simple (GLMakie)
│   └── animation_double.jl   Animacion del pendulo doble (GLMakie)
├── notebooks/          Notebooks de exploracion (Pluto)
├── figures/            Graficas y animaciones generadas
├── docs/               Documentacion del proyecto
│   ├── *.md                  Documentos de investigacion (teoria)
│   ├── resumen_tecnico/      Informe tecnico (LaTeX + PDF) + make_report_figs.jl
│   └── resumen_ejecutivo/    Resumen ejecutivo (LaTeX + PDF, maximo 5 paginas)
└── _archive/           Versiones anteriores (draft_v1, draft_v2, draft_v3)
```

Las dos configuraciones comparten el modulo `Controller` (LQR, Ackermann,
Riccati) y las funciones de analisis de `Linearization`, que son genericas en
la dimension del estado.

## Uso

```bash
# Instalar dependencias (una sola vez)
julia setup.jl

# Pendulo simple (Configuracion I, estado en R^4)
julia main_simple.jl

# Pendulo doble (Configuracion II, estado en R^6)
julia main_double.jl
```

O de forma interactiva:

```bash
julia --project=.
julia> include("main_simple.jl")   # o include("main_double.jl")
```

## Ejecución interactiva con Pluto (paso a paso)

Los notebooks de `notebooks/` permiten correr todo el análisis de forma
interactiva: se mueve un slider (una masa, una longitud, un peso del LQR) y
Pluto recalcula automáticamente las matrices, los eigenvalores, la ganancia `K`,
las gráficas y la animación. Este es el proceso completo.

### Paso 0. Requisito (una sola vez)

Instalar las dependencias del proyecto. Desde la carpeta del proyecto:

```bash
julia setup.jl
```

(Alternativa equivalente: `julia --project=. -e "using Pkg; Pkg.instantiate()"`.)

### Paso 1. Abrir Julia en la carpeta del proyecto

```bash
cd ruta/al/Proyecto
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
primera vez puede tardar un poco mientras precompila. Deja esta terminal
abierta: es el servidor; si la cierras, se cierra Pluto.

### Paso 3. Abrir un notebook

En la pantalla de inicio de Pluto, en el campo **"Open a notebook"**, escribe o
pega la ruta del notebook que quieras y pulsa **Open**:

- `notebooks/01_exploracion_simple.jl` — péndulo simple (estado en R^4)
- `notebooks/02_exploracion_doble.jl` — péndulo doble (estado en R^6)

La primera celda activa el proyecto (`Pkg.activate("..")`) y la segunda carga
los módulos de `src/`. La primera apertura precompila CairoMakie y
DifferentialEquations (puede tardar varios minutos); las siguientes son rápidas.

### Paso 4. Interactuar con los sliders

Al cambiar cualquier slider, todas las celdas que dependen de él se recalculan
solas. Controles disponibles:

| Notebook | Sliders de parámetros | Sliders de control | Condición inicial |
|---|---|---|---|
| `01_exploracion_simple.jl` | `M`, `m`, `Lbar`, `g`, `b` | `Q11`, `Q33`, `R` | `theta0` |
| `02_exploracion_doble.jl` | `M`, `m1`, `m2`, `L1`, `L2`, `g` | `Q` (pos, `theta1`, `theta2`), `R` | `theta1_0`, `theta2_0` |

Verás reaccionar en vivo: el eigenvalor inestable, la ganancia `K`, los polos de
lazo cerrado, las gráficas de respuesta y la animación (con su propio slider de
tiempo).

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

- **Simple:** espectro de lazo abierto $\{+4.21,\,0,\,-0.077,\,-4.23\}$ (1 modo
  inestable); LQR $K=(-3.16,-4.69,-45.39,-10.93)$; Ackermann con polos
  $\{-1,-2,-3,-4\}$ da $K=(-1.75,-3.75,-39.01,-9.60)$.
- **Doble:** espectro de lazo abierto $\{+8.57,\,+4.09,\,0,\,0,\,-4.09,\,-8.57\}$
  (2 modos inestables); rangos de controlabilidad y observabilidad $6/6$;
  LQR $K=(3.16,\,5.82,\,-191.55,\,-10.99,\,228.32,\,36.14)$.

Ambas configuraciones reproducen los valores del informe (`docs/resumen_tecnico/`).

## Nota sobre el entorno

`Manifest.toml` proviene de la versión v3. Tras clonar o tras unificar el
entorno, conviene ejecutar una vez:

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.resolve()
```
