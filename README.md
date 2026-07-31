# Péndulo invertido sobre un carro — Álgebra Lineal Aplicada

Análisis y control del péndulo invertido sobre un carro mediante herramientas de
álgebra lineal: representación en espacio de estados, análisis espectral,
controlabilidad y observabilidad de Kalman, asignación de polos (Ackermann) y
regulación lineal cuadrática (LQR). El proyecto estudia **tres configuraciones**
—un péndulo simple, uno doble y uno triple— que comparten la misma maquinaria
genérica de análisis y control, y las trata como tres instancias de una sola
familia de $N$ eslabones.

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
├── main_triple.jl                Pipeline del pendulo TRIPLE (Configuracion III)
├── README.md
├── src/
│   ├── model_simple.jl           Pendulo simple: parametros, EOM y lazo cerrado
│   ├── model_double.jl           Pendulo doble: parametros, EOM y lazo cerrado
│   ├── model_nlink.jl            Modelo GENERICO de N eslabones (contiene a I, II y III)
│   ├── model_triple.jl           Pendulo triple: capa delgada sobre ModelNLink
│   ├── linearization.jl          Linealizacion (simple, doble, generica en N), espectro, Kalman
│   ├── controller.jl             LQR (Riccati via Hamiltoniano) y Ackermann (genericos)
│   ├── metrics.jl                Margen PBH, condicionamiento, elipsoide, ZOH (genericos)
│   ├── observer.jl               Observador de Luenberger por dualidad, separacion
│   ├── sweep.jl                  Barridos y fronteras de operabilidad (genericos)
│   ├── animation_simple.jl       Animacion del pendulo simple (GLMakie)
│   ├── animation_double.jl       Animacion del pendulo doble (GLMakie)
│   └── animation_triple.jl       Animacion de N eslabones (GLMakie)
├── test/
│   └── runtests.jl               Pruebas de regresion
├── notebooks/                    Exploradores interactivos (Pluto)
│   ├── 01_exploracion_simple.jl
│   └── 02_exploracion_doble.jl                  
└── docs/
    ├── guia_maestra.md             Documento maestro (teoria + codigo + defensa)
    ├── segunda_entrega.md          Plan de trabajo de la segunda entrega
    ├── resumen_tecnico/
    │   ├── entrega_1/              Primera entrega (Configuraciones I y II)
    │   │   ├── resumen_tecnico.tex
    │   │   ├── resumen_tecnico.pdf
    │   │   ├── make_report_figs.jl
    │   │   └── figs/
    │   └── entrega_2/              Segunda entrega (Configuracion III y analisis)
    │       ├── resumen_tecnico.tex
    │       ├── resumen_tecnico.pdf
    │       ├── make_report_figs.jl
    │       └── figs/
    ├── resumen_ejecutivo/
    │   ├── entrega_1/
    │   │   ├── resumen_ejecutivo.tex
    │   │   └── resumen_ejecutivo.pdf
    │   └── entrega_2/
    │       ├── resumen_ejecutivo.tex
    │       └── resumen_ejecutivo.pdf
    ├── segunda_entrega/            Atlas de operabilidad
    │   ├── make_atlas_figs.jl
    │   └── figs/                   atlas_a.png ... atlas_d.png
    └── presentacion/               Solo primera entrega (la segunda no la lleva)
        ├── presentacion.tex        Diapositivas (Beamer, 16:9, 20 minutos)
        ├── presentacion.pdf
        ├── make_slide_figs.jl
        └── figs/
```

Cada entrega tiene su propio informe técnico y su resumen ejecutivo, con sus
figuras. El informe de la segunda entrega toma además los cuatro mapas del atlas
de `docs/segunda_entrega/figs/` mediante `\graphicspath`, para no duplicar ni el
cálculo ni los archivos.

Las tres configuraciones comparten los módulos `Controller` (LQR, Ackermann,
Riccati), `Metrics` (margen PBH, condicionamiento, elipsoide de no saturación) y
las funciones de análisis de `Linearization`, todos **genéricos en la dimensión
del estado**: el mismo código analiza el sistema de $\mathbb{R}^4$, el de
$\mathbb{R}^6$ y el de $\mathbb{R}^8$ sin cambios.

`ModelNLink` va un paso más allá: es un modelo genérico de $N$ eslabones que
**contiene a las tres configuraciones como casos particulares** ($N=1$ con barra
uniforme, $N=2$ con masas puntuales, $N=3$ con barras uniformes). Las
Configuraciones I y II conservan su implementación propia e independiente, de
modo que la coincidencia entre ambas es una prueba de regresión genuina y no un
acto de fe.


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

# Pendulo triple (Configuracion III, estado en R^8)
julia main_triple.jl
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
- **Triple:** espectro de lazo abierto
  $\{+18.06,\ +9.77,\ +4.38,\ 0,\ -0.06,\ -4.40,\ -9.78,\ -18.06\}$ (3 modos
  inestables); rangos $8/8$;
  LQR $K=(-3.16,\,-6.17,\,-317.43,\,-4.37,\,911.95,\,37.73,\,-667.37,\,-61.42)$
  con $\|K\|_2 = 1176$.

Las tres configuraciones reproducen los valores del informe técnico
(`docs/resumen_tecnico/`) y están blindadas por `test/runtests.jl`.

### Comparativa de las tres configuraciones

Con los parámetros por defecto de cada una:

| Cantidad | Simple | Doble | Triple |
|---|---|---|---|
| $n$ | 4 | 6 | 8 |
| Modos inestables | 1 | 2 | 3 |
| $\lambda_{\max}$ [s$^{-1}$] | 4.21 | 8.57 | 18.06 |
| $\tau_{\text{caída}}$ [ms] | 237.5 | 116.7 | 55.4 |
| $\operatorname{cond}\mathcal{C}$ | $3.6\times10^{1}$ | $1.6\times10^{4}$ | $2.2\times10^{8}$ |
| $\hat\mu$ (margen PBH) | $1.72\times10^{-2}$ | $2.56\times10^{-3}$ | $4.92\times10^{-4}$ |
| $\|K\|_2$ | 47.0 | 300.5 | 1176.0 |
| $\theta_0$ garantizado sin saturar ($u_{\max}=50$ N) | $55.1°$ | $13.4°$ | $8.7°$ |
| $h_{\max}$ (ZOH) [ms] | 197.4 | 77.3 | 32.0 |

El rango de Kalman vale $n$ en las tres, pero $\operatorname{cond}\mathcal{C}$
crece unos cuatro órdenes de magnitud por eslabón: es un teorema correcto que
deja de ser un buen algoritmo. El margen PBH, en cambio, decae suavemente y sí
sirve para comparar.

## Regenerar las figuras del informe

Las figuras de respuesta temporal que aparecen en el informe técnico se generan
con un script aparte (usa CairoMakie, salida estática):

```bash
julia --project=. docs/resumen_tecnico/entrega_1/make_report_figs.jl
```

Esto reescribe `docs/resumen_tecnico/entrega_1/figs/` y reporta las métricas
(tiempo de asentamiento, esfuerzo de control pico) que se citan en la discusión.

Las de la segunda entrega, análogamente:

```bash
julia --project=. -t auto docs/resumen_tecnico/entrega_2/make_report_figs.jl
```

Además de las figuras, este último imprime la tabla comparativa de las tres
configuraciones, las fronteras prácticas, la tabla de subconjuntos de sensores y
la comparación de regiones de atracción con y sin observador.

Las figuras de las diapositivas (solo primera entrega) se regeneran así:

```bash
julia --project=. docs/presentacion/make_slide_figs.jl
```

Para compilar cualquiera de los cuatro documentos, desde su propia carpeta:

```bash
latexmk -pdf resumen_tecnico.tex
```

## Observador de estado y dualidad

El control $u = -K\mathbf{x}$ usa el estado **completo**, que en un sistema real
no se mide: nadie pone un sensor de velocidad angular en cada articulación. Como
$\operatorname{rank}\mathcal{O} = n$, el estado se puede reconstruir con un
observador de Luenberger

$$\dot{\hat{\mathbf{x}}} = A\hat{\mathbf{x}} + Bu + L(\mathbf{y} - C\hat{\mathbf{x}}),
\qquad \dot{\mathbf{e}} = (A - LC)\,\mathbf{e}$$

y por la **dualidad de Kalman** la ganancia sale de la misma rutina que diseñó
$K$, sin escribir ningún algoritmo nuevo: $L = \operatorname{lqr}(A^\top, C^\top, Q_o, R_o)^\top$.

El **principio de separación** sale de dos líneas de álgebra lineal: en
coordenadas $(\mathbf{x}, \mathbf{e})$ la matriz de lazo cerrado es triangular
por bloques,

$$\begin{pmatrix}A - BK & BK\\ 0 & A - LC\end{pmatrix}
\quad\Longrightarrow\quad
\sigma = \sigma(A-BK)\ \cup\ \sigma(A-LC)$$

así que controlador y observador se diseñan por separado. Verificado
numéricamente a $3\times10^{-12}$ en las tres configuraciones.

**¿Qué sensores hacen falta?** De los $2^4-1 = 15$ subconjuntos de
$\{x, \theta_1, \theta_2, \theta_3\}$ en la Configuración III:

- Todo subconjunto que **incluya la posición del carro** es observable ($8/8$).
- Ninguno que la excluya lo es: se quedan en $7/8$. Midiendo solo ángulos, la
  posición absoluta del carro es inobservable, porque la dinámica angular es
  invariante a trasladar el carro.
- El juego **mínimo** es un único sensor: la posición del carro. Basta para
  reconstruir los ocho estados — pero con margen PBH $1.5\times10^{-6}$, dos
  órdenes por debajo del juego completo ($1.6\times10^{-4}$). Observable no es
  lo mismo que practicable.

## Atlas de operabilidad

Barre el espacio de parámetros para identificar dónde el controlador funciona
bien, dónde apenas y dónde ya es imposible:

```bash
julia --project=. -t auto docs/segunda_entrega/make_atlas_figs.jl
```

El flag `-t auto` importa: los mapas B y D integran el modelo no lineal en cada
punto de la malla y se paralelizan sobre las filas. Los resultados se cachean en
`docs/segunda_entrega/data/*.csv`; para forzar el recómputo, borrar esa carpeta
o exportar `ATLAS_FORCE=1`.

| Mapa | Ejes | Qué muestra |
|---|---|---|
| A | $(\ell_3/\ell_1,\ m_3/m_1)$ | Región estructuralmente sana (margen PBH) |
| B | $(u_{\max},\ \theta_0)$ | **La frontera de operabilidad práctica** |
| C | $(M,\ m_{\text{total}})$ | El óptimo interior de la masa del carro |
| D | $(R,\ \theta_0)$ | Compromiso agresividad-robustez |

Conviene distinguir **dos tipos de imposibilidad**:

- **Estructural:** el par $(A,B)$ deja de ser controlable, o se acerca tanto que
  ningún $K$ razonable sirve. Es una propiedad *del sistema*. Se mide con el
  margen PBH.
- **Práctica:** el sistema es controlable, pero el actuador satura, el riel se
  acaba o el muestreo es muy lento. Es una propiedad *de la implementación*. Se
  mide simulando el modelo no lineal y buscando la frontera por bisección.

Cuál de las tres restricciones prácticas está activa **cambia con la
configuración**:

| Configuración | Cota elipsoidal $c^\star$ | Frontera real | Restricción activa |
|---|---|---|---|
| Simple | $0.962$ rad | $0.580$ rad | **Riel** (1.5 m) |
| Doble | $0.234$ rad | $0.345$ rad | Saturación |
| Triple | $0.152$ rad | $0.230$ rad | Saturación |

## Pruebas de regresión

```bash
julia --project=. test/runtests.jl
```

Verifican de forma reproducible lo que antes solo se había comprobado a mano:

- El modelo genérico de $N$ eslabones (`ModelNLink`) reproduce las
  Configuraciones I y II, que son implementaciones **independientes** de la
  misma física.
- Los jacobianos **analíticos** escritos a mano en `linearization.jl` coinciden
  con la diferenciación automática de las ecuaciones no lineales.
- Los espectros y las ganancias $K$ publicados en los cuatro documentos.
- `solve_care` contra `MatrixEquations.arec`, el residuo de la CARE,
  $P = P^\top \succ 0$ y Cayley-Hamilton.
- Los valores de referencia del módulo `Metrics`.

> No se usa `Pkg.test()`: el repositorio es un **entorno** de Julia, no un
> paquete (`Project.toml` no declara `name`, `uuid` ni `version`, y `Pkg.test()`
> exige las tres). Ejecutar el archivo directamente da el mismo resultado sin
> reestructurar el proyecto.

## Nota sobre el entorno

`Manifest.toml` fija las versiones exactas para reproducibilidad. Tras clonar el
repositorio conviene ejecutar una vez:

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
```
