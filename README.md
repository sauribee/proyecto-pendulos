# Péndulo invertido sobre un carro — Álgebra Lineal Aplicada

Análisis y control del péndulo invertido sobre un carro mediante herramientas de
álgebra lineal: representación en espacio de estados, análisis espectral,
controlabilidad y observabilidad de Kalman, asignación de polos (Ackermann),
regulación lineal cuadrática (LQR), márgenes de controlabilidad por valores
singulares y estimación de estado por dualidad.

El proyecto estudia **tres configuraciones** —un péndulo simple, uno doble y uno
triple— que comparten la misma maquinaria genérica de análisis y control, y las
trata como tres instancias de una sola familia de $N$ eslabones. La tesis que lo
recorre: *el álgebra lineal escala sin cambios conceptuales al crecer el número
de eslabones, pero el condicionamiento no.*

Universidad Nacional de Colombia, Sede Medellín — Facultad de Ciencias.
Autores: Mateo Bedoya Rojas, Camilo Alejandro Patiño Osorio, Santiago Uribe Echavarría.

---

## Contenido

- [Las dos entregas](#las-dos-entregas)
- [Estructura del repositorio](#estructura-del-repositorio)
- [Requisitos](#requisitos)
- [Uso rápido](#uso-rápido)
- [Exploración interactiva con Pluto](#exploración-interactiva-con-pluto-paso-a-paso)
- [Resultados](#resultados)
- [Cuando el rango deja de servir: el margen PBH](#cuando-el-rango-deja-de-servir-el-margen-pbh)
- [Observador de estado y dualidad](#observador-de-estado-y-dualidad)
- [Atlas de operabilidad](#atlas-de-operabilidad)
- [Documentos](#documentos)
- [Regenerar figuras y compilar los PDFs](#regenerar-figuras-y-compilar-los-pdfs)
- [Pruebas de regresión](#pruebas-de-regresión)
- [Nota sobre el entorno](#nota-sobre-el-entorno)

---

## Las dos entregas

| | **Primera entrega** | **Segunda entrega** |
|---|---|---|
| Configuraciones | I (simple, $\mathbb{R}^4$) y II (doble, $\mathbb{R}^6$) | III (triple, $\mathbb{R}^8$) y la familia genérica de $N$ eslabones |
| Contenido | Euler-Lagrange, linealización, espectro, Kalman, LQR vía Riccati, Ackermann, simulación no lineal | Margen PBH, condicionamiento, elipsoide de no saturación, discretización ZOH, atlas de operabilidad, observador de Luenberger, pruebas de regresión |
| Entregables | Informe técnico, resumen ejecutivo, presentación | Informe técnico, resumen ejecutivo |

La segunda entrega ataca los cinco huecos que dejaba la primera: la
observabilidad se analizaba y no se usaba; no había noción de límite de validez;
el criterio de rango es binario y frágil, y no se decía; no había pruebas
automatizadas; y con solo dos configuraciones la afirmación "el método escala"
era anecdótica.

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
│   ├── observer.jl               Observador de Luenberger por dualidad, separacion, sensores
│   ├── sweep.jl                  Barridos y fronteras de operabilidad (genericos)
│   ├── animation_simple.jl       Animacion del pendulo simple (GLMakie)
│   ├── animation_double.jl       Animacion del pendulo doble (GLMakie)
│   └── animation_triple.jl       Animacion de N eslabones (GLMakie)
├── test/
│   └── runtests.jl               618 pruebas de regresion en ocho conjuntos
├── notebooks/                    Exploradores interactivos (Pluto)
│   ├── 01_exploracion_simple.jl  Configuracion I   + metricas, observador y barridos
│   ├── 02_exploracion_doble.jl   Configuracion II  + metricas, observador y barridos
│   └── 03_exploracion_triple.jl  Configuracion III + metricas, observador y barridos
├── figures/                      Salidas de los pipelines (ignorada por git)
└── docs/
    ├── Entrega_1/                  PRIMERA ENTREGA (Configuraciones I y II)
    │   ├── informe_tecnico/
    │   │   ├── informe_tecnico.tex
    │   │   ├── informe_tecnico.pdf
    │   │   ├── make_report_figs.jl
    │   │   └── figs/
    │   ├── resumen_ejecutivo/
    │   │   ├── resumen_ejecutivo.tex
    │   │   └── resumen_ejecutivo.pdf
    │   └── presentacion/
    │       ├── presentacion.tex
    │       ├── presentacion.pdf
    │       ├── make_slide_figs.jl
    │       └── figs/
    └── Entrega_2/                  SEGUNDA ENTREGA (Configuracion III y analisis)
        ├── informe_tecnico/
        │   ├── informe_tecnico.tex
        │   ├── informe_tecnico.pdf
        │   ├── make_report_figs.jl
        │   └── figs/
        ├── resumen_ejecutivo/
        │   ├── resumen_ejecutivo.tex
        │   └── resumen_ejecutivo.pdf
        └── atlas/                  Atlas de operabilidad
            ├── make_atlas_figs.jl
            ├── figs/               atlas_a.png ... atlas_d.png
            └── data/               cache CSV de los barridos (ignorado por git)
```

Cada entrega es autocontenida: su informe técnico, su resumen ejecutivo y sus
figuras viven bajo la misma carpeta. Dentro de una entrega, los documentos se
referencian entre sí con rutas cortas (`../informe_tecnico/figs/`,
`../atlas/figs/`) mediante `\graphicspath`, de modo que ninguna figura se
duplica y ningún cálculo se repite.

### El diseño del código

Las tres configuraciones comparten los módulos `Controller` (LQR, Ackermann,
Riccati), `Metrics` (margen PBH, condicionamiento, elipsoide de no saturación,
ZOH), `Sweep` (barridos y fronteras) y `Observer` (dualidad, separación,
sensores), todos **genéricos en la dimensión del estado**: el mismo código
analiza el sistema de $\mathbb{R}^4$, el de $\mathbb{R}^6$ y el de
$\mathbb{R}^8$ sin cambios. Agregar la Configuración III no obligó a tocar una
sola línea de `Controller`.

`ModelNLink` va un paso más allá: es un modelo genérico de $N$ eslabones que
**contiene a las tres configuraciones como casos particulares** ($N=1$ con barra
uniforme, $N=2$ con masas puntuales, $N=3$ con barras uniformes). Las
Configuraciones I y II conservan su implementación propia e independiente, de
modo que la coincidencia entre ambas es una prueba de regresión genuina y no un
acto de fe: la discrepancia máxima medida es $5.3\times10^{-15}$ para $N=1$ y
**exactamente cero** para $N=2$.

## Requisitos

- [Julia](https://julialang.org/) (probado con la versión del `Manifest.toml`).
- Las dependencias se declaran en `Project.toml` y se fijan en `Manifest.toml`.
- Para compilar los PDFs: una distribución de LaTeX con `latexmk` (MiKTeX o TeX Live).

Paquetes principales: `DifferentialEquations` (solver `Tsit5`),
`CairoMakie` / `GLMakie` (gráficas y animación), `MatrixEquations` y
`ForwardDiff` (verificación de la Riccati y de los jacobianos, en las pruebas) y
`Pluto` / `PlutoUI` (notebooks interactivos). Todo lo declarado se usa: cada
dependencia de `Project.toml` aparece en algún archivo del repositorio.

Las rutinas de álgebra lineal centrales —linealización, Riccati, Ackermann,
margen PBH, discretización ZOH— **se programaron a mano desde las
definiciones**; las librerías especializadas se usan solo como verificación.

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
julia> include("main_triple.jl")   # o main_simple.jl, o main_double.jl
```

Cada pipeline ejecuta el flujo completo: define parámetros, simula la respuesta
libre (sin control), linealiza y analiza (eigenvalores, controlabilidad,
observabilidad), diseña el o los controladores, simula el lazo cerrado, genera
las gráficas comparativas en `figures/` y produce la animación. El de la
Configuración III añade las métricas de viabilidad y el observador de Luenberger.

El **límite del actuador** no es el mismo en los tres, porque la exigencia crece
con la complejidad: 50 N en el simple, 100 N en el doble y 150 N en el triple.
Las métricas que se comparan entre configuraciones (el ángulo garantizado sin
saturar, las fronteras prácticas) se calculan todas con $u_{\max} = 50$ N, para
que la comparación sea directa.

## Exploración interactiva con Pluto (paso a paso)

Los notebooks de `notebooks/` permiten correr todo el análisis de forma
interactiva: se mueve un slider (una masa, una longitud, un peso del LQR) y
Pluto recalcula automáticamente las matrices, los eigenvalores, la ganancia $K$,
las métricas, el observador, las gráficas y la animación.

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
- `notebooks/03_exploracion_triple.jl` — péndulo triple (estado en $\mathbb{R}^8$)

La primera celda activa el proyecto y la segunda carga los módulos de `src/`. La
primera apertura precompila CairoMakie y DifferentialEquations (puede tardar
varios minutos); las siguientes son rápidas.

### Paso 4. Interactuar con los sliders

Al cambiar cualquier slider, todas las celdas dependientes se recalculan solas.

| Notebook | Sliders de parámetros | Sliders de control | Condición inicial |
|---|---|---|---|
| `01_exploracion_simple.jl` | `M`, `m`, `Lbar`, `g`, `b` | `Q11`, `Q33`, `R` | `theta0` |
| `02_exploracion_doble.jl` | `M`, `m1`, `m2`, `L1`, `L2`, `g` | `Q` (pos, `theta1`, `theta2`), `R` | `theta1_0`, `theta2_0` |
| `03_exploracion_triple.jl` | `M`, `m1`–`m3`, `l1`–`l3`, `g`, `b` | `Q` (pos, ángulos), `R` | `theta0` (los tres) |

Los tres llevan una **segunda parte** con todo lo que agregó la segunda entrega,
reaccionando a los mismos sliders y con controles propios (`u_max`, `h`, `Qo`,
el parámetro a barrer, la longitud del riel):

- **margen PBH** de control y observación, condicionamiento de Kalman y el modo
  que hace de cuello de botella;
- **elipsoide de no saturación**: $c^\star$ y el ángulo garantizado;
- **periodo de muestreo máximo**, con los polos discretos entrando y saliendo del
  círculo unitario al mover $h$;
- **observador de Luenberger** por dualidad: polos de $A-LC$, principio de
  separación verificado en vivo, tabla de subconjuntos de sensores y simulación
  con el estimador arrancando en cero;
- **barridos 1D** de parámetros y la frontera práctica por bisección sobre el
  modelo no lineal.

### Paso 5. Exportar la animación (opcional)

Al final de la primera parte de cada notebook hay una casilla (checkbox). Al
marcarla se genera el archivo de animación en `figures/`:

- simple: `04_comparacion_libre_vs_lqr.gif` (y `05_..._.mp4` si hay `ffmpeg`)
- doble: `08_doble_exploracion.gif`
- triple: `12_triple_exploracion.gif`

Desmárcala para no regenerar el archivo en cada cambio de slider.

### Paso 6. Cerrar

Guarda el notebook (Pluto guarda solo el `.jl`), cierra la pestaña del navegador
y detén el servidor con `Ctrl-C` en la terminal de Julia.

## Resultados

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
  con $\lVert K\rVert_2 = 1176$.

Sobre el modelo **no lineal** en lazo cerrado, desde la condición inicial nominal
de cada configuración:

| | Simple (LQR) | Simple (Ackermann) | Doble | Triple |
|---|---|---|---|---|
| $t_s$ [s] | 1.45 | 1.54 | 1.75–1.77 | 1.84–1.88 |
| $\max\lvert u\rvert$ [N] | 6.8 | 5.9 | 4.2 | 7.29 |
| Excursión del carro [m] | 0.33 | — | 0.43 | 0.47 |

Las tres configuraciones reproducen los valores de los informes técnicos
(`docs/Entrega_*/informe_tecnico/`) y están blindadas por `test/runtests.jl`.

### Comparativa de las tres configuraciones

Con los parámetros por defecto de cada una:

| Cantidad | Simple | Doble | Triple |
|---|---|---|---|
| $n$ | 4 | 6 | 8 |
| Modos inestables | 1 | 2 | 3 |
| $\lambda_{\max}$ [s$^{-1}$] | 4.21 | 8.57 | 18.06 |
| $\tau_{\text{caída}}$ [ms] | 237.5 | 116.7 | 55.4 |
| $\mathrm{cond}\,\mathcal{C}$ | $3.6\times10^{1}$ | $1.6\times10^{4}$ | $2.2\times10^{8}$ |
| $\hat\mu$ (margen PBH) | $1.72\times10^{-2}$ | $2.56\times10^{-3}$ | $4.92\times10^{-4}$ |
| $\lVert K\rVert_2$ | 47.0 | 300.5 | 1176.0 |
| $\theta_0$ garantizado sin saturar ($u_{\max}=50$ N) | $55.1°$ | $13.4°$ | $8.7°$ |
| $h_{\max}$ (ZOH) [ms] | 197.4 | 77.3 | 32.0 |

El triple exige muestrear por encima de **31 Hz**, y el margen relativo respecto
al tiempo de caída se estrecha al crecer la complejidad
($h_{\max}\lambda_{\max} = 0.83 \to 0.66 \to 0.58$): no basta con escalar la
frecuencia proporcionalmente al modo dominante.

## Cuando el rango deja de servir: el margen PBH

Éste es el resultado central de la segunda entrega. La matriz de Kalman
$\mathcal{C}=[B\ \ AB\ \cdots\ A^{n-1}B]$ es una **base de Krylov**, notoriamente
mal condicionada. Medido sobre una familia normalizada (longitud total 1 m, masa
total 0.6 kg repartidas en $N$ barras uniformes, $M=1$ kg, $b=0.1$):

| $N$ | $n$ | $\lambda_{\max}$ [s$^{-1}$] | $\hat\mu$ (PBH) | $\mathrm{cond}\,\mathcal{C}$ | $\lVert K\rVert_2$ |
|---|---|---|---|---|---|
| 1 | 4  | 4.51  | $1.31\times10^{-2}$ | $4.8\times10^{1}$  | 52.9 |
| 2 | 6  | 10.59 | $1.79\times10^{-3}$ | $6.2\times10^{4}$  | 248.7 |
| 3 | 8  | 18.06 | $4.92\times10^{-4}$ | $2.2\times10^{8}$  | 1176.0 |
| 4 | 10 | 26.54 | $1.90\times10^{-4}$ | $1.7\times10^{12}$ | 5266.8 |
| 5 | 12 | 35.64 | $9.03\times10^{-5}$ | $3.3\times10^{16}$ | 23221.9 |

> Ojo: esta familia **no** coincide con los parámetros por defecto de las
> Configuraciones I y II, que tienen otras masas. Es una serie aparte,
> construida para que la comparación entre valores de $N$ sea limpia.

Con $N=5$ el condicionamiento alcanza $1/\varepsilon_{\text{máq}}$: el rango
numérico de $\mathcal{C}$ es **ruido**, y sin embargo el sistema *sigue siendo*
controlable en aritmética exacta. La herramienta que lo reemplaza es la versión
cuantitativa del test de Popov–Belevitch–Hautus,

$$\mu_{\text{PBH}}=\min_{\lambda\in\sigma(A)}\ \sigma_{\min}\big[\,A-\lambda I \ \big|\ B\,\big],
\qquad \hat\mu=\frac{\mu_{\text{PBH}}}{\lVert[\,A\ \ B\,]\rVert_2}$$

que por el **teorema de Eckart–Young** es exactamente la distancia en norma 2 al
conjunto de matrices de rango deficiente: mide *cuánto habría que perturbar el
sistema para volverlo incontrolable*. A diferencia de
$\mathrm{cond}\,\mathcal{C}$, que estalla casi quince órdenes de magnitud,
el margen PBH decae poco más de dos órdenes —un factor 145— y sigue siendo
utilizable para comparar y para barrer.

## Observador de estado y dualidad

El control $u = -K\mathbf{x}$ usa el estado **completo**, que en un sistema real
no se mide: nadie pone un sensor de velocidad angular en cada articulación. Como
$\mathrm{rank}\,\mathcal{O} = n$, el estado se puede reconstruir con un
observador de Luenberger

$$\dot{\hat{\mathbf{x}}} = A\hat{\mathbf{x}} + Bu + L(\mathbf{y} - C\hat{\mathbf{x}}),
\qquad \dot{\mathbf{e}} = (A - LC)\,\mathbf{e}$$

y por la **dualidad de Kalman** la ganancia sale de la misma rutina que diseñó
$K$, sin escribir ningún algoritmo nuevo: $L = \mathrm{lqr}(A^\top, C^\top, Q_o, R_o)^\top$.

El **principio de separación** sale de dos líneas de álgebra lineal: en
coordenadas $(\mathbf{x}, \mathbf{e})$ la matriz de lazo cerrado es triangular
por bloques,

$$\begin{pmatrix}A - BK & BK\\ 0 & A - LC\end{pmatrix}
\quad\Longrightarrow\quad
\sigma = \sigma(A-BK)\ \cup\ \sigma(A-LC)$$

así que controlador y observador se diseñan por separado. Verificado
numéricamente en las tres configuraciones: el error máximo al emparejar
espectros es del orden de $10^{-14}$ (simple) y $10^{-13}$ (doble), y de
$3.4\times10^{-12}$ en el triple, que es el peor de los tres y el que cita el
informe.

**¿Qué sensores hacen falta?** De los $2^4-1 = 15$ subconjuntos de
$\{x, \theta_1, \theta_2, \theta_3\}$ en la Configuración III:

- Todo subconjunto que **incluya la posición del carro** es observable ($8/8$).
  Son los $2^3 = 8$ que la contienen.
- Ninguno que la excluya lo es: los $7$ restantes se quedan en $7/8$. Midiendo
  solo ángulos, la posición absoluta del carro es inobservable, porque la
  dinámica angular es invariante a trasladar el carro. Falta exactamente una
  dirección.
- El juego **mínimo** es un único sensor: la posición del carro. Basta para
  reconstruir los ocho estados — pero con margen PBH $1.5\times10^{-6}$, dos
  órdenes por debajo del juego completo ($1.6\times10^{-4}$). Observable no es
  lo mismo que practicable.

**El observador no sale gratis.** El principio de separación es un teorema sobre
el sistema **lineal**. Sobre el modelo no lineal la región de atracción se encoge
entre la mitad y la cuarta parte, y por un motivo contrario al que cabría
esperar: el estimador arranca en cero, así que el mando arranca en cero también y
durante el transitorio el controlador aplica fuerza *de menos*, no de más —desde
la misma condición inicial su pico llega a ser menor que el del lazo con estado
completo: 25.8 N frente a 31.7 N desviando un ángulo 0.10 rad—. Lo que crece es
la **excursión**: mientras la estimación converge, la planta cae gobernada por un
mando que todavía no la describe, y el ángulo máximo del transitorio se
multiplica por un factor de entre 1.1 y 2.8. Es eso lo que acaba sacándola del
régimen donde la linealización vale:

| Ángulos desviados | Estado completo [rad] | Con observador [rad] | Ratio |
|---|---|---|---|
| 1 | 0.2605 | 0.1452 | 56 % |
| 2 | 0.1419 | 0.0450 | 32 % |
| 3 | 0.3018 | 0.0743 | 25 % |

La condición inicial nominal del informe (los tres eslabones a 0.10 rad) está
*dentro* de la región con estado completo y *fuera* de la del lazo con
observador.

## Atlas de operabilidad

Barre el espacio de parámetros para identificar dónde el controlador funciona
bien, dónde apenas y dónde ya es imposible:

```bash
julia --project=. -t auto docs/Entrega_2/atlas/make_atlas_figs.jl
```

El flag `-t auto` importa: los mapas B y D integran el modelo no lineal en cada
punto de la malla y se paralelizan sobre las filas. Los resultados se cachean en
`docs/Entrega_2/atlas/data/*.csv`; para forzar el recómputo, borrar esa carpeta
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

Y dentro de la práctica, **cuál restricción está activa cambia con la
configuración**. La columna "sin riel" repite el cálculo eliminando la
restricción de carrera; es contra ésa, y no contra la nominal, que hay que
comparar la cota elipsoidal:

| Configuración | Cota $c^\star$ [rad] | Frontera real [rad] | Sin riel [rad] | Restricción activa |
|---|---|---|---|---|
| Simple | 0.962 | 0.580 | 1.063 | **Riel** (1.5 m) |
| Doble | 0.234 | 0.345 | 0.345 | Saturación |
| Triple | 0.152 | 0.230 | 0.230 | Saturación |

Otros dos hallazgos del atlas: el margen PBH tiene un **máximo interior** en la
masa del carro (alrededor de 0.13 kg), y el **eslabón superior corto es el
peligroso** —con $\ell_3 = 2$ cm el modo dominante sube a 42.9 s$^{-1}$ y el
margen se hunde a $3.4\times10^{-5}$—. Nótese además que $c^\star$ y $\hat\mu$
apuntan en direcciones opuestas al barrer $M$: no existe *la* mejor
configuración sin decir antes qué se optimiza.

## Documentos

| Documento | Ruta | Páginas |
|---|---|---|
| Informe técnico, 1.ª entrega | `docs/Entrega_1/informe_tecnico/` | 20 |
| Resumen ejecutivo, 1.ª entrega | `docs/Entrega_1/resumen_ejecutivo/` | 3 |
| Presentación (Beamer 16:9, 20 min) | `docs/Entrega_1/presentacion/` | 35 |
| Informe técnico, 2.ª entrega | `docs/Entrega_2/informe_tecnico/` | 19 |
| Resumen ejecutivo, 2.ª entrega | `docs/Entrega_2/resumen_ejecutivo/` | 4 |

Cada entrega se documenta por completo en su informe técnico. El resumen
ejecutivo de cada una sirve de guía de lectura y remite, sección por sección, al
informe correspondiente.

## Regenerar figuras y compilar los PDFs

Las figuras de los informes **no** se dibujan a mano: las produce un script que
reutiliza los módulos de `src/` e imprime por consola las métricas que se citan
en el texto.

```bash
# Figuras del informe de la primera entrega
julia --project=. docs/Entrega_1/informe_tecnico/make_report_figs.jl

# Figuras del informe de la segunda entrega
julia --project=. -t auto docs/Entrega_2/informe_tecnico/make_report_figs.jl

# Mapas del atlas de operabilidad
julia --project=. -t auto docs/Entrega_2/atlas/make_atlas_figs.jl

# Figuras de las diapositivas (solo primera entrega)
julia --project=. docs/Entrega_1/presentacion/make_slide_figs.jl
```

El script de la segunda entrega imprime además la tabla comparativa de las tres
configuraciones, las fronteras prácticas, la tabla de subconjuntos de sensores y
la comparación de regiones de atracción con y sin observador.

Cada documento se compila **desde su propia carpeta**, porque las rutas de las
figuras son relativas. Las `.png` deben existir antes de compilar:

```bash
cd docs/Entrega_1/informe_tecnico    && latexmk -pdf informe_tecnico.tex
cd ../resumen_ejecutivo              && latexmk -pdf resumen_ejecutivo.tex
cd ../presentacion                   && latexmk -pdf presentacion.tex
cd ../../Entrega_2/informe_tecnico   && latexmk -pdf informe_tecnico.tex
cd ../resumen_ejecutivo              && latexmk -pdf resumen_ejecutivo.tex
```

## Pruebas de regresión

```bash
julia --project=. -t auto test/runtests.jl
```

Son **618 pruebas** en ocho conjuntos (unos 45 s con `-t auto` tras
precompilar). Verifican de forma reproducible lo que antes solo se había
comprobado a mano:

1. **Modelos** — el modelo genérico de $N$ eslabones reproduce las
   Configuraciones I y II, que son implementaciones **independientes** de la
   misma física; $\mathbf{M}(q)$ es simétrica definida positiva en
   configuraciones angulares aleatorias; la fricción articular es disipativa.
2. **Linealización** — los jacobianos **analíticos** escritos a mano coinciden
   con la diferenciación automática de las ecuaciones no lineales, tanto en $A$
   como en $B$; y el espectro tiene la estructura $\pm\lambda$ más cero doble
   predicha, para $N=1,\dots,5$.
3. **Valores publicados** — los espectros y las ganancias $K$ de los documentos
   de las dos entregas.
4. **Controller** — `solve_care` contra `MatrixEquations.arec`, el residuo de la
   CARE, $P = P^\top \succ 0$ y Cayley-Hamilton.
5. **Métricas** — margen PBH, condicionamiento, $c^\star$, ángulo garantizado,
   discretización ZOH exacta y criterio de Schur.
6. **Degradación con $N$** — la tabla de arriba, con el decaimiento monótono del
   margen PBH frente al crecimiento del condicionamiento.
7. **Atlas de operabilidad** — periodo de muestreo máximo, el óptimo interior de
   la masa del carro, la geometría derivada al barrer longitudes, y que la cota
   elipsoidal queda dentro de la frontera de saturación.
8. **Observador y dualidad** — la dualidad exacta, el principio de separación,
   la partición de subconjuntos de sensores y el encogimiento de la región de
   atracción al estimar en vez de medir.

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
