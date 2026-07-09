# Guía maestra del proyecto: péndulo invertido sobre un carro

Documento maestro de referencia. Reúne, relaciona y explica en profundidad
**toda** la implementación del repositorio: la teoría de álgebra lineal, el
modelado físico, el código en Julia, los resultados numéricos y los cuatro
artefactos de documentación (informe técnico, resumen ejecutivo, presentación
Beamer y notebooks). Su objetivo es doble: servir de mapa para entender el
proyecto de principio a fin y ser la herramienta de preparación de la
presentación oral.

- **Curso:** Álgebra Lineal Aplicada. Universidad Nacional de Colombia, Sede Medellín, Facultad de Ciencias.
- **Autores:** Mateo Bedoya Rojas, Camilo Alejandro Patiño Osorio, Santiago Uribe Echavarría.
- **Lenguaje:** Julia (análisis, control, simulación, animación); LaTeX (documentos); Pluto (exploración interactiva).

---

## Índice

1. [Cómo usar este documento](#1-cómo-usar-este-documento)
2. [Panorama general: la gran idea](#2-panorama-general-la-gran-idea)
3. [Arquitectura del repositorio y flujo de datos](#3-arquitectura-del-repositorio-y-flujo-de-datos)
4. [Fundamento físico: modelado por Euler-Lagrange](#4-fundamento-físico-modelado-por-euler-lagrange)
5. [Linealización y espacio de estados](#5-linealización-y-espacio-de-estados)
6. [Las herramientas de álgebra lineal](#6-las-herramientas-de-álgebra-lineal)
7. [Diseño de controladores](#7-diseño-de-controladores)
8. [De la teoría al código: mapeo concepto a implementación](#8-de-la-teoría-al-código-mapeo-concepto-a-implementación)
9. [Resultados numéricos e interpretación](#9-resultados-numéricos-e-interpretación)
10. [Los artefactos de documentación y cómo se relacionan](#10-los-artefactos-de-documentación-y-cómo-se-relacionan)
11. [Guía para la presentación (20 minutos)](#11-guía-para-la-presentación-20-minutos)
12. [Defensa: preguntas frecuentes y puntos delicados](#12-defensa-preguntas-frecuentes-y-puntos-delicados)
13. [Glosario de símbolos y checklist final](#13-glosario-de-símbolos-y-checklist-final)
14. [Guía de ejecución paso a paso](#14-guía-de-ejecución-paso-a-paso)

---

## 1. Cómo usar este documento

Hay tres formas de leerlo según lo que se necesite:

- **Para entender el proyecto entero:** leer en orden las secciones 2 a 9. Es la
  cadena lógica completa (problema, modelo, linealización, análisis, control,
  código, resultados).
- **Para preparar la presentación:** leer la sección 2 (mensaje central), luego
  saltar a las secciones 11 y 12 (guion, tiempos, preguntas de defensa), y usar
  las secciones 4 a 7 como material de consulta para dominar cada concepto.
- **Para navegar el código:** la sección 3 (arquitectura) y la sección 8 (mapeo
  teoría-código) son el punto de entrada.

Convención de notación en todo el documento: el estado es $\mathbf{x}$, la
entrada de control (fuerza sobre el carro) es $u$, y las matrices del espacio de
estados son $A, B, C, D$. El ángulo $\theta$ se mide **desde la vertical
superior**: $\theta = 0$ es el péndulo erguido (equilibrio inestable).

---

## 2. Panorama general: la gran idea

### 2.1 El problema en una frase

Un carro se desliza sobre un riel horizontal y sobre él se articula un péndulo.
La posición de interés es la **vertical superior**, con el péndulo apuntando
hacia arriba. Esa posición es un equilibrio, pero **inestable**: la mínima
perturbación hace que el péndulo caiga. El objetivo de control es aplicar una
fuerza horizontal $u(t)$ **únicamente sobre el carro** para mantener el péndulo
erguido.

### 2.2 Por qué es interesante

Tres rasgos, presentes también en problemas reales (estabilización de cohetes,
robots bípedos, vehículos autobalanceados):

| Rasgo | Significado | Consecuencia técnica |
|---|---|---|
| **Subactuado** | Más grados de libertad que actuadores | $B$ tiene una sola columna; una entrada escalar debe gobernar todo el estado |
| **Inestable** | En lazo abierto no se sostiene | $A$ tiene al menos un eigenvalor con $\operatorname{Re}(\lambda) > 0$ |
| **No lineal** | Aparecen $\sin\theta$, $\cos\theta$ | No es de la forma $\dot{\mathbf{x}} = A\mathbf{x} + B u$; hay que linealizar |

### 2.3 El mensaje central del proyecto

> **Casi todas las preguntas relevantes sobre el sistema se responden con álgebra lineal.**

Esta es la tesis que hay que transmitir en la presentación. En concreto:

- ¿Cómo evolucionan las pequeñas desviaciones del equilibrio? Lo dice la
  **exponencial matricial** $e^{At}$.
- ¿El sistema es estable o inestable? Lo deciden los **eigenvalores** de $A$.
- ¿Podemos, con la sola fuerza del carro, influir en todos los modos? Es una
  pregunta de **rango** (controlabilidad).
- ¿Podemos reconstruir el estado completo con pocas mediciones? También es de
  **rango** (observabilidad).
- ¿Cómo diseñar el control? Consiste en **reubicar los eigenvalores** de una
  matriz mediante realimentación.

### 2.4 Dos configuraciones de complejidad creciente

El proyecto estudia dos versiones del mismo sistema, donde la segunda se obtiene
de la primera aumentando la complejidad:

| | **Configuración I: simple** | **Configuración II: doble** |
|---|---|---|
| Descripción | Una barra rígida uniforme articulada al carro | Dos eslabones en serie, masas puntuales |
| Coordenadas generalizadas | $q = (x, \theta)$, $d = 2$ | $q = (x, \theta_1, \theta_2)$, $d = 3$ |
| Estado | $\mathbf{x} \in \mathbb{R}^4$ | $\mathbf{x} \in \mathbb{R}^6$ |
| Modos inestables | 1 | 2 |
| Fricción | Sí, coeficiente $b$ | No |
| Inercia | Barra uniforme, $I = \tfrac{1}{12}m(2L)^2$ | Masas puntuales (inercia concentrada) |
| Controladores | LQR y Ackermann | LQR |

La idea pedagógica del doble es **mostrar que el mismo razonamiento algebraico
escala** (linealizar, estudiar el espectro, verificar controlabilidad, diseñar la
realimentación) sin cambios conceptuales cuando crecen los grados de libertad:
solo cambian el tamaño de las matrices y el esfuerzo de control.

---

## 3. Arquitectura del repositorio y flujo de datos

### 3.1 Árbol de archivos comentado

```
proyecto-pendulos/
├── Project.toml              Entorno Julia unico de todo el proyecto (dependencias)
├── Manifest.toml             Versiones exactas resueltas (reproducibilidad)
├── setup.jl                  Instala e instancia las dependencias
├── main_simple.jl            Pipeline completo del pendulo SIMPLE (Config. I)
├── main_double.jl            Pipeline completo del pendulo DOBLE (Config. II)
├── README.md                 Instrucciones de uso
├── src/
│   ├── model_simple.jl       Modulo Model:        parametros, EOM no lineales, lazo cerrado (simple)
│   ├── model_double.jl       Modulo ModelDouble:  parametros, EOM no lineales, lazo cerrado (doble)
│   ├── linearization.jl      Modulo Linearization: A,B,C,D, espectro, Kalman (generico en n)
│   ├── controller.jl         Modulo Controller:   LQR (Riccati/Hamiltoniano) y Ackermann (genericos)
│   ├── animation_simple.jl   Modulo Animation:       animacion GLMakie del simple
│   └── animation_double.jl   Modulo AnimationDouble: animacion GLMakie del doble
├── notebooks/
│   ├── 01_exploracion_simple.jl   Explorador interactivo Pluto (R^4)
│   └── 02_exploracion_doble.jl    Explorador interactivo Pluto (R^6)
└── docs/
    ├── guia_maestra.md            Este documento
    ├── resumen_ejecutivo/         Resumen ejecutivo (LaTeX + PDF), max 5 paginas
    ├── resumen_tecnico/           Informe tecnico (LaTeX + PDF) + make_report_figs.jl + figs/
    └── presentacion/              Diapositivas Beamer 16:9 + make_slide_figs.jl + figs/
```

### 3.2 El principio de diseño: física separada de álgebra

La decisión arquitectónica más importante es **separar la física (distinta en
cada configuración) de los algoritmos de álgebra lineal (idénticos)**:

- La **física** vive en `model_simple.jl` y `model_double.jl`. Cada uno tiene sus
  propios parámetros, sus ecuaciones de movimiento no lineales y su función de
  lazo cerrado. Son módulos distintos porque la física es distinta.
- Los **algoritmos** de `linearization.jl` (excepto las dos funciones de
  linealización, que sí son específicas) y **todo** `controller.jl` son
  **genéricos en la dimensión del estado $n$**. La estructura `StateSpaceModel`
  es la interfaz común: una vez construida, el mismo código de controlabilidad,
  observabilidad, LQR y Ackermann opera igual sobre el sistema $4\times4$ del
  simple y el $6\times6$ del doble.

En consecuencia, `design_lqr`, `solve_care`, `design_pole_placement`,
`controllability_matrix`, `observability_matrix` y `print_analysis` se escribieron
**una sola vez** y el doble los reutiliza sin cambios. Este es un argumento fuerte
para la presentación: el álgebra lineal es indiferente a que el estado viva en
$\mathbb{R}^4$ o en $\mathbb{R}^6$.

### 3.3 Flujo de datos de cada pipeline

Ambos `main_*.jl` reproducen el mismo orden lógico, que es también el orden del
informe técnico:

```
model_*.jl          linearization.jl              controller.jl            simulacion no lineal
(fisica f(x,u))  →  (A,B,C,D + espectro + Kalman) → (ganancia K por LQR/Acker) → (verificacion en lazo cerrado) → figuras
```

Concretamente, `main_simple.jl` ejecuta ocho pasos: (1) parámetros, (2)
simulación libre sin control, (3) linealización y análisis, (4) diseño LQR, (5)
diseño Ackermann, (6) simulación en lazo cerrado, (7) gráficas comparativas, (8)
animación. `main_double.jl` sigue siete pasos análogos (sin el paso Ackermann
separado, porque para $n=6$ elegir polos a mano no es evidente y se prefiere el
LQR).

### 3.4 Dependencias (Project.toml)

| Paquete | Uso en el proyecto |
|---|---|
| `DifferentialEquations` | Integración numérica de las EOM no lineales (solver `Tsit5`, Runge-Kutta explícito de orden 5) |
| `LinearAlgebra` | `eigen`, `rank`, `inv`, `cholesky`, productos matriciales (el núcleo del proyecto) |
| `CairoMakie` | Gráficas estáticas para informes y notebooks (salida en CPU) |
| `GLMakie` | Animaciones interactivas del carro-péndulo (los `main_*.jl`) |
| `ControlSystems`, `MatrixEquations` | Apoyo y verificación (Riccati, Lyapunov) |
| `Symbolics` | Derivaciones simbólicas de apoyo |
| `Pluto`, `PlutoUI` | Notebooks interactivos con sliders |
| `Printf` | Salida formateada por consola |

Nota clave: las rutinas de álgebra lineal centrales (linealización, Riccati,
Ackermann) **se programaron a mano desde las definiciones**, no se delegaron a
librerías. Las librerías especializadas se usan como verificación. Esto es
deliberado: hace explícito el papel de cada concepto del curso.

---

## 4. Fundamento físico: modelado por Euler-Lagrange

Esta sección corresponde a la Sección 2 del informe técnico y a los frames
"Coordenadas generalizadas", "El principio variacional", "Ecuaciones de
Euler-Lagrange" y las dos configuraciones de la presentación.

### 4.1 Los espacios del problema

Antes de escribir una ecuación conviene precisar en qué conjunto vive cada objeto.

- **Coordenadas generalizadas** $q = (q_1, \dots, q_d)$: un conjunto *mínimo* de
  números independientes que determinan por completo la configuración una vez
  impuestas las restricciones (barras rígidas, carro que no se levanta del riel).
  El entero $d$ es el número de **grados de libertad**.
  - Simple: $q = (x, \theta)$, $d = 2$.
  - Doble: $q = (x, \theta_1, \theta_2)$, $d = 3$.

- **Espacio de configuración** $M$: el conjunto de configuraciones admisibles.
  Es una variedad diferenciable de dimensión $d$. La posición del carro aporta un
  factor $\mathbb{R}$; el ángulo vive en la circunferencia $S^1$ (porque $\theta$
  y $\theta + 2\pi$ son la misma configuración). Así:
  - Simple: $M = \mathbb{R} \times S^1$.
  - Doble: $M = \mathbb{R} \times S^1 \times S^1$.
  - Cerca del equilibrio (ángulos pequeños) se identifica localmente $S^1$ con
    $\mathbb{R}$, y se trata $M$ como $\mathbb{R}^2$ o $\mathbb{R}^3$.

- **Espacio de estados**: el **fibrado tangente** $\mathrm{T}M$, cuyos puntos son
  los pares posición-velocidad $(q, \dot{q})$, de dimensión $2d$. La razón de que
  el estado sea $(q, \dot{q})$ y no solo $q$: las leyes de la mecánica son
  ecuaciones de **segundo orden**, y para predecir el futuro no basta la posición,
  hace falta también la velocidad. Dos configuraciones iguales con velocidades
  distintas tienen futuros distintos.
  - Simple: estado en $\mathbb{R}^4$.
  - Doble: estado en $\mathbb{R}^6$.

- **Espacio de caminos** $\mathcal{P}$: el conjunto (de dimensión infinita) de
  trayectorias suaves con extremos fijos $q_a, q_b$. Sobre él se define el
  **funcional de acción**
  $$ S[\gamma] = \int_{t_1}^{t_2} L\big(\gamma(t), \dot\gamma(t)\big)\, dt. $$

### 4.2 El lagrangiano y el principio de Hamilton

El **lagrangiano** es la función
$$ L: \mathrm{T}M \to \mathbb{R}, \qquad L(q, \dot q) = T(q, \dot q) - V(q), $$
donde $T$ es la energía cinética (depende de posiciones y velocidades) y $V$ la
energía potencial (depende solo de posiciones).

El **principio de Hamilton** afirma que las trayectorias físicas hacen
**estacionaria** la acción ($\delta S = 0$): son los puntos críticos del
funcional $S$. Es el análogo, en dimensión infinita, de buscar dónde se anula la
derivada de una función.

Imponer $\delta S = 0$ conduce a las **ecuaciones de Euler-Lagrange**, para cada
coordenada $i = 1, \dots, d$:
$$ \frac{d}{dt}\!\left(\frac{\partial L}{\partial \dot q_i}\right) - \frac{\partial L}{\partial q_i} = Q_i, $$
donde $Q_i$ es la **fuerza generalizada no conservativa** asociada a $q_i$ (en
nuestro caso: la fuerza de control $u$ y la fricción $-b\dot x$). Cuando todas las
fuerzas derivan de un potencial, $Q_i = 0$ y se recupera la forma usual.

**El plan de trabajo** (idéntico para ambas configuraciones):

1. Elegir las coordenadas generalizadas $q$.
2. Escribir las posiciones de cada masa y derivarlas para obtener las velocidades.
3. Formar $T$ y $V$, y con ellas $L = T - V$.
4. Aplicar Euler-Lagrange para obtener las ecuaciones no lineales.
5. Linealizar y reescribir en el espacio de estados.

### 4.3 Configuración I: péndulo simple (deducción completa)

**Coordenadas y posiciones.** Con $q = (x, \theta)$ y $\theta$ desde la vertical
superior, el centro de masa de la barra (a distancia $L$ del pivote) está en
$$ \mathbf{r}_{\mathrm{cm}} = (x + L\sin\theta,\ L\cos\theta). $$

**Velocidades.** Derivando respecto del tiempo,
$$ \dot{\mathbf{r}}_{\mathrm{cm}} = (\dot x + L\dot\theta\cos\theta,\ -L\dot\theta\sin\theta), $$
cuyo módulo al cuadrado es $\|\dot{\mathbf{r}}_{\mathrm{cm}}\|^2 = \dot x^2 + 2L\dot x\dot\theta\cos\theta + L^2\dot\theta^2$.

**Energías y lagrangiano.** La energía cinética suma traslación del carro,
traslación del CM de la barra y rotación de la barra ($\tfrac12 I\dot\theta^2$):
$$ T = \tfrac12 M\dot x^2 + \tfrac12 m\big(\dot x^2 + 2L\dot x\dot\theta\cos\theta + L^2\dot\theta^2\big) + \tfrac12 I\dot\theta^2, $$
$$ V = mgL\cos\theta, \qquad L = T - V. $$

**Ecuaciones de movimiento.** Aplicando Euler-Lagrange a $q_1 = x$ con
$Q_1 = u - b\dot x$ y a $q_2 = \theta$ con $Q_2 = 0$:
$$
\begin{aligned}
(M+m)\ddot x + mL\ddot\theta\cos\theta - mL\dot\theta^2\sin\theta + b\dot x &= u, \\
(I+mL^2)\ddot\theta + mL\ddot x\cos\theta - mgL\sin\theta &= 0.
\end{aligned}
$$

En forma de matriz de masa (la que aparece en el código):
$$
\underbrace{\begin{pmatrix} M+m & mL\cos\theta \\ mL\cos\theta & I+mL^2 \end{pmatrix}}_{\mathbf{M}(q)}
\begin{pmatrix} \ddot x \\ \ddot\theta \end{pmatrix} =
\begin{pmatrix} u - b\dot x + mL\dot\theta^2\sin\theta \\ mgL\sin\theta \end{pmatrix}.
$$
Su determinante evaluado en $\theta = 0$ es
$$ D_0 = (M+m)(I+mL^2) - (mL)^2, $$
la cantidad clave que reaparece en todas las entradas de $A$ y $B$.

> **La firma de la inestabilidad.** El término $-mgL\sin\theta$ es el que vuelve
> inestable el equilibrio superior: al linealizar produce un coeficiente
> *positivo* que originará un eigenvalor real positivo. Si el péndulo colgara
> hacia abajo, ese término aparecería con signo opuesto y el equilibrio sería
> estable. **La diferencia entre "colgar" y "estar invertido" es, literalmente,
> un signo.**

### 4.4 Configuración II: péndulo doble (deducción completa)

**Coordenadas y posiciones.** Con $q = (x, \theta_1, \theta_2)$, masas puntuales
$m_1$ (articulación intermedia) y $m_2$ (extremo), barras de masa despreciable:
$$
\mathbf{r}_1 = (x + L_1\sin\theta_1,\ L_1\cos\theta_1), \qquad
\mathbf{r}_2 = \mathbf{r}_1 + (L_2\sin\theta_2,\ L_2\cos\theta_2).
$$

**Energías.**
$$ T = \tfrac12 M\dot x^2 + \tfrac12 m_1\|\dot{\mathbf{r}}_1\|^2 + \tfrac12 m_2\|\dot{\mathbf{r}}_2\|^2, $$
$$ V = (m_1+m_2)gL_1\cos\theta_1 + m_2 gL_2\cos\theta_2. $$

**Ecuaciones de movimiento.** Aplicando Euler-Lagrange a las tres coordenadas y
agrupando, se llega a la **forma matricial estándar de la mecánica**, la llamada
*ecuación del manipulador*:
$$ \mathbf{M}(q)\,\ddot q + \mathbf{C}(q, \dot q)\,\dot q + \mathbf{G}(q) = \mathbf{F}\,u, $$
con $\mathbf{F} = (1, 0, 0)^\top$, **matriz de masa**
$$
\mathbf{M}(q) = \begin{pmatrix}
M + m_1 + m_2 & L_1(m_1+m_2)c_1 & L_2 m_2 c_2 \\
L_1(m_1+m_2)c_1 & L_1^2(m_1+m_2) & L_1 L_2 m_2 c_{12} \\
L_2 m_2 c_2 & L_1 L_2 m_2 c_{12} & L_2^2 m_2
\end{pmatrix},
$$
donde $c_1 = \cos\theta_1$, $c_2 = \cos\theta_2$, $c_{12} = \cos(\theta_1 - \theta_2)$;
**vector de gravedad**
$$ \mathbf{G}(q) = \big(0,\ -(m_1+m_2)gL_1\sin\theta_1,\ -m_2 gL_2\sin\theta_2\big)^\top; $$
y $\mathbf{C}(q,\dot q)\dot q$ los términos centrífugos y de Coriolis.

**Propiedad clave para el código:** $\mathbf{M}(q)$ es **simétrica y definida
positiva** (es la Hessiana de la energía cinética). Por eso es invertible y el
sistema $\mathbf{M}(q)\ddot q = \text{rhs}$ se resuelve para $\ddot q$ mediante
la **factorización de Cholesky** $\mathbf{M} = LL^\top$ (más barata y estable que
invertir o que una LU genérica).

**Verificación por reducción.** Tomando $m_2 = 0$ en el doble (e $I = b = 0$ en el
simple), ambas derivaciones coinciden: el doble se reduce al simple. Esta
consistencia es una primera validación de las ecuaciones.

---

## 5. Linealización y espacio de estados

Corresponde a la Sección 2.4-2.6 del informe y a los frames "Linealización" y
"Espacio de estados" de la presentación. Esta sección es deliberadamente
detallada: explica **de dónde sale cada entrada** de las matrices $A$ y $B$, por
qué aparece la cantidad $D_0$, y por qué el procedimiento funciona igual para el
simple ($\mathbb{R}^4$) y el doble ($\mathbb{R}^6$). Es, en muchos sentidos, el
puente entre la física de la sección 4 y el álgebra lineal de la sección 6.

### 5.1 El objetivo: llegar a la forma de espacio de estados

Toda la teoría de control lineal opera sobre la estructura
$$ \dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u}, \qquad \mathbf{y} = C\mathbf{x} + D\mathbf{u}. $$

Cada matriz tiene un nombre y un papel:

| Matriz | Nombre | Papel |
|---|---|---|
| $A \in \mathbb{R}^{n\times n}$ | Matriz de dinámica | La transformación lineal que rige la evolución de las desviaciones |
| $B \in \mathbb{R}^{n\times m}$ | Matriz de entrada | Cómo entra la fuerza de control ($m = 1$: una sola columna) |
| $C \in \mathbb{R}^{p\times n}$ | Matriz de salida | Qué se mide con los sensores |
| $D \in \mathbb{R}^{p\times m}$ | Transmisión directa | Efecto instantáneo de la entrada sobre la salida; aquí $D = \mathbf{0}$ |

El problema es que las ecuaciones de movimiento de la sección 4 **no** tienen esa
forma: contienen $\sin\theta$, $\cos\theta$ y productos de velocidades. Son de la
forma general $\dot{\mathbf{x}} = f(\mathbf{x}, u)$ con $f$ **no lineal**. El
objetivo de esta sección es fabricar $A$ y $B$ a partir de $f$.

### 5.2 Qué significa linealizar

Linealizar es reemplazar el campo vectorial no lineal $f(\mathbf{x}, u)$ por su
**mejor aproximación lineal** en torno a un punto de operación. La justificación
es doble:

- **Física:** el objetivo de control es *mantener* el péndulo cerca del
  equilibrio erguido. Si el sistema permanece cerca de ese punto, las desviaciones
  son pequeñas y una descripción de primer orden basta.
- **Matemática:** una función suave se parece a su plano tangente cerca del punto;
  el error del truncamiento es de segundo orden, es decir, del tamaño del
  *cuadrado* de la desviación.

El primer paso, entonces, es escribir las EOM como un sistema de primer orden
$\dot{\mathbf{x}} = f(\mathbf{x}, u)$. Para el simple, con
$\mathbf{x} = (x, \dot x, \theta, \dot\theta)$, esto significa despejar las
aceleraciones de la sección 4.3:
$$
f(\mathbf{x}, u) = \begin{pmatrix}
\dot x \\[2pt]
\ddot x(\mathbf{x}, u) \\[2pt]
\dot\theta \\[2pt]
\ddot\theta(\mathbf{x}, u)
\end{pmatrix},
\qquad
\begin{pmatrix} \ddot x \\ \ddot\theta \end{pmatrix}
= \mathbf{M}(q)^{-1}\,\mathbf{h}(\mathbf{x}, u),
$$
donde, recordando la forma de matriz de masa de la sección 4.3,
$$
\mathbf{M}(q) = \begin{pmatrix} M+m & mL\cos\theta \\ mL\cos\theta & I+mL^2 \end{pmatrix},
\qquad
\mathbf{h}(\mathbf{x}, u) = \begin{pmatrix} u - b\dot x + mL\dot\theta^2\sin\theta \\ mgL\sin\theta \end{pmatrix}.
$$
Aquí $\mathbf{h}$ agrupa **todo lo que no es la aceleración**: la fuerza de
control, la fricción, el término centrífugo y la gravedad. Es el lado derecho
(*right-hand side*).

### 5.3 El punto de equilibrio

Un **equilibrio** es un estado $\mathbf{x}^\star$ en el que, sin control, el
sistema no cambia: $f(\mathbf{x}^\star, 0) = \mathbf{0}$. Físicamente, el péndulo
erguido y todo en reposo:
$$ \mathbf{x}^\star = (x^\star,\ 0,\ 0,\ 0), \qquad u^\star = 0. $$
Comprobación directa: en ese estado $\dot x = 0$ y $\dot\theta = 0$ (las filas 1 y
3 de $f$ se anulan), y $\mathbf{h}(\mathbf{x}^\star, 0) = \mathbf{0}$ porque
$u = 0$, $\dot x = 0$, $\dot\theta = 0$ y $\sin 0 = 0$ (las filas 2 y 4 se anulan
al multiplicar por $\mathbf{M}^{-1}$). **Que $\mathbf{h}$ se anule en el
equilibrio será la clave de todo lo que sigue** (subsección 5.6).

Observación: la posición del carro $x^\star$ es **libre** (cualquier posición es
un equilibrio; el sistema es invariante ante traslaciones del carro). Por eso $x$
no aparecerá en ninguna otra ecuación y la **primera columna de $A$ es
idénticamente cero**. Tomamos $x^\star = 0$ por comodidad.

### 5.4 El desarrollo de Taylor y los jacobianos $A$, $B$

El desarrollo de Taylor de primer orden de $f$ alrededor de
$(\mathbf{x}^\star, 0)$ es
$$
f(\mathbf{x}, u) \approx
\underbrace{f(\mathbf{x}^\star, 0)}_{=\,\mathbf{0}}
+ \underbrace{\left.\frac{\partial f}{\partial \mathbf{x}}\right|_{(\mathbf{x}^\star,0)}}_{=:\,A}
(\mathbf{x} - \mathbf{x}^\star)
+ \underbrace{\left.\frac{\partial f}{\partial u}\right|_{(\mathbf{x}^\star,0)}}_{=:\,B}
(u - 0).
$$
Definiendo la desviación $\mathbf{x} \leftarrow \mathbf{x} - \mathbf{x}^\star$ (y
conservando el símbolo para aligerar la notación), el sistema queda, a primer
orden, en la forma buscada $\dot{\mathbf{x}} = A\mathbf{x} + Bu$. Las matrices son
**jacobianos**: $A$ es la matriz de todas las derivadas parciales de cada
componente de $f$ respecto de cada componente del estado; $B$, respecto de la
entrada.

Para las funciones trigonométricas, linealizar equivale a las aproximaciones de
ángulo pequeño, que son exactamente el truncamiento de la serie de Taylor:
$$ \sin\theta = \theta - \tfrac{\theta^3}{6} + \cdots \approx \theta, \qquad
   \cos\theta = 1 - \tfrac{\theta^2}{2} + \cdots \approx 1. $$
Además, todo producto de dos cantidades pequeñas es de segundo orden y se
descarta. Por ejemplo, si $\theta \sim 0.1$, entonces $\theta^2 \sim 0.01$ es un
orden de magnitud menor; y términos como $\dot\theta^2\theta$ o $\dot x\dot\theta$
son productos de pequeños por pequeños, despreciables a primer orden.

> **Consistencia local:** el modelo lineal es fiel solo mientras los ángulos sean
> pequeños. Pero ese es justamente el régimen en que opera un controlador
> estabilizante: una vez erguido, el péndulo permanece cerca de $\theta = 0$. (Y
> los resultados muestran que funciona incluso bastante lejos, hasta $\approx 23°$:
> ver sección 9.5.)

### 5.5 La estructura de bloques: de dónde salen los "unos"

Antes de calcular nada, la sola forma de $f$ dicta la mitad de $A$. El estado se
ordena como (posición, velocidad) para cada grado de libertad, y $f$ dice que la
derivada de cada posición **es** su velocidad:
$$ \frac{d}{dt}x = \dot x, \qquad \frac{d}{dt}\theta = \dot\theta. $$
Estas son ecuaciones ya lineales y exactas (no hay nada que aproximar). Su
jacobiano produce los **unos** en las posiciones $A_{1,2}$ y $A_{3,4}$ del simple
(y $A_{1,2}, A_{3,4}, A_{5,6}$ del doble). Así, $A$ tiene una estructura de
bloques en la que las **filas impares** (las de posición) son triviales y toda la
física no trivial vive en las **filas pares** (las de aceleración). Este patrón
—cada par (posición, velocidad) acoplado por un doble integrador
$\left(\begin{smallmatrix} 0 & 1 \\ \ast & 0 \end{smallmatrix}\right)$— es
genérico de los sistemas mecánicos y reaparece en el análisis de Jordan de la
sección 6.3.

Queda entonces por calcular solo el jacobiano de las aceleraciones,
$\partial \ddot q / \partial \mathbf{x}$ y $\partial \ddot q / \partial u$.

### 5.6 El truco decisivo: la matriz de masa se evalúa en el equilibrio

Las aceleraciones son $\ddot q = \mathbf{M}(q)^{-1}\mathbf{h}(\mathbf{x}, u)$. A
primera vista, derivar esto respecto del estado parece complicado, porque **tanto**
$\mathbf{M}^{-1}$ **como** $\mathbf{h}$ dependen del estado. Por la regla del
producto, la derivada respecto de una componente $x_j$ es
$$
\frac{\partial \ddot q}{\partial x_j}
= \underbrace{\frac{\partial \mathbf{M}^{-1}}{\partial x_j}\,\mathbf{h}}_{\text{(I)}}
+ \underbrace{\mathbf{M}^{-1}\,\frac{\partial \mathbf{h}}{\partial x_j}}_{\text{(II)}}.
$$
Aquí está el punto elegante: **evaluado en el equilibrio, el término (I) se
anula**, porque $\mathbf{h}(\mathbf{x}^\star, 0) = \mathbf{0}$ (subsección 5.3).
La variación de $\mathbf{M}^{-1}$ multiplica a algo que vale cero, así que **no
contribuye al jacobiano**. Solo sobrevive el término (II), en el que
$\mathbf{M}^{-1}$ aparece **evaluada en el equilibrio**, es decir, como una matriz
**constante**:
$$
\mathbf{M}_0 := \mathbf{M}(q)\big|_{\theta=0}
= \begin{pmatrix} M+m & mL \\ mL & I+mL^2 \end{pmatrix}
\qquad (\text{usando } \cos 0 = 1).
$$

En consecuencia, el bloque no trivial de $A$ y todo $B$ se obtienen con la receta
compacta
$$
\boxed{\;
A_{\text{acel}} = \mathbf{M}_0^{-1}\,
\left.\frac{\partial \mathbf{h}}{\partial \mathbf{x}}\right|_{\text{eq}},
\qquad
B_{\text{acel}} = \mathbf{M}_0^{-1}\,
\left.\frac{\partial \mathbf{h}}{\partial u}\right|_{\text{eq}}.
\;}
$$
Este es el motivo por el que en el informe se linealiza "el sistema $2\times2$"
evaluando la matriz de masa en $\theta = 0$: no es una aproximación adicional,
es una **consecuencia exacta** de que el lado derecho se anula en el equilibrio.

### 5.7 De dónde sale $D_0$

$D_0$ es, simplemente, el **determinante de la matriz de masa evaluada en el
equilibrio**:
$$
D_0 := \det \mathbf{M}_0
= (M+m)(I+mL^2) - (mL)^2.
$$
Aparece en todas las entradas de $A$ y $B$ porque invertir una matriz $2\times2$
introduce su determinante en el denominador:
$$
\mathbf{M}_0^{-1} = \frac{1}{D_0}
\begin{pmatrix} I+mL^2 & -mL \\ -mL & M+m \end{pmatrix}.
$$
Dos hechos importantes sobre $D_0$:

- **Es estrictamente positivo.** $\mathbf{M}_0$ es simétrica y **definida
  positiva** (es la Hessiana de la energía cinética, que es una forma cuadrática
  positiva en las velocidades). Por el criterio de los menores principales, su
  determinante $D_0 > 0$. Esto garantiza que $\mathbf{M}_0$ es invertible y que
  las aceleraciones están bien definidas.
- **Valor numérico** (parámetros de la sección 9.1): $M+m = 1.3$,
  $I+mL^2 = 0.025 + 0.3(0.5)^2 = 0.1$, $(mL)^2 = (0.15)^2 = 0.0225$, de modo que
  $$ D_0 = 1.3 \times 0.1 - 0.0225 = 0.1075. $$

### 5.8 Derivación completa del péndulo simple, paso a paso

Ahora ensamblamos $A$ y $B$ ejecutando la receta de 5.6.

**Paso 1: jacobiano del lado derecho.** Derivamos
$\mathbf{h} = (u - b\dot x + mL\dot\theta^2\sin\theta,\ mgL\sin\theta)^\top$
respecto de $(x, \dot x, \theta, \dot\theta)$ y evaluamos en el equilibrio
(donde $\theta = 0$, $\dot\theta = 0$):

- $\partial h_1/\partial \dot x = -b$.
- $\partial h_1/\partial \theta = mL\dot\theta^2\cos\theta \to 0$ (porque
  $\dot\theta^2 = 0$: el término centrífugo, cuadrático en velocidad, **no
  contribuye**).
- $\partial h_1/\partial \dot\theta = 2mL\dot\theta\sin\theta \to 0$.
- $\partial h_2/\partial \theta = mgL\cos\theta \to mgL$ (aquí está la gravedad).
- El resto de derivadas son cero.

Es decir,
$$
\left.\frac{\partial \mathbf{h}}{\partial \mathbf{x}}\right|_{\text{eq}}
= \begin{pmatrix} 0 & -b & 0 & 0 \\ 0 & 0 & mgL & 0 \end{pmatrix},
\qquad
\left.\frac{\partial \mathbf{h}}{\partial u}\right|_{\text{eq}}
= \begin{pmatrix} 1 \\ 0 \end{pmatrix}.
$$

**Paso 2: multiplicar por $\mathbf{M}_0^{-1}$.** Para el bloque de aceleraciones
de $A$,
$$
A_{\text{acel}}
= \frac{1}{D_0}
\begin{pmatrix} I+mL^2 & -mL \\ -mL & M+m \end{pmatrix}
\begin{pmatrix} 0 & -b & 0 & 0 \\ 0 & 0 & mgL & 0 \end{pmatrix}
= \frac{1}{D_0}
\begin{pmatrix}
0 & -b(I+mL^2) & -m^2gL^2 & 0 \\
0 & bmL & (M+m)mgL & 0
\end{pmatrix}.
$$
La primera fila es $\ddot x$ (fila 2 de $A$) y la segunda es $\ddot\theta$ (fila 4
de $A$). Para $B$,
$$
B_{\text{acel}} = \frac{1}{D_0}
\begin{pmatrix} I+mL^2 & -mL \\ -mL & M+m \end{pmatrix}
\begin{pmatrix} 1 \\ 0 \end{pmatrix}
= \frac{1}{D_0}\begin{pmatrix} I+mL^2 \\ -mL \end{pmatrix}.
$$

**Paso 3: insertar en la estructura de bloques** (los unos de 5.5). El resultado
es exactamente lo que aparece en `linearize_system`:
$$
A = \begin{pmatrix}
0 & 1 & 0 & 0 \\
0 & -\dfrac{b(I+mL^2)}{D_0} & -\dfrac{m^2gL^2}{D_0} & 0 \\
0 & 0 & 0 & 1 \\
0 & \dfrac{bmL}{D_0} & \dfrac{(M+m)mgL}{D_0} & 0
\end{pmatrix},
\quad
B = \begin{pmatrix} 0 \\ \dfrac{I+mL^2}{D_0} \\ 0 \\ -\dfrac{mL}{D_0} \end{pmatrix}.
$$

**Verificación numérica** (con $D_0 = 0.1075$): $A_{2,2} = -0.093$,
$A_{2,3} = -2.053$, $A_{4,2} = 0.140$, $A_{4,3} = 17.795$, $B_2 = 0.930$,
$B_4 = -1.395$. Coinciden con las matrices del Apéndice B del informe.

### 5.9 La firma de la inestabilidad (el signo)

La entrada crítica es $A_{4,3} = +(M+m)mgL/D_0 > 0$. Es la traducción algebraica
de la observación física de la sección 4.3: **el signo positivo del coeficiente
gravitatorio es lo que vuelve inestable el equilibrio superior**. Aislando el
bloque angular (sección 6.2), este coeficiente positivo produce eigenvalores
$\pm\sqrt{\kappa}$ con $\kappa > 0$, uno de ellos real y positivo. Si el péndulo
colgara, el término gravitatorio cambiaría de signo, $A_{4,3}$ sería negativo, y
los eigenvalores serían imaginarios puros (oscilación estable). Al implementar, es
un punto fácil de equivocar; se verifica comparando el espectro obtenido con el
esperado (un eigenvalor real positivo).

### 5.10 Derivación del péndulo doble (el mismo método en $\mathbb{R}^6$)

El doble se linealiza **exactamente igual**, solo que ahora la matriz de masa es
$3\times3$ y el equilibrio es $\theta_1 = \theta_2 = 0$. Partiendo de la ecuación
del manipulador de la sección 4.4,
$$
\mathbf{M}(q)\ddot q + \mathbf{C}(q,\dot q)\dot q + \mathbf{G}(q) = \mathbf{F}u
\quad\Longrightarrow\quad
\ddot q = \mathbf{M}(q)^{-1}\big[\mathbf{F}u - \mathbf{C}(q,\dot q)\dot q - \mathbf{G}(q)\big].
$$
De nuevo el lado derecho se anula en el equilibrio ($u = 0$; $\mathbf{G}(0) = 0$
porque $\sin 0 = 0$; y $\mathbf{C}\dot q = 0$ porque $\dot q = 0$), así que vale
el mismo truco de 5.6: se evalúa $\mathbf{M}$ en el equilibrio y solo se linealiza
el lado derecho.

**El término de Coriolis desaparece.** Un detalle que conviene resaltar:
$\mathbf{C}(q,\dot q)\dot q$ es **cuadrático en las velocidades** (contiene
$\dot\theta_i^2$ y productos $\dot\theta_i\dot\theta_j$). Su jacobiano respecto de
$\dot q$ trae siempre un factor $\dot q$, que se anula en el equilibrio; y su
jacobiano respecto de $q$ también trae factores $\dot q$. Por tanto, **Coriolis no
contribuye en absoluto a la linealización**. Solo sobreviven la gravedad (vía
$-\partial \mathbf{G}/\partial q$) y la fuerza de control (vía $\mathbf{F}$).

**Los ingredientes.** La matriz de masa en el equilibrio ($c_1 = c_2 = c_{12} = 1$):
$$
\mathbf{M}_0 = \begin{pmatrix}
M+m_1+m_2 & L_1(m_1+m_2) & L_2 m_2 \\
L_1(m_1+m_2) & L_1^2(m_1+m_2) & L_1 L_2 m_2 \\
L_2 m_2 & L_1 L_2 m_2 & L_2^2 m_2
\end{pmatrix},
$$
y el jacobiano de la gravedad, con $\mathbf{G} = (0,\ -(m_1+m_2)gL_1\sin\theta_1,\ -m_2 gL_2\sin\theta_2)^\top$:
$$
\left.-\frac{\partial \mathbf{G}}{\partial q}\right|_{\text{eq}}
= \begin{pmatrix}
0 & 0 & 0 \\
0 & (m_1+m_2)gL_1 & 0 \\
0 & 0 & m_2 gL_2
\end{pmatrix}
\qquad (\text{usando } \cos 0 = 1).
$$

**El resultado.** Las columnas de $\theta_1$ y $\theta_2$ del bloque de
aceleraciones de $A$, y la columna $B$, son
$$
A_{\text{acel}} = \mathbf{M}_0^{-1}\left(-\frac{\partial \mathbf{G}}{\partial q}\right),
\qquad
B_{\text{acel}} = \mathbf{M}_0^{-1}\,\mathbf{F}, \quad \mathbf{F} = (1,0,0)^\top.
$$
Es decir, **$B$ es literalmente la primera columna de $\mathbf{M}_0^{-1}$**. Vale
la pena verificarla, porque sale sorprendentemente limpia. El código da
$B_{\text{acel}} = (1/M,\ -1/(ML_1),\ 0)^\top$; comprobando que
$\mathbf{M}_0\,(1/M, -1/(ML_1), 0)^\top = (1, 0, 0)^\top$:
$$
\begin{aligned}
\text{fila 1:}\ & \tfrac{M+m_1+m_2}{M} - \tfrac{L_1(m_1+m_2)}{ML_1} = \tfrac{M+m_1+m_2-(m_1+m_2)}{M} = 1, \\
\text{fila 2:}\ & \tfrac{L_1(m_1+m_2)}{M} - \tfrac{L_1^2(m_1+m_2)}{ML_1} = 0, \\
\text{fila 3:}\ & \tfrac{L_2 m_2}{M} - \tfrac{L_1 L_2 m_2}{ML_1} = 0. \quad\checkmark
\end{aligned}
$$
Lectura física de $B$: la fuerza sobre el carro lo acelera ($1/M$), induce una
rotación en el eslabón inferior ($-1/(ML_1)$), pero **no toca directamente el
eslabón superior a primer orden** ($B_6 = 0$); a este solo llega el control de
forma indirecta, a través del acoplamiento en $A$.

El bloque completo de $A$ (haciendo el mismo producto para las columnas de
gravedad) es el que está escrito entrada por entrada en `linearize_system_double`:
$$
A = \begin{pmatrix}
0 & 1 & 0 & 0 & 0 & 0 \\
0 & 0 & -\dfrac{g(m_1+m_2)}{M} & 0 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 & 0 \\
0 & 0 & \dfrac{g(M+m_1)(m_1+m_2)}{ML_1 m_1} & 0 & -\dfrac{g m_2}{L_1 m_1} & 0 \\
0 & 0 & 0 & 0 & 0 & 1 \\
0 & 0 & -\dfrac{g(m_1+m_2)}{L_2 m_1} & 0 & \dfrac{g(m_1+m_2)}{L_2 m_1} & 0
\end{pmatrix},
\quad
B = \Big(0,\ \tfrac1M,\ 0,\ -\tfrac{1}{ML_1},\ 0,\ 0\Big)^\top.
$$

**¿Y el "$D_0$" del doble?** No aparece un único denominador común como en el
simple. La razón es que $\mathbf{M}_0^{-1} = \frac{1}{\det\mathbf{M}_0}\operatorname{adj}(\mathbf{M}_0)$:
al multiplicar por la matriz de gravedad (que es diagonal y comparte factores con
la adjunta), parte del $\det\mathbf{M}_0$ se **cancela** con los cofactores,
dejando denominadores más simples ($M$, $ML_1 m_1$, $L_2 m_1$) en lugar de un
$\det\mathbf{M}_0$ único. El análogo conceptual de $D_0$ sigue siendo
$\det\mathbf{M}_0 > 0$ (por la misma razón de definición positiva); solo que en el
doble ya no queda visible tras la simplificación. Los coeficientes gravitatorios
**positivos** $A_{4,3}$ y $A_{6,5}$ son, otra vez, la firma de la inestabilidad,
ahora de **dos** modos.

### 5.11 Las matrices de salida $C$ y $D$

$C$ y $D$ no salen de linealizar la dinámica, sino de **decidir qué se mide**. En
este proyecto se miden posiciones (del carro y los ángulos), no velocidades,
porque son las magnitudes accesibles con sensores directos. Para el simple:
$$
C = \begin{pmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \end{pmatrix}, \qquad D = \mathbf{0}
\qquad (p = 2:\ \text{mide } x \text{ y } \theta).
$$
Para el doble:
$$
C = \begin{pmatrix} 1&0&0&0&0&0 \\ 0&0&1&0&0&0 \\ 0&0&0&0&1&0 \end{pmatrix}, \qquad D = \mathbf{0}
\qquad (p = 3:\ \text{mide } x,\ \theta_1,\ \theta_2).
$$
$D = \mathbf{0}$ porque la fuerza de control no aparece **instantáneamente** en las
mediciones: afecta a las posiciones solo de forma diferida, a través de la
dinámica. Que las velocidades no se midan pero puedan **reconstruirse** es
precisamente lo que garantiza la observabilidad (sección 6.6).

### 5.12 Verificaciones cruzadas

El proyecto valida estas derivaciones de tres formas independientes:

1. **Reducción de configuraciones.** Tomando $m_2 = 0$ en el doble (e $I = b = 0$
   en el simple), ambas linealizaciones coinciden: el doble se reduce al simple.
2. **Analítico contra numérico.** En `linearization.jl` las matrices se escriben
   entrada por entrada con la forma analítica del jacobiano; el informe reporta
   que coinciden con el jacobiano **numérico** de las EOM no lineales hasta
   $\sim 10^{-11}$.
3. **Espectro esperado.** El signo del término gravitatorio se verifica
   comprobando que el espectro tiene el número correcto de eigenvalores reales
   positivos (uno en el simple, dos en el doble; sección 9.2).

---

## 6. Las herramientas de álgebra lineal

Este es el corazón teórico del proyecto (Sección 3 del informe; sección
"Herramientas de álgebra lineal" de la presentación). Cada herramienta responde
una pregunta concreta sobre el sistema.

### 6.1 La exponencial matricial: el operador solución

En una dimensión, $\dot x = ax$ tiene solución $x(t) = e^{at}x(0)$. La
generalización a $n$ dimensiones reemplaza el número $a$ por la matriz $A$:
$$ e^{At} = \sum_{k=0}^{\infty} \frac{(At)^k}{k!}, $$
serie convergente para todo $t$ y toda $A$. La solución del sistema es
$$ \mathbf{x}(t) = e^{At}\mathbf{x}_0 + \int_0^t e^{A(t-s)}B\,\mathbf{u}(s)\,ds. $$

$e^{At}$ es el **operador solución** (o matriz de transición de estados): aplicado
al estado inicial, devuelve el estado en el instante $t$. Es literalmente el flujo
del sistema, y toda la información dinámica del sistema libre está en él.

**El puente con los eigenvalores.** Si $A$ es diagonalizable, $A = P\Lambda P^{-1}$,
entonces
$$ e^{At} = P\,\operatorname{diag}(e^{\lambda_1 t}, \dots, e^{\lambda_n t})\,P^{-1}, $$
y el cálculo se reduce a exponenciales escalares de los eigenvalores. De aquí la
conclusión central:

> $\mathbf{x}(t) \to \mathbf{0}$ cuando $t \to \infty$ **si y solo si**
> $\operatorname{Re}(\lambda_i) < 0$ para todo eigenvalor $\lambda_i$ de $A$. El
> signo de la parte real de los eigenvalores decide la estabilidad.

### 6.2 Eigenvalores, polinomio característico y Cayley-Hamilton

Un escalar $\lambda$ es eigenvalor de $A$ si $(A - \lambda I)\mathbf{v} = \mathbf{0}$
con $\mathbf{v} \neq \mathbf{0}$, lo que ocurre exactamente cuando
$$ p(\lambda) = \det(A - \lambda I) = 0. $$
El **polinomio característico** $p(\lambda)$ tiene grado $n$; sus raíces son los
eigenvalores, y su conjunto es el **espectro** $\sigma(A)$.

**Ejemplo revelador (bloque angular del simple, sin fricción).** Aislando el
sub-bloque $(\theta, \dot\theta)$ con $\ddot\theta = \kappa\,\theta$,
$\kappa = (M+m)mgL/D_0 > 0$:
$$ A_\theta = \begin{pmatrix} 0 & 1 \\ \kappa & 0 \end{pmatrix}, \qquad
p(\lambda) = \lambda^2 - \kappa \ \Longrightarrow\ \lambda = \pm\sqrt{\kappa}. $$
Como $\kappa > 0$, uno de los eigenvalores es **real y positivo**: el equilibrio
es inestable. Con los parámetros del informe, $\sqrt{\kappa} \approx 4.22$, muy
próximo al eigenvalor exacto $+4.21$ del sistema completo (la pequeña diferencia
la introducen la fricción y el acoplamiento). **Este ejemplo es oro para la
presentación:** en dos líneas muestra de dónde sale la inestabilidad.

**Teorema de Cayley-Hamilton.** Toda matriz satisface su propio polinomio
característico: $p(A) = 0$. La consecuencia práctica es que
$$ A^n = -\big(c_{n-1}A^{n-1} + \cdots + c_1 A + c_0 I\big), $$
es decir, $A^n$ (y toda potencia superior) es combinación lineal de
$I, A, \dots, A^{n-1}$. **Las potencias "altas" de $A$ no aportan nada nuevo.**
Esta es la razón de fondo de que las matrices de controlabilidad y observabilidad
solo necesiten potencias hasta $A^{n-1}$.

### 6.3 Polinomio minimal, descomposición primaria y forma de Jordan

La fórmula limpia $e^{At} = P\,\operatorname{diag}(\dots)\,P^{-1}$ supone que $A$
es diagonalizable. ¿Cuándo lo es?

- **Polinomio minimal** $m(\lambda)$: el polinomio mónico de grado mínimo con
  $m(A) = 0$. Divide a todo polinomio que anule a $A$ (en particular al
  característico) y sus raíces son exactamente los eigenvalores.
- **Descomposición primaria:** si $m(\lambda) = \prod_i p_i(\lambda)$ con los
  $p_i$ coprimos, el espacio se descompone en suma directa de subespacios
  $A$-invariantes.
- **Corolario (diagonalización):** sobre $\mathbb{C}$,
  $m(\lambda) = \prod_i (\lambda - \alpha_i)^{d_i}$. Si todos los $d_i = 1$
  (raíces simples), $A$ es **diagonalizable**. Los factores con $d_i > 1$ dan los
  bloques de la **forma canónica de Jordan**.

**Por qué esto no es decorativo: el doble ejercita Jordan de verdad.** Un punto
sofisticado y muy defendible en la presentación:

- Las matrices de **lazo cerrado** $A - BK$ de ambas configuraciones, y la $A$ de
  **lazo abierto del simple**, tienen eigenvalores distintos: son diagonalizables
  y la forma de Jordan se reduce al caso diagonal.
- En cambio, la $A$ de lazo abierto del **doble** tiene un eigenvalor
  $\lambda = 0$ de multiplicidad algebraica **dos** pero geométrica **uno**: es
  **defectuosa**. El bloque responsable es el **doble integrador del carro**,
  $\left(\begin{smallmatrix} 0 & 1 \\ 0 & 0 \end{smallmatrix}\right)$, un bloque
  de Jordan $2\times2$ genuino. Su exponencial arrastra un término secular
  $t\,e^{0\cdot t} = t$, que es físicamente la **deriva del carro a velocidad
  constante** ante una perturbación de velocidad: $x(t) = x_0 + \dot x_0\, t$.
- La fricción del simple ($b \neq 0$) parte ese doble integrador en dos
  eigenvalores distintos $\{0, -0.077\}$, y por eso su $A$ **sí** es
  diagonalizable. El doble, sin fricción, conserva el bloque de Jordan.

Conclusión: la maquinaria de Jordan no es una red de seguridad; el doble la
**usa**, y es la que **justifica** que la diagonalización empleada en los demás
casos sea lícita.

### 6.4 Estabilidad: Hurwitz y Lyapunov

- **Matriz de Hurwitz:** $A$ es de Hurwitz (o estable) si todos sus eigenvalores
  tienen parte real estrictamente negativa. El **criterio espectral** dice que el
  origen es asintóticamente estable si y solo si $A$ es de Hurwitz.
- **Lyapunov para sistemas lineales:** $A$ es de Hurwitz si y solo si, para
  cualquier $Q \succ 0$ simétrica, la ecuación de Lyapunov
  $$ A^\top P + PA = -Q $$
  tiene una única solución simétrica $P \succ 0$. Esta caracterización algebraica
  (que no requiere calcular eigenvalores) reaparece en el diseño óptimo del LQR.

### 6.5 Controlabilidad: el criterio de rango de Kalman

Para estabilizar un sistema inestable necesitamos poder **influir**, mediante la
entrada, en todos sus modos.

- **Definición:** el par $(A, B)$ es completamente controlable si para
  cualesquiera estados $\mathbf{x}_0, \mathbf{x}_f$ existe un tiempo finito y una
  entrada que lleva el sistema de uno a otro.
- **Matriz de controlabilidad:**
  $$ \mathcal{C} = \begin{pmatrix} B & AB & A^2B & \cdots & A^{n-1}B \end{pmatrix}. $$
  (Se trunca en $A^{n-1}B$ por Cayley-Hamilton: agregar más columnas no aumenta
  el espacio generado.)
- **Criterio de rango de Kalman:** $(A, B)$ es completamente controlable si y solo
  si $\operatorname{rank}(\mathcal{C}) = n$.
- **Test alternativo PBH (Popov-Belevitch-Hautus):** equivalentemente,
  $\operatorname{rank}[A - \lambda I \mid B] = n$ para todo $\lambda$. Falla en un
  eigenvalor $\lambda_i$ justo cuando existe un eigenvector izquierdo ortogonal a
  las columnas de $B$: un modo que la entrada "no toca". Esta lectura modo a modo
  es útil para diagnosticar qué grado de libertad quedaría sin control.

La importancia para el diseño es directa: la controlabilidad es **exactamente** la
condición que permite reubicar a voluntad los eigenvalores mediante
realimentación de estado.

### 6.6 Observabilidad y dualidad

En la práctica no se miden todos los estados, solo la salida $\mathbf{y} = C\mathbf{x}$
(aquí: posición y ángulos, **no** velocidades). La observabilidad pregunta si esas
mediciones bastan para reconstruir el estado completo.

- **Matriz de observabilidad:**
  $$ \mathcal{O} = \begin{pmatrix} C \\ CA \\ CA^2 \\ \vdots \\ CA^{n-1} \end{pmatrix}. $$
- **Criterio:** $(C, A)$ es completamente observable si y solo si
  $\operatorname{rank}(\mathcal{O}) = n$.
- **Dualidad:** $(A, B)$ es controlable si y solo si $(B^\top, A^\top)$ es
  observable, y viceversa. Refleja la simetría algebraica entre $\mathcal{C}$ y
  $\mathcal{O}^\top$: todo teorema sobre controlabilidad se traduce
  automáticamente en uno sobre observabilidad.

Verificar $\operatorname{rank}(\mathcal{O}) = n$ confirma que, midiendo solo
posiciones y ángulos, se pueden reconstruir también las velocidades, lo que
legitima el uso de la realimentación de estado $\mathbf{u} = -K\mathbf{x}$.

---

## 7. Diseño de controladores

Corresponde a la Sección 3.7-3.8 del informe y a la sección "Diseño de
controladores" de la presentación.

### 7.1 Realimentación de estado y asignación de polos (Ackermann)

La ley de control por realimentación de estado es $\mathbf{u} = -K\mathbf{x}$.
Sustituida en el sistema, transforma la dinámica en
$$ \dot{\mathbf{x}} = (A - BK)\mathbf{x}, $$
estable si y solo si $A - BK$ es de Hurwitz. La pregunta de diseño es qué
espectros $\sigma(A - BK)$ podemos **imponer** eligiendo $K$.

- **Teorema de asignación arbitraria de polos:** si $(A, B)$ es completamente
  controlable, para cualquier conjunto $\{\mu_1, \dots, \mu_n\}$ cerrado bajo
  conjugación existe $K$ con $\sigma(A - BK) = \{\mu_1, \dots, \mu_n\}$. Para una
  entrada ($m = 1$), $K$ es **única**. Esta es la razón por la que la
  controlabilidad importa tanto: garantiza que podemos colocar los polos donde
  queramos, en particular en el semiplano izquierdo.

- **Fórmula de Ackermann** (para $m = 1$), con polinomio deseado
  $p_d(\lambda) = \prod_i (\lambda - \mu_i)$:
  $$ K = \mathbf{e}_n^\top\, \mathcal{C}^{-1}\, p_d(A), \qquad \mathbf{e}_n^\top = (0, \dots, 0, 1). $$
  Requiere $\mathcal{C}$ invertible (controlabilidad) y usa Cayley-Hamilton
  ($p_d(A)$ es una suma finita de potencias). Es elegante y exacta, pero **exige
  elegir a mano los polos**, lo que para orden alto (el doble, $n = 6$) no es
  evidente.

### 7.2 Control óptimo: el regulador lineal cuadrático (LQR)

En lugar de prescribir los polos, el LQR los obtiene minimizando un compromiso
entre rapidez de regulación y esfuerzo de control.

- **Problema LQR:** dadas matrices de peso $Q \succeq 0$ (penaliza desviaciones
  del estado) y $R \succ 0$ (penaliza el esfuerzo de control), minimizar
  $$ J = \int_0^{\infty} \big(\mathbf{x}^\top Q\, \mathbf{x} + \mathbf{u}^\top R\, \mathbf{u}\big)\, dt. $$

- **Solución del LQR:** si $(A, B)$ es controlable, el control óptimo es
  $\mathbf{u} = -K\mathbf{x}$ con
  $$ K = R^{-1}B^\top P, $$
  donde $P \succ 0$ es la solución estabilizante de la **ecuación algebraica de
  Riccati (CARE)**
  $$ A^\top P + PA - PBR^{-1}B^\top P + Q = 0. $$
  El lazo cerrado $A - BK$ resulta automáticamente de Hurwitz.

### 7.3 Cómo se resuelve la CARE: la matriz hamiltoniana

La CARE es una ecuación matricial **cuadrática** en $P$. Su solución estabilizante
se obtiene del álgebra lineal de la **matriz hamiltoniana**
$$ H = \begin{pmatrix} A & -BR^{-1}B^\top \\ -Q & -A^\top \end{pmatrix} \in \mathbb{R}^{2n\times 2n}. $$

El algoritmo (tal como está en `solve_care`):

1. Calcular eigenvalores y eigenvectores de $H$.
2. Seleccionar los $n$ asociados a $\operatorname{Re}(\lambda) < 0$ (el
   **subespacio invariante estable**).
3. Partir esos eigenvectores en dos bloques $V_1, V_2 \in \mathbb{R}^{n\times n}$
   y formar $P = V_2 V_1^{-1}$.

**Por qué hay exactamente $n$ eigenvalores estables.** Escribiendo
$\mathbf{x}^\top Q\,\mathbf{x} = \|Q^{1/2}\mathbf{x}\|^2$ (con $Q^{1/2}$ la raíz
cuadrada de $Q$), el costo "observa" el estado a través de la salida ficticia
$\mathbf{z} = Q^{1/2}\mathbf{x}$. Que $(A, B)$ sea controlable y $(Q^{1/2}, A)$
observable (ambos verificados aquí con Kalman) garantiza que ningún modo inestable
queda fuera del alcance del control ni oculto para el costo, de modo que $H$ tiene
exactamente $n$ eigenvalores estables y ninguno sobre el eje imaginario. **Una vez
más, todo se reduce a eigenvalores, eigenvectores y subespacios invariantes.**

### 7.4 LQR frente a Ackermann: cuándo cada uno

| | **Ackermann** | **LQR** |
|---|---|---|
| Qué se especifica | Los polos deseados, a mano | Los pesos $Q$, $R$ |
| Resultado | Polos exactos prescritos | Polos óptimos (compromiso costo-esfuerzo) |
| Escala a $n$ alto | Difícil (elegir 6 polos no es evidente) | Natural (misma elección de pesos) |
| Herramientas de álgebra lineal | $\mathcal{C}^{-1}$, Cayley-Hamilton | Hamiltoniano, subespacio invariante estable |
| Uso en el proyecto | Simple (comparación) | Simple y doble |

---

## 8. De la teoría al código: mapeo concepto a implementación

Esta sección conecta cada concepto teórico con su realización en Julia. Es el
puente entre las secciones 4-7 y los archivos de `src/`.

### 8.1 Tabla maestra concepto → función → ubicación

| Concepto (sección) | Función Julia | Archivo | Notas |
|---|---|---|---|
| Parámetros físicos simple | `SystemParams`, `default_params` | `model_simple.jl` | Barra uniforme: $I = \tfrac{1}{12}m(2L)^2$ |
| EOM no lineales simple | `nonlinear_eom!` | `model_simple.jl` | Resuelve el sistema $2\times2$ con $D$ analítico |
| Lazo cerrado simple | `closed_loop_eom!` | `model_simple.jl` | $u = -K\mathbf{x}$, con saturación opcional |
| Parámetros físicos doble | `SystemParamsDouble`, `default_params_double` | `model_double.jl` | Masas puntuales, sin fricción |
| EOM no lineales doble | `nonlinear_eom_double!` | `model_double.jl` | $\mathbf{M}(q)\ddot q = \text{rhs}$ vía Cholesky |
| Lazo cerrado doble | `closed_loop_eom_double!` | `model_double.jl` | Reutiliza la física, pasa $u$ ya calculado |
| Linealización $A,B,C,D$ simple | `linearize_system` | `linearization.jl` | Jacobiano analítico entrada por entrada |
| Linealización $A,B,C,D$ doble | `linearize_system_double` | `linearization.jl` | Jacobiano analítico; coincide con el numérico a $10^{-11}$ |
| Espectro $\sigma(A)$ | `eigen` dentro de las anteriores | `linearization.jl` | Guardado en `StateSpaceModel` |
| Matriz de controlabilidad $\mathcal{C}$ | `controllability_matrix`, `check_controllability` | `linearization.jl` | $[B\ AB\ \dots\ A^{n-1}B]$, `rank == n` |
| Matriz de observabilidad $\mathcal{O}$ | `observability_matrix`, `check_observability` | `linearization.jl` | $[C; CA; \dots; CA^{n-1}]$, `rank == n` |
| CARE vía Hamiltoniano | `solve_care` | `controller.jl` | Subespacio invariante estable, $P = V_2 V_1^{-1}$ |
| LQR | `design_lqr` | `controller.jl` | $K = R^{-1}B^\top P$ |
| Ackermann | `design_pole_placement` | `controller.jl` | $K = \mathbf{e}_n^\top \mathcal{C}^{-1} p_d(A)$ |
| Animación simple | `animate_pendulum` | `animation_simple.jl` | GLMakie, observables reactivos |
| Animación doble | `animate_pendulum_double` | `animation_double.jl` | Dos eslabones; misma interfaz `save_animation` |

### 8.2 Los tres extractos que hay que saber explicar

Estos tres fragmentos condensan la relación teoría-código y aparecen en la
presentación (frame "Los dos pseudocódigos centrales" y "De $\mathbf{M}(q)\ddot q$
a las aceleraciones").

**a) CARE por la matriz hamiltoniana** (`controller.jl`):

```julia
function solve_care(A, B, Q, R)
    n = size(A, 1)
    S = B * inv(R) * B'
    H = [A -S; -Q -A']                        # matriz hamiltoniana 2n x 2n
    eig = eigen(H)
    stable = findall(real.(eig.values) .< 0)  # subespacio invariante estable
    V = eig.vectors[:, stable]
    V1 = V[1:n, :]; V2 = V[n+1:2n, :]
    P = real.(V2 / V1)                         # P = V2 * inv(V1)
    return (P + P') / 2                         # simetrizar
end
```
Cada línea es un concepto: construir $H$, eigendescomposición, seleccionar el
subespacio invariante estable (los $n$ eigenvalores con $\operatorname{Re} < 0$),
formar $P = V_2 V_1^{-1}$, simetrizar para limpiar errores numéricos.

**b) Ackermann** (`controller.jl`):

```julia
function design_pole_placement(A, B, desired_poles)
    n = size(A, 1)
    C_ctrl = hcat([A^k * B for k in 0:n-1]...)      # matriz de controlabilidad
    rank(C_ctrl) < n && error("No controlable")     # criterio de Kalman
    phi_A = prod(A - p*I for p in desired_poles)    # p_d(A), Cayley-Hamilton
    en = zeros(1, n); en[1, n] = 1.0                # e_n'
    K = en * inv(C_ctrl) * real.(phi_A)             # formula de Ackermann
end
```
El polinomio deseado se evalúa en $A$ con aritmética compleja (para admitir pares
conjugados) y el resultado se toma real.

**c) EOM del doble por Cholesky** (`model_double.jl`):

```julia
Mq = [ M+m1+m2       (m1+m2)*L1*c1  m2*L2*c2;
       (m1+m2)*L1*c1 (m1+m2)*L1^2   m2*L1*L2*c12;
       m2*L2*c2      m2*L1*L2*c12   m2*L2^2 ]
rhs = [ F + (m1+m2)*L1*s1*w1^2 + m2*L2*s2*w2^2,
        (m1+m2)*g*L1*s1 - m2*L1*L2*s12*w2^2,
        m2*g*L2*s2 + m2*L1*L2*s12*w1^2 ]
qdd = cholesky(Symmetric(Mq)) \ rhs   # aceleraciones via Cholesky
```
Punto sutil que conviene tener claro: escribir $\mathbf{M}(q)\ddot q = \text{rhs}$
**no** es una factorización, es *despejar*. Como $\ddot q$ entra linealmente, se
pasan al lado derecho todos los términos sin aceleraciones. La factorización
aparece al **resolver**: `cholesky(Symmetric(Mq))` factoriza
$\mathbf{M} = LL^\top$ (existe exactamente para matrices simétricas definidas
positivas) y `\` resuelve dos sistemas triangulares, sin calcular $\mathbf{M}^{-1}$.

### 8.3 Detalles de implementación que suelen preguntarse

- **Convención de signo (crítica):** el término gravitatorio usa el signo del
  equilibrio *superior*, esto es $A_{2,3} = -m^2gL^2/D_0$ y
  $A_{4,3} = +(M+m)mgL/D_0$ en el simple. El signo opuesto describiría el péndulo
  colgando (estable). Se verifica comparando el espectro obtenido con el esperado
  (un eigenvalor real positivo para el simple, dos para el doble).
- **Solver numérico:** `Tsit5()` (Tsitouras 5/4), un Runge-Kutta explícito de
  orden 5. Las EOM se integran con `saveat` fino (0.005-0.01 s).
- **Saturación del actuador:** el lazo cerrado admite `saturate` (clamp de la
  fuerza): $\pm50$ N en el simple, $\pm100$ N en el doble. Las fuerzas reales pico
  quedan muy por debajo (ver 9.3).
- **Genericidad:** `print_controller_summary` y las rutinas de análisis generan
  etiquetas $K_1, \dots, K_n$ automáticamente, funcionando para cualquier $n$.

---

## 9. Resultados numéricos e interpretación

Corresponde a la Sección 5 del informe y a los frames de resultados de la
presentación. **Todos estos valores se generan corriendo el código, no se
transcriben a mano.**

### 9.1 Parámetros usados

- **Simple:** barra uniforme de longitud total $2L = 1$ m, $M = 1.0$ kg,
  $m = 0.3$ kg, $L = 0.5$ m, $I = \tfrac{1}{12}m(2L)^2 = 0.025$ kg·m$^2$,
  $b = 0.1$ N·s/m, $g = 9.81$ m/s$^2$. Con estos valores, $D_0 = 0.1075$.
- **Doble:** $M = 1.0$ kg, $m_1 = m_2 = 0.3$ kg, $L_1 = L_2 = 0.5$ m,
  $g = 9.81$ m/s$^2$, sin fricción.

### 9.2 Lazo abierto: espectros y rangos

| Configuración | Eigenvalores de $A$ | Diagnóstico |
|---|---|---|
| Simple ($n = 4$) | $+4.21,\ 0,\ -0.077,\ -4.23$ | Inestable (1 modo) |
| Doble ($n = 6$) | $+8.57,\ +4.09,\ 0,\ 0,\ -4.09,\ -8.57$ | Inestable (2 modos) |

- El eigenvalor **positivo** ($+4.21$ en el simple, $+8.57$ y $+4.09$ en el doble)
  es la firma de la inestabilidad predicha por el signo del término gravitatorio.
- Los eigenvalores **nulos** corresponden al modo de traslación libre del carro
  (no hay fuerza recuperadora sobre su posición). En el doble, el $\lambda = 0$
  doble es el bloque de Jordan del que se habló en 6.3.
- **Rangos de Kalman:** en ambos casos $\operatorname{rank}(\mathcal{C}) = n$ y
  $\operatorname{rank}(\mathcal{O}) = n$ ($4/4$ para el simple, $6/6$ para el
  doble). Los dos sistemas son **completamente controlables y observables**, así
  que su espectro puede reubicarse por realimentación.

### 9.3 Lazo cerrado: ganancias y polos

| Método | Ganancia $K$ y polos de lazo cerrado |
|---|---|
| Simple — LQR, $Q = \operatorname{diag}(1,0,10,0)$, $R = 0.1$ | $K = (-3.16, -4.69, -45.39, -10.93)$;  $\sigma = \{-4.48\pm1.56i,\ -1.01\pm0.95i\}$ |
| Simple — Ackermann, polos $\{-1,-2,-3,-4\}$ | $K = (-1.75, -3.75, -39.01, -9.60)$;  $\sigma = \{-1, -2, -3, -4\}$ |
| Doble — LQR, $Q = \operatorname{diag}(1,0,10,0,10,0)$, $R = 0.1$ | $K = (3.16, 5.82, -191.55, -10.99, 228.32, 36.14)$;  $\sigma = \{-8.69\pm1.00i,\ -4.31\pm1.74i,\ -0.89\pm0.82i\}$ |

Interpretación:

- En todos los casos $A - BK$ es de **Hurwitz**: el controlador estabiliza el
  equilibrio erguido.
- Las dos vías del simple producen ganancias distintas **porque persiguen
  objetivos distintos** (LQR minimiza el costo; Ackermann fuerza polos
  prescritos), pero ambas estabilizan. La diferencia es de criterio de diseño, no
  de capacidad.
- Las magnitudes mucho mayores de $K$ en el doble reflejan que estabilizar **dos**
  modos inestables, uno de ellos rápido ($\lambda = +8.57$), exige un esfuerzo de
  control considerablemente mayor.

### 9.4 Métricas de respuesta (verificación sobre el modelo no lineal)

Se integra el modelo **no lineal** en lazo cerrado con $\mathbf{u} = -K\mathbf{x}$.
El tiempo de asentamiento $t_s$ usa el umbral absoluto $|\theta| < 0.02$ rad
($\approx 1.1°$).

| Caso | $t_s$ [s] | $\max|u|$ [N] | Excursión del carro |
|---|---|---|---|
| Simple — LQR | 1.45 | 6.8 | $\approx 0.33$ m |
| Simple — Ackermann | 1.54 | 5.9 | — |
| Doble — LQR | 1.77 | 4.2 | $\approx 0.43$ m |

- El LQR del simple ubica polos complejos y asienta algo más rápido, con respuesta
  levemente sub-amortiguada. Ackermann, con polos reales, da respuesta sin
  oscilación y esfuerzo pico ligeramente menor.
- En el doble, pese a que las ganancias son un orden de magnitud mayores, para
  perturbaciones pequeñas la fuerza pico se mantiene moderada (4.2 N); el costo se
  traslada a una **mayor excursión del carro** (0.43 m frente a 0.33 m del
  simple), que debe maniobrar más para equilibrar dos eslabones acoplados.

### 9.5 Robustez frente a la condición inicial (solo en la presentación)

El frame "Distintos casos" muestra el LQR del simple desde ángulos iniciales
crecientes. Aunque se diseñó sobre el modelo **lineal**, estabiliza el modelo no
lineal incluso desde $\theta_0 = 0.40$ rad ($\approx 23°$), muy fuera de la
vecindad de linealización:

| $\theta_0$ | $t_s$ [s] | $\max|u|$ [N] |
|---|---|---|
| 0.10 rad | 1.30 | 4.5 |
| 0.40 rad ($\approx 23°$) | 3.49 | 18.2 |

En el plano de fase, todas las trayectorias convergen al origen. Este es un
resultado fuerte para cerrar la presentación: **el diseño lineal es robusto más
allá de donde la teoría lo garantiza estrictamente.**

### 9.6 Lectura geométrica (modos propios)

Con $A_c = A - BK$ diagonalizable,
$$ e^{A_c t} = \sum_{i=1}^{n} e^{\mu_i t}\, \mathbf{v}_i \mathbf{w}_i^\top, $$
con $\mathbf{v}_i, \mathbf{w}_i$ los eigenvectores derecho e izquierdo. Cada
sumando es la proyección de la trayectoria sobre el $i$-ésimo modo propio y decae
a la tasa $|\operatorname{Re}(\mu_i)|$. Los polos más a la izquierda son los modos
que se extinguen más rápido. Así, la descomposición espectral de $A_c$ explica la
forma con que las trayectorias vuelven al origen.

---

## 10. Los artefactos de documentación y cómo se relacionan

El proyecto tiene cuatro artefactos, con niveles de detalle decrecientes-crecientes
y audiencias distintas. Entender su relación es clave para saber de dónde sacar
cada cosa al preparar la presentación.

### 10.1 Mapa de artefactos

| Artefacto | Ubicación | Audiencia | Rol |
|---|---|---|---|
| **Informe técnico** | `docs/resumen_tecnico/` | Evaluador que quiere el detalle completo | Fuente de verdad: todas las definiciones, teoremas, demostraciones y derivaciones |
| **Resumen ejecutivo** | `docs/resumen_ejecutivo/` | Lectura rápida (máx. 5 páginas) | Mapa del proyecto; remite al técnico en cada sección |
| **Presentación** | `docs/presentacion/` | Exposición oral (20 min) | Selección visual de lo esencial |
| **Notebooks Pluto** | `notebooks/` | Exploración interactiva | Demostración en vivo; recalculan todo al mover sliders |
| **Esta guía** | `docs/guia_maestra.md` | Preparación integral | Relaciona todo lo anterior con el código |

### 10.2 El informe técnico (la fuente de verdad)

Estructura (≈ 1660 líneas de LaTeX, formato artículo tipo tesis con Palatino):

1. **El problema** (subactuado, inestable, no lineal; las dos configuraciones).
2. **Formalización:** espacios del problema, Euler-Lagrange, deducción del simple
   y del doble, linealización, espacio de estados.
3. **Herramientas de álgebra lineal:** exponencial matricial, eigenvalores y
   Cayley-Hamilton, polinomio minimal y Jordan, Hurwitz y Lyapunov,
   controlabilidad, observabilidad, asignación de polos, LQR.
4. **Metodología:** el flujo de cómputo en Julia (modelo → linealización →
   controlador → simulación).
5. **Resultados:** lazo abierto, lazo cerrado, respuesta temporal, lectura
   geométrica.
6. **Conclusiones.**
7. **Apéndices:** A (derivación de los lagrangianos), B (verificación numérica con
   las matrices explícitas), C (implementación: árbol de archivos y extractos de
   código).

Usa entornos de teorema numerados (`definicion`, `teorema`, `proposicion`, etc.) y
cajas `clave` para los resultados centrales. Las figuras (`simple_respuesta.png`,
`simple_lqr_vs_acker.png`, `doble_respuesta.png`) las genera
`make_report_figs.jl`.

### 10.3 El resumen ejecutivo

Máximo 5 páginas, mismo formato visual que el técnico. Su propósito explícito es
**servir de guía para recorrer el proyecto**: en cada sección remite al informe
técnico con la sección exacta. Contiene una tabla única de resultados clave (la
`tab:exec`, que condensa $n$, modos inestables, rangos, $t_s$, $\max|u|$ de ambas
configuraciones) y una sola figura (la respuesta del simple). Es el documento a
leer primero si alguien llega nuevo al proyecto.

### 10.4 La presentación Beamer (frame por frame)

Formato 16:9, 11pt, Palatino, tema sobrio construido sobre una paleta propia
(marino `#1B4079`, petróleo `#4D7C8A`, gris verde `#7F9C96`, salvia `#8FAD88`,
limón `#C7DB94`). Duración objetivo: 20 minutos. Estructura por secciones:

| Sección | Frames | Contenido |
|---|---|---|
| Portada + contenido | 2 | Título, autores, índice |
| **El problema** | 2 | Los tres rasgos; las dos configuraciones (con diagramas TikZ) |
| **Formalización y modelado** | 8 | Coordenadas y espacios, principio variacional, Euler-Lagrange, config. I, config. II, $\mathbf{M}(q)\ddot q = \text{rhs}$ (con código), linealización, espacio de estados (simple y doble) |
| **Herramientas de álgebra lineal** | 7 | Exponencial, eigenvalores y polinomio característico, Cayley-Hamilton, polinomio minimal, Jordan (el doble la ejercita), Hurwitz/Lyapunov, controlabilidad, observabilidad/dualidad |
| **Diseño de controladores** | 4 | Asignación de polos y Ackermann, LQR, Riccati/Hamiltoniano, los dos pseudocódigos |
| **Implementación y resultados** | 5 | Arquitectura del código, espectros y rangos (con mapa de polos), ganancias y polos, respuesta del simple, robustez (casos $\theta_0$), respuesta del doble |
| **Conclusiones** | 1 | Cierre: el mensaje central |
| Referencias + cierre | 2 | Bibliografía; "¡Gracias!" |

Las figuras (`slides_polos.png`, `slides_simple_respuesta.png`,
`slides_simple_casos.png`, `slides_doble_respuesta.png`) las genera
`make_slide_figs.jl` con la paleta de la presentación.

### 10.5 Los notebooks Pluto (demostración en vivo)

Dos notebooks reactivos (`01_exploracion_simple.jl`, `02_exploracion_doble.jl`):
al mover un slider, **todas las celdas dependientes se recalculan solas** (matrices
$A, B$, eigenvalores, ganancia $K$, polos de lazo cerrado, gráficas y animación).

| Notebook | Sliders de parámetros | Sliders de control | Condición inicial |
|---|---|---|---|
| Simple | $M$, $m$, $L_{bar}$, $g$, $b$ | $Q_{11}$, $Q_{33}$, $R$ | $\theta_0$ |
| Doble | $M$, $m_1$, $m_2$, $L_1$, $L_2$, $g$ | $Q$ (pos, $\theta_1$, $\theta_2$), $R$ | $\theta_{1,0}$, $\theta_{2,0}$ |

El notebook simple incluye, además del análisis, visualizaciones ricas: panel de
las cuatro variables de estado, señal de control con banda de saturación, retrato
de fase $(\theta, \dot\theta)$, balance energético (cinética, potencial, total),
mapa de eigenvalores lazo abierto vs. cerrado, y una animación mecánica del
carro-péndulo con exportación a GIF/MP4. Son ideales para una **demo en vivo**: por
ejemplo, subir $\theta_0$ y mostrar cómo el LQR sigue estabilizando, o poner
$g = 1.62$ (Luna) y discutir si es más fácil de controlar.

### 10.6 Los scripts de figuras (reproducibilidad)

- `make_report_figs.jl`: genera las figuras del informe técnico (CairoMakie,
  salida estática) y **reporta por consola las métricas** ($t_s$, $\max|u|$,
  $\max|x|$) que se citan en la discusión. Se ejecuta con
  `julia --project=. docs/resumen_tecnico/make_report_figs.jl`.
- `make_slide_figs.jl`: genera las figuras de la presentación con la paleta
  Beamer, incluyendo el mapa de polos y el barrido de condiciones iniciales.

Ambos **reutilizan los módulos de `src/`**: no duplican física ni algoritmos.
Definen su propia función `settling_time` con banda absoluta de 0.02 rad.

---

## 11. Guía para la presentación (20 minutos)

### 11.1 El arco narrativo

La presentación cuenta una historia con una tesis clara. El arco es:

1. **Enganche (problema):** "Un péndulo que se cae, y una sola fuerza para
   sostenerlo. Parece un problema de física, pero es un problema de álgebra
   lineal."
2. **Construcción (modelado y linealización):** de la física no lineal a
   $\dot{\mathbf{x}} = A\mathbf{x} + Bu$. Aquí la clave es la *firma de la
   inestabilidad* (el signo).
3. **Núcleo (herramientas):** el espectro diagnostica; el rango decide si podemos
   controlar y observar; Cayley-Hamilton justifica por qué basta hasta $A^{n-1}$.
4. **Resolución (control):** reubicar eigenvalores con Ackermann y LQR; la Riccati
   se resuelve con un subespacio invariante.
5. **Prueba (resultados):** los controladores estabilizan el modelo *no lineal*,
   y el doble muestra que todo escala.
6. **Cierre (tesis):** un puñado de conceptos de álgebra lineal basta para
   modelar, analizar y controlar un sistema dinámico no trivial.

### 11.2 Reparto de tiempo sugerido (20 min)

| Bloque | Tiempo | Frames aprox. |
|---|---|---|
| El problema y las dos configuraciones | 2.5 min | 2 |
| Modelado (Euler-Lagrange, ambas config.) | 4 min | 5 |
| Linealización y espacio de estados | 2.5 min | 3 |
| Herramientas (espectro, Cayley-Hamilton, Jordan, Kalman) | 4.5 min | 7 |
| Control (Ackermann, LQR, Riccati, código) | 3.5 min | 4 |
| Resultados (espectros, respuestas, robustez, doble) | 2.5 min | 5 |
| Conclusiones | 0.5 min | 1 |

Ritmo objetivo: unas 40-45 s por frame de contenido. Si el tiempo aprieta, los
frames comprimibles son "Polinomio minimal y descomposición primaria" (se puede
fusionar con Jordan) y uno de los dos frames de espacio de estados (mostrar solo
el simple y mencionar que el doble es análogo pero $6\times6$).

### 11.3 Reparto entre los tres expositores (sugerencia)

Como son tres autores, una división natural por afinidad con el arco:

- **Expositor 1 (planteamiento):** el problema, las dos configuraciones, el
  modelado por Euler-Lagrange.
- **Expositor 2 (teoría):** linealización, espacio de estados y las herramientas
  de álgebra lineal (el bloque más conceptual).
- **Expositor 3 (control y resultados):** Ackermann, LQR/Riccati, la arquitectura
  del código y todos los resultados, incluida la robustez.

Las conclusiones puede cerrarlas quien abrió, para dar simetría.

### 11.4 Los cinco momentos que hay que clavar

Si algo tiene que salir perfecto, es esto:

1. **La firma de la inestabilidad (el signo).** "El signo positivo de $A_{4,3}$ es
   la diferencia entre colgar y estar invertido." Es memorable y conecta física
   con álgebra.
2. **El ejemplo del bloque angular $\lambda = \pm\sqrt{\kappa}$.** En dos líneas se
   ve de dónde sale el eigenvalor positivo. $\sqrt{\kappa} \approx 4.22$ casi
   coincide con el $+4.21$ real.
3. **Por qué el doble ejercita Jordan.** El $\lambda = 0$ defectuoso, el doble
   integrador del carro, el término secular $t$ (deriva a velocidad constante).
   Distingue una exposición buena de una excelente.
4. **La CARE se resuelve con un subespacio invariante estable.** "Resolver una
   ecuación cuadrática matricial se reduce a eigenvalores y eigenvectores del
   Hamiltoniano." Es el clímax del "todo es álgebra lineal".
5. **La verificación sobre el modelo no lineal.** El controlador se diseñó sobre
   el modelo lineal, pero se prueba sobre el no lineal, y funciona incluso desde
   $23°$. Cierra el círculo.

### 11.5 Frases de transición listas para usar

- Problema → modelado: *"Para responder estas preguntas con álgebra lineal,
  primero necesitamos una matriz. Vamos a obtenerla."*
- Modelado → linealización: *"Estas ecuaciones tienen senos y cosenos; no son
  lineales. Pero solo nos importa la vecindad del equilibrio, y ahí Taylor nos
  basta."*
- Espacio de estados → herramientas: *"Ya tenemos la matriz $A$. Ahora la
  interrogamos: ¿qué nos dice sobre el sistema?"*
- Herramientas → control: *"Sabemos que es inestable pero controlable. Controlable
  significa que podemos mover los eigenvalores. Hagámoslo."*
- Control → resultados: *"Diseñamos sobre el modelo lineal. La prueba de fuego es
  el modelo no lineal."*

---

## 12. Defensa: preguntas frecuentes y puntos delicados

Preguntas plausibles de un evaluador, con respuestas preparadas.

**P: ¿Por qué linealizan si el sistema es no lineal? ¿No pierden validez?**
R: Linealizamos porque el objetivo es *mantener* el péndulo cerca del equilibrio,
y ahí las desviaciones son pequeñas. El modelo lineal es fiel localmente. Lo
importante es que **verificamos el controlador sobre el modelo no lineal
completo**, y estabiliza incluso desde $23°$, muy fuera de la vecindad de
linealización (sección 9.5). Así que la linealización es una herramienta de
*diseño*, no una limitación del resultado.

**P: ¿Por qué el término gravitatorio hace inestable el sistema?**
R: Con $\theta$ medido desde la vertical superior, el término $-mgL\sin\theta$ al
linealizar da un coeficiente positivo en $A_{4,3}$. El bloque angular es entonces
$\left(\begin{smallmatrix}0&1\\\kappa&0\end{smallmatrix}\right)$ con $\kappa > 0$,
cuyos eigenvalores son $\pm\sqrt{\kappa}$: uno positivo. Si el péndulo colgara, el
signo se invertiría y los eigenvalores serían imaginarios puros (oscilación
estable). Es literalmente un signo.

**P: ¿Qué significan los eigenvalores nulos del lazo abierto?**
R: El modo de traslación libre del carro: no hay fuerza recuperadora sobre su
posición, así que una perturbación de posición no genera dinámica que la corrija.
En el doble, el $\lambda = 0$ es doble y defectuoso (bloque de Jordan): corresponde
al doble integrador del carro, cuya respuesta a una perturbación de velocidad es
una deriva a velocidad constante, $x(t) = x_0 + \dot x_0 t$.

**P: ¿Por qué la matriz de controlabilidad solo llega hasta $A^{n-1}B$?**
R: Por Cayley-Hamilton. $A^n$ es combinación lineal de $I, A, \dots, A^{n-1}$,
luego $A^n B$ ya está en el espacio generado por $B, AB, \dots, A^{n-1}B$. Agregar
más columnas no aumenta el rango.

**P: ¿Diferencia entre LQR y Ackermann? ¿Cuál es mejor?**
R: No hay uno mejor; persiguen objetivos distintos. Ackermann fuerza polos que
uno elige a mano (control exacto de la ubicación, pero hay que saber elegir).
El LQR minimiza un costo cuadrático $\int (\mathbf{x}^\top Q\mathbf{x} +
u^\top R u)\,dt$ y **deriva** los polos óptimos de esa elección de pesos. El LQR
escala mejor: para el doble ($n = 6$) elegir 6 polos a mano no es evidente, pero
elegir dos pesos sí. Por eso el doble usa solo LQR.

**P: ¿Cómo se resuelve la ecuación de Riccati, que es no lineal?**
R: Es cuadrática en $P$, pero su solución estabilizante se obtiene por álgebra
lineal: se forma la matriz hamiltoniana $H$ ($2n \times 2n$), se calculan sus
eigenvectores, se toma el subespacio invariante estable (los $n$ eigenvalores con
parte real negativa) y se forma $P = V_2 V_1^{-1}$. Que haya exactamente $n$
estables lo garantizan la controlabilidad de $(A,B)$ y la observabilidad de
$(Q^{1/2}, A)$.

**P: ¿Por qué usan Cholesky y no simplemente `inv(Mq)`?**
R: Porque $\mathbf{M}(q)$ es simétrica y definida positiva, y para ese caso
Cholesky ($\mathbf{M} = LL^\top$) es más rápida y numéricamente más estable que
una inversión o una LU genérica. Además, no necesitamos $\mathbf{M}^{-1}$
explícita: solo resolver $\mathbf{M}\ddot q = \text{rhs}$, que Cholesky hace con
dos sustituciones triangulares.

**P: ¿Miden todas las variables de estado?**
R: No. $C$ mide solo posición del carro y ángulos, no las velocidades. Por eso
importa la observabilidad: verificamos $\operatorname{rank}(\mathcal{O}) = n$, lo
que garantiza que las velocidades se pueden **reconstruir** a partir de las
mediciones. Esto legitima usar realimentación de estado completo (en una
implementación real, con un observador/estimador).

**P: ¿Las ganancias del doble son enormes ($228$, $-191$). ¿Es un problema?**
R: Reflejan que hay **dos** modos inestables, uno rápido ($\lambda = +8.57$), que
exigen más autoridad de control. Pero para perturbaciones pequeñas la fuerza pico
real es moderada (4.2 N); el costo se traslada a que el carro se desplace más
(0.43 m). Las ganancias grandes multiplican estados pequeños.

**Puntos delicados a no equivocar:**

- El ángulo se mide desde la vertical **superior** ($\theta = 0$ arriba). No es la
  convención habitual de péndulo colgante.
- El doble usa **masas puntuales sin fricción**; el simple, **barra uniforme con
  inercia y fricción $b$**. No son "el mismo sistema con un eslabón más": difieren
  en las hipótesis físicas.
- $D_0$ es el determinante de la matriz de masa **evaluado en $\theta = 0$**, no
  en general.
- La convención de ángulos del doble: cada $\theta_i$ se mide desde **su propia**
  vertical (no relativo al eslabón anterior). Esto es lo que da las matrices
  linealizadas limpias.

---

## 13. Glosario de símbolos y checklist final

### 13.1 Glosario

| Símbolo | Significado |
|---|---|
| $q$, $\dot q$ | Coordenadas generalizadas y sus velocidades |
| $d$ | Número de grados de libertad ($d = 2$ simple, $d = 3$ doble) |
| $n$ | Dimensión del estado ($n = 2d$: 4 simple, 6 doble) |
| $\mathbf{x}$ | Vector de estado |
| $u$ | Entrada de control (fuerza horizontal sobre el carro) |
| $\mathbf{y}$ | Salida (lo que se mide) |
| $M$ | Espacio de configuración (variedad) / masa del carro (según contexto) |
| $\mathrm{T}M$ | Fibrado tangente (espacio de estados) |
| $L = T - V$ | Lagrangiano (energía cinética menos potencial) |
| $S[\gamma]$ | Funcional de acción |
| $A, B, C, D$ | Matrices del espacio de estados |
| $D_0$ | Determinante de la matriz de masa del simple en $\theta = 0$ |
| $\mathbf{M}(q)$ | Matriz de masa (simétrica definida positiva) |
| $\sigma(A)$ | Espectro (conjunto de eigenvalores) de $A$ |
| $e^{At}$ | Exponencial matricial (operador solución) |
| $\mathcal{C}$, $\mathcal{O}$ | Matrices de controlabilidad y observabilidad |
| $K$ | Ganancia de realimentación ($u = -K\mathbf{x}$) |
| $Q$, $R$ | Pesos del LQR (estado y esfuerzo) |
| $P$ | Solución de la ecuación de Riccati |
| $H$ | Matriz hamiltoniana ($2n \times 2n$) |
| $\kappa$ | Coeficiente del bloque angular, $(M+m)mgL/D_0$ |
| $t_s$ | Tiempo de asentamiento ($|\theta| < 0.02$ rad) |

### 13.2 Checklist antes de la presentación

- [ ] Compilar el PDF de la presentación y verificar que las cuatro figuras
      aparecen (regenerarlas con `make_slide_figs.jl` si hace falta).
- [ ] Tener a mano los valores clave de memoria: espectros ($+4.21$; $+8.57$),
      rangos ($4/4$, $6/6$), $t_s$ ($1.45$ s, $1.77$ s).
- [ ] Ensayar los cinco momentos que hay que clavar (sección 11.4).
- [ ] Preparar (opcional) una demo en vivo del notebook: subir $\theta_0$ o poner
      gravedad lunar y comentar.
- [ ] Repasar la sección 12 (defensa) con los tres expositores.
- [ ] Confirmar la división de bloques y los tiempos (secciones 11.2 y 11.3).
- [ ] Tener claro el mensaje de cierre: *un puñado de conceptos de álgebra lineal
      es el lenguaje natural para modelar, analizar y controlar un sistema
      dinámico.*

### 13.3 Para profundizar

- **Detalle teórico completo:** informe técnico (`docs/resumen_tecnico/resumen_tecnico.pdf`),
  con las demostraciones y los dos apéndices de derivación y verificación.
- **Panorama rápido:** resumen ejecutivo (`docs/resumen_ejecutivo/resumen_ejecutivo.pdf`).
- **Código ejecutable:** `main_simple.jl`, `main_double.jl` y los módulos de `src/`.
- **Exploración interactiva:** notebooks de `notebooks/` (ver instrucciones en el
  `README.md`).
- **Bibliografía base:** Ogata (*Modern Control Engineering*), Chen (*Linear System
  Theory and Design*), Hoffman-Kunze y Strang (*Linear Algebra*), Hirsch-Smale-Devaney
  (*Differential Equations, Dynamical Systems*), Golub-Van Loan (*Matrix
  Computations*).

---

## 14. Guía de ejecución paso a paso

Esta sección explica, de forma concreta y en orden, cómo correr todo el código
del proyecto: desde instalar las dependencias hasta reproducir las figuras y los
números de la sección 9. Los comandos están pensados para **Windows con
PowerShell** (el entorno de trabajo), pero funcionan igual en cualquier terminal
donde `julia` esté en el `PATH`. Las rutas con `/` funcionan dentro de Julia
incluso en Windows.

Supuesto de partida: el repositorio está en
`C:\Users\surib\Proyectos\proyecto-pendulos`. Todos los comandos se ejecutan
**desde esa carpeta** salvo que se indique lo contrario.

### 14.0 Qué genera cada script (mapa rápido)

| Quiero... | Ejecuto | Produce |
|---|---|---|
| Instalar dependencias (1 sola vez) | `julia setup.jl` | Entorno listo (precompilado) |
| Analizar y controlar el simple | `julia main_simple.jl` | Consola + `figures/01`, `02`, `03_animacion_lqr.mp4` |
| Analizar y controlar el doble | `julia main_double.jl` | Consola + `figures/06`, `07_doble_animacion_lqr.mp4` |
| Figuras del informe técnico | `julia --project=. docs/resumen_tecnico/make_report_figs.jl` | `docs/resumen_tecnico/figs/*.png` + métricas en consola |
| Figuras de la presentación | `julia --project=. docs/presentacion/make_slide_figs.jl` | `docs/presentacion/figs/*.png` + métricas en consola |
| Exploración interactiva | `julia --project=.` y luego `import Pluto; Pluto.run()` | Notebook en el navegador |

La carpeta `figures/` se crea sola la primera vez (no hay que crearla a mano).

### 14.1 Paso 0: requisito previo (instalar Julia)

Se necesita [Julia](https://julialang.org/) instalado y accesible desde la
terminal. Para comprobarlo:

```powershell
julia --version
```

Debe imprimir una versión (por ejemplo `julia version 1.10.x`). Si el comando no
se reconoce, hay que instalar Julia y asegurarse de que su carpeta `bin` esté en
el `PATH` del sistema.

### 14.2 Paso 1: instalar las dependencias (una sola vez)

Desde la carpeta del proyecto:

```powershell
cd C:\Users\surib\Proyectos\proyecto-pendulos
julia setup.jl
```

Qué hace `setup.jl`: activa el entorno local (el de `Project.toml`), instala los
ocho paquetes (`DifferentialEquations`, `CairoMakie`, `GLMakie`,
`ControlSystems`, `MatrixEquations`, `Symbolics`, `Pluto`, `PlutoUI`) y los
**precompila**.

Advertencia importante sobre el tiempo: esta primera vez **puede tardar varios
minutos** (incluso 10-20 según la máquina y la conexión), porque descarga y
precompila paquetes pesados. Es normal y solo ocurre una vez. Las ejecuciones
posteriores son rápidas.

Alternativa equivalente sin usar `setup.jl`:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

`Pkg.instantiate()` reproduce **exactamente** las versiones fijadas en
`Manifest.toml`, lo que garantiza un entorno idéntico al de los autores.

### 14.3 Paso 2: ejecutar el pipeline del péndulo simple

```powershell
julia main_simple.jl
```

Este script no necesita `--project=.` porque activa el entorno internamente
(`Pkg.activate(@__DIR__)`). Ejecuta los ocho pasos del pipeline e imprime el
análisis por consola. Al terminar habrá generado, en `figures/`:

- `01_comparativa_lqr.png` — cuatro paneles: ángulo, posición, fuerza de control
  y velocidades (libre vs. LQR).
- `02_comparativa_ackermann.png` — ángulo y posición, LQR vs. Ackermann.
- `03_animacion_lqr.mp4` — animación del carro-péndulo estabilizado por LQR.

Qué esperar en la consola (debe coincidir con la sección 9):

- Espectro de lazo abierto $\{+4.21,\ 0,\ -0.077,\ -4.23\}$, diagnóstico
  **INESTABLE**.
- Controlabilidad y observabilidad: `rank = 4 / 4` → **CONTROLABLE** y
  **OBSERVABLE**.
- Ganancia LQR $K = (-3.16, -4.69, -45.39, -10.93)$ y polos de lazo cerrado todos
  estables.
- Ganancia Ackermann $K = (-1.75, -3.75, -39.01, -9.60)$ con polos exactos
  $\{-1, -2, -3, -4\}$.

Nota: `main_simple.jl` usa **GLMakie** para la animación, que abre una ventana
gráfica. Requiere un entorno con soporte de pantalla (en un Windows de escritorio
normal funciona sin problema). Si se ejecuta en una máquina sin display, ver
14.9.

### 14.4 Paso 3: ejecutar el pipeline del péndulo doble

```powershell
julia main_double.jl
```

Análogo al anterior, para la Configuración II. Genera en `figures/`:

- `06_doble_comparativa_lqr.png` — ángulos $\theta_1$ y $\theta_2$, posición y
  fuerza (libre vs. LQR).
- `07_doble_animacion_lqr.mp4` — animación de los dos eslabones estabilizados.

Qué esperar en consola:

- Espectro de lazo abierto $\{+8.57,\ +4.09,\ 0,\ 0,\ -4.09,\ -8.57\}$, dos modos
  inestables.
- `rank = 6 / 6` en controlabilidad y observabilidad.
- Ganancia LQR $K = (3.16, 5.82, -191.55, -10.99, 228.32, 36.14)$.

### 14.5 Paso 4 (opcional): uso interactivo desde el REPL

En lugar de correr el script de un tirón, se puede trabajar dentro de una sesión
de Julia, lo que evita reprecompilar entre ejecuciones sucesivas:

```powershell
julia --project=.
```

Y dentro del prompt `julia>`:

```julia
include("main_simple.jl")   # corre todo el pipeline del simple
include("main_double.jl")   # o el del doble
```

Aquí el flag `--project=.` sí importa: activa el entorno del proyecto para que
Julia encuentre todos los paquetes. Ventaja del modo interactivo: tras el primer
`include`, las variables (`ss`, `lqr_result`, `sol_lqr`, etc.) quedan disponibles
para inspeccionarlas a mano, por ejemplo:

```julia
ss.eigenvalues          # espectro de A
lqr_result.K            # ganancia LQR
check_controllability(ss)
```

### 14.6 Paso 5: exploración interactiva con Pluto (paso a paso)

Los notebooks recalculan **todo** el análisis al mover un slider. Son ideales
para una demo en vivo durante la presentación.

1. **Abrir Julia en la carpeta del proyecto** (el `--project=.` es imprescindible
   para que el notebook encuentre los paquetes):

   ```powershell
   cd C:\Users\surib\Proyectos\proyecto-pendulos
   julia --project=.
   ```

2. **Lanzar el servidor de Pluto** desde el prompt `julia>`:

   ```julia
   import Pluto
   Pluto.run()
   ```

   Esto abre Pluto en el navegador (normalmente `http://localhost:1234`). La
   primera vez puede tardar mientras precompila. **Deja esa terminal abierta:** es
   el servidor; si la cierras, Pluto se detiene.

3. **Abrir un notebook.** En la pantalla de inicio de Pluto, en el campo *"Open a
   notebook"*, pega la ruta y pulsa **Open**:
   - `notebooks/01_exploracion_simple.jl` — péndulo simple ($\mathbb{R}^4$).
   - `notebooks/02_exploracion_doble.jl` — péndulo doble ($\mathbb{R}^6$).

   La primera apertura precompila CairoMakie y DifferentialEquations (puede tardar
   varios minutos); las siguientes son rápidas.

4. **Interactuar con los sliders.** Al cambiar cualquiera, todas las celdas
   dependientes se recalculan solas. Verás reaccionar en vivo el eigenvalor
   inestable, la ganancia $K$, los polos de lazo cerrado, las gráficas y la
   animación.

   | Notebook | Sliders de parámetros | Sliders de control | Condición inicial |
   |---|---|---|---|
   | Simple | `M`, `m`, `Lbar`, `g`, `b` | `Q11`, `Q33`, `R` | `theta0` |
   | Doble | `M`, `m1`, `m2`, `L1`, `L2`, `g` | `Q` (pos, `theta1`, `theta2`), `R` | `theta1_0`, `theta2_0` |

5. **Exportar la animación (opcional).** Al final de cada notebook hay una casilla
   (checkbox). Al marcarla se genera el archivo en `figures/`:
   - Simple: `04_comparacion_libre_vs_lqr.gif` (y `05_..._.mp4` si hay `ffmpeg`).
   - Doble: `08_doble_exploracion.gif`.

   Desmárcala después para no regenerar el archivo en cada cambio de slider.

6. **Cerrar.** Guarda el notebook (Pluto guarda solo el `.jl`), cierra la pestaña
   del navegador y detén el servidor con `Ctrl-C` en la terminal de Julia.

Ideas de demo para la presentación (ver también 11.5): subir `theta0` a 0.40 rad
y mostrar que el LQR sigue estabilizando; poner `g = 1.62` (gravedad lunar) y
discutir si es más fácil de controlar; aumentar `Q33` y ver cómo crece `K_theta`.

### 14.7 Paso 6: regenerar las figuras de los documentos

Las figuras de respuesta temporal que aparecen en el informe y en las
diapositivas **no** se dibujan a mano: las produce un script que reutiliza los
módulos de `src/`. Además, cada script **imprime por consola las métricas**
($t_s$, $\max|u|$, $\max|x|$) que se citan en los textos.

Figuras del informe técnico:

```powershell
julia --project=. docs/resumen_tecnico/make_report_figs.jl
```

Reescribe `docs/resumen_tecnico/figs/` con `simple_respuesta.png`,
`simple_lqr_vs_acker.png` y `doble_respuesta.png`, y reporta las métricas
(Simple LQR: $t_s = 1.45$ s, $\max|u| = 6.8$ N; etc.).

Figuras de la presentación:

```powershell
julia --project=. docs/presentacion/make_slide_figs.jl
```

Reescribe `docs/presentacion/figs/` con `slides_polos.png`,
`slides_simple_respuesta.png`, `slides_simple_casos.png` y
`slides_doble_respuesta.png` (con la paleta de color de las diapositivas).

Estos scripts usan **CairoMakie** (salida estática, sin ventana), por lo que
funcionan también en máquinas sin entorno gráfico.

### 14.8 Paso 7: compilar los documentos LaTeX

Requiere una distribución de LaTeX instalada (por ejemplo MiKTeX o TeX Live) con
`pdflatex` y `latexmk`. **Importante:** las figuras `.png` deben existir antes de
compilar (paso 14.7), o LaTeX fallará al no encontrarlas.

```powershell
# Informe tecnico
cd C:\Users\surib\Proyectos\proyecto-pendulos\docs\resumen_tecnico
latexmk -pdf resumen_tecnico.tex

# Resumen ejecutivo (reutiliza las figuras del informe tecnico)
cd C:\Users\surib\Proyectos\proyecto-pendulos\docs\resumen_ejecutivo
latexmk -pdf resumen_ejecutivo.tex

# Presentacion
cd C:\Users\surib\Proyectos\proyecto-pendulos\docs\presentacion
latexmk -pdf presentacion.tex
```

Si no se dispone de `latexmk`, sirve `pdflatex nombre.tex` ejecutado dos veces
(la segunda pasada resuelve las referencias cruzadas). El resumen ejecutivo
apunta a `../resumen_tecnico/figs/`, así que basta con haber generado las figuras
del informe.

### 14.9 Orden recomendado de principio a fin

Si se parte de cero y se quiere reproducir todo el proyecto, este es el orden
mínimo:

```powershell
cd C:\Users\surib\Proyectos\proyecto-pendulos

# 1. Instalar dependencias (una sola vez, tarda varios minutos)
julia setup.jl

# 2. Correr los dos pipelines (analisis, control, figuras y animaciones)
julia main_simple.jl
julia main_double.jl

# 3. Regenerar las figuras de los documentos y ver las metricas
julia --project=. docs/resumen_tecnico/make_report_figs.jl
julia --project=. docs/presentacion/make_slide_figs.jl

# 4. (Opcional) Compilar los PDFs
cd docs\resumen_tecnico;  latexmk -pdf resumen_tecnico.tex
cd ..\resumen_ejecutivo;  latexmk -pdf resumen_ejecutivo.tex
cd ..\presentacion;       latexmk -pdf presentacion.tex
```

### 14.10 Solución de problemas comunes

| Síntoma | Causa probable | Solución |
|---|---|---|
| `Package X not found` o `not installed` | El entorno no está instanciado | Ejecutar `julia setup.jl`, o `julia --project=. -e "using Pkg; Pkg.instantiate()"` |
| La primera ejecución tarda muchísimo | Precompilación de paquetes pesados | Es normal y ocurre una sola vez; las siguientes son rápidas |
| Error de GLMakie / no abre ventana | Máquina sin entorno gráfico (servidor, SSH) | Usar los scripts `make_*_figs.jl` (CairoMakie, estáticos) en vez de `main_*.jl`, o correr los notebooks |
| No se genera el `.mp4` del notebook | Falta `ffmpeg` en el sistema | Instalar `ffmpeg`, o usar la exportación a `.gif` (no requiere `ffmpeg`) |
| `julia` no se reconoce | Julia no está en el `PATH` | Reinstalar Julia marcando "Add to PATH", o añadir su carpeta `bin` al `PATH` |
| LaTeX falla por figura faltante | Las `.png` no se han generado | Correr antes `make_report_figs.jl` y `make_slide_figs.jl` |
| Los notebooks no encuentran los paquetes | Julia abierto sin `--project=.` | Cerrar y reabrir con `julia --project=.` antes de `Pluto.run()` |
| Números distintos a los de la sección 9 | Se cambiaron parámetros o pesos | Revisar `default_params`, `default_params_double` y las matrices `Q`, `R` de los `main_*.jl` |

---

*Documento maestro del proyecto. Toda la información numérica proviene de ejecutar
el código de `src/` y los pipelines `main_*.jl`; toda la información teórica está
desarrollada en detalle en el informe técnico. Ante cualquier discrepancia, la
fuente de verdad es el código en ejecución y el informe técnico.*
