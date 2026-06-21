# Conceptos Teóricos: Física, Ecuaciones Diferenciales y Cálculo de Variaciones

## Fundamentos para el Modelado del Péndulo Invertido sobre Carro

---

## Tabla de Contenidos

1. [Mecánica Newtoniana](#1-mecánica-newtoniana)
2. [Mecánica Lagrangiana](#2-mecánica-lagrangiana)
3. [Cálculo de Variaciones](#3-cálculo-de-variaciones)
4. [Ecuaciones Diferenciales Ordinarias](#4-ecuaciones-diferenciales-ordinarias)
5. [Sistemas de EDOs Lineales y Representación en Espacio de Estados](#5-sistemas-de-edos-lineales-y-representación-en-espacio-de-estados)
6. [Estabilidad y Retratos de Fase](#6-estabilidad-y-retratos-de-fase)
7. [Linealización de Sistemas No Lineales](#7-linealización-de-sistemas-no-lineales)
8. [Exponencial Matricial y Solución de Sistemas Lineales](#8-exponencial-matricial-y-solución-de-sistemas-lineales)
9. [Sistemas Forzados y Resonancia](#9-sistemas-forzados-y-resonancia)
10. [Energía, Amortiguamiento y Disipación](#10-energía-amortiguamiento-y-disipación)
11. [Síntesis: De la Física al Control del Péndulo Invertido](#11-síntesis-de-la-física-al-control-del-péndulo-invertido)
12. [Referencias](#12-referencias)

---

## 1. Mecánica Newtoniana

### 1.1. Las Leyes de Newton

El punto de partida de toda la mecánica clásica son las tres leyes de Newton, formuladas en los *Principia Mathematica* (1687). Para el proyecto del péndulo invertido, la segunda ley es la más relevante.

**Segunda Ley de Newton.** Para una partícula de masa $m$ sujeta a una fuerza neta $\mathbf{F}$, la ecuación de movimiento es:

$$\mathbf{F} = m\,\mathbf{a} = m\,\frac{d^2\mathbf{r}}{dt^2}$$

donde $\mathbf{r}(t)$ es el vector de posición y $\mathbf{a}$ la aceleración. En el contexto de sistemas mecánicos con múltiples cuerpos, esta ley se aplica a cada componente por separado. Como señala Olver (§10.5, Ec. 10.62) [2], para un sistema de $n$ masas, la segunda ley conduce directamente a:

$$\text{Fuerza} = \text{Masa} \times \text{Aceleración}$$

que en forma vectorial para una cadena de masas y resortes se escribe:

$$M\,\frac{d^2\mathbf{u}}{dt^2} = \mathbf{f}(t) - K\,\mathbf{u}$$

donde $M = \text{diag}(m_1, \ldots, m_n)$ es la matriz diagonal de masas (positiva definida) y $K = A^T C A$ es la **matriz de rigidez** construida a partir de la matriz de incidencia reducida $A$ y la matriz diagonal de constantes elásticas $C$.

**Referencia:** Olver & Shakiban, §10.5, Ec. (10.62–10.63) [2].

### 1.2. Diagrama de Cuerpo Libre del Sistema Carro-Péndulo

Para el sistema del péndulo invertido sobre carro, se identifican dos cuerpos:

- **El carro** (masa $M$): se mueve horizontalmente a lo largo de un riel, sujeto a una fuerza de control $u(t)$ y a las fuerzas de reacción del péndulo en el pivote.
- **El péndulo** (masa $m$, longitud $l$): unido al carro por un pivote sin fricción. Su posición angular $\theta$ se mide respecto a la vertical superior (posición de equilibrio inestable).

Las coordenadas del centro de masa del péndulo (modelado como masa puntual en el extremo de una varilla sin masa, o como varilla uniforme) son:

$$x_p = x + l\sin\theta, \qquad y_p = l\cos\theta$$

donde $x$ es la posición horizontal del carro. La velocidad del péndulo se obtiene diferenciando:

$$\dot{x}_p = \dot{x} + l\dot{\theta}\cos\theta, \qquad \dot{y}_p = -l\dot{\theta}\sin\theta$$

**Referencia:** Ogata, Cap. 3, Ejemplos 3-5 y 3-6 [1].

### 1.3. Aplicación de la Segunda Ley al Carro y al Péndulo

Aplicando $F = ma$ por separado a cada cuerpo y resolviendo las fuerzas de reacción internas del pivote, se llega a las **ecuaciones de movimiento no lineales acopladas**. Ogata (§3, Ec. 3-20, 3-21) presenta el caso simplificado de masa puntual ($I = 0$) [1]:

$$Ml\ddot{\theta} = (M + m)g\theta - u$$
$$M\ddot{x} = u - mg\theta$$

Estas son las ecuaciones ya linealizadas (válidas para $\theta$ pequeño). Las ecuaciones no lineales completas, antes de linealizar, son:

$$(M + m)\ddot{x} + ml\ddot{\theta}\cos\theta - ml\dot{\theta}^2\sin\theta = u$$
$$ml^2\ddot{\theta} + ml\ddot{x}\cos\theta - mgl\sin\theta = 0$$

**Referencia:** Ogata, Cap. 3, Ej. 3-5 y 3-6, Ecuaciones (3-20) a (3-23) [1].

### 1.4. Momento de Inercia y Modelo de Varilla Uniforme

Para un modelo más realista donde el péndulo es una varilla uniforme de masa $m$ y longitud $2l$ (con centro de masa a distancia $l$ del pivote), el momento de inercia respecto al pivote es:

$$I_{\text{pivote}} = I_{cm} + ml^2 = \frac{1}{3}m(2l)^2 = \frac{4}{3}ml^2$$

Para una varilla de longitud $l$ con pivote en un extremo:

$$I_{\text{pivote}} = \frac{1}{3}ml^2$$

El efecto del momento de inercia es modificar los coeficientes en las ecuaciones de movimiento. En el caso más general, la ecuación angular queda:

$$(I + ml^2)\ddot{\theta} + ml\ddot{x}\cos\theta - mgl\sin\theta = 0$$

Cuando $I = 0$ (masa puntual), se recuperan las ecuaciones simplificadas de Ogata. El momento de inercia introduce un acoplamiento más fuerte entre las dinámicas traslacional y rotacional.

**Referencia:** Ogata, Cap. 3, Ej. 3-5 [1]; cualquier texto de mecánica clásica, por ejemplo Goldstein, *Classical Mechanics*, Cap. 5 [6].

### 1.5. Coordenadas Generalizadas y Grados de Libertad

El sistema carro-péndulo tiene **dos grados de libertad**: la posición del carro $x(t)$ y el ángulo del péndulo $\theta(t)$. Estas son las **coordenadas generalizadas** del sistema, un concepto que será fundamental al pasar a la formulación lagrangiana.

Un grado de libertad corresponde a una coordenada independiente necesaria para especificar completamente la configuración del sistema. Para un sistema con $N$ partículas en el espacio tridimensional sujeto a $k$ restricciones holonómicas, el número de grados de libertad es:

$$n = 3N - k$$

En nuestro caso, el carro está restringido a moverse en una dimensión, y el péndulo está fijado al carro por un pivote, resultando en $n = 2$.

**Referencia:** Goldstein, *Classical Mechanics*, §1.3 [6]; Landau & Lifshitz, *Mechanics*, §1 [7].

---

## 2. Mecánica Lagrangiana

### 2.1. Motivación: ¿Por Qué el Formalismo Lagrangiano?

La mecánica newtoniana requiere identificar **todas** las fuerzas que actúan sobre cada cuerpo, incluyendo las fuerzas de reacción internas (las del pivote, en nuestro caso). Esto se vuelve engorroso para sistemas con restricciones.

El formalismo lagrangiano ofrece ventajas fundamentales:

1. **Trabaja directamente con coordenadas generalizadas**, eliminando automáticamente las fuerzas de restricción.
2. **Es invariante bajo cambios de coordenadas**: se puede elegir el sistema de coordenadas más conveniente.
3. **Se basa en cantidades escalares** (energía cinética y potencial), no vectoriales (fuerzas), lo que simplifica el planteamiento.
4. **Se conecta directamente con el cálculo de variaciones** a través del principio de Hamilton.

**Referencia:** Landau & Lifshitz, *Mechanics*, §1–2 [7]; Goldstein, *Classical Mechanics*, Cap. 1 [6].

### 2.2. Energía Cinética y Energía Potencial

**Energía cinética** del sistema completo. Para el carro-péndulo (modelo masa puntual):

$$T = \frac{1}{2}M\dot{x}^2 + \frac{1}{2}m\left(\dot{x}_p^2 + \dot{y}_p^2\right)$$

Sustituyendo las expresiones de $\dot{x}_p$ y $\dot{y}_p$ de la sección 1.2:

$$T = \frac{1}{2}M\dot{x}^2 + \frac{1}{2}m\left[(\dot{x} + l\dot{\theta}\cos\theta)^2 + (l\dot{\theta}\sin\theta)^2\right]$$

Expandiendo y simplificando con $\cos^2\theta + \sin^2\theta = 1$:

$$T = \frac{1}{2}(M + m)\dot{x}^2 + ml\dot{x}\dot{\theta}\cos\theta + \frac{1}{2}ml^2\dot{\theta}^2$$

Si el péndulo es una varilla uniforme con momento de inercia $I$ respecto a su centro de masa, se añade el término rotacional:

$$T = \frac{1}{2}(M + m)\dot{x}^2 + ml\dot{x}\dot{\theta}\cos\theta + \frac{1}{2}(I + ml^2)\dot{\theta}^2$$

**Energía potencial.** Tomando como referencia la posición del pivote:

$$V = mgl\cos\theta$$

Nótese que para el péndulo **invertido**, el equilibrio está en $\theta = 0$ (péndulo arriba). En esta posición, $V = mgl$ es un **máximo** de la energía potencial, lo que confirma que es un equilibrio inestable desde el punto de vista energético.

**Referencia:** Ogata, Cap. 3, Ej. 3-5 [1]; Goldstein, Cap. 1 [6].

### 2.3. El Lagrangiano

El **lagrangiano** del sistema se define como la diferencia entre la energía cinética y la energía potencial:

$$\mathcal{L}(q, \dot{q}, t) = T - V$$

donde $q = (x, \theta)^T$ son las coordenadas generalizadas y $\dot{q} = (\dot{x}, \dot{\theta})^T$ las velocidades generalizadas.

Para el péndulo invertido (modelo masa puntual):

$$\mathcal{L} = \frac{1}{2}(M+m)\dot{x}^2 + ml\dot{x}\dot{\theta}\cos\theta + \frac{1}{2}ml^2\dot{\theta}^2 - mgl\cos\theta$$

El lagrangiano contiene **toda la información dinámica** del sistema. Las ecuaciones de movimiento se obtienen a partir de él mediante las ecuaciones de Euler-Lagrange.

**Referencia:** Landau & Lifshitz, §1–2 [7]; Goldstein, §1.4 [6].

### 2.4. Las Ecuaciones de Euler-Lagrange

Para un sistema con coordenadas generalizadas $q_i$ ($i = 1, \ldots, n$) y fuerza generalizada no conservativa $Q_i$ asociada a la coordenada $q_i$, las ecuaciones de Euler-Lagrange son:

$$\frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{q}_i}\right) - \frac{\partial \mathcal{L}}{\partial q_i} = Q_i, \qquad i = 1, \ldots, n$$

Cuando todas las fuerzas son conservativas (derivables de un potencial), el lado derecho es cero. En nuestro caso, la fuerza de control $u(t)$ actúa sobre el carro, de modo que:

$$Q_x = u(t), \qquad Q_\theta = 0$$

**Derivación para la coordenada $x$:**

$$\frac{\partial \mathcal{L}}{\partial \dot{x}} = (M+m)\dot{x} + ml\dot{\theta}\cos\theta$$

$$\frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{x}}\right) = (M+m)\ddot{x} + ml\ddot{\theta}\cos\theta - ml\dot{\theta}^2\sin\theta$$

$$\frac{\partial \mathcal{L}}{\partial x} = 0$$

Por tanto:

$$(M+m)\ddot{x} + ml\ddot{\theta}\cos\theta - ml\dot{\theta}^2\sin\theta = u$$

**Derivación para la coordenada $\theta$:**

$$\frac{\partial \mathcal{L}}{\partial \dot{\theta}} = ml\dot{x}\cos\theta + ml^2\dot{\theta}$$

$$\frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{\theta}}\right) = ml\ddot{x}\cos\theta - ml\dot{x}\dot{\theta}\sin\theta + ml^2\ddot{\theta}$$

$$\frac{\partial \mathcal{L}}{\partial \theta} = -ml\dot{x}\dot{\theta}\sin\theta + mgl\sin\theta$$

Sustituyendo y simplificando (los términos $-ml\dot{x}\dot{\theta}\sin\theta$ se cancelan):

$$ml^2\ddot{\theta} + ml\ddot{x}\cos\theta - mgl\sin\theta = 0$$

Estas son exactamente las ecuaciones de movimiento no lineales del sistema, obtenidas sin necesidad de analizar fuerzas de reacción internas.

**Referencia:** Goldstein, §1.4 [6]; Ogata, Cap. 3 [1]; Landau & Lifshitz, §2 [7].

### 2.5. Fuerzas Generalizadas y Fuerzas No Conservativas

No todas las fuerzas derivan de un potencial. La fuerza de control $u(t)$ y la fricción son ejemplos de fuerzas no conservativas. Para incluirlas, se calcula el **trabajo virtual**:

$$\delta W = \sum_i Q_i \, \delta q_i$$

donde $Q_i$ es la fuerza generalizada conjugada a $q_i$. Para la fuerza horizontal $u$ aplicada al carro:

$$\delta W = u\,\delta x + 0 \cdot \delta\theta$$

Si hay fricción viscosa $b\dot{x}$ oponiéndose al movimiento del carro:

$$Q_x = u - b\dot{x}$$

Y si hay fricción en el pivote con coeficiente $\beta$:

$$Q_\theta = -\beta\dot{\theta}$$

**Referencia:** Goldstein, §1.5 [6].

### 2.6. Momentos Generalizados y el Hamiltoniano

El **momento generalizado** conjugado a la coordenada $q_i$ se define como:

$$p_i = \frac{\partial \mathcal{L}}{\partial \dot{q}_i}$$

Para nuestro sistema:

$$p_x = (M+m)\dot{x} + ml\dot{\theta}\cos\theta$$
$$p_\theta = ml\dot{x}\cos\theta + ml^2\dot{\theta}$$

El **hamiltoniano** se define mediante la transformada de Legendre:

$$\mathcal{H} = \sum_i p_i \dot{q}_i - \mathcal{L} = T + V$$

Para sistemas conservativos, $\mathcal{H}$ representa la **energía total** del sistema y se conserva a lo largo de las trayectorias. Este hecho es fundamental para analizar la estabilidad del péndulo invertido: en ausencia de control, la energía total se conserva y el sistema oscila (péndulo normal) o diverge (péndulo invertido).

**Referencia:** Goldstein, §8.1 [6]; Landau & Lifshitz, §6–7 [7].

---

## 3. Cálculo de Variaciones

### 3.1. Motivación y Contexto Histórico

El cálculo de variaciones es la rama del análisis matemático que estudia la **optimización de funcionales** — funciones cuyo dominio es un espacio de funciones. Fue desarrollado por Euler y Lagrange en el siglo XVIII, originalmente para resolver el problema de la **braquistócrona** (la curva de descenso más rápido) planteado por Johann Bernoulli en 1696.

La conexión con la mecánica es profunda: el **principio de Hamilton** establece que las trayectorias físicas de un sistema mecánico son aquellas que hacen estacionaria (no necesariamente mínima) una cierta funcional — la **acción**. Las ecuaciones de Euler-Lagrange surgen como la condición necesaria de estacionariedad.

**Referencia:** Gelfand & Fomin, *Calculus of Variations*, Cap. 1 [8]; Goldstein, §2.1 [6].

### 3.2. Funcionales

Un **funcional** es una correspondencia que asigna un número real a cada función de un cierto espacio. Formalmente, si $\mathcal{C}$ es un espacio de funciones suficientemente regulares, un funcional es una aplicación:

$$J: \mathcal{C} \to \mathbb{R}$$

El ejemplo más importante para nosotros es la **funcional de acción**:

$$J[q] = \int_{t_1}^{t_2} \mathcal{L}(q(t), \dot{q}(t), t)\, dt$$

donde $\mathcal{L}$ es el lagrangiano y $q(t)$ es una trayectoria que satisface condiciones de frontera fijas: $q(t_1) = q_1$, $q(t_2) = q_2$.

A diferencia del cálculo ordinario donde se buscan los extremos de funciones $f: \mathbb{R}^n \to \mathbb{R}$, aquí se buscan las funciones $q(t)$ que hacen estacionaria a $J[q]$.

**Referencia:** Gelfand & Fomin, §1.1–1.3 [8].

### 3.3. La Primera Variación y la Condición de Estacionariedad

Sea $q^*(t)$ una trayectoria candidata y $\eta(t)$ una variación admisible, es decir, una función suficientemente suave con $\eta(t_1) = \eta(t_2) = 0$. La **primera variación** de $J$ en la dirección $\eta$ es:

$$\delta J = \left.\frac{d}{d\epsilon}\right|_{\epsilon=0} J[q^* + \epsilon\,\eta]$$

Calculando explícitamente:

$$J[q^* + \epsilon\,\eta] = \int_{t_1}^{t_2} \mathcal{L}(q^* + \epsilon\eta, \dot{q}^* + \epsilon\dot{\eta}, t)\, dt$$

Diferenciando bajo el signo integral y evaluando en $\epsilon = 0$:

$$\delta J = \int_{t_1}^{t_2} \left[\frac{\partial \mathcal{L}}{\partial q}\eta + \frac{\partial \mathcal{L}}{\partial \dot{q}}\dot{\eta}\right] dt$$

Integrando por partes el segundo término (usando que $\eta$ se anula en los extremos):

$$\delta J = \int_{t_1}^{t_2} \left[\frac{\partial \mathcal{L}}{\partial q} - \frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{q}}\right]\eta(t)\, dt$$

La condición $\delta J = 0$ para toda variación admisible $\eta$ implica, por el **lema fundamental del cálculo de variaciones**, que el integrando debe anularse idénticamente.

**Referencia:** Gelfand & Fomin, §2–3 [8].

### 3.4. Deducción de la Ecuación de Euler-Lagrange desde el Cálculo de Variaciones

De la condición de estacionariedad $\delta J = 0$ se obtiene la **ecuación de Euler-Lagrange**:

$$\frac{\partial \mathcal{L}}{\partial q} - \frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{q}}\right) = 0$$

Para un sistema con $n$ coordenadas generalizadas $q_1, \ldots, q_n$, esta ecuación se escribe para cada coordenada:

$$\frac{\partial \mathcal{L}}{\partial q_i} - \frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{q}_i}\right) = 0, \qquad i = 1, \ldots, n$$

Es un resultado profundo: las ecuaciones de movimiento de la mecánica clásica son **condiciones necesarias de optimalidad** para la funcional de acción. Esto no es una suposición adicional a las leyes de Newton; se puede demostrar que ambas formulaciones son equivalentes para sistemas holonómicos.

**El lema fundamental del cálculo de variaciones** (también llamado lema de Du Bois-Reymond) establece que si $f$ es continua y $\int_a^b f(t)\eta(t)\,dt = 0$ para toda función suave $\eta$ con $\eta(a) = \eta(b) = 0$, entonces $f(t) = 0$ para todo $t \in [a, b]$. Este lema es el paso lógico clave que convierte la condición integral en una ecuación diferencial puntual.

**Referencia:** Gelfand & Fomin, §3 [8]; Goldstein, §2.2 [6].

### 3.5. El Principio de Hamilton (Principio de Mínima Acción)

El **principio de Hamilton** afirma que, de todas las trayectorias posibles $q(t)$ que conectan la configuración $q(t_1)$ con $q(t_2)$, la trayectoria físicamente realizada es aquella que hace estacionaria la **acción**:

$$S[q] = \int_{t_1}^{t_2} \mathcal{L}(q, \dot{q}, t)\, dt$$

Es decir:

$$\delta S = 0$$

Este principio unifica toda la mecánica clásica. Es más general que las leyes de Newton porque:

1. **Es independiente del sistema de coordenadas** — funciona igual en coordenadas cartesianas, polares, o cualesquiera generalizadas.
2. **Se extiende naturalmente a campos** — la formulación lagrangiana de la electrodinámica, la relatividad general y la mecánica cuántica parte de este mismo principio.
3. **Proporciona un marco natural para las restricciones** — las restricciones holonómicas se incorporan automáticamente al elegir coordenadas generalizadas adecuadas.

**Nota importante:** El término "mínima acción" es históricamente incorrecto en general. La trayectoria física hace la acción *estacionaria* (su primera variación es cero), pero no necesariamente mínima. Puede ser un mínimo, un máximo, o un punto de silla. En la mayoría de los casos mecánicos simples (trayectorias suficientemente cortas), sí es un mínimo.

**Referencia:** Goldstein, §2.1 [6]; Landau & Lifshitz, §2 [7]; Gelfand & Fomin, §4 [8].

### 3.6. Condiciones de Contorno y Variaciones con Extremos Libres

En la derivación estándar, se exige que $\eta(t_1) = \eta(t_2) = 0$ (extremos fijos). Pero en muchos problemas, uno o ambos extremos son libres. Al integrar por partes sin suponer que $\eta$ se anula en los extremos, aparecen **términos de frontera**:

$$\delta J = \left[\frac{\partial \mathcal{L}}{\partial \dot{q}}\eta\right]_{t_1}^{t_2} + \int_{t_1}^{t_2}\left[\frac{\partial \mathcal{L}}{\partial q} - \frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{q}}\right]\eta\,dt$$

Si un extremo es libre, la condición de estacionariedad impone la **condición natural de frontera** (o condición de transversalidad):

$$\left.\frac{\partial \mathcal{L}}{\partial \dot{q}}\right|_{t = t_{\text{libre}}} = 0$$

Estas condiciones son relevantes en problemas de control óptimo, donde la trayectoria no está completamente prescrita en los extremos.

**Referencia:** Gelfand & Fomin, §5 [8].

### 3.7. Restricciones y Multiplicadores de Lagrange en el Cálculo de Variaciones

Cuando el sistema está sujeto a restricciones de la forma $g(q, \dot{q}, t) = 0$, se modifica la funcional introduciendo un multiplicador de Lagrange $\lambda(t)$:

$$J^*[q, \lambda] = \int_{t_1}^{t_2}\left[\mathcal{L}(q, \dot{q}, t) + \lambda(t)\,g(q, \dot{q}, t)\right]dt$$

La condición de estacionariedad de $J^*$ con respecto a $q$ y a $\lambda$ simultáneamente genera las ecuaciones de Euler-Lagrange modificadas más la restricción como ecuación adicional.

En el contexto del péndulo invertido, la restricción de que la masa del péndulo está unida al carro por una varilla rígida de longitud $l$ es una **restricción holonómica** que se elimina al elegir $(x, \theta)$ como coordenadas generalizadas. No hace falta usar multiplicadores en este caso.

**Referencia:** Gelfand & Fomin, §12 [8]; Goldstein, §2.4 [6].

### 3.8. Segunda Variación y Condiciones Suficientes

La **segunda variación** determina si un extremal es un mínimo, máximo o punto de silla:

$$\delta^2 J = \left.\frac{d^2}{d\epsilon^2}\right|_{\epsilon=0} J[q^* + \epsilon\eta]$$

La condición $\delta^2 J > 0$ para toda variación admisible no nula es necesaria para un mínimo. Esta condición se relaciona con la **ecuación accesoria de Jacobi** y con la existencia de **puntos conjugados** a lo largo del extremal.

En el contexto de control, la segunda variación está conectada con las condiciones de suficiencia del **principio del máximo de Pontryagin** y con la formulación LQR (regulador cuadrático lineal), donde la funcional a minimizar es cuadrática en el estado y el control.

**Referencia:** Gelfand & Fomin, §24–26 [8].

### 3.9. Conexión con el Control Óptimo

El control óptimo puede verse como una extensión del cálculo de variaciones donde la funcional a optimizar depende también de una función de control $u(t)$:

$$J[x, u] = \int_{t_0}^{t_f} L(x(t), u(t), t)\, dt$$

sujeto a la restricción dinámica $\dot{x} = f(x, u, t)$.

En el problema LQR del péndulo invertido, la funcional de costo es:

$$J = \int_0^\infty \left(x^T Q\, x + u^T R\, u\right) dt$$

donde $Q \geq 0$ y $R > 0$ son matrices de peso. La ecuación de Euler-Lagrange de este problema variacional conduce, tras manipulaciones algebraicas, a la **ecuación algebraica de Riccati**:

$$A^T P + PA - PBR^{-1}B^T P + Q = 0$$

cuya solución $P$ determina la ganancia óptima $K = R^{-1}B^T P$. Esta es la conexión directa entre el cálculo de variaciones y el diseño de controladores modernos.

**Referencia:** Ogata, Cap. 10, Ej. 10-5 [1]; Kirk, *Optimal Control Theory*, Cap. 4 [9].

---

## 4. Ecuaciones Diferenciales Ordinarias

### 4.1. Clasificación y Conceptos Fundamentales

Una **ecuación diferencial ordinaria** (EDO) de orden $n$ es una ecuación de la forma:

$$F\left(t, u, \frac{du}{dt}, \frac{d^2u}{dt^2}, \ldots, \frac{d^nu}{dt^n}\right) = 0$$

Las ecuaciones del péndulo invertido son EDOs de segundo orden no lineales acopladas. Los conceptos clave son:

- **Orden**: el de la derivada más alta presente. Nuestro sistema es de orden 2.
- **Linealidad**: una EDO es lineal si es lineal en la función incógnita y sus derivadas. Las ecuaciones del péndulo son no lineales (contienen $\sin\theta$, $\cos\theta$, $\dot{\theta}^2\sin\theta$).
- **Autonomía**: una EDO es autónoma si no depende explícitamente de $t$. Nuestro sistema es autónomo cuando $u$ es constante.

**Referencia:** Olver & Shakiban, §7.4 [2]; Boyce & DiPrima, *Elementary Differential Equations*, Cap. 2–3 [10].

### 4.2. EDOs Lineales de Segundo Orden con Coeficientes Constantes

La forma general es:

$$a\ddot{u} + b\dot{u} + cu = f(t)$$

El **ansatz exponencial** $u = e^{\lambda t}$ conduce a la **ecuación característica**:

$$a\lambda^2 + b\lambda + c = 0$$

Las raíces $\lambda_{1,2}$ determinan el comportamiento cualitativo de las soluciones:

- **Raíces reales distintas** $\lambda_1 \neq \lambda_2$: solución general $u(t) = c_1 e^{\lambda_1 t} + c_2 e^{\lambda_2 t}$.
- **Raíz doble** $\lambda_1 = \lambda_2 = \lambda$: solución general $u(t) = (c_1 + c_2 t)e^{\lambda t}$.
- **Raíces complejas conjugadas** $\lambda = \mu \pm i\nu$: solución general $u(t) = e^{\mu t}(c_1\cos\nu t + c_2\sin\nu t)$.

Olver (§7.4, Ec. 7.50–7.53) desarrolla este procedimiento en detalle para el operador $L = D^2 - 2D - 3$, obteniendo las soluciones $e^{3x}$ y $e^{-x}$ a partir de las raíces de la ecuación característica $\lambda^2 - 2\lambda - 3 = 0$ [2].

**Referencia:** Olver & Shakiban, §7.4, Teorema 7.34 [2].

### 4.3. El Problema de Valor Inicial y Existencia y Unicidad

Un **problema de valor inicial** (PVI) consiste en una EDO junto con condiciones iniciales:

$$\frac{d\mathbf{u}}{dt} = \mathbf{f}(\mathbf{u}, t), \qquad \mathbf{u}(t_0) = \mathbf{u}_0$$

El **teorema de Picard-Lindelöf** (o de existencia y unicidad) garantiza que si $\mathbf{f}$ es continua y Lipschitz-continua en $\mathbf{u}$, entonces existe una única solución local. Para sistemas lineales, la solución es global (existe para todo $t$).

Para el sistema del péndulo invertido, el PVI completo requiere especificar cuatro condiciones iniciales: $x(0)$, $\dot{x}(0)$, $\theta(0)$, $\dot{\theta}(0)$. Esto concuerda con el hecho de que el sistema es equivalente a cuatro EDOs de primer orden.

**Referencia:** Boyce & DiPrima, §2.4 [10]; Olver & Shakiban, Teorema 7.34 [2].

### 4.4. Reducción de Orden: Conversión a Sistema de Primer Orden

Toda EDO de orden $n$ puede convertirse en un sistema de $n$ EDOs de primer orden mediante la introducción de nuevas variables. Para una EDO de segundo orden $\ddot{u} = g(u, \dot{u}, t)$, definimos:

$$u_1 = u, \qquad u_2 = \dot{u}$$

El sistema equivalente es:

$$\dot{u}_1 = u_2, \qquad \dot{u}_2 = g(u_1, u_2, t)$$

Para el péndulo invertido, con $\mathbf{x} = (\theta, \dot{\theta}, x, \dot{x})^T$, las dos EDOs de segundo orden se convierten en cuatro de primer orden. Como señala Olver (§10.1), esta es la **construcción del plano de fase** [2], donde el sistema de segundo orden se reduce a primer orden al costo de duplicar la dimensión.

**Referencia:** Olver & Shakiban, §10.1, Ec. (10.8) [2]; Ogata, §2-4 [1].

### 4.5. Superposición y Espacio de Soluciones

Para un sistema lineal homogéneo $L[\mathbf{u}] = 0$, el **principio de superposición** (Teorema 7.30 de Olver [2]) garantiza que el espacio de soluciones es un **subespacio vectorial**: toda combinación lineal de soluciones es otra solución.

$$L[c_1\mathbf{u}_1 + c_2\mathbf{u}_2] = c_1 L[\mathbf{u}_1] + c_2 L[\mathbf{u}_2] = 0$$

Para una EDO lineal de orden $n$, el espacio de soluciones del sistema homogéneo tiene dimensión exactamente $n$ (Teorema 7.34 de Olver) [2]. Por tanto, una vez halladas $n$ soluciones linealmente independientes, toda solución es combinación lineal de ellas.

Para el sistema **no homogéneo** $L[\mathbf{u}] = \mathbf{f}(t)$, la solución general tiene la forma:

$$\mathbf{u}(t) = \mathbf{u}_p(t) + \mathbf{z}(t)$$

donde $\mathbf{u}_p$ es una solución particular y $\mathbf{z}$ es la solución general del sistema homogéneo asociado.

**Referencia:** Olver & Shakiban, §7.4, Teorema 7.30 y Teorema 7.34 [2].

---

## 5. Sistemas de EDOs Lineales y Representación en Espacio de Estados

### 5.1. Sistemas de Primer Orden

Un sistema lineal autónomo de primer orden tiene la forma:

$$\frac{d\mathbf{u}}{dt} = A\,\mathbf{u}$$

donde $A$ es una matriz $n \times n$ constante y $\mathbf{u}(t) \in \mathbb{R}^n$. El ansatz $\mathbf{u}(t) = e^{\lambda t}\mathbf{v}$ conduce a la ecuación de eigenvalores:

$$A\mathbf{v} = \lambda\mathbf{v}$$

Cada eigenpar $(\lambda_i, \mathbf{v}_i)$ produce una solución $\mathbf{u}_i(t) = e^{\lambda_i t}\mathbf{v}_i$, y la solución general (cuando $A$ es diagonalizable) es:

$$\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{v}_1 + c_2 e^{\lambda_2 t}\mathbf{v}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{v}_n$$

**Referencia:** Olver & Shakiban, §10.1 y §10.2 [2].

### 5.2. La Representación en Espacio de Estados

La **representación en espacio de estados** es la forma estándar para sistemas de control. Ogata (§2-4) la define como [1]:

$$\dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u}$$
$$\mathbf{y} = C\mathbf{x} + D\mathbf{u}$$

donde:
- $\mathbf{x} \in \mathbb{R}^n$ es el **vector de estado** (contiene toda la información necesaria para predecir el comportamiento futuro del sistema).
- $\mathbf{u} \in \mathbb{R}^m$ es el **vector de entrada** (señales de control).
- $\mathbf{y} \in \mathbb{R}^p$ es el **vector de salida** (variables medidas).
- $A \in \mathbb{R}^{n \times n}$ es la **matriz de estado** (dinámica interna).
- $B \in \mathbb{R}^{n \times m}$ es la **matriz de entrada** (cómo el control afecta al estado).
- $C \in \mathbb{R}^{p \times n}$ es la **matriz de salida** (qué se mide).
- $D \in \mathbb{R}^{p \times m}$ es la **matriz de transmisión directa** (usualmente cero).

**Referencia:** Ogata, §2-4 [1].

### 5.3. Espacio de Estados del Péndulo Invertido

Definiendo el vector de estado $\mathbf{x} = (\theta,\; \dot{\theta},\; x,\; \dot{x})^T$ y linealizando alrededor de $\theta = 0$ (las aproximaciones $\sin\theta \approx \theta$, $\cos\theta \approx 1$, $\dot{\theta}^2 \approx 0$ son válidas para ángulos pequeños), Ogata (Cap. 10, Ej. 10-5) obtiene las matrices numéricas para $M = 2$ kg, $m = 0.1$ kg, $l = 0.5$ m [1]:

$$A = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 20.601 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ -0.4905 & 0 & 0 & 0 \end{pmatrix}, \quad B = \begin{pmatrix} 0 \\ -1 \\ 0 \\ 0.5 \end{pmatrix}, \quad C = \begin{pmatrix} 0 & 0 & 1 & 0 \end{pmatrix}$$

La presencia del elemento $a_{21} = 20.601 > 0$ en la matriz $A$ refleja la inestabilidad inherente del péndulo invertido: un desplazamiento angular produce una aceleración angular en la misma dirección.

**Referencia:** Ogata, Cap. 10, Ej. 10-5 [1].

### 5.4. La Función de Transferencia

La relación entre la representación en espacio de estados y la función de transferencia se obtiene aplicando la transformada de Laplace (asumiendo condiciones iniciales nulas). Ogata (§2-4) demuestra que [1]:

$$G(s) = C(sI - A)^{-1}B + D$$

La función de transferencia contiene la misma información dinámica que la representación en espacio de estados, pero en el dominio de la frecuencia. Los **polos** de $G(s)$ coinciden con los eigenvalores de $A$ (bajo condiciones de controlabilidad y observabilidad), y determinan la estabilidad del sistema.

**Referencia:** Ogata, §2-4 [1].

---

## 6. Estabilidad y Retratos de Fase

### 6.1. Definición Formal de Estabilidad

La noción de estabilidad se formaliza mediante las definiciones de Lyapunov. Para el sistema $\dot{\mathbf{u}} = A\mathbf{u}$, el equilibrio $\mathbf{u} = 0$ es:

- **Estable**: si para todo $\epsilon > 0$ existe $\delta > 0$ tal que $\|\mathbf{u}(t_0)\| < \delta$ implica $\|\mathbf{u}(t)\| < \epsilon$ para todo $t > t_0$.
- **Asintóticamente estable**: si es estable y además $\lim_{t \to \infty} \mathbf{u}(t) = 0$.
- **Inestable**: si no es estable.

**Referencia:** Olver & Shakiban, §10.2 [2].

### 6.2. Criterio de Estabilidad por Eigenvalores

El resultado fundamental (Teorema 10.16 de Olver [2]) establece:

**Teorema (Estabilidad de sistemas lineales).** El sistema $\dot{\mathbf{u}} = A\mathbf{u}$ tiene equilibrio asintóticamente estable si y solo si **todos** los eigenvalores de $A$ tienen parte real negativa: $\text{Re}(\lambda_i) < 0$ para todo $i$.

Si al menos un eigenvalor tiene parte real positiva, el equilibrio es inestable.

Esto se demuestra observando que las soluciones contienen funciones de la forma $t^k e^{\mu t}\cos\nu t$ y $t^k e^{\mu t}\sin\nu t$, donde $\lambda = \mu + i\nu$ es un eigenvalor. Estas funciones tienden a cero cuando $\mu < 0$ (Lema 10.15 de Olver) [2], y divergen cuando $\mu > 0$.

**Para el péndulo invertido sin control**: la matriz $A$ tiene un eigenvalor con parte real positiva ($\lambda \approx +4.54$), confirmando la inestabilidad del equilibrio vertical superior. El controlador (LQR o por asignación de polos) modifica los eigenvalores de lazo cerrado $A - BK$ para que todos tengan parte real negativa.

**Referencia:** Olver & Shakiban, §10.2, Teorema 10.16 y Teorema 10.19 [2].

### 6.3. Clasificación de Retratos de Fase en 2D

Para un sistema planar $\dot{\mathbf{u}} = A\mathbf{u}$ con $A \in \mathbb{R}^{2 \times 2}$, la clasificación completa del retrato de fase depende del **traza** $\tau = \text{tr}(A)$, el **determinante** $\delta = \det(A)$, y el **discriminante** $\Delta = \tau^2 - 4\delta$ (Olver, §10.3, Proposición 10.22) [2]:

**Eigenvalores reales distintos** ($\Delta > 0$):
- **Nodo estable**: $\delta > 0$, $\tau < 0$ (ambos $\lambda < 0$). Soluciones convergen al origen.
- **Punto silla**: $\delta < 0$ (eigenvalores de signo opuesto). Equilibrio inestable.
- **Nodo inestable**: $\delta > 0$, $\tau > 0$ (ambos $\lambda > 0$). Soluciones divergen.

**Eigenvalores complejos conjugados** ($\Delta < 0$):
- **Foco estable**: $\tau < 0$. Soluciones espiralan hacia el origen.
- **Centro**: $\tau = 0$. Órbitas periódicas elípticas (estable, no asintóticamente).
- **Foco inestable**: $\tau > 0$. Soluciones espiralan hacia afuera.

**Eigenvalor doble** ($\Delta = 0$):
- **Nodo impropio estable/inestable**: según el signo de $\tau$.
- **Estrella estable/inestable**: cuando $A = \lambda I$.

La Figura 10.4 de Olver [2] muestra las regiones del plano $(\tau, \delta)$ correspondientes a cada tipo de retrato de fase. La frontera entre estabilidad e inestabilidad es la región $\tau = 0$, $\delta > 0$ (centros).

**Referencia:** Olver & Shakiban, §10.3, Figura 10.3 y 10.4, Proposición 10.22 [2].

### 6.4. Estabilidad Estructural

Un sistema es **estructuralmente estable** si su comportamiento cualitativo no cambia bajo perturbaciones pequeñas de la matriz $A$. Como explica Olver (§10.3) [2]:

- Los sistemas asintóticamente estables y los inestables son estructuralmente estables.
- Los centros y nodos impropios **no** son estructuralmente estables: una perturbación infinitesimal puede convertir un centro en un foco estable o inestable.

Esto es relevante para el péndulo invertido: el punto de equilibrio invertido es un punto silla en el plano de fase $(\theta, \dot{\theta})$, que es estructuralmente estable como tipo de equilibrio (sigue siendo silla bajo perturbaciones). La tarea del controlador no es cambiar la naturaleza del equilibrio sin control, sino crear un nuevo sistema dinámico (en lazo cerrado) donde el equilibrio sea un nodo o foco estable.

**Referencia:** Olver & Shakiban, §10.3 [2].

---

## 7. Linealización de Sistemas No Lineales

### 7.1. El Proceso de Linealización

Las ecuaciones del péndulo invertido son **no lineales**. Para aplicar las herramientas del álgebra lineal y la teoría de control lineal, es necesario linealizar alrededor de un punto de equilibrio.

Dado el sistema no lineal:

$$\dot{\mathbf{x}} = \mathbf{f}(\mathbf{x}, \mathbf{u})$$

con un punto de equilibrio $(\mathbf{x}_0, \mathbf{u}_0)$ tal que $\mathbf{f}(\mathbf{x}_0, \mathbf{u}_0) = 0$, la linealización consiste en una expansión de Taylor de primer orden:

$$\dot{\mathbf{x}} \approx \mathbf{f}(\mathbf{x}_0, \mathbf{u}_0) + \left.\frac{\partial \mathbf{f}}{\partial \mathbf{x}}\right|_{(\mathbf{x}_0, \mathbf{u}_0)} (\mathbf{x} - \mathbf{x}_0) + \left.\frac{\partial \mathbf{f}}{\partial \mathbf{u}}\right|_{(\mathbf{x}_0, \mathbf{u}_0)} (\mathbf{u} - \mathbf{u}_0)$$

Definiendo las desviaciones $\Delta\mathbf{x} = \mathbf{x} - \mathbf{x}_0$ y $\Delta\mathbf{u} = \mathbf{u} - \mathbf{u}_0$, y las matrices jacobianas:

$$A = \left.\frac{\partial \mathbf{f}}{\partial \mathbf{x}}\right|_{(\mathbf{x}_0, \mathbf{u}_0)}, \qquad B = \left.\frac{\partial \mathbf{f}}{\partial \mathbf{u}}\right|_{(\mathbf{x}_0, \mathbf{u}_0)}$$

el sistema linealizado es:

$$\Delta\dot{\mathbf{x}} = A\,\Delta\mathbf{x} + B\,\Delta\mathbf{u}$$

Ogata (§2-7) presenta este procedimiento como el método estándar para obtener modelos lineales de sistemas físicos no lineales [1].

**Referencia:** Ogata, §2-7 [1]; Olver & Shakiban, §10.4, p. 605 [2].

### 7.2. La Matriz Jacobiana

La **matriz jacobiana** $J = \partial\mathbf{f}/\partial\mathbf{x}$ es la generalización multivariable de la derivada. Para $\mathbf{f}: \mathbb{R}^n \to \mathbb{R}^n$, es la matriz $n \times n$ cuyo elemento $(i, j)$ es:

$$(J)_{ij} = \frac{\partial f_i}{\partial x_j}$$

Olver describe la jacobiana como "la linealización de $\mathbf{f}(\mathbf{u})$ en $\mathbf{u}_0$" y la denota $A = \mathbf{f}'(\mathbf{u}_0) = (\partial f_i/\partial u_j)$ (§10.4, p. 605) [2].

Para el péndulo invertido, la jacobiana evaluada en $(\theta, \dot{\theta}, x, \dot{x}) = (0, 0, 0, 0)$ se calcula diferenciando las ecuaciones no lineales. Las aproximaciones clave son:

$$\left.\frac{\partial}{\partial\theta}(\sin\theta)\right|_{\theta=0} = \cos(0) = 1, \qquad \left.\frac{\partial}{\partial\theta}(\cos\theta)\right|_{\theta=0} = -\sin(0) = 0$$

$$\left.\frac{\partial}{\partial\dot{\theta}}(\dot{\theta}^2\sin\theta)\right|_{(\theta,\dot{\theta})=(0,0)} = 0$$

Estas derivadas producen los coeficientes de la matriz $A$ linealizada.

**Referencia:** Olver & Shakiban, §10.4, p. 605 [2]; Ogata, §2-7 [1].

### 7.3. El Teorema de la Variedad Central (Center Manifold Theorem)

El **Teorema de la Variedad Central**, descrito por Olver (§10.4, p. 605) [2], extiende los resultados de estabilidad lineal al caso no lineal. Establece que cerca de un punto de equilibrio $\mathbf{u}_0$ del sistema no lineal $\dot{\mathbf{u}} = \mathbf{f}(\mathbf{u})$, la dinámica se organiza en tres **variedades invariantes**:

- **Variedad estable** $W^s$: tangente al subespacio estable de la linealización. Las soluciones sobre $W^s$ convergen al equilibrio exponencialmente.
- **Variedad inestable** $W^u$: tangente al subespacio inestable. Las soluciones se alejan del equilibrio.
- **Variedad central** $W^c$: tangente al subespacio central. La dinámica sobre $W^c$ depende de los términos no lineales y no puede determinarse solo con la linealización.

La implicación práctica para el péndulo invertido es que la **linealización es válida** para determinar la estabilidad cuando ningún eigenvalor tiene parte real cero. Dado que el equilibrio invertido tiene eigenvalores con partes reales estrictamente positivas y negativas (es un punto silla en la dinámica libre), la linealización captura correctamente la inestabilidad.

**Referencia:** Olver & Shakiban, §10.4, p. 605 [2]; Hirsch & Smale, *Differential Equations, Dynamical Systems, and Linear Algebra*, Cap. 9 [11].

### 7.4. Aproximación de Ángulo Pequeño

La linealización del péndulo utiliza las aproximaciones para $\theta$ pequeño:

$$\sin\theta \approx \theta, \qquad \cos\theta \approx 1, \qquad \dot{\theta}^2\sin\theta \approx 0$$

La primera es la **aproximación de Taylor de primer orden**: $\sin\theta = \theta - \theta^3/6 + \cdots \approx \theta$. La segunda descarta el término $\theta^2/2$ en $\cos\theta = 1 - \theta^2/2 + \cdots$. La tercera descarta productos de términos pequeños.

Estas aproximaciones transforman las ecuaciones no lineales del péndulo en un sistema lineal con coeficientes constantes, permitiendo la aplicación directa de todas las herramientas de álgebra lineal y teoría de control.

**Cuantificación del error:** Para $|\theta| < 0.1$ rad ($\approx 5.7°$), el error relativo de $\sin\theta \approx \theta$ es menor al 0.17%. Para $|\theta| < 0.3$ rad ($\approx 17°$), el error es menor al 1.5%. Esta es la región de validez del modelo linealizado.

**Referencia:** Ogata, Cap. 3, Ej. 3-6 [1].

---

## 8. Exponencial Matricial y Solución de Sistemas Lineales

### 8.1. Definición de la Exponencial Matricial

La **exponencial matricial** $e^{tA}$ es la generalización natural de la exponencial escalar al caso matricial. Se define mediante la serie de potencias:

$$e^{tA} = I + tA + \frac{t^2}{2!}A^2 + \frac{t^3}{3!}A^3 + \cdots = \sum_{k=0}^{\infty} \frac{t^k}{k!}A^k$$

Esta serie converge para toda matriz $A$ y todo $t$ real. Olver (§10.4) la introduce como la herramienta central para resolver sistemas lineales [2].

**Propiedades fundamentales:**

1. $e^{0 \cdot A} = I$ (la identidad).
2. $\frac{d}{dt}e^{tA} = Ae^{tA} = e^{tA}A$.
3. $e^{(s+t)A} = e^{sA}e^{tA}$.
4. $(e^{tA})^{-1} = e^{-tA}$.
5. Si $A$ es diagonalizable con $A = PDP^{-1}$, entonces $e^{tA} = Pe^{tD}P^{-1}$, donde $e^{tD} = \text{diag}(e^{\lambda_1 t}, \ldots, e^{\lambda_n t})$.

**Referencia:** Olver & Shakiban, §10.4 [2].

### 8.2. Solución del Sistema Homogéneo

La solución del PVI $\dot{\mathbf{u}} = A\mathbf{u}$, $\mathbf{u}(t_0) = \mathbf{u}_0$ es:

$$\mathbf{u}(t) = e^{(t-t_0)A}\mathbf{u}_0$$

La exponencial matricial $e^{tA}$ actúa como la **matriz fundamental** del sistema: sus columnas forman un conjunto completo de soluciones linealmente independientes.

Cuando $A$ es diagonalizable con eigenvalores $\lambda_1, \ldots, \lambda_n$ y eigenvectores $\mathbf{v}_1, \ldots, \mathbf{v}_n$, la solución se descompone en **modos**:

$$\mathbf{u}(t) = c_1 e^{\lambda_1 t}\mathbf{v}_1 + c_2 e^{\lambda_2 t}\mathbf{v}_2 + \cdots + c_n e^{\lambda_n t}\mathbf{v}_n$$

Cada término $e^{\lambda_i t}\mathbf{v}_i$ es una **eigensolución** o modo normal del sistema. Los coeficientes $c_i$ se determinan por las condiciones iniciales.

**Referencia:** Olver & Shakiban, §10.4, Ec. (10.38–10.39) [2].

### 8.3. Solución del Sistema Inhomogéneo: Variación de Parámetros

Para el sistema forzado $\dot{\mathbf{u}} = A\mathbf{u} + \mathbf{f}(t)$, Olver (§10.4, Ec. 10.57–10.58) presenta el método de **variación de parámetros** [2]:

Se propone $\mathbf{u}(t) = e^{tA}\mathbf{v}(t)$, lo que conduce a:

$$\frac{d\mathbf{v}}{dt} = e^{-tA}\mathbf{f}(t)$$

Integrando:

$$\mathbf{u}(t) = e^{(t-t_0)A}\mathbf{u}_0 + \int_{t_0}^{t} e^{(t-s)A}\mathbf{f}(s)\,ds$$

El primer término es la **respuesta libre** (dependiente de las condiciones iniciales) y la integral es la **respuesta forzada** (debida a la entrada externa). Esta fórmula es la base de la **integral de convolución** en la teoría de sistemas lineales.

**Referencia:** Olver & Shakiban, §10.4, Ec. (10.56–10.58) [2].

---

## 9. Sistemas Forzados y Resonancia

### 9.1. Forzamiento Periódico de Sistemas de Segundo Orden

Olver (§10.6) analiza en profundidad el forzamiento periódico de sistemas masa-resorte, que es el modelo fundamental para entender la dinámica del péndulo [2]. Para un sistema escalar sin amortiguamiento:

$$m\ddot{u} + ku = \alpha\cos\gamma t$$

la solución depende de la relación entre la **frecuencia natural** $\omega = \sqrt{k/m}$ y la **frecuencia de forzamiento** $\gamma$. Cuando $\gamma \neq \omega$, existe una solución particular (Olver, Ec. 10.100):

$$u_p(t) = \frac{\alpha}{m(\omega^2 - \gamma^2)}\cos\gamma t$$

La solución general es **cuasi-periódica**, combinación de la vibración natural y la respuesta al forzamiento (Olver, Ec. 10.101) [2].

### 9.2. Resonancia

Cuando $\gamma = \omega$ (forzamiento a la frecuencia natural), el ansatz trigonométrico estándar falla. La solución involucra un término que **crece linealmente** con el tiempo (Olver, Ec. 10.105) [2]:

$$u(t) = \frac{\alpha}{2m\omega}t\sin\omega t + r\cos(\omega t - \delta)$$

Este crecimiento ilimitado de la amplitud es el fenómeno de **resonancia**. La amplitud crece como $t$, lo que eventualmente lleva al sistema fuera del régimen lineal.

### 9.3. Efecto del Amortiguamiento

Con amortiguamiento viscoso $\beta\dot{u}$, la ecuación forzada se convierte en:

$$m\ddot{u} + \beta\dot{u} + ku = \alpha\cos\gamma t$$

La solución particular persistente tiene la forma (Olver, Ec. 10.107) [2]:

$$u_p(t) = \frac{\alpha}{\sqrt{m^2(\omega^2 - \gamma^2)^2 + \beta^2\gamma^2}}\cos(\gamma t - \varepsilon)$$

donde $\varepsilon$ es un desfase. La amplitud es máxima cuando $\gamma = \omega$ (resonancia), donde vale $\alpha/(\beta\omega)$. El amortiguamiento limita la respuesta resonante a un valor finito, pero potencialmente grande si $\beta$ es pequeño.

La solución general (Olver, Ec. 10.109) [2] incluye un **transitorio** $re^{-\mu t}\cos(\nu t - \delta)$ que decae exponencialmente. A largo plazo, solo persiste la respuesta forzada.

### 9.4. Batimientos (Beats)

Cuando $\gamma \approx \omega$ pero $\gamma \neq \omega$, la solución presenta **batimientos** (Olver, Ec. 10.103) [2]: oscilaciones rápidas cuya amplitud varía lentamente con frecuencia $|\omega - \gamma|/2$. La solución puede escribirse como:

$$u(t) = \frac{2\alpha}{m(\omega^2 - \gamma^2)}\sin\left(\frac{\omega + \gamma}{2}t\right)\sin\left(\frac{\omega - \gamma}{2}t\right)$$

El primer factor es una oscilación rápida (frecuencia promedio), modulada por una envolvente lenta (segundo factor).

### 9.5. Resonancia en Sistemas Multidimensionales

Para un sistema matricial forzado (Olver, Ec. 10.113) [2]:

$$M\ddot{\mathbf{u}} + K\mathbf{u} = \cos(\gamma t)\,\mathbf{a}$$

la resonancia ocurre cuando $\gamma^2$ coincide con un eigenvalor generalizado del par $(K, M)$. El ansatz $\mathbf{u}(t) = \cos(\gamma t)\,\mathbf{w}$ conduce al sistema algebraico $(K - \gamma^2 M)\mathbf{w} = \mathbf{a}$ (Olver, Ec. 10.115) [2]. Si $\gamma^2$ es un eigenvalor, el sistema tiene solución si y solo si $\mathbf{a}$ es ortogonal al eigenespacio correspondiente (por la alternativa de Fredholm).

**Referencia para toda esta sección:** Olver & Shakiban, §10.6, Ecuaciones (10.95–10.115) [2].

---

## 10. Energía, Amortiguamiento y Disipación

### 10.1. Sistemas Conservativos y la Conservación de Energía

Para un sistema conservativo (sin fricción ni control), la energía total $E = T + V$ se conserva:

$$\frac{dE}{dt} = 0$$

Para el péndulo invertido libre (sin control, sin fricción):

$$E = \frac{1}{2}(M+m)\dot{x}^2 + ml\dot{x}\dot{\theta}\cos\theta + \frac{1}{2}ml^2\dot{\theta}^2 + mgl\cos\theta$$

La conservación de $E$ implica que las trayectorias en el espacio de fases se mueven sobre **superficies de energía constante**. Para un péndulo simple (sin carro), esto produce las conocidas curvas de nivel del retrato de fase: libraciones (oscilaciones alrededor del equilibrio inferior) y rotaciones (el péndulo da vueltas completas), separadas por las **separatrices** que pasan por el punto silla.

### 10.2. Sistemas Disipativos

El amortiguamiento (fricción viscosa $b\dot{x}$ en el carro, $\beta\dot{\theta}$ en el pivote) hace que la energía decrezca:

$$\frac{dE}{dt} = -b\dot{x}^2 - \beta\dot{\theta}^2 \leq 0$$

La función $E$ actúa como una **función de Lyapunov**: su decrecimiento a lo largo de las trayectorias garantiza la estabilidad del equilibrio inferior del péndulo (donde $E$ tiene un mínimo).

### 10.3. Sistemas Hamiltonianos

Los **sistemas hamiltonianos** son una clase especial de sistemas conservativos donde las ecuaciones de movimiento se escriben en forma canónica:

$$\dot{q}_i = \frac{\partial \mathcal{H}}{\partial p_i}, \qquad \dot{p}_i = -\frac{\partial \mathcal{H}}{\partial q_i}$$

Esta formulación preserva la estructura simpléctica del espacio de fases, lo que tiene consecuencias profundas: los flujos hamiltonianos conservan volumen en el espacio de fases (teorema de Liouville), y los centros del sistema linealizado corresponden a verdaderos centros del sistema no lineal (no se convierten en focos bajo perturbación hamiltoniana).

El péndulo invertido sin control es un sistema hamiltoniano. Al agregar control y/o amortiguamiento, deja de serlo, lo que permite la estabilización asintótica (imposible en sistemas hamiltonianos, donde los centros no pueden convertirse en nodos estables).

**Referencia:** Olver & Shakiban, §10.3, Ec. (10.25) [2]; Goldstein, Cap. 8 [6].

### 10.4. Oscilaciones Amortiguadas: Tipos de Amortiguamiento

Para el sistema escalar $m\ddot{u} + \beta\dot{u} + ku = 0$, la ecuación característica es:

$$m\lambda^2 + \beta\lambda + k = 0 \implies \lambda = \frac{-\beta \pm \sqrt{\beta^2 - 4mk}}{2m}$$

Dependiendo del discriminante $\beta^2 - 4mk$:

- **Subamortiguado** ($\beta^2 < 4mk$): eigenvalores complejos conjugados con parte real negativa. Solución: $u(t) = re^{-\mu t}\cos(\nu t - \delta)$. Oscilaciones decrecientes.
- **Críticamente amortiguado** ($\beta^2 = 4mk$): eigenvalor doble real negativo. Retorno al equilibrio sin oscilación, lo más rápido posible.
- **Sobreamortiguado** ($\beta^2 > 4mk$): dos eigenvalores reales negativos distintos. Retorno lento al equilibrio sin oscilación.

En el contexto del control del péndulo invertido, el diseñador elige los polos del sistema en lazo cerrado (eigenvalores de $A - BK$) de modo que el sistema sea **subamortiguado** (respuesta rápida con oscilaciones leves) o **críticamente amortiguado** (sin sobreimpulso), según los requisitos de rendimiento.

**Referencia:** Olver & Shakiban, §10.5 [2]; Ogata, Cap. 5 [1].

---

## 11. Síntesis: De la Física al Control del Péndulo Invertido

### 11.1. El Camino Completo

El desarrollo del modelo del péndulo invertido sigue una cadena lógica precisa que conecta todas las disciplinas cubiertas en este documento:

1. **Física (Mecánica Newtoniana o Lagrangiana)** $\longrightarrow$ Se obtienen las ecuaciones de movimiento no lineales del sistema carro-péndulo. La formulación lagrangiana (basada en el cálculo de variaciones) es más elegante y evita las fuerzas de reacción.

2. **Cálculo de variaciones** $\longrightarrow$ Proporciona la justificación teórica profunda de las ecuaciones de Euler-Lagrange como condiciones de estacionariedad del funcional de acción. Además, el control óptimo LQR es en sí mismo un problema variacional.

3. **Linealización (cálculo multivariable)** $\longrightarrow$ La expansión de Taylor/jacobiana transforma el sistema no lineal en un sistema lineal con coeficientes constantes, válido cerca del equilibrio.

4. **Ecuaciones diferenciales** $\longrightarrow$ La solución del sistema linealizado se expresa en términos de eigenvalores y eigenvectores. La exponencial matricial $e^{tA}$ es la solución fundamental.

5. **Álgebra lineal** $\longrightarrow$ Los eigenvalores determinan la estabilidad, la controlabilidad y la observabilidad del sistema. El diseño del controlador (asignación de polos, LQR) es un problema algebraico.

### 11.2. Dos Formulaciones, Una Realidad

El proyecto compara dos caminos hacia las ecuaciones de movimiento:

| Aspecto | Newton ($F = ma$) | Lagrange ($\delta S = 0$) |
|---------|-------------------|---------------------------|
| Cantidades fundamentales | Fuerzas (vectores) | Energías (escalares) |
| Tratamiento de restricciones | Fuerzas de reacción explícitas | Eliminadas por coord. generalizadas |
| Invariancia de coordenadas | No | Sí |
| Herramienta matemática | Análisis vectorial | Cálculo de variaciones |
| Conexión con control óptimo | Indirecta | Directa (LQR $\leftrightarrow$ Euler-Lagrange) |
| Resultado final | Idénticas ecuaciones de movimiento | Idénticas ecuaciones de movimiento |

Ambos caminos producen las mismas ecuaciones diferenciales, pero el formalismo lagrangiano ofrece una conexión más profunda con los principios variacionales que fundamentan el control óptimo.

### 11.3. Diagrama Conceptual

```
CÁLCULO DE VARIACIONES
       │
       ▼
  Principio de Hamilton (δS = 0)
       │
       ▼
  Ecuaciones de Euler-Lagrange ◄──── Lagrangiano L = T - V
       │                                      ▲
       ▼                                      │
  EDOs no lineales acopladas ◄──── También derivables por Newton (F=ma)
       │
       ▼
  Linealización (Jacobiana, Taylor)
       │
       ▼
  Sistema lineal: ẋ = Ax + Bu ◄──── ÁLGEBRA LINEAL
       │                              (eigenvalores, rango,
       ▼                               espacios fundamentales)
  Análisis de estabilidad
  (eigenvalores de A en semiplano derecho → inestable)
       │
       ▼
  Diseño de controlador ◄──── Control óptimo LQR
  (K tal que A-BK estable)        (otro problema variacional:
       │                           min ∫(x'Qx + u'Ru)dt)
       ▼                                │
  Sistema estabilizado               Ecuación de Riccati
  (todos los eigenvalores              (A'P + PA - PBR⁻¹B'P + Q = 0)
   en semiplano izquierdo)
```

---

## 12. Referencias

[1] **Ogata, K.** *Ingeniería de Control Moderna*, 5ª edición. Pearson Educación, 2010. — Capítulos 2 (espacio de estados, función de transferencia, linealización §2-7), 3 (modelado del péndulo invertido, Ej. 3-5 y 3-6), 9 (controlabilidad §9-6, observabilidad §9-7), 10 (diseño de controladores, asignación de polos, Ej. 10-5).

[2] **Olver, P. J. & Shakiban, C.** *Applied Linear Algebra*, 2ª edición. Springer, 2018. — Capítulo 6 (equilibrio: masas-resortes, redes eléctricas, estructuras), Capítulo 7 §7.4 (operadores diferenciales lineales, ecuaciones características, Teorema 7.34), Capítulo 8 (eigenvalores, eigenvectores, SVD, PCA), Capítulo 10 (dinámica: §10.1 técnicas de solución, §10.2 estabilidad Teorema 10.16, §10.3 clasificación de retratos de fase 2D, §10.4 exponencial matricial y teorema de la variedad central, §10.5 dinámica de estructuras y modos normales, §10.6 forzamiento y resonancia).

[3] **Lauwens, B. & Downey, A. B.** *Introducción a la Programación en Julia* (traducción de ThinkJulia). — Referencia para implementación computacional.

[4] **Strogatz, S. H.** *Nonlinear Dynamics and Chaos*, 2ª edición. Westview Press, 2015. — Referencia complementaria sobre sistemas no lineales, bifurcaciones y retratos de fase.

[5] **Hirsch, M. W., Smale, S. & Devaney, R. L.** *Differential Equations, Dynamical Systems, and an Introduction to Chaos*, 3ª edición. Academic Press, 2012. — Referencia complementaria sobre el teorema de la variedad central y estabilidad de sistemas no lineales.

[6] **Goldstein, H., Poole, C. & Safko, J.** *Classical Mechanics*, 3ª edición. Addison-Wesley, 2002. — Referencia fundamental sobre mecánica lagrangiana (Cap. 1–2), principio de Hamilton (Cap. 2), mecánica hamiltoniana (Cap. 8), y momento de inercia (Cap. 5).

[7] **Landau, L. D. & Lifshitz, E. M.** *Mechanics* (Course of Theoretical Physics, Vol. 1), 3ª edición. Butterworth-Heinemann, 1976. — Referencia clásica sobre el principio de mínima acción (§1–2), coordenadas generalizadas y simetrías.

[8] **Gelfand, I. M. & Fomin, S. V.** *Calculus of Variations*. Dover Publications, 2000 (reimpresión del original de 1963). — Referencia fundamental sobre funcionales (Cap. 1), ecuación de Euler-Lagrange (Cap. 2–3), condiciones de frontera (Cap. 5), restricciones y multiplicadores de Lagrange (Cap. 12), segunda variación y condiciones suficientes (Cap. 24–26).

[9] **Kirk, D. E.** *Optimal Control Theory: An Introduction*. Dover Publications, 2004. — Referencia sobre control óptimo, principio del máximo de Pontryagin, LQR y su conexión con el cálculo de variaciones.

[10] **Boyce, W. E. & DiPrima, R. C.** *Elementary Differential Equations and Boundary Value Problems*, 11ª edición. Wiley, 2017. — Referencia complementaria sobre EDOs lineales, existencia y unicidad, y métodos de solución.

[11] **Hirsch, M. W. & Smale, S.** *Differential Equations, Dynamical Systems, and Linear Algebra*. Academic Press, 1974. — Referenciado por Olver como [41] en su bibliografía. Tratamiento riguroso de la relación entre sistemas dinámicos y álgebra lineal.

---

*Documento preparado como parte del proyecto de investigación sobre el péndulo invertido para el curso de Álgebra Lineal Aplicada. Los contenidos de las secciones 6, 8 y 9 se basan extensamente en Olver & Shakiban [2], Capítulo 10. Los contenidos de la sección 5 integran material de Ogata [1] y Olver [2]. Las secciones 2 y 3 se fundamentan en Goldstein [6], Landau & Lifshitz [7] y Gelfand & Fomin [8]. Todas las referencias a secciones y ecuaciones específicas de los textos han sido verificadas contra los documentos del proyecto.*
