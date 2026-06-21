# El Péndulo Invertido: Modelado, Control y Álgebra Lineal Aplicada

---

## Introducción

El péndulo invertido sobre un carro es uno de los sistemas dinámicos más estudiados en ingeniería y física aplicada. Este sistema —una barra rígida balanceándose en posición vertical sobre un carro móvil— es inherentemente inestable, y constituye el banco de pruebas canónico donde convergen la mecánica lagrangiana, la teoría de control moderna y el álgebra lineal.

Su relevancia radica en que captura, en un modelo relativamente simple, los desafíos fundamentales de estabilizar sistemas no lineales subactuados. Aparece como modelo subyacente en cohetes, vehículos autobalanceados (como el Segway), robots bípedos y en el control postural humano.

Este informe desarrolla la derivación completa de las ecuaciones de movimiento mediante Euler-Lagrange, su linealización en espacio de estados, el diseño de controladores y las conexiones profundas con conceptos de álgebra lineal. La implementación computacional se realizará en Julia.

---

## 1. Anatomía del Sistema: Carro, Péndulo y Fuerzas

El sistema carro-péndulo invertido consta de tres elementos fundamentales:

- Un **carro de masa** $M$ que se desplaza sobre un riel horizontal restringido al eje $x$, impulsado por una fuerza externa $F$ (la entrada de control).
- Una **barra rígida de masa** $m$ **y longitud** $L$, articulada al carro mediante una junta de revolución, cuyo centro de masa se ubica a una distancia $l = L/2$ del pivote.
- El sistema posee **dos grados de libertad** (posición del carro $x$ y ángulo del péndulo $\theta$) pero solo **una entrada de control** ($F$), lo que lo clasifica como un sistema *subactuado*.

### 1.1 Variables de estado

Las cuatro variables que describen completamente la configuración dinámica del sistema son:

| Variable       | Descripción                              |
|:--------------:|:-----------------------------------------|
| $x$            | Posición horizontal del carro            |
| $\dot{x}$      | Velocidad lineal del carro               |
| $\theta$       | Ángulo del péndulo respecto a la vertical ($\theta = 0$ es la posición invertida, hacia arriba) |
| $\dot{\theta}$  | Velocidad angular del péndulo            |

### 1.2 Parámetros del sistema

Los parámetros físicos del modelo son:

| Parámetro | Descripción                                    | Valor típico      |
|:---------:|:-----------------------------------------------|:-------------------|
| $M$       | Masa del carro                                 | $0.5$ kg           |
| $m$       | Masa del péndulo                               | $0.2$ kg           |
| $L$       | Longitud total del péndulo                     | $0.6$ m            |
| $l$       | Distancia del pivote al centro de masa ($L/2$) | $0.3$ m            |
| $g$       | Aceleración gravitacional                      | $9.81$ m/s²        |
| $b$       | Coeficiente de fricción viscosa del carro      | $0.1$ N·s/m        |
| $I_{cm}$  | Momento de inercia respecto al centro de masa  | $mL^2/12$          |
| $I_{piv}$ | Momento de inercia respecto al pivote          | $mL^2/3$           |

### 1.3 Diagrama de cuerpo libre

En el diagrama de cuerpo libre:

- El **carro** experimenta: la fuerza aplicada $F$, la fricción $-b\dot{x}$, su peso $Mg$ (compensado por la normal del riel), y las fuerzas de reacción del pivote $P_x$ (horizontal) y $P_y$ (vertical).
- El **péndulo** experimenta: su peso $mg$ en el centro de masa, y las reacciones $-P_x$ y $-P_y$ en el punto de articulación (por la tercera ley de Newton).

### 1.4 El equilibrio inestable

El equilibrio que interesa estabilizar es el **equilibrio inestable superior** ($\theta = 0$, péndulo apuntando hacia arriba). Cualquier perturbación infinitesimal causa que el péndulo caiga exponencialmente, y es precisamente este comportamiento el que el controlador debe contrarrestar.

---

## 2. Derivación Mediante Ecuaciones de Euler-Lagrange

El formalismo lagrangiano ofrece una derivación elegante y sistemática basada en energías, evitando el cálculo explícito de fuerzas de restricción. Se presentan dos versiones: el modelo de masa puntual (más simple, frecuente en textos de control) y el modelo de péndulo físico (barra uniforme, más riguroso).

### 2.1 Coordenadas generalizadas y posición del centro de masa

Se eligen las coordenadas generalizadas $q = (x, \theta)$.

**Modelo de masa puntual** (masa $m$ concentrada a distancia $\ell$ del pivote):

$$x_p = x - \ell \sin\theta, \qquad y_p = \ell \cos\theta$$

**Modelo de barra uniforme** (centro de masa a $L/2$ del pivote):

$$x_p = x + \frac{L}{2}\sin\theta, \qquad y_p = \frac{L}{2}\cos\theta$$

### 2.2 Energía cinética

**Modelo de masa puntual.** Las velocidades del péndulo son:

$$\dot{x}_p = \dot{x} - \ell\dot{\theta}\cos\theta, \qquad \dot{y}_p = -\ell\dot{\theta}\sin\theta$$

La energía cinética total del sistema es:

$$T = \frac{1}{2}M\dot{x}^2 + \frac{1}{2}m\left(\dot{x}_p^2 + \dot{y}_p^2\right)$$

Expandiendo:

$$T = \frac{1}{2}(M+m)\dot{x}^2 - m\ell\,\dot{x}\,\dot{\theta}\cos\theta + \frac{1}{2}m\ell^2\dot{\theta}^2$$

**Modelo de barra uniforme.** Incluyendo la energía cinética rotacional $\frac{1}{2}I_{cm}\dot{\theta}^2$:

$$T = \frac{1}{2}(M+m)\dot{x}^2 + \frac{1}{2}m\dot{x}L\dot{\theta}\cos\theta + \frac{mL^2}{6}\dot{\theta}^2$$

### 2.3 Energía potencial

Tomando como referencia la altura del pivote:

- **Masa puntual:** $\quad V = mg\ell\cos\theta$
- **Barra uniforme:** $\quad V = mg\dfrac{L}{2}\cos\theta$

### 2.4 El lagrangiano

Para el modelo de masa puntual:

$$\mathcal{L} = T - V = \frac{1}{2}(M+m)\dot{x}^2 - m\ell\,\dot{x}\,\dot{\theta}\cos\theta + \frac{1}{2}m\ell^2\dot{\theta}^2 - mg\ell\cos\theta$$

Las fuerzas generalizadas (no conservativas) son: $Q_x = F$ y $Q_\theta = 0$ (la junta no está actuada).

### 2.5 Ecuación de Euler-Lagrange para $q_1 = x$

Se aplica:

$$\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{x}} - \frac{\partial \mathcal{L}}{\partial x} = F$$

Calculando las derivadas parciales:

$$\frac{\partial \mathcal{L}}{\partial \dot{x}} = (M+m)\dot{x} - m\ell\dot{\theta}\cos\theta$$

$$\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{x}} = (M+m)\ddot{x} - m\ell\ddot{\theta}\cos\theta + m\ell\dot{\theta}^2\sin\theta$$

$$\frac{\partial \mathcal{L}}{\partial x} = 0$$

**Ecuación 1 — Dinámica del carro:**

$$(M+m)\ddot{x} - m\ell\cos\theta\;\ddot{\theta} + m\ell\sin\theta\;\dot{\theta}^2 = F$$

### 2.6 Ecuación de Euler-Lagrange para $q_2 = \theta$

Se aplica:

$$\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{\theta}} - \frac{\partial \mathcal{L}}{\partial \theta} = 0$$

Calculando:

$$\frac{\partial \mathcal{L}}{\partial \dot{\theta}} = -m\ell\dot{x}\cos\theta + m\ell^2\dot{\theta}$$

$$\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{\theta}} = -m\ell\ddot{x}\cos\theta + m\ell\dot{x}\dot{\theta}\sin\theta + m\ell^2\ddot{\theta}$$

$$\frac{\partial \mathcal{L}}{\partial \theta} = m\ell\dot{x}\dot{\theta}\sin\theta + mg\ell\sin\theta$$

Los términos $m\ell\dot{x}\dot{\theta}\sin\theta$ se cancelan, resultando en:

**Ecuación 2 — Dinámica del péndulo:**

$$\ell\,\ddot{\theta} - \ddot{x}\cos\theta - g\sin\theta = 0$$

### 2.7 De las ecuaciones escalares a la forma matricial: la ecuación del manipulador

Las Ecuaciones 1 y 2 forman un sistema de dos ecuaciones con dos incógnitas ($\ddot{x}$ y $\ddot{\theta}$). El paso a forma matricial no es un truco notacional: es la aplicación directa de álgebra lineal para organizar un sistema de ecuaciones según el tipo de término que contiene. La idea es clasificar cada término de ambas ecuaciones en una de tres categorías y agrupar cada categoría en su propia estructura matricial.

#### Paso previo: homogeneizar dimensiones

Para que ambas ecuaciones tengan las mismas unidades (fuerza, es decir N = kg·m/s²) y el sistema matricial sea dimensionalmente consistente, se multiplica la Ecuación 2 por $m\ell$:

$$m\ell^2\,\ddot{\theta} - m\ell\cos\theta\;\ddot{x} - mg\ell\sin\theta = 0$$

Reordenando para que $\ddot{x}$ aparezca primero:

$$-m\ell\cos\theta\;\ddot{x} + m\ell^2\;\ddot{\theta} - mg\ell\sin\theta = 0$$

#### Clasificación de términos

Ahora se tienen las dos ecuaciones en forma compatible:

**Ecuación 1:** $(M+m)\ddot{x} - m\ell\cos\theta\;\ddot{\theta} + m\ell\sin\theta\;\dot{\theta}^2 = F$

**Ecuación 2 (escalada):** $-m\ell\cos\theta\;\ddot{x} + m\ell^2\;\ddot{\theta} - mg\ell\sin\theta = 0$

Cada término pertenece a exactamente una de estas tres categorías:

**Términos de aceleración** (contienen $\ddot{x}$ o $\ddot{\theta}$) — son las incógnitas a despejar:

- En Ec.1: $(M+m)\ddot{x}$ y $-m\ell\cos\theta\;\ddot{\theta}$
- En Ec.2: $-m\ell\cos\theta\;\ddot{x}$ y $m\ell^2\;\ddot{\theta}$

**Términos de velocidad** (contienen $\dot{x}$ o $\dot{\theta}$ pero no aceleraciones) — fuerzas centrípetas y de Coriolis:

- En Ec.1: $m\ell\sin\theta\;\dot{\theta}^2$
- En Ec.2: ninguno (en el modelo sin fricción)

**Términos de posición** (dependen solo de $\theta$) — fuerzas gravitacionales:

- En Ec.1: ninguno
- En Ec.2: $-mg\ell\sin\theta$

Y al lado derecho del signo igual están las **fuerzas externas**: $F$ y $0$.

#### Construcción de cada matriz

**La matriz de masa $H(q)$.** Se agrupan los coeficientes que multiplican a las aceleraciones. La primera fila corresponde a la Ecuación 1, la segunda a la Ecuación 2:

$$H(q) = \begin{bmatrix} (M+m) & -m\ell\cos\theta \\ -m\ell\cos\theta & m\ell^2 \end{bmatrix}$$

Verificación directa: la multiplicación $H(q)\ddot{q}$ reproduce exactamente los términos de aceleración:

$$H(q)\begin{bmatrix} \ddot{x} \\ \ddot{\theta} \end{bmatrix} = \begin{bmatrix} (M+m)\ddot{x} + (-m\ell\cos\theta)\ddot{\theta} \\ (-m\ell\cos\theta)\ddot{x} + (m\ell^2)\ddot{\theta} \end{bmatrix}$$

que son, efectivamente, los términos de aceleración de cada ecuación.

**La matriz de Coriolis y fuerzas centrípetas $C(q, \dot{q})$.** El término de velocidad en la Ecuación 1 es $m\ell\sin\theta\;\dot{\theta}^2$, que se puede factorizar como $(m\ell\dot{\theta}\sin\theta) \cdot \dot{\theta}$. Esto se organiza como una matriz que multiplica al vector de velocidades:

$$C(q, \dot{q})\dot{q} = \begin{bmatrix} 0 & m\ell\dot{\theta}\sin\theta \\ 0 & 0 \end{bmatrix} \begin{bmatrix} \dot{x} \\ \dot{\theta} \end{bmatrix} = \begin{bmatrix} m\ell\sin\theta\;\dot{\theta}^2 \\ 0 \end{bmatrix}$$

El término $m\ell\sin\theta\;\dot{\theta}^2$ es la **fuerza centrípeta** que el péndulo en rotación ejerce sobre el carro. Físicamente, cuando el péndulo gira rápido ($\dot{\theta}$ grande), "jala" al carro lateralmente.

**El vector gravitacional $G(q)$.** Los términos que dependen solo de la configuración (no de velocidades ni aceleraciones):

$$G(q) = \begin{bmatrix} 0 \\ -mg\ell\sin\theta \end{bmatrix}$$

El primer componente es cero porque la gravedad no actúa directamente en la dirección horizontal del carro. El segundo componente es el torque gravitacional sobre el péndulo.

**El vector de fuerzas generalizadas.** Contiene las fuerzas y torques externos aplicados a cada coordenada generalizada:

$$\tau = \begin{bmatrix} F \\ 0 \end{bmatrix}$$

La fuerza $F$ actúa sobre $x$ (el carro), y no hay torque externo aplicado directamente sobre $\theta$ (la junta es pasiva).

#### La ecuación del manipulador completa

Juntando las cuatro piezas:

$$H(q)\,\ddot{q} + C(q, \dot{q})\,\dot{q} + G(q) = \tau$$

Escrita explícitamente:

$$\begin{bmatrix} M+m & -m\ell\cos\theta \\ -m\ell\cos\theta & m\ell^2 \end{bmatrix} \begin{bmatrix} \ddot{x} \\ \ddot{\theta} \end{bmatrix} + \begin{bmatrix} 0 & m\ell\dot{\theta}\sin\theta \\ 0 & 0 \end{bmatrix}\begin{bmatrix} \dot{x} \\ \dot{\theta} \end{bmatrix} + \begin{bmatrix} 0 \\ -mg\ell\sin\theta \end{bmatrix} = \begin{bmatrix} F \\ 0 \end{bmatrix}$$

Esta forma se denomina "ecuación del manipulador" porque es la forma estándar en robótica para describir la dinámica de cualquier cadena articulada, desde un péndulo simple hasta un brazo robótico de $n$ grados de libertad.

#### ¿Por qué $H(q)$ se llama "matriz de masa"?

El nombre proviene de la analogía directa con la segunda ley de Newton. Para una partícula, $F = ma$ relaciona fuerza con aceleración a través de la masa $m$ (un escalar). En un sistema con múltiples coordenadas generalizadas, esa relación se generaliza a:

$$\tau = H(q)\,\ddot{q} + \ldots$$

donde $H(q)$ es la **versión matricial de la masa**: conecta fuerzas generalizadas con aceleraciones generalizadas. Así como $m$ cuantifica la resistencia de un cuerpo a ser acelerado, $H(q)$ cuantifica la resistencia del sistema completo a ser acelerado en cada dirección del espacio de configuración.

#### Propiedades de $H(q)$ y su significado desde el álgebra lineal

La matriz de masa tiene tres propiedades fundamentales, todas consecuencia de la física y expresadas en lenguaje de álgebra lineal:

**Simetría: $H(q) = H(q)^T$.** Esto no es coincidencia — viene del hecho de que la energía cinética es una **forma cuadrática** en las velocidades generalizadas:

$$T = \frac{1}{2}\dot{q}^T H(q)\, \dot{q}$$

Las formas cuadráticas siempre producen matrices simétricas (o pueden simetrizarse sin cambiar el valor de $T$). En nuestro caso, el elemento fuera de la diagonal $-m\ell\cos\theta$ aparece tanto en la posición $(1,2)$ como en la $(2,1)$, reflejando que el acoplamiento inercial entre el carro y el péndulo es recíproco.

**Definitud positiva.** Para cualquier $\dot{q} \neq 0$, se cumple $\dot{q}^T H(q)\, \dot{q} > 0$. Físicamente, esto dice que si algo se mueve ($\dot{q} \neq 0$), entonces la energía cinética es estrictamente positiva. No existe ningún movimiento con energía cinética cero o negativa.

**Invertibilidad.** Es consecuencia directa de la definitud positiva (toda matriz definida positiva tiene determinante positivo, por lo tanto es invertible). Esto garantiza que siempre se pueden despejar las aceleraciones:

$$\ddot{q} = H(q)^{-1}\left(\tau - C(q, \dot{q})\,\dot{q} - G(q)\right)$$

Esta es la forma que se le pasa directamente a un solver numérico como `DifferentialEquations.jl` para simular el sistema.

#### Dependencia de la configuración

La notación $H(q)$ enfatiza que la matriz **cambia con la configuración** del sistema. El elemento de acoplamiento $-m\ell\cos\theta$ varía con el ángulo del péndulo:

- Cuando $\theta = 0$ (péndulo vertical): $\cos\theta = 1$, el acoplamiento es **máximo**. Mover el carro produce la máxima influencia sobre la aceleración angular.
- Cuando $\theta = \pi/2$ (péndulo horizontal): $\cos\theta = 0$, el acoplamiento se **anula**. En esta configuración instantánea, la aceleración del carro no afecta la aceleración angular del péndulo.

Esta dependencia es una característica intrínseca de los sistemas no lineales, y es precisamente lo que desaparece al linealizar (donde $\cos\theta \approx 1$ se fija como constante).

#### ¿Para qué sirve esta forma matricial?

La forma $H\ddot{q} + C\dot{q} + G = \tau$ no es solo estética. Tiene tres ventajas concretas para el proyecto:

**Para simular el sistema** (resolver numéricamente): como $H$ es invertible, se escribe $\ddot{q} = H^{-1}(\tau - C\dot{q} - G)$ y se pasa directamente a `DifferentialEquations.jl`.

**Para linealizar** (pasar a espacio de estados): la estructura matricial facilita calcular los jacobianos $\partial f / \partial \mathbf{x}$ y $\partial f / \partial u$ que producen las matrices $A$ y $B$ de la Sección 3.

**Para análisis de propiedades** (álgebra lineal pura): la simetría y definitud positiva de $H$ garantizan existencia y unicidad de solución, y permiten definir métricas de energía en el espacio de configuración.

### 2.8 Solución explícita para las aceleraciones

Dado que $H(q)$ es invertible, se pueden despejar $\ddot{x}$ y $\ddot{\theta}$ algebraicamente. Resolviendo el sistema $2 \times 2$:

$$\ddot{x} = \frac{F + m\ell\sin\theta\left(\ell\dot{\theta}^2 + g\cos\theta\right)}{M + m\sin^2\theta}$$

$$\ddot{\theta} = \frac{-F\cos\theta - m\ell\dot{\theta}^2\cos\theta\sin\theta - (M+m)g\sin\theta}{\ell(M + m\sin^2\theta)}$$

El denominador $M + m\sin^2\theta$ es siempre positivo (nunca se anula), lo cual es consistente con la invertibilidad de $H(q)$.

### 2.9 Modelo de barra uniforme con fricción

Para el modelo más realista (barra uniforme, fricción $b$):

$$(M+m)\ddot{x} + \frac{mL}{2}\left(\ddot{\theta}\cos\theta - \dot{\theta}^2\sin\theta\right) + b\dot{x} = F$$

$$\frac{mL^2}{3}\ddot{\theta} + \frac{mL}{2}\ddot{x}\cos\theta - mg\frac{L}{2}\sin\theta = 0$$

---

## 3. Representación en Espacio de Estados

### 3.1 ¿Qué es el espacio de estados y por qué se necesita?

Las ecuaciones obtenidas en la Sección 2 son **no lineales** y de **segundo orden** — contienen productos trigonométricos ($\sin\theta\cos\theta$), productos de variables ($\dot{\theta}^2\sin\theta$), y segundas derivadas ($\ddot{x}$, $\ddot{\theta}$). Para aplicar las herramientas de la teoría de control lineal (eigenvalores, controlabilidad, LQR), necesitamos transformarlas en un formato estándar.

La **representación en espacio de estados** es ese formato estándar. La idea es reformular cualquier sistema de EDOs de orden $n$ como un sistema equivalente de EDOs de **primer orden**, y luego (si se desea aplicar control lineal) linealizar alrededor de un punto de equilibrio. El resultado es un sistema matricial compacto que se presta completamente al análisis con herramientas de álgebra lineal.

#### Reducción a primer orden

El primer paso es puramente mecánico. Cualquier EDO de segundo orden se puede convertir en dos EDOs de primer orden introduciendo las velocidades como nuevas variables. Para el péndulo invertido, se definen las cuatro **variables de estado**:

$$x_1 = x, \qquad x_2 = \dot{x}, \qquad x_3 = \theta, \qquad x_4 = \dot{\theta}$$

Entonces las derivadas temporales de estas variables son:

$$\dot{x}_1 = x_2 \quad \text{(definición)}$$

$$\dot{x}_2 = \ddot{x} = f_1(x_1, x_2, x_3, x_4, u) \quad \text{(Ecuación 1 despejada)}$$

$$\dot{x}_3 = x_4 \quad \text{(definición)}$$

$$\dot{x}_4 = \ddot{\theta} = f_2(x_1, x_2, x_3, x_4, u) \quad \text{(Ecuación 2 despejada)}$$

donde $f_1$ y $f_2$ son las expresiones de la Sección 2.8. Agrupando en un vector $\mathbf{x} = [x_1, x_2, x_3, x_4]^T$, el sistema completo se escribe como:

$$\dot{\mathbf{x}} = f(\mathbf{x}, u)$$

Este es el **modelo no lineal en espacio de estados**. Es un sistema de 4 ecuaciones diferenciales de primer orden, completamente equivalente a las 2 ecuaciones de segundo orden originales. La función $f: \mathbb{R}^4 \times \mathbb{R} \to \mathbb{R}^4$ es no lineal (contiene senos, cosenos, productos).

#### Linealización: del modelo no lineal al modelo lineal

El paso crucial es la **linealización** alrededor del punto de equilibrio $\mathbf{x}_0 = [0, 0, 0, 0]^T$, $u_0 = 0$. La herramienta matemática es la expansión de Taylor de primer orden de $f$:

$$f(\mathbf{x}, u) \approx f(\mathbf{x}_0, u_0) + \underbrace{\frac{\partial f}{\partial \mathbf{x}}\bigg|_{(\mathbf{x}_0, u_0)}}_{A} (\mathbf{x} - \mathbf{x}_0) + \underbrace{\frac{\partial f}{\partial u}\bigg|_{(\mathbf{x}_0, u_0)}}_{B} (u - u_0)$$

Como $f(\mathbf{x}_0, u_0) = 0$ (el equilibrio es un punto fijo), y definiendo las desviaciones $\delta\mathbf{x} = \mathbf{x} - \mathbf{x}_0$, $\delta u = u - u_0$, se obtiene:

$$\dot{\delta\mathbf{x}} = A\,\delta\mathbf{x} + B\,\delta u$$

Las matrices $A$ y $B$ son los **jacobianos** de $f$ evaluados en el equilibrio. Es decir, la matriz $A$ es la matriz $4 \times 4$ de derivadas parciales:

$$A_{ij} = \frac{\partial f_i}{\partial x_j}\bigg|_{eq}$$

y el vector $B$ es:

$$B_i = \frac{\partial f_i}{\partial u}\bigg|_{eq}$$

El significado físico de la linealización es que, para perturbaciones **pequeñas** alrededor del equilibrio, el sistema no lineal se comporta *aproximadamente* como un sistema lineal. La calidad de esta aproximación se degrada a medida que $\theta$ se aleja de cero.

### 3.2 Las cuatro matrices del espacio de estados

El modelo linealizado completo tiene la forma:

$$\dot{\mathbf{x}} = A\mathbf{x} + Bu$$

$$\mathbf{y} = C\mathbf{x} + Du$$

donde cada matriz tiene un rol específico:

**Matriz de dinámica $A$ ($4 \times 4$).** Captura cómo el estado evoluciona por sí mismo, sin intervención del controlador. Si $u = 0$, el sistema evoluciona según $\dot{\mathbf{x}} = A\mathbf{x}$, cuya solución es $\mathbf{x}(t) = e^{At}\mathbf{x}(0)$. Los **eigenvalores de $A$** determinan completamente si el sistema es estable (todos con parte real negativa), inestable (al menos uno con parte real positiva), o marginalmente estable (parte real cero).

**Matriz de entrada $B$ ($4 \times 1$).** Mapea la acción de control $u$ a tasas de cambio del estado. Describe *cómo* la fuerza aplicada al carro afecta la aceleración de cada variable. La dimensión de $B$ refleja que tenemos 4 estados pero solo 1 entrada de control.

**Matriz de salida $C$ ($p \times 4$, donde $p$ es el número de mediciones).** Selecciona qué variables de estado son directamente medibles. En la práctica, se miden posiciones ($x$, $\theta$) con encoders o potenciómetros, pero no las velocidades. $C$ formaliza qué información está disponible para el controlador.

**Matriz de transmisión directa $D$ ($p \times 1$).** Captura el efecto instantáneo de la entrada sobre la salida. Para sistemas mecánicos como el péndulo invertido, $D = 0$ porque una fuerza produce aceleración (efecto de segundo orden), no un cambio instantáneo en posición o ángulo.

### 3.3 Obtención de $A$ y $B$ para el péndulo invertido

Para obtener las matrices $A$ y $B$ concretas, se linealizan las ecuaciones de la Sección 2.9 (modelo con momento de inercia $I$ y fricción $b$) aplicando la aproximación de ángulo pequeño:

$$\sin\theta \approx \theta, \qquad \cos\theta \approx 1, \qquad \dot{\theta}^2 \approx 0$$

Las ecuaciones linealizadas son:

$$(M+m)\ddot{x} + ml\,\ddot{\theta} + b\dot{x} = u$$

$$(I + ml^2)\ddot{\theta} + ml\,\ddot{x} - mgl\,\theta = 0$$

Resolviendo simultáneamente para $\ddot{x}$ y $\ddot{\theta}$, con el denominador $p = I(M+m) + Mml^2$:

$$\ddot{x} = \frac{-(I+ml^2)\,b\,\dot{x} + m^2gl^2\,\theta + (I+ml^2)\,u}{p}$$

$$\ddot{\theta} = \frac{-ml\,b\,\dot{x} + mgl(M+m)\,\theta + ml\,u}{p}$$

Estas expresiones, junto con las identidades $\dot{x}_1 = x_2$ y $\dot{x}_3 = x_4$, dan directamente las matrices:

$$A = \begin{bmatrix} 0 & 1 & 0 & 0 \\ 0 & \dfrac{-(I+ml^2)b}{p} & \dfrac{m^2gl^2}{p} & 0 \\ 0 & 0 & 0 & 1 \\ 0 & \dfrac{-mlb}{p} & \dfrac{mgl(M+m)}{p} & 0 \end{bmatrix}$$

$$B = \begin{bmatrix} 0 \\ \dfrac{I+ml^2}{p} \\ 0 \\ \dfrac{ml}{p} \end{bmatrix}$$

$$C = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \end{bmatrix}, \qquad D = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$

### 3.4 Valores numéricos e interpretación

Con los parámetros estándar ($M = 0.5$ kg, $m = 0.2$ kg, $b = 0.1$ N·s/m, $l = 0.3$ m, $I = 0.006$ kg·m², $g = 9.8$ m/s²), se obtiene $p = 0.0132$ y:

$$A = \begin{bmatrix} 0 & 1 & 0 & 0 \\ 0 & -0.1818 & 2.6727 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & -0.4545 & 31.1818 & 0 \end{bmatrix}, \qquad B = \begin{bmatrix} 0 \\ 1.8182 \\ 0 \\ 4.5455 \end{bmatrix}$$

**Lectura fila por fila de $A$:**

- **Fila 1** $[0, 1, 0, 0]$: $\dot{x}_1 = x_2$, es decir, la derivada de la posición es la velocidad (trivial, por definición).
- **Fila 2** $[0, -0.18, 2.67, 0]$: $\dot{x}_2 = -0.18\,\dot{x} + 2.67\,\theta$. La aceleración del carro depende de la fricción (frena) y del ángulo (el peso del péndulo inclinado empuja al carro).
- **Fila 3** $[0, 0, 0, 1]$: $\dot{x}_3 = x_4$, la derivada del ángulo es la velocidad angular (trivial).
- **Fila 4** $[0, -0.45, 31.18, 0]$: $\dot{x}_4 = -0.45\,\dot{x} + 31.18\,\theta$. La aceleración angular depende fuertemente del ángulo (el $31.18$ es el **término dominante de inestabilidad gravitacional**: cualquier desviación angular se amplifica violentamente).

**Lectura de $B$:**

- Las componentes $B_1 = 0$ y $B_3 = 0$ indican que la fuerza no cambia instantáneamente las posiciones.
- $B_2 = 1.82$ y $B_4 = 4.55$ indican que la fuerza produce aceleraciones. El péndulo responde con mayor ganancia ($4.55$ vs $1.82$) porque tiene menor inercia efectiva que el carro.

### 3.5 Controlabilidad — Criterio de Kalman

El sistema $(A, B)$ es completamente controlable si y solo si la **matriz de controlabilidad**:

$$\mathcal{C} = \begin{bmatrix} B & AB & A^2B & A^3B \end{bmatrix}$$

tiene **rango** $n = 4$ (igual a la dimensión del estado).

La intuición detrás de este criterio es que $B$ dice qué direcciones del espacio de estados se pueden afectar directamente con $u$; $AB$ dice qué nuevas direcciones se alcanzan después de un paso de evolución dinámica; $A^2B$ después de dos pasos; y así sucesivamente. Si después de $n-1$ pasos las columnas generan todo $\mathbb{R}^n$, entonces toda dirección es alcanzable.

Para el péndulo invertido, $\text{rank}(\mathcal{C}) = 4$: **el sistema es controlable**. Físicamente, aunque solo se empuja al carro, el acoplamiento dinámico permite controlar indirectamente el ángulo del péndulo.

### 3.6 Observabilidad — Criterio de Kalman

El sistema $(A, C)$ es completamente observable si la **matriz de observabilidad**:

$$\mathcal{O} = \begin{bmatrix} C \\ CA \\ CA^2 \\ CA^3 \end{bmatrix}$$

tiene **rango 4**. Para el péndulo invertido con $C$ midiendo posición del carro y ángulo, $\text{rank}(\mathcal{O}) = 4$: **el sistema es observable**. Incluso midiendo solo una de las dos posiciones ($x$ o $\theta$), el sistema permanece observable: las velocidades se pueden reconstruir a partir del historial de mediciones gracias al acoplamiento dinámico.

### 3.7 Diseño de controladores

#### Realimentación de estado completo

La ley de control $u = -K\mathbf{x}$ transforma el sistema en lazo cerrado con dinámica gobernada por:

$$A_{cl} = A - BK$$

La estabilidad depende de que todos los **eigenvalores de $A - BK$ tengan parte real negativa**.

#### Asignación de polos (Pole Placement)

Se eligen ubicaciones deseadas para los polos en lazo cerrado $s_1, s_2, s_3, s_4$ (todos con parte real negativa) y se calcula $K$ tal que:

$$\det(sI - A + BK) = (s - s_1)(s - s_2)(s - s_3)(s - s_4)$$

El método requiere que $(A, B)$ sea controlable. Se puede usar la fórmula de Ackermann o la igualación directa de coeficientes del polinomio característico.

#### Regulador Cuadrático Lineal (LQR)

El LQR minimiza la función de costo:

$$J = \int_0^\infty \left(\mathbf{x}^T Q\, \mathbf{x} + u^T R\, u\right) dt$$

donde:

- $Q$ (semidefinida positiva, $4 \times 4$): penaliza desviaciones del estado.
- $R$ (definida positiva, escalar para una entrada): penaliza el esfuerzo de control.

La ganancia óptima es:

$$K = R^{-1}B^TP$$

donde $P$ es la solución de la **ecuación algebraica de Riccati continua (CARE)**:

$$A^TP + PA - PBR^{-1}B^TP + Q = 0$$

**Ventajas del LQR sobre la asignación de polos:**

- Método sistemático (sin necesidad de adivinar ubicaciones de polos).
- Garantiza estabilidad cuando el sistema es controlable.
- Ofrece márgenes de robustez garantizados (margen de ganancia infinito, margen de fase $\geq 60°$ en SISO).

#### Observador de Luenberger

Cuando no se miden todas las variables de estado (típicamente las velocidades $\dot{x}$ y $\dot{\theta}$), se diseña un observador:

$$\dot{\hat{\mathbf{x}}} = A\hat{\mathbf{x}} + Bu + L(\mathbf{y} - C\hat{\mathbf{x}})$$

El error de estimación evoluciona según $\dot{e} = (A - LC)e$; los eigenvalores de $(A - LC)$ determinan la convergencia. La regla práctica es colocar los polos del observador **4 a 10 veces más rápidos** que los del controlador. El **principio de separación** garantiza que controlador y observador pueden diseñarse independientemente.

---

## 4. Del Lagrangiano al Controlador: Dos Enfoques que se Complementan

La relación entre el enfoque de Euler-Lagrange y el de teoría de control forma un **pipeline de dos etapas**, no una dicotomía.

### 4.1 Etapa de modelado (Euler-Lagrange)

El enfoque lagrangiano pertenece a la etapa de **modelado**: produce las ecuaciones diferenciales no lineales que describen la planta física. Es un enfoque **basado en energía y libre de coordenadas de restricción** que trabaja naturalmente con coordenadas generalizadas y evita calcular fuerzas internas en las juntas. Es especialmente ventajoso para sistemas multi-cuerpo complejos (péndulos dobles, triples, manipuladores robóticos).

### 4.2 Etapa de diseño (Teoría de Control)

El enfoque de control pertenece a la etapa de **diseño**: toma las ecuaciones obtenidas, las linealiza y construye la ley de retroalimentación que estabiliza el sistema. Se centra en la representación matricial y ofrece herramientas poderosas para análisis de estabilidad, diseño de controladores y verificación de propiedades estructurales (controlabilidad, observabilidad).

### 4.3 El puente: linealización

El puente entre ambos enfoques es la **linealización por expansión de Taylor**. Las ecuaciones no lineales obtenidas vía Euler-Lagrange se evalúan alrededor del equilibrio inestable, aplicando las aproximaciones de ángulo pequeño. Este paso produce directamente las matrices $A$ y $B$ del espacio de estados. Sin las ecuaciones del modelo, no hay nada que linealizar ni controlar.

### 4.4 Tabla comparativa

| Aspecto                  | Euler-Lagrange                                  | Teoría de Control                              |
|:------------------------:|:------------------------------------------------|:-----------------------------------------------|
| **Objetivo**             | Obtener las ecuaciones de movimiento            | Diseñar la ley de retroalimentación            |
| **Base conceptual**      | Energía (cinética y potencial)                  | Variables de estado y matrices                 |
| **Resultado**            | EDOs no lineales acopladas                      | Ley de control $u = -K\mathbf{x}$              |
| **Maneja no linealidades** | Sí, naturalmente                               | Requiere linealización previa                  |
| **Fuerzas internas**     | No las necesita (se eliminan automáticamente)   | No son relevantes (trabaja con el modelo ya derivado) |
| **Escalabilidad**        | Excelente para sistemas multi-cuerpo            | Puede volverse complejo con muchos estados     |
| **Herramientas clave**   | Lagrangiano, coordenadas generalizadas          | Eigenvalores, rango, Riccati                   |
| **¿Cuándo usarlo?**      | Siempre que se necesite derivar la dinámica     | Siempre que se necesite estabilizar o regular   |

### 4.5 Unificación: Lagrangianos Controlados

Existe una unificación profunda de ambos enfoques: el método de **Lagrangianos Controlados** (Bloch, Leonard y Marsden, Caltech, 2000) diseña controladores modificando el propio lagrangiano, de modo que la dinámica en lazo cerrado preserva la estructura lagrangiana. Este método demuestra que la estabilización puede entenderse como la creación efectiva de un pozo de energía mediante retroalimentación.

---

## 5. El Álgebra Lineal como Lenguaje Unificador

Los conceptos de álgebra lineal no son meras herramientas auxiliares: son el **lenguaje fundamental** que conecta modelado, análisis y diseño.

### 5.1 Eigenvalores y estabilidad

Los eigenvalores de la matriz $A$ determinan completamente la estabilidad en lazo abierto. Para el ejemplo numérico, los eigenvalores de $A$ son aproximadamente:

$$\lambda \approx \{0,\; -0.14,\; +5.57,\; -5.60\}$$

El eigenvalor **positivo** $\lambda \approx +5.57$ confirma la inestabilidad: las perturbaciones crecen como $e^{5.57t}$. Después de aplicar retroalimentación $u = -K\mathbf{x}$, los eigenvalores de $A - BK$ se ubican en el semiplano izquierdo complejo, garantizando estabilidad asintótica.

### 5.2 Eigenvectores y modos del sistema

Cada par eigenvalor-eigenvector $(\lambda_i, \mathbf{v}_i)$ corresponde a un **modo natural** del sistema:

- Eigenvalores reales negativos: modos de decaimiento exponencial.
- Eigenvalores positivos: modos de crecimiento (inestabilidad).
- Pares complejos conjugados $a \pm j\omega$: modos oscilatorios.

Los eigenvectores definen las *formas modales*, es decir, qué estados participan en cada modo y con qué amplitud relativa. El modo inestable del péndulo involucra principalmente el ángulo $\theta$ creciendo mientras el carro se desplaza.

### 5.3 Rango de matrices y propiedades estructurales

- **Controlabilidad:** exige que las columnas de $\mathcal{C} = [B,\; AB,\; A^2B,\; A^3B]$ generen todo $\mathbb{R}^4$ (rango completo).
- **Observabilidad:** exige lo análogo para las filas de $\mathcal{O}$.

Estas son condiciones puramente de álgebra lineal: el rango de una matriz determina si el espacio de estados es alcanzable (controlabilidad) o distinguible (observabilidad).

### 5.4 La ecuación de Riccati como problema matricial

El LQR requiere resolver:

$$A^TP + PA - PBR^{-1}B^TP + Q = 0$$

Esta es una ecuación matricial cuadrática. Su solución $P$ es una matriz simétrica definida positiva que codifica el compromiso óptimo entre regulación del estado y esfuerzo de control. La ganancia $K = R^{-1}B^TP$ modifica directamente la estructura de eigenvalores, "reflejando" los eigenvalores inestables al semiplano estable.

### 5.5 La matriz de masa del lagrangiano

La matriz $H(q)$ en la ecuación del manipulador es simétrica y definida positiva, lo que garantiza su invertibilidad. La energía cinética $T = \frac{1}{2}\dot{q}^T H(q)\dot{q}$ define una **forma bilineal** en el espacio tangente, y las propiedades espectrales de $H$ (sus eigenvalores son las inercias principales) determinan las escalas naturales de tiempo del sistema mecánico.

---

## 6. Implementación Computacional en Julia

El ecosistema de Julia ofrece un flujo de trabajo completo para todas las etapas del proyecto: modelado, linealización, diseño de controladores y simulación/animación.

### 6.1 Stack de paquetes

| Paquete                        | Rol en el proyecto                                      |
|:-------------------------------|:--------------------------------------------------------|
| `DifferentialEquations.jl`     | Integración numérica de las EDOs no lineales (solvers `Tsit5()`, `Rodas4()`, etc.) |
| `ForwardDiff.jl`               | Linealización automática vía diferenciación automática (cálculo de jacobianos $A$, $B$) |
| `ControlSystems.jl`            | Funciones de control: `lqr`, `place`, `ctrb`, `obsv`, `are`, `ss` |
| `RobustAndOptimalControl.jl`   | Tutorial completo del carro-péndulo, herramientas avanzadas |
| `MatrixEquations.jl`           | Resolución directa de la ecuación de Riccati             |
| `Symbolics.jl`                 | Derivación simbólica de las ecuaciones de movimiento     |
| `GLMakie.jl`                   | Visualización y animación del péndulo invertido          |
| `Pluto.jl`                     | Notebooks interactivos para documentación del proyecto   |

### 6.2 Ejemplo: linealización automática

```julia
using ForwardDiff, ControlSystemsBase

# Definir la dinámica no lineal (Sección 2.8)
function cartpole(x, u)
    # x = [posición, velocidad, ángulo, vel_angular]
    # u = fuerza aplicada al carro
    # Retorna dx/dt
    ...
end

# Punto de equilibrio
x0 = [0.0, 0.0, 0.0, 0.0]
u0 = [0.0]

# Jacobianos automáticos → matrices A y B
Ac = ForwardDiff.jacobian(x -> cartpole(x, u0), x0)
Bc = ForwardDiff.jacobian(u -> cartpole(x0, u), u0)

# Modelo en espacio de estados
Cc = [1 0 0 0; 0 0 1 0]
sys = ss(Ac, Bc, Cc, 0)
```

### 6.3 Ejemplo: diseño LQR y verificación

```julia
using ControlSystems, LinearAlgebra

# Verificar controlabilidad
C_ctrl = ctrb(Ac, Bc)
println("Rango de controlabilidad: ", rank(C_ctrl))  # Debe ser 4

# Diseñar controlador LQR
Q = diagm([1.0, 0, 10.0, 0])  # Penalizar posición y ángulo
R = [1.0]                       # Costo de control
K = lqr(sys, Q, R)

# Eigenvalores en lazo cerrado (todos deben tener Re < 0)
A_cl = Ac - Bc * K
println("Eigenvalores lazo cerrado: ", eigvals(A_cl))
```

---

## 7. Conclusión

El péndulo invertido concentra, en un sistema de apenas cuatro variables de estado, una densidad extraordinaria de conceptos de álgebra lineal aplicada.

La secuencia conceptual completa del proyecto es:

**Mecánica lagrangiana** $\rightarrow$ **Ecuaciones no lineales** $\rightarrow$ **Ecuación del manipulador** $H\ddot{q} + C\dot{q} + G = \tau$ $\rightarrow$ **Reducción a primer orden** $\rightarrow$ **Linealización (jacobiano)** $\rightarrow$ **Espacio de estados** $(A,B,C,D)$ $\rightarrow$ **Análisis espectral (eigenvalores)** $\rightarrow$ **Controlabilidad/observabilidad (rango)** $\rightarrow$ **Diseño del controlador (Riccati / polos)** $\rightarrow$ **Estabilización**

Cada paso es, en su esencia, una operación de álgebra lineal. El hecho de que un sistema con un eigenvalor positivo ($\lambda \approx +5.57$) pueda estabilizarse completamente mediante una simple ley $u = -K\mathbf{x}$ — moviendo todos los eigenvalores al semiplano izquierdo — es quizá la demostración más tangible del poder del álgebra lineal: la retroalimentación transforma la estructura espectral de un operador lineal, convirtiendo inestabilidad en estabilidad.

---

## Referencias

1. Ogata, K. *Modern Control Engineering*, 5th ed. Prentice Hall, 2010.
2. Franklin, G.F., Powell, J.D. y Emami-Naeini, A. *Feedback Control of Dynamic Systems*, 8th ed. Pearson.
3. Åström, K.J. y Murray, R.M. *Feedback Systems: An Introduction for Scientists and Engineers*. Princeton University Press. Disponible en: [https://www.cds.caltech.edu/~murray/amwiki](https://www.cds.caltech.edu/~murray/amwiki)
4. Goldstein, H., Poole, C. y Safko, J. *Classical Mechanics*, 3rd ed. Addison-Wesley, 2002.
5. Boubaker, O. "The Inverted Pendulum Benchmark in Nonlinear Control Theory: A Survey." *Int. J. Advanced Robotic Systems*, 10(5), 2013.
6. Lundberg, K.H. y Barton, T.W. "History of Inverted-Pendulum Systems." IFAC, 2009.
7. Bloch, A.M., Leonard, N.E. y Marsden, J.E. "Controlled Lagrangians and the Stabilization of Mechanical Systems." *IEEE Trans. Automatic Control*, 45(12), 2000.
8. Control Tutorials for MATLAB and Simulink — University of Michigan: [https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum](https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum)
9. RobustAndOptimalControl.jl — Cart-Pole Tutorial: [https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/cartpole/](https://juliacontrol.github.io/RobustAndOptimalControl.jl/dev/cartpole/)
10. Rackauckas, C. y Nie, Q. "DifferentialEquations.jl — A Performant and Feature-Rich Ecosystem for Solving Differential Equations in Julia." *JORS*, 5(1), 2017.
11. Tedrake, R. *Underactuated Robotics: Algorithms for Walking, Running, Swimming, Flying, and Manipulation*. MIT. Disponible en: [https://underactuated.mit.edu](https://underactuated.mit.edu)
12. IIT Kharagpur Virtual Labs — "Controllability and Observability of Inverted Pendulum on Cart": [http://vlabs.iitkgp.ac.in/dctrl/Exp6/theory.html](http://vlabs.iitkgp.ac.in/dctrl/Exp6/theory.html)
13. Narayan, A. "Cartpole Dynamics and Control": [https://ashwinnarayan.com/post/cartpole-dynamics/](https://ashwinnarayan.com/post/cartpole-dynamics/)
14. Piedrahita, D. "Cart-Pole Control": [https://danielpiedrahita.wordpress.com/portfolio/cart-pole-control/](https://danielpiedrahita.wordpress.com/portfolio/cart-pole-control/)
