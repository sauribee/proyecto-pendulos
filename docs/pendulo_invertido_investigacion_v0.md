# El Péndulo Invertido: Modelado, Control y Álgebra Lineal Aplicada

---

## Introducción

El péndulo invertido sobre un carro es uno de los sistemas dinámicos más estudiados en ingeniería y física aplicada. Este sistema —una barra rígida balanceándose en posición vertical sobre un carro móvil— es inherentemente inestable, y constituye el banco de pruebas canónico donde convergen la mecánica lagrangiana, la teoría de control moderna y el álgebra lineal.

Su relevancia radica en que captura, en un modelo relativamente simple, los desafíos fundamentales de estabilizar sistemas no lineales subactuados. Aparece como modelo subyacente en cohetes, vehículos autobalanceados (como el Segway), robots bípedos y en el control postural humano.

Este informe desarrolla la derivación completa de las ecuaciones de movimiento mediante Euler-Lagrange, su linealización en espacio de estados, el diseño de controladores y las conexiones profundas con conceptos de álgebra lineal.

---

## 1. Anatomía del Sistema: Carro, Péndulo y Fuerzas

El sistema carro-péndulo invertido consta de tres elementos fundamentales:

- Un **carro de masa** $M$ que se desplaza sobre un riel horizontal restringido al eje $x$, impulsado por una fuerza externa $F$ (la entrada de control).
- Una **barra rígida de masa** $m$ **y longitud** $L$, articulada al carro mediante una junta de revolución, cuyo centro de masa se ubica a una distancia $l = L/2$ del pivote.
- El sistema posee **dos grados de libertad** (posición del carro $x$ y ángulo del péndulo $\theta$) pero solo **una entrada de control** ($F$), lo que lo clasifica como un sistema *subactuado*.

### 1.1 Variables de estado

Las cuatro variables que describen completamente la configuración dinámica del sistema son:

| Variable       | Descripción                                                                              |
|:--------------:|:-----------------------------------------------------------------------------------------|
| $x$            | Posición horizontal del carro                                                            |
| $\dot{x}$      | Velocidad lineal del carro                                                               |
| $\theta$       | Ángulo del péndulo respecto a la vertical ($\theta = 0$ es la posición invertida)        |
| $\dot{\theta}$ | Velocidad angular del péndulo                                                            |

### 1.2 Parámetros del sistema

Los parámetros físicos del modelo son:

| Parámetro | Descripción                                    | Valor típico       |
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

### 2.7 Forma matricial (ecuación del manipulador)

Las dos ecuaciones se escriben en la forma $H(q)\ddot{q} + C(q,\dot{q})\dot{q} + G(q) = \tau$:

$$\begin{bmatrix} M+m & -m\ell\cos\theta \\ -m\ell\cos\theta & m\ell^2 \end{bmatrix} \begin{bmatrix} \ddot{x} \\ \ddot{\theta} \end{bmatrix} + \begin{bmatrix} 0 & m\ell\dot{\theta}\sin\theta \\ 0 & 0 \end{bmatrix}\begin{bmatrix} \dot{x} \\ \dot{\theta} \end{bmatrix} + \begin{bmatrix} 0 \\ -mg\ell\sin\theta \end{bmatrix} = \begin{bmatrix} F \\ 0 \end{bmatrix}$$

La **matriz de masa** $H(q)$ es simétrica y definida positiva (siempre invertible), lo que garantiza que las aceleraciones $\ddot{x}$ y $\ddot{\theta}$ pueden resolverse explícitamente.

### 2.8 Solución explícita para las aceleraciones

Resolviendo el sistema algebraicamente:

$$\ddot{x} = \frac{F + m\ell\sin\theta\left(\ell\dot{\theta}^2 + g\cos\theta\right)}{M + m\sin^2\theta}$$

$$\ddot{\theta} = \frac{-F\cos\theta - m\ell\dot{\theta}^2\cos\theta\sin\theta - (M+m)g\sin\theta}{\ell(M + m\sin^2\theta)}$$

### 2.9 Modelo de barra uniforme con fricción

Para el modelo más realista (barra uniforme, fricción $b$):

$$(M+m)\ddot{x} + \frac{mL}{2}\left(\ddot{\theta}\cos\theta - \dot{\theta}^2\sin\theta\right) + b\dot{x} = F$$

$$\frac{mL^2}{3}\ddot{\theta} + \frac{mL}{2}\ddot{x}\cos\theta - mg\frac{L}{2}\sin\theta = 0$$

---

## 3. Representación en Espacio de Estados y Diseño de Controladores

### 3.1 Linealización alrededor del equilibrio inestable

Para diseñar controladores lineales, se linealizan las ecuaciones no lineales alrededor del punto de equilibrio superior ($x=0$, $\dot{x}=0$, $\theta=0$, $\dot{\theta}=0$) mediante la **aproximación de ángulo pequeño**:

$$\sin\theta \approx \theta, \qquad \cos\theta \approx 1, \qquad \dot{\theta}^2 \approx 0$$

Sustituyendo en las ecuaciones del modelo con momento de inercia $I$:

$$(M+m)\ddot{x} + ml\,\ddot{\theta} + b\dot{x} = u$$

$$(I + ml^2)\ddot{\theta} + ml\,\ddot{x} - mgl\,\theta = 0$$

Resolviendo simultáneamente para $\ddot{x}$ y $\ddot{\theta}$, con el denominador:

$$p = I(M+m) + Mml^2$$

se obtiene:

$$\ddot{x} = \frac{-(I+ml^2)\,b\,\dot{x} + m^2gl^2\,\theta + (I+ml^2)\,u}{p}$$

$$\ddot{\theta} = \frac{-ml\,b\,\dot{x} + mgl(M+m)\,\theta + ml\,u}{p}$$

### 3.2 Matrices del espacio de estados

Definiendo el vector de estado $\mathbf{x} = [x,\; \dot{x},\; \theta,\; \dot{\theta}]^T$, el sistema linealizado adopta la forma canónica:

$$\dot{\mathbf{x}} = A\mathbf{x} + Bu, \qquad \mathbf{y} = C\mathbf{x} + Du$$

con las matrices:

$$A = \begin{bmatrix} 0 & 1 & 0 & 0 \\ 0 & \dfrac{-(I+ml^2)b}{p} & \dfrac{m^2gl^2}{p} & 0 \\ 0 & 0 & 0 & 1 \\ 0 & \dfrac{-mlb}{p} & \dfrac{mgl(M+m)}{p} & 0 \end{bmatrix}$$

$$B = \begin{bmatrix} 0 \\ \dfrac{I+ml^2}{p} \\ 0 \\ \dfrac{ml}{p} \end{bmatrix}$$

$$C = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \end{bmatrix}, \qquad D = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$

### 3.3 Valores numéricos

Con los parámetros estándar ($M = 0.5$ kg, $m = 0.2$ kg, $b = 0.1$ N·s/m, $l = 0.3$ m, $I = 0.006$ kg·m², $g = 9.8$ m/s²), se obtiene $p = 0.0132$ y:

$$A = \begin{bmatrix} 0 & 1 & 0 & 0 \\ 0 & -0.1818 & 2.6727 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & -0.4545 & 31.1818 & 0 \end{bmatrix}, \qquad B = \begin{bmatrix} 0 \\ 1.8182 \\ 0 \\ 4.5455 \end{bmatrix}$$

**Interpretación física de las matrices:**

- **Matriz $A$:** codifica la dinámica interna. El elemento $A_{4,3} \approx 31.18$ es el término dominante de inestabilidad gravitacional: una desviación angular produce aceleración angular creciente.
- **Matriz $B$:** mapea la fuerza de control a aceleraciones. El péndulo responde con mayor ganancia ($4.55$ vs $1.82$) por su menor inercia efectiva.
- **Matriz $C$:** selecciona las variables medidas (posición del carro y ángulo).
- **Matriz $D = 0$:** la fuerza no afecta instantáneamente las posiciones (solo aceleraciones).

### 3.4 Controlabilidad — Criterio de Kalman

El sistema $(A, B)$ es completamente controlable si y solo si la **matriz de controlabilidad**:

$$\mathcal{C} = \begin{bmatrix} B & AB & A^2B & A^3B \end{bmatrix}$$

tiene **rango** $n = 4$ (igual a la dimensión del estado).

Para el péndulo invertido, $\text{rank}(\mathcal{C}) = 4$: **el sistema es controlable**. Físicamente, esto significa que aplicando secuencias apropiadas de fuerza al carro, es posible llevar las cuatro variables de estado desde cualquier condición inicial a cualquier estado deseado en tiempo finito, a pesar de que solo se controla directamente el carro. El acoplamiento dinámico entre carro y péndulo permite el control indirecto del ángulo.

### 3.5 Observabilidad — Criterio de Kalman

El sistema $(A, C)$ es completamente observable si la **matriz de observabilidad**:

$$\mathcal{O} = \begin{bmatrix} C \\ CA \\ CA^2 \\ CA^3 \end{bmatrix}$$

tiene **rango 4**. Para el péndulo invertido con $C$ midiendo posición del carro y ángulo, $\text{rank}(\mathcal{O}) = 4$: **el sistema es observable**. Incluso midiendo solo una de las dos posiciones ($x$ o $\theta$), el sistema permanece observable: las velocidades se pueden reconstruir a partir del historial de mediciones.

### 3.6 Diseño de controladores

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

| Aspecto                    | Euler-Lagrange                                  | Teoría de Control                                     |
|:--------------------------:|:-----------------------------------------------:|:-----------------------------------------------------:|
| **Objetivo**               | Obtener las ecuaciones de movimiento            | Diseñar la ley de retroalimentación                   |
| **Base conceptual**        | Energía (cinética y potencial)                  | Variables de estado y matrices                        |
| **Resultado**              | EDOs no lineales acopladas                      | Ley de control $u = -K\mathbf{x}$                     |
| **Maneja no linealidades** | Sí, naturalmente                                | Requiere linealización previa                         |
| **Fuerzas internas**       | No las necesita (se eliminan automáticamente)   | No son relevantes (trabaja con el modelo ya derivado) |
| **Escalabilidad**          | Excelente para sistemas multi-cuerpo            | Puede volverse complejo con muchos estados            |
| **Herramientas clave**     | Lagrangiano, coordenadas generalizadas          | Eigenvalores, rango, Riccati                          |
| **¿Cuándo usarlo?**        | Siempre que se necesite derivar la dinámica     | Siempre que se necesite estabilizar o regular         |  

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

### 5.5 Linealización como operación de álgebra lineal

El acto de linealizar es calcular el **jacobiano** (matriz de derivadas parciales) de la función no lineal evaluado en el punto de equilibrio:

$$A = \frac{\partial f}{\partial \mathbf{x}}\bigg|_{eq}, \qquad B = \frac{\partial f}{\partial u}\bigg|_{eq}$$

### 5.6 La matriz de masa del lagrangiano

La matriz $H(q)$ en la ecuación del manipulador es simétrica y definida positiva, lo que garantiza su invertibilidad: una propiedad de álgebra lineal esencial para que las aceleraciones estén bien definidas.

---

## 6. Implementación Computacional

### 6.1 Julia

El ecosistema de Julia ofrece un flujo de trabajo completo y moderno:

- **`DifferentialEquations.jl`** (SciML): solvers de alto rendimiento para las ecuaciones no lineales (métodos `Tsit5()`, `Rodas4()`, etc.).
- **`ForwardDiff.jl`**: linealización automática mediante diferenciación automática, calculando los jacobianos $A$ y $B$ numéricamente sin derivación manual.
- **`ControlSystems.jl`**: funciones directas como `lqr(sys, Q, R)`, `place(A, B, poles)`, `ctrb(A, B)`, `obsv(A, C)`, y `are()` para la ecuación de Riccati.
- **`RobustAndOptimalControl.jl`**: incluye un tutorial completo del carro-péndulo.
- **`Multibody.jl`**: modelado basado en componentes físicos, integrado con `ModelingToolkit.jl`.

Ejemplo de linealización automática en Julia:

```julia
using ForwardDiff, ControlSystemsBase

# Linealización automática en el punto de equilibrio
Ac = ForwardDiff.jacobian(x -> cartpole(x, u0), x0)
Bc = ForwardDiff.jacobian(u -> cartpole(x0, u), u0)
sys = ss(Ac, Bc, Cc, 0)  # Modelo en espacio de estados
```

### 6.2 MATLAB/Simulink

El recurso más completo es el *Control Tutorials for MATLAB and Simulink* de la Universidad de Michigan, que cubre el ciclo completo desde modelado hasta control digital con código ejecutable.

### 6.3 Python

- **`python-control`:** replica la funcionalidad de MATLAB (`control.lqr`, `control.ss`, `control.ctrb`).
- **`scipy.integrate.solve_ivp`:** resuelve las ODEs no lineales.
- **`scipy.linalg.solve_continuous_are`:** resuelve la ecuación de Riccati.
- **`SymPy`:** permite derivaciones simbólicas de las ecuaciones de movimiento.

---

## 7. Conclusión

El péndulo invertido concentra, en un sistema de apenas cuatro variables de estado, una densidad extraordinaria de conceptos de álgebra lineal aplicada.

La secuencia conceptual completa es:

**Mecánica lagrangiana** $\rightarrow$ **Ecuaciones no lineales** $\rightarrow$ **Linealización (jacobiano)** $\rightarrow$ **Espacio de estados $(A,B,C,D)$** $\rightarrow$ **Análisis espectral (eigenvalores)** $\rightarrow$ **Verificación de controlabilidad/observabilidad (rango)** $\rightarrow$ **Diseño del controlador (Riccati / asignación de polos)** $\rightarrow$ **Estabilización**

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
