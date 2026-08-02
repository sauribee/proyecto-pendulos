# Segunda entrega: plan de extensión del proyecto

Documento de planeación de la segunda entrega. Propone una **tercera
configuración** (péndulo triple sobre carro, estado en $\mathbb{R}^8$) y un
catálogo de **módulos independientes** para robustecer lo ya implementado, con
el énfasis puesto en una funcionalidad nueva: un **atlas de operabilidad** que
barre el espacio de parámetros para identificar dónde el controlador funciona
bien, dónde funciona a duras penas y dónde ya es imposible controlar el sistema.

- **Curso:** Álgebra Lineal Aplicada. Universidad Nacional de Colombia, Sede Medellín, Facultad de Ciencias.
- **Autores:** Mateo Bedoya Rojas, Camilo Alejandro Patiño Osorio, Santiago Uribe Echavarría.
- **Base:** primera entrega (Configuraciones I y II), documentada en `docs/guia_maestra.md`.

---

## Índice

0. [Cómo usar este documento](#0-cómo-usar-este-documento)
1. [Punto de partida: qué hay y qué falta](#1-punto-de-partida-qué-hay-y-qué-falta)
2. [Configuración III: el péndulo triple sobre carro](#2-configuración-iii-el-péndulo-triple-sobre-carro)
3. [El atlas de operabilidad: barridos y fronteras](#3-el-atlas-de-operabilidad-barridos-y-fronteras)
4. [Catálogo de módulos para robustecer](#4-catálogo-de-módulos-para-robustecer)
5. [Rutas recomendadas](#5-rutas-recomendadas)
6. [Prompts listos para usar](#6-prompts-listos-para-usar)
7. [Impacto en los documentos LaTeX](#7-impacto-en-los-documentos-latex)
8. [Checklist de la segunda entrega](#8-checklist-de-la-segunda-entrega)

---

## 0. Cómo usar este documento

El documento está escrito como un **menú, no como una secuencia obligatoria**.
Cada bloque de trabajo es un *módulo* independiente identificado como `M1`, `M2`,
etc. De cada módulo se da:

- **Qué es** y **por qué vale la pena** (el argumento de álgebra lineal que aporta).
- **La matemática** que hay detrás, desarrollada, no citada.
- **La API propuesta** (nombres de archivos, módulos y funciones) coherente con lo existente.
- **Esfuerzo estimado** y **dependencias** con otros módulos.
- **El prompt** listo para copiar y pegar.

La sección 5 propone tres rutas ya armadas (mínima, media, ambiciosa) por si se
prefiere no escoger módulo por módulo.

### 0.1 Los números de este documento están verificados

Todas las cifras de la Configuración III que aparecen aquí (espectro, rangos,
márgenes, ganancia $K$, polos de lazo cerrado) **se calcularon** con un
prototipo independiente que reproduce, con el mismo constructor genérico, los
resultados publicados de las Configuraciones I y II:

| Validación | Espectro obtenido | Coincide con `README.md` |
|---|---|---|
| Config I | $\{+4.2105,\ 0,\ -0.0769,\ -4.2266\}$ | Sí |
| Config I, $K$ | $(-3.16,\ -4.69,\ -45.39,\ -10.93)$ | Sí |
| Config II | $\{+8.5726,\ +4.0941,\ 0,\ 0,\ -4.0941,\ -8.5726\}$ | Sí |
| Config II, $K$ | $(3.16,\ 5.82,\ -191.55,\ -10.99,\ 228.32,\ 36.14)$ | Sí |

Es decir: las fórmulas generales de la sección 2 **contienen como casos
particulares** a los dos modelos ya implementados. Los valores de la
Configuración III que se dan en 2.7 son, por tanto, **metas de verificación**:
si la implementación no los reproduce, hay un error.

Lo que **no** está verificado y queda explícitamente marcado como *a medir*: todo
lo que requiere integrar el modelo no lineal (tiempos de asentamiento, fuerza
pico, ángulo máximo recuperable). Para eso hace falta `DifferentialEquations`, y
esos números son precisamente el producto de los módulos que se proponen.

### 0.2 Convenciones del proyecto que hay que respetar

Estas son las convenciones **observadas en el código actual**, no preferencias
generales. Cualquier código nuevo debe seguirlas:

| Aspecto | Convención vigente |
|---|---|
| Ángulos | $\theta_j$ medidos desde la **vertical superior**; $\theta_j = 0$ es el equilibrio erguido |
| Orden del estado | **Intercalado**: $[x,\ \dot x,\ \theta_1,\ \omega_1,\ \theta_2,\ \omega_2,\ \dots]$ |
| Parámetros | `Base.@kwdef struct SystemParamsX` con unidades SI en comentario de línea |
| Constructor | `default_params_x()` que devuelve los valores del informe |
| Dinámica | `nonlinear_eom_x!(dx, x, p, t)` con `p.params` y `p.F` (o `p.controller`) |
| Lazo cerrado | `closed_loop_eom_x!` que **reutiliza** la EOM no lineal, con `p.K` y `p.saturate` opcional |
| Forma mecánica | Armar $\mathbf{M}(q)$ explícita y resolver `cholesky(Symmetric(Mq)) \ rhs` |
| Linealización | `linearize_system_x(params) -> StateSpaceModel`, jacobiano **analítico** |
| Genericidad | `Controller` y los análisis de `Linearization` **no** dependen de $n$ |
| Comentarios | En **español sin tildes ni eñes**, sin emojis ni símbolos decorativos |
| Encabezados | Banner `# ===...===` con título, descripción y variables de estado |
| Salidas | Figuras a `figures/` (ignorada por git), numeradas `01_`, `02_`, ... |
| Scripts de doc | `Pkg.activate(PROJ_ROOT)` + `include` de `src/`; nunca duplican física |
| Notación en LaTeX y markdown | Todo símbolo matemático en LaTeX, sin excepción |
| Commits | En español, `tipo: descripcion` (feat, fix, refactor, docs, test) |

> **Nota sobre los comentarios.** El código actual usa comentarios en español sin
> tildes (`# Matriz de masa M(q), con q = [pos, theta1, theta2]`). Se mantiene esa
> convención por coherencia interna del repositorio, aunque difiera de la
> preferencia global de comentar en inglés. Cambiarla ahora obligaría a reescribir
> los seis módulos existentes.

---

## 1. Punto de partida: qué hay y qué falta

### 1.1 Lo que ya está resuelto

La primera entrega cubre, para dos configuraciones, la cadena completa:

$$\text{Euler-Lagrange} \rightarrow \mathbf{M}(q)\ddot q = \text{rhs} \rightarrow \text{linealización} \rightarrow (A,B,C,D) \rightarrow \sigma(A),\ \operatorname{rank}\mathcal{C},\ \operatorname{rank}\mathcal{O} \rightarrow K \rightarrow \text{simulación no lineal}$$

con dos vías de diseño (LQR vía Riccati resuelta por el Hamiltoniano, y
Ackermann), animaciones, notebooks reactivos y cuatro artefactos de
documentación. El diseño del código es sólido: la física está separada del
álgebra, y `Controller` es genérico en la dimensión del estado.

### 1.2 Las cinco grietas

Un lector crítico del informe puede señalar estos huecos. Cada uno es una
oportunidad de extensión, y están ordenados por lo llamativo que resulta el
hueco:

1. **La observabilidad se analiza pero no se usa.** Se calcula
   $\operatorname{rank}\mathcal{O} = n$ en ambas configuraciones y ahí termina.
   No hay observador, y sin embargo el control implementado es $u = -K\mathbf{x}$
   con **estado completo**, que en un sistema real no se mide. Es la grieta más
   visible: se demuestra que se *puede* reconstruir el estado y no se reconstruye.

2. **No hay noción de límite de validez.** Se reporta que el LQR estabiliza desde
   $\theta_0 = 0.40$ rad, pero no *hasta dónde* llega, ni con qué actuador, ni
   para qué parámetros. La frase "el diseño lineal es robusto más allá de donde la
   teoría lo garantiza" es cierta pero no está cuantificada.

3. **El criterio de rango es frágil y no se dice.** `rank(C_ctrl)` sobre la matriz
   de Kalman es una decisión binaria tomada con una tolerancia numérica implícita.
   Como se muestra en 3.3, para el triple $\operatorname{cond}(\mathcal{C}) \approx
   2\times 10^{8}$ y para cinco eslabones $\approx 3.3\times 10^{16}$: en ese punto el
   rango numérico **deja de significar algo**. Hay que sustituirlo (o acompañarlo)
   por una medida continua.

4. **No hay pruebas automatizadas.** Los jacobianos analíticos de
   `linearize_system` y `linearize_system_double` están escritos a mano. Coinciden
   con el modelo no lineal (lo verifiqué), pero nada en el repositorio lo comprueba
   de forma reproducible.

5. **Solo dos puntos en el eje de complejidad.** El mensaje "el mismo razonamiento
   algebraico escala" se sostiene sobre dos casos. Con tres, la tendencia deja de
   ser anecdótica y se vuelve una curva.

La propuesta de esta segunda entrega ataca las cinco: la Configuración III cubre
la grieta 5, el atlas de operabilidad cubre las grietas 2 y 3, y los módulos del
catálogo cubren las grietas 1 y 4.

---

## 2. Configuración III: el péndulo triple sobre carro

### 2.1 Por qué el triple y no otra cosa

Se consideraron cuatro candidatos:

| Candidato | A favor | En contra |
|---|---|---|
| **Péndulo triple** | Escala la misma maquinaria a $\mathbb{R}^8$; 3 modos inestables; la degradación numérica se vuelve visible y **medible**; permite la comparación limpia I-II-III | Ninguno serio: es el paso natural |
| Péndulo de Furuta (rotacional) | Geometría distinta, atractivo visual | Sigue en $\mathbb{R}^4$; obliga a rehacer toda la cinemática sin ganar dimensión |
| Acrobot / Pendubot | Subactuación más severa | Sin carro: rompe la continuidad con I y II |
| Doble con barras flexibles | Muy realista | Introduce EDPs o modos modales: se sale del alcance del curso |

Gana el triple, y por una razón que va más allá de "es más grande": es la
configuración donde **el criterio de rango se rompe**. Con
$\operatorname{cond}(\mathcal{C}) \approx 2\times 10^{8}$, decir
"$\operatorname{rank}(\mathcal{C}) = 8$, luego es controlable" es formalmente
correcto pero prácticamente vacío. Eso obliga a introducir el **test de
Popov-Belevitch-Hautus** y la **distancia a la incontrolabilidad**, que son
álgebra lineal de primera categoría (valores singulares, perturbación de
matrices) y que no aparecen en la primera entrega. La Configuración III no es
más de lo mismo: es la excusa para una herramienta nueva.

### 2.2 Coordenadas, espacios y convenciones

Tres eslabones articulados en serie sobre el carro. Coordenadas generalizadas:

$$q = (x,\ \theta_1,\ \theta_2,\ \theta_3) \in \mathbb{R}^4, \qquad d = 4 \text{ grados de libertad}$$

- $x$: posición horizontal del carro.
- $\theta_j$: ángulo del eslabón $j$ medido **desde la vertical superior**, igual
  que en I y II. Son ángulos **absolutos** (respecto a la vertical), no relativos
  al eslabón anterior. Esta es la misma convención de `model_double.jl` y hay que
  conservarla.

El espacio de estados es el fibrado tangente $T\mathcal{Q} \cong \mathbb{R}^8$,
con el orden intercalado del proyecto:

$$\mathbf{x} = (x,\ \dot x,\ \theta_1,\ \omega_1,\ \theta_2,\ \omega_2,\ \theta_3,\ \omega_3)^\top \in \mathbb{R}^8$$

La entrada sigue siendo escalar: $u \in \mathbb{R}$, la fuerza horizontal sobre
el carro. **Cuatro grados de libertad, un actuador**: la subactuación pasa de
$2:1$ (Config. I) a $3:1$ (Config. II) a $4:1$.

Los eslabones se modelan como **barras uniformes** con inercia, no como masas
puntuales. Esto es deliberado: unifica el tratamiento de la Configuración I
(barra con $I = \tfrac{1}{12}m\ell^2$) con la estructura multi-eslabón de la
Configuración II (que usa masas puntuales), y hace de la Configuración III la
generalización que las contiene a ambas.

Parámetros de cada eslabón $j \in \{1,2,3\}$:

| Símbolo | Significado |
|---|---|
| $m_j$ | masa del eslabón [kg] |
| $\ell_j$ | longitud del eslabón [m] |
| $a_j$ | distancia de la articulación inferior al centro de masa [m] |
| $I_j$ | momento de inercia respecto al centro de masa [kg·m$^2$] |

Para una barra uniforme: $a_j = \ell_j/2$ y $I_j = \tfrac{1}{12}m_j\ell_j^2$.
Para una masa puntual en el extremo: $a_j = \ell_j$ y $I_j = 0$. Un solo modelo
sirve para los dos casos, que es justo lo que permite recuperar I y II.

### 2.3 Energías y lagrangiano

Posición del centro de masa del eslabón $j$ (acumulando los eslabones inferiores):

$$x_{c_j} = x + \sum_{i<j}\ell_i\sin\theta_i + a_j\sin\theta_j, \qquad
  y_{c_j} = \sum_{i<j}\ell_i\cos\theta_i + a_j\cos\theta_j$$

Derivando y agrupando, la energía cinética es

$$T = \tfrac12 M\dot x^2 + \sum_{j}\left[\tfrac12 m_j\left(\dot x_{c_j}^2 + \dot y_{c_j}^2\right) + \tfrac12 I_j\dot\theta_j^2\right]$$

y la potencial, $V = g\sum_j m_j y_{c_j}$.

Al expandir aparecen tres agrupaciones de parámetros que **organizan todo el
modelo**. Definimos, para $j = 1,\dots,N$:

$$\boxed{\ \beta_j = m_j a_j + \ell_j\!\!\sum_{i>j}m_i\ },\qquad
  \boxed{\ J_j = I_j + m_j a_j^2 + \ell_j^2\!\!\sum_{i>j}m_i\ },\qquad
  \boxed{\ \Gamma_{jk} = \ell_{\min(j,k)}\,\beta_{\max(j,k)}\ }$$

Interpretación física:

- $\beta_j$ es el **momento estático** del eslabón $j$ y de todo lo que carga
  encima: la masa propia por su brazo, más la longitud completa del eslabón por
  la masa de los que están más arriba.
- $J_j$ es el **momento de inercia efectivo** del eslabón $j$ respecto a su
  articulación, incluyendo la carga de los superiores (Steiner más el resto).
- $\Gamma_{jk}$ es el **acoplamiento inercial** entre los eslabones $j$ y $k$.

Con estas definiciones la potencial se escribe de forma notablemente compacta:

$$V = g\sum_{j=1}^{N}\beta_j\cos\theta_j$$

Vale la pena detenerse aquí: **el mismo $\beta_j$ que acopla el carro con el
eslabón $j$ en la energía cinética es el que pesa la gravedad sobre ese
eslabón**. Esa coincidencia no es casual (ambas provienen del mismo momento
estático) y es la razón estructural de que la inestabilidad y la
controlabilidad estén ligadas: el canal por el que la gravedad desestabiliza es
el mismo por el que el carro puede actuar.

### 2.4 Las ecuaciones de movimiento

Aplicando $\frac{d}{dt}\frac{\partial L}{\partial\dot q_i} - \frac{\partial L}{\partial q_i} = Q_i$ con $L = T - V$:

**Ecuación del carro:**

$$\left(M + \sum_{i=1}^{N}m_i\right)\ddot x + \sum_{j=1}^{N}\beta_j\left(\cos\theta_j\,\ddot\theta_j - \sin\theta_j\,\dot\theta_j^2\right) + b\,\dot x = u$$

**Ecuación del eslabón $j$:**

$$\beta_j\cos\theta_j\,\ddot x + J_j\ddot\theta_j + \sum_{k\neq j}\Gamma_{jk}\left[\cos(\theta_j-\theta_k)\,\ddot\theta_k + \sin(\theta_j-\theta_k)\,\dot\theta_k^2\right] = g\beta_j\sin\theta_j - \tau_j^{\text{fric}}$$

En forma mecánica estándar $\mathbf{M}(q)\ddot q + \mathbf{C}(q,\dot q)\dot q + \mathbf{G}(q) = \mathbf{F}u$, la matriz de masa es

$$\mathbf{M}(q) = \begin{pmatrix}
M + \sum_i m_i & \beta_1\cos\theta_1 & \cdots & \beta_N\cos\theta_N \\
\beta_1\cos\theta_1 & J_1 & \cdots & \Gamma_{1N}\cos(\theta_1-\theta_N) \\
\vdots & \vdots & \ddots & \vdots \\
\beta_N\cos\theta_N & \Gamma_{N1}\cos(\theta_N-\theta_1) & \cdots & J_N
\end{pmatrix}$$

simétrica y definida positiva (es la Hessiana de $T$ respecto a $\dot q$), lo que
justifica seguir usando `cholesky(Symmetric(Mq)) \ rhs` como en
`model_double.jl`.

El **signo positivo** de $g\beta_j\sin\theta_j$ en el lado derecho es la firma de
la inestabilidad, exactamente igual que el $+(M+m)mgL/D_0$ de la Configuración I:
para $\theta_j$ pequeño el término es $+g\beta_j\theta_j$, una realimentación
positiva.

**Fricción.** Se propone dejarla parametrizada y **desactivada por defecto** en
las articulaciones ($c_j = 0$), conservando solo la fricción viscosa del carro
$b\dot x$ como en la Configuración I. Si se activa, lo físicamente correcto es
un par **relativo** entre eslabones contiguos:

$$\tau_j^{\text{fric}} = c_j(\dot\theta_j - \dot\theta_{j-1}) - c_{j+1}(\dot\theta_{j+1} - \dot\theta_j), \qquad \dot\theta_0 := 0$$

que produce una matriz de amortiguamiento tridiagonal simétrica.

### 2.5 La formulación genérica en $N$ y por qué contiene a I y II

Las fórmulas de 2.3 y 2.4 están escritas para $N$ eslabones. Instanciándolas:

**$N = 1$, barra uniforme.** Con la notación de la Configuración I, $L$ es la
distancia al centro de masa, de modo que la longitud completa de la barra es
$\ell_1 = 2L$; entonces $a_1 = \ell_1/2 = L$ e
$I_1 = \tfrac{1}{12}m\ell_1^2 = \tfrac{1}{12}m(2L)^2$, y no hay eslabones superiores:

$$\beta_1 = mL, \qquad J_1 = I + mL^2$$

y la ecuación del carro queda $(M+m)\ddot x + mL(\cos\theta\,\ddot\theta - \sin\theta\,\dot\theta^2) + b\dot x = u$, con
$mL\cos\theta\,\ddot x + (I+mL^2)\ddot\theta = mgL\sin\theta$. **Es exactamente
`model_simple.jl`**, línea por línea.

**$N = 2$, masas puntuales** ($a_j = \ell_j$, $I_j = 0$):

$$\beta_1 = m_1L_1 + L_1m_2 = (m_1+m_2)L_1, \quad \beta_2 = m_2L_2$$
$$J_1 = m_1L_1^2 + L_1^2m_2 = (m_1+m_2)L_1^2, \quad J_2 = m_2L_2^2, \quad \Gamma_{12} = m_2L_1L_2$$

que reproduce **término a término** la `Mq` y el `rhs` de `model_double.jl`.

Esto tiene una consecuencia práctica de primer orden: **un módulo genérico en $N$
puede validarse contra los dos modelos ya entregados**. No es un refactor que
haya que creer; es un refactor que se puede *demostrar* correcto con dos pruebas
numéricas. Y esa demostración es material de informe: muestra que las tres
configuraciones no son tres problemas, sino tres puntos de una misma familia.

### 2.6 Linealización

Evaluando en $\theta_j = 0$, $\dot q = 0$, $u = 0$: $\cos\to 1$, $\sin\theta_j\to\theta_j$,
$\cos(\theta_j-\theta_k)\to 1$, y los términos cuadráticos en $\dot\theta$
desaparecen (son de segundo orden). Queda un sistema lineal con **matriz de masa
constante**, que es el "truco decisivo" ya documentado en la guía maestra
(sección 5.6):

$$\mathbf{M}_0\,\ddot q = \mathbf{G}_0\,q - \mathbf{F}_0\,\dot q + \mathbf{e}_1 u$$

con $\mathbf{M}_0 = \mathbf{M}(0)$ constante,
$\mathbf{G}_0 = g\operatorname{diag}(0,\ \beta_1,\ \dots,\ \beta_N)$ y
$\mathbf{F}_0$ la matriz de fricción. Despejando, en orden **por bloques**
$(q,\dot q)$:

$$A_{\text{bloques}} = \begin{pmatrix}\mathbf{0} & \mathbf{I}\\[2pt] \mathbf{M}_0^{-1}\mathbf{G}_0 & -\mathbf{M}_0^{-1}\mathbf{F}_0\end{pmatrix},
\qquad
B_{\text{bloques}} = \begin{pmatrix}\mathbf{0}\\[2pt] \mathbf{M}_0^{-1}\mathbf{e}_1\end{pmatrix}$$

El proyecto usa orden **intercalado**, que es esta misma matriz conjugada por una
permutación $\Pi$: $A = \Pi A_{\text{bloques}}\Pi^\top$. El espectro no cambia
(matrices semejantes), pero el código debe ensamblar directamente la versión
intercalada para ser coherente con `linearize_system_double`.

**Predicción estructural del espectro.** Si no hay fricción ($\mathbf{F}_0 = 0$),
buscar soluciones $q(t) = ve^{\lambda t}$ da el problema generalizado

$$\mathbf{G}_0\,v = \lambda^2\,\mathbf{M}_0\,v$$

de modo que los eigenvalores de $A$ son las **raíces cuadradas** de los
eigenvalores de $\mathbf{M}_0^{-1}\mathbf{G}_0$, y por tanto vienen en pares
$\pm\lambda$. Como $\mathbf{G}_0$ tiene una fila y una columna nulas (el carro no
tiene fuerza recuperadora), $\mathbf{M}_0^{-1}\mathbf{G}_0$ es singular y aparece
$\lambda = 0$ doble, en un **bloque de Jordan** $2\times 2$. Para $N$ eslabones:

$$\sigma(A) = \{\pm\lambda_1,\ \pm\lambda_2,\ \dots,\ \pm\lambda_N,\ 0,\ 0\}, \qquad N \text{ modos inestables}$$

Esto se cumple exactamente en la Configuración II
($\{\pm 4.09, \pm 8.57, 0, 0\}$) y sirve como **verificación a priori** de la
implementación del triple. La fricción del carro rompe ligeramente la simetría y
desdobla el cero: en la Configuración III con $b = 0.1$ se obtiene
$\{0,\ -0.0625\}$ en lugar de $\{0, 0\}$, y cada par $\pm\lambda_j$ queda
desbalanceado en la tercera cifra. Es exactamente lo que ya pasa en la
Configuración I, donde el par $\{+4.21,\ -4.23\}$ está desbalanceado por el mismo
motivo.

### 2.7 Parámetros propuestos y resultados esperados

**Parámetros por defecto de la Configuración III.** Se eligen para que la
comparación con la Configuración II sea *limpia*: misma masa total de péndulo
(0.6 kg) y misma longitud total (1 m), redistribuidas en tres barras uniformes.

| Parámetro | Valor | Justificación |
|---|---|---|
| $M$ | 1.0 kg | Igual que I y II |
| $m_1 = m_2 = m_3$ | 0.2 kg | Total 0.6 kg, igual que la Config. II |
| $\ell_1 = \ell_2 = \ell_3$ | $1/3$ m | Total 1 m, igual que I y II |
| $a_j$ | $\ell_j/2 = 1/6$ m | Barra uniforme |
| $I_j$ | $\tfrac{1}{12}m_j\ell_j^2 = 1.852\times10^{-3}$ kg·m$^2$ | Barra uniforme |
| $g$ | 9.81 m/s$^2$ | Igual que I y II |
| $b$ | 0.1 N·s/m | Igual que la Config. I |
| $c_j$ | 0 | Sin fricción articular por defecto |

Con ellos, $\beta = (0.16667,\ 0.1,\ 0.03333)$ y $J = (0.05185,\ 0.02963,\ 0.00741)$.

**Metas de verificación (calculadas, deben reproducirse):**

Espectro de lazo abierto:

$$\sigma(A) = \{+18.0613,\ +9.7740,\ +4.3777,\ 0,\ -0.0625,\ -4.3990,\ -9.7813,\ -18.0647\}$$

**Tres modos inestables**, con la estructura casi antisimétrica predicha en 2.6.
El modo dominante $\lambda_{\max} = 18.06$ s$^{-1}$ implica una constante de
tiempo de caída de $1/\lambda_{\max} = 55$ ms: el sistema se desploma cuatro veces
más rápido que la Configuración I.

Rangos y márgenes:

| Cantidad | Valor |
|---|---|
| $\operatorname{rank}\mathcal{C}$ | $8/8$ |
| $\operatorname{rank}\mathcal{O}$ (midiendo $x,\theta_1,\theta_2,\theta_3$) | $8/8$ |
| $\operatorname{cond}(\mathcal{C})$ | $2.157\times 10^{8}$ |
| Margen PBH normalizado (ver 3.2) | $4.916\times 10^{-4}$ |

LQR con $Q = \operatorname{diag}(1,0,10,0,10,0,10,0)$ y $R = 0.1$ (los mismos
pesos de las otras dos configuraciones, extendidos):

$$K = (-3.16,\ -6.17,\ -317.43,\ -4.37,\ 911.95,\ 37.73,\ -667.37,\ -61.42), \qquad \|K\|_2 = 1176.0$$

$$\sigma(A - BK) = \{-0.84\pm0.79i,\ -4.57\pm2.32i,\ -10.02\pm1.73i,\ -18.14\pm0.97i\}$$

Nótese que el polo más rápido de lazo cerrado ($-18.14$) es prácticamente el
espejo del modo inestable más rápido de lazo abierto ($+18.06$): el LQR, con
$R$ pequeño, tiende a **reflejar** los polos inestables al semiplano izquierdo.
Es una propiedad conocida del regulador barato y se puede señalar en el informe.

Matrices para contraste directo (redondeadas a 4 cifras):

$$\mathbf{M}_0 = \begin{pmatrix}
1.6 & 0.16667 & 0.1 & 0.03333\\
0.16667 & 0.05185 & 0.03333 & 0.01111\\
0.1 & 0.03333 & 0.02963 & 0.01111\\
0.03333 & 0.01111 & 0.01111 & 0.00741
\end{pmatrix},
\qquad
B = \begin{pmatrix}0\\ 0.9455\\ 0\\ -3.6\\ 0\\ 0.9818\\ 0\\ -0.3273\end{pmatrix}$$

Las filas no triviales de $A$ (orden intercalado, columnas $[x,\dot x,\theta_1,\omega_1,\theta_2,\omega_2,\theta_3,\omega_3]$):

| Fila | $\dot x$ | $\theta_1$ | $\theta_2$ | $\theta_3$ |
|---|---|---|---|---|
| $\ddot x$ (fila 2) | $-0.0945$ | $-5.886$ | $+0.9632$ | $-0.1070$ |
| $\dot\omega_1$ (fila 4) | $+0.3600$ | $+141.264$ | $-95.3532$ | $+10.5948$ |
| $\dot\omega_2$ (fila 6) | $-0.0982$ | $-158.922$ | $+194.559$ | $-51.0477$ |
| $\dot\omega_3$ (fila 8) | $+0.0327$ | $+52.974$ | $-153.143$ | $+105.306$ |

**Lo que queda por medir** (requiere integrar el modelo no lineal): tiempo de
asentamiento, fuerza pico, excursión del carro y ángulo máximo recuperable. La
expectativa razonada, extrapolando de I a II: $t_s$ del orden de 2 s, fuerza pico
mayor que los 4.2 N del doble por la mayor exigencia de ancho de banda, y una
región de atracción **notablemente más pequeña**. Sobre esto último sí hay una
cota rigurosa, y está en 3.2.

### 2.8 Puntos delicados y riesgos

Vale la pena tenerlos presentes antes de empezar:

1. **La rigidez numérica sube.** Con $\lambda_{\max} = 18.06$ y polos de lazo
   cerrado en $-18.14$, el sistema es más rígido. `Tsit5` sigue sirviendo, pero
   conviene bajar `saveat` a 0.002 y vigilar `reltol`. Si aparecen problemas,
   `Rodas5` (implícito) es la alternativa; **no cambiar el solver por defecto** sin
   comprobar que `Tsit5` falla, para no romper la comparabilidad con I y II.

2. **La saturación se vuelve determinante.** Con $\|K\| = 1176$, una desviación de
   $\theta_2 = 0.05$ rad ya pide $911.95\times0.05 \approx 46$ N solo por ese
   término. El límite de 50 N usado en la Configuración I es **insuficiente** para
   el triple. Se propone $u_{\max} = 150$ N por defecto, y que el valor sea un
   parámetro del barrido, no una constante escondida.

3. **La animación necesita más cuidado.** Tres eslabones de $1/3$ m se ven
   pequeños; conviene escalar el dibujo y ajustar los límites del eje, siguiendo
   `link_positions` de `animation_double.jl`, que ya generaliza sin problema.

4. **No confiar en `rank`.** Con $\operatorname{cond}(\mathcal{C}) \approx 2\times10^8$,
   `rank(C_ctrl)` con la tolerancia por defecto todavía da 8, pero es el último
   caso en que se puede confiar. Ver la sección 3.

5. **El LQR con $R = 0.1$ pide mucho ancho de banda.** Si en el barrido se quiere
   un controlador más suave, subir $R$ y documentarlo. Mantener $R = 0.1$ para el
   caso nominal, por comparabilidad.

### 2.9 Qué archivos se crean y cuáles se tocan

Propuesta de estructura, con la decisión arquitectónica explicada abajo:

```
src/
  model_nlink.jl          NUEVO   Modelo generico de N eslabones (modulo ModelNLink)
  model_triple.jl         NUEVO   Instancia N=3 (modulo ModelTriple), API del proyecto
  linearization.jl        TOCADO  Se agrega linearize_system_triple (y _nlink)
  animation_triple.jl     NUEVO   Animacion de 3 eslabones (modulo AnimationTriple)
  controller.jl           INTACTO Ya es generico en n
  model_simple.jl         INTACTO
  model_double.jl         INTACTO
main_triple.jl            NUEVO   Pipeline completo de la Configuracion III
test/runtests.jl          NUEVO   Verifica que ModelNLink reproduce I, II y III
```

**Decisión arquitectónica recomendada: genérico *aditivo*, no refactor
destructivo.** Es decir: crear `model_nlink.jl` como módulo nuevo y **no tocar**
`model_simple.jl` ni `model_double.jl`. Razones:

- La primera entrega ya está calificada y compilada; cualquier cambio en esos dos
  archivos obliga a revalidar los cuatro documentos.
- Al dejar los tres modelos coexistiendo, la comparación `ModelNLink(N=1)` contra
  `Model` se convierte en una **prueba de regresión genuina**: dos
  implementaciones independientes que deben coincidir. Si se refactoriza
  destructivamente, esa prueba desaparece.
- `model_triple.jl` queda como una capa delgada sobre `ModelNLink` que expone la
  API del proyecto (`SystemParamsTriple`, `default_params_triple`,
  `nonlinear_eom_triple!`, ...), de modo que `main_triple.jl` se lea igual que
  `main_double.jl`.

La alternativa conservadora (escribir `model_triple.jl` autónomo, calcado de
`model_double.jl` con una $\mathbf{M}(q)$ de $4\times4$ escrita a mano) es
perfectamente válida y **más rápida**, a costa de duplicar física y perder la
prueba cruzada. Si el tiempo aprieta, es una retirada honorable.

---

## 3. El atlas de operabilidad: barridos y fronteras

Esta es la funcionalidad central que se pide: **barrer el espacio de parámetros
para decir dónde el controlador sirve y dónde ya no**.

### 3.1 Las tres preguntas

1. **¿Dónde funciona bien?** ¿Existe una configuración de parámetros que maximice
   el margen de control? (Adelanto: **sí**, y no está en los extremos.)
2. **¿Dónde funciona a duras penas?** ¿Cómo se degrada la capacidad de control al
   alejarse de ese óptimo?
3. **¿Dónde es imposible?** ¿Cuál es la frontera —en masa, en longitud, en fuerza
   de actuador, en ángulo inicial, en período de muestreo— más allá de la cual no
   hay controlador que sirva?

Conviene distinguir **dos tipos de imposibilidad**, porque el informe debe
separarlas:

- **Imposibilidad estructural:** el par $(A,B)$ deja de ser controlable (o se
  acerca tanto que ningún $K$ razonable sirve). Es una propiedad **del sistema**,
  independiente del controlador. Se mide con el margen PBH.
- **Imposibilidad práctica:** el sistema es controlable, pero el actuador se
  satura, el riel se acaba o el período de muestreo es muy largo. Es una
  propiedad **de la implementación**. Se mide con simulación no lineal y
  bisección.

### 3.2 Las métricas, una por una

#### (a) Modo inestable dominante

$$\lambda_{\max} = \max_i \operatorname{Re}\lambda_i(A), \qquad \tau_{\text{caída}} = 1/\lambda_{\max}$$

Es la escala de tiempo del problema: el controlador debe actuar en una fracción
de $\tau_{\text{caída}}$. Barato de calcular y muy interpretable.

#### (b) Rango de Kalman, y por qué no basta

$\operatorname{rank}\mathcal{C} = n$ es una condición **binaria** evaluada con una
tolerancia numérica. El problema es que $\mathcal{C} = [B, AB, \dots, A^{n-1}B]$
es una base de Krylov, notoriamente mal condicionada. Los números medidos sobre
la **familia normalizada** (longitud total 1 m, masa total 0.6 kg repartidas en
$N$ barras uniformes, $M=1$ kg, $b=0.1$), que permite comparar $N$ contra $N$:

| $N$ | $n$ | $\operatorname{cond}(\mathcal{C})$ |
|---|---|---|
| 1 | 4 | $4.8\times10^{1}$ |
| 2 | 6 | $6.2\times10^{4}$ |
| **3** | **8** | $2.2\times10^{8}$ |
| 4 | 10 | $1.7\times10^{12}$ |
| 5 | 12 | $3.3\times10^{16}$ |

Con $N=5$ el condicionamiento alcanza $1/\varepsilon_{\text{máq}}$: el rango
numérico es **ruido**. Y sin embargo el sistema *sigue siendo* controlable en
aritmética exacta. Esta tabla, por sí sola, es un resultado publicable en el
informe: muestra experimentalmente por qué el criterio de rango, siendo un
teorema correcto, es un mal algoritmo.

#### (c) Margen de Popov-Belevitch-Hautus (la métrica principal)

El test PBH dice: $(A,B)$ es controlable si y solo si

$$\operatorname{rank}\big[\,A - \lambda I \ \big|\ B\,\big] = n \quad \text{para todo } \lambda\in\mathbb{C}$$

y basta comprobarlo en $\lambda \in \sigma(A)$. La versión **cuantitativa** susti-
tuye el rango por el menor valor singular:

$$\mu_{\text{PBH}} = \min_{\lambda\in\sigma(A)}\ \sigma_{\min}\big[\,A-\lambda I\ \big|\ B\,\big]$$

Esto ya no es binario: mide **cuánto habría que perturbar** el sistema para
volverlo incontrolable (por el teorema de Eckart-Young, $\sigma_{\min}$ es
exactamente la distancia en norma 2 al conjunto de matrices de rango deficiente).
Para que sea comparable entre configuraciones de distinta escala, se normaliza:

$$\hat\mu = \frac{\mu_{\text{PBH}}}{\big\|[\,A\ \ B\,]\big\|_2}$$

Valores medidos (misma familia: longitud total 1 m, masa total 0.6 kg):

| $N$ | $\lambda_{\max}$ | $\tau$ [ms] | $\hat\mu$ | $\operatorname{cond}\mathcal{C}$ | $\|K\|$ |
|---|---|---|---|---|---|
| 1 | 4.51 | 222 | $1.31\times10^{-2}$ | $4.8\times10^{1}$ | 52.9 |
| 2 | 10.59 | 94 | $1.79\times10^{-3}$ | $6.2\times10^{4}$ | 248.7 |
| **3** | **18.06** | **55** | $4.92\times10^{-4}$ | $2.2\times10^{8}$ | 1176.0 |
| 4 | 26.54 | 38 | $1.90\times10^{-4}$ | $1.7\times10^{12}$ | 5266.8 |
| 5 | 35.64 | 28 | $9.03\times10^{-5}$ | $3.3\times10^{16}$ | 23221.9 |

A diferencia de $\operatorname{cond}(\mathcal{C})$, que estalla, $\hat\mu$ decae
suavemente: es una métrica **utilizable** para hacer mapas de calor.

> **Nota técnica.** $\mu_{\text{PBH}}$ evaluado solo sobre $\sigma(A)$ es una cota
> superior de la verdadera *distancia a la incontrolabilidad* de Eising,
> $\mu_E = \min_{s\in\mathbb{C}}\sigma_{\min}[A-sI \mid B]$, cuyo mínimo puede
> caer fuera del espectro. Para el barrido, la versión sobre el espectro es
> suficiente y cuesta $n$ descomposiciones SVD. Si se quiere $\mu_E$, se refina
> con una malla en $\mathbb{C}$ alrededor de los eigenvalores; se sugiere como
> extensión opcional, no como requisito.

#### (d) Costo óptimo y esfuerzo

$$J^\star(\mathbf{x}_0) = \mathbf{x}_0^\top P\,\mathbf{x}_0, \qquad \|K\|_2$$

donde $P$ es la solución de la CARE que ya calcula `solve_care`. $J^\star$ es el
costo mínimo alcanzable desde $\mathbf{x}_0$: una medida directa y sin ambigüedad
de "qué tan difícil es". $\|K\|$ mide la amplificación del ruido de sensores.

> **Advertencia:** el gramiano de controlabilidad de horizonte infinito
> $W_c = \int_0^\infty e^{At}BB^\top e^{A^\top t}dt$ **no existe** aquí, porque $A$
> es inestable y la integral diverge. Es un error fácil de cometer. Si se quiere
> un gramiano, hay que usar horizonte finito $W_c(T)$ (calculable con el truco de
> Van Loan sobre una matriz $2n\times2n$) o el gramiano de **lazo cerrado**, con
> $A-BK$, que sí es Hurwitz. Recomendación: usar $J^\star$ y evitar el asunto.

#### (e) El elipsoide de no saturación (cota rigurosa y barata)

Esta es, en mi opinión, la pieza más elegante que puede aportar la segunda
entrega, porque es **álgebra lineal pura que produce una garantía dura**.

Con $V(\mathbf{x}) = \mathbf{x}^\top P\mathbf{x}$ y el control $u = -K\mathbf{x}$,
el mayor conjunto de nivel $\mathcal{E}_c = \{\mathbf{x} : \mathbf{x}^\top P\mathbf{x}\le c\}$
que cabe dentro de la franja de no saturación $\{|K\mathbf{x}|\le u_{\max}\}$ se
obtiene resolviendo

$$\max_{\mathbf{x}^\top P\mathbf{x}\le c}|K\mathbf{x}| = \sqrt{c\;K P^{-1}K^\top} \;\le\; u_{\max}
\qquad\Longrightarrow\qquad
\boxed{\ c^\star = \frac{u_{\max}^2}{K P^{-1}K^\top}\ }$$

(el máximo de una forma lineal sobre un elipsoide es la norma dual, un ejercicio
de Cauchy-Schwarz). Dentro de $\mathcal{E}_{c^\star}$ el control **nunca satura**
y, para el sistema lineal, $\dot V = -\mathbf{x}^\top(Q + K^\top RK)\mathbf{x} < 0$,
de modo que $\mathcal{E}_{c^\star}$ es **invariante** y todas las trayectorias que
empiezan ahí convergen.

Además, como $K = R^{-1}B^\top P$, la expresión se simplifica a
$KP^{-1}K^\top = R^{-1}B^\top P B R^{-1}$, que para $R$ escalar es
$B^\top P B / R^2$: no hace falta invertir $P$.

Particularizando a la condición inicial estándar del proyecto,
$\mathbf{x}_0 = \theta_0\,\mathbf{e}_3$ (solo el primer eslabón desviado), la
condición $\mathbf{x}_0\in\mathcal{E}_{c^\star}$ da

$$|\theta_0| \;\le\; \sqrt{c^\star / P_{33}}$$

Con $u_{\max} = 50$ N y los pesos nominales, los valores calculados son:

| Configuración | $P_{33}$ | $c^\star$ | $\theta_0$ garantizado sin saturar |
|---|---|---|---|
| Simple | 24.80 | 22.94 | $0.962$ rad $\approx 55.1°$ |
| Doble | 163.82 | 8.99 | $0.234$ rad $\approx 13.4°$ |
| **Triple** | **161.97** | **3.73** | $\mathbf{0.152}$ **rad** $\approx \mathbf{8.7°}$ |

Una progresión $55° \to 13° \to 9°$ que cuantifica exactamente lo que en la
primera entrega solo se podía decir de forma cualitativa.

> **Honestidad sobre el alcance.** $\mathcal{E}_{c^\star}$ garantiza no saturación
> y convergencia **del sistema lineal**. Para el sistema **no lineal** hace falta,
> además, acotar el residuo de la linealización y exigir $\dot V<0$ con él
> incluido; el resultado es un $c$ menor. Por eso la cota elipsoidal se
> complementa con la bisección numérica del punto (f), que da la estimación
> práctica. Presentar $c^\star$ como si fuera la región de atracción no lineal
> sería un error, y probablemente una pregunta de defensa.

#### (f) Fronteras por bisección (la prueba de la verdad)

Tres barridos por bisección sobre el **modelo no lineal con saturación**, que es
donde se decide si el controlador realmente sirve:

**Criterio de éxito** (hay que fijarlo y documentarlo, porque toda la frontera
depende de él):

$$\text{éxito} \iff \Big(\max_j|\theta_j(T)| < 0.02\ \text{rad}\Big) \ \wedge\ \Big(\max_t |x(t)| \le x_{\text{riel}}\Big) \ \wedge\ \Big(\max_{t,j}|\theta_j(t)| < \pi/2\Big)$$

con $T = 10$ s y $x_{\text{riel}} = 1.5$ m. El umbral de 0.02 rad es el mismo de
`settling_time` en `make_report_figs.jl`: hay que reutilizarlo, no inventar otro.

1. **Ángulo máximo recuperable** $\theta_{\max}(p)$: bisección sobre $\theta_0$.
   Responde "¿desde qué tan inclinado puede recuperarse?".
2. **Actuador mínimo** $u_{\min}(p)$: bisección sobre $u_{\max}$ con $\theta_0$
   fijo. Responde "¿qué motor hay que comprar?".
3. **Período de muestreo máximo** $h_{\max}(p)$: bisección sobre $h$, discretizando
   con retenedor de orden cero. Responde "¿a qué frecuencia hay que muestrear?".

El tercero es especialmente bonito porque es **puro exponencial matricial**, y
conecta directamente con la sección 6.1 de la guía maestra. Con
$A_d = e^{Ah}$ y $B_d = \left(\int_0^h e^{A\tau}d\tau\right)B$, ambos obtenidos de
una sola exponencial de la matriz aumentada $\begin{pmatrix}A & B\\ 0 & 0\end{pmatrix}$
(truco de Van Loan), se pide que $A_d - B_dK$ sea **Schur** ($|\lambda_i|<1$),
que es el análogo discreto de Hurwitz. Valores medidos con los **parámetros por
defecto** de cada configuración:

| Configuración | $h_{\max}$ [ms] | $1/\lambda_{\max}$ [ms] | $h_{\max}\cdot\lambda_{\max}$ |
|---|---|---|---|
| Simple | 197.4 | 237.5 | 0.83 |
| Doble | 77.3 | 116.7 | 0.66 |
| **Triple** | **32.0** | **55.4** | **0.58** |

Es decir: el triple exige muestrear a más de **31 Hz**, y el margen relativo
respecto al tiempo de caída se estrecha ($0.83 \to 0.66 \to 0.58$). Un resultado
concreto, cuantitativo y con lectura de ingeniería.

### 3.3 Evidencia preliminar: los barridos ya dan resultados

Estos barridos ya los corrí en el prototipo. Confirman que **hay algo que
encontrar**, que es lo que había que verificar antes de invertir en el módulo.

**Masa del eslabón superior $m_3$** (con $m_1 = m_2 = 0.2$ kg, $\ell_j = 1/3$ m):

| $m_3$ [kg] | $\lambda_{\max}$ | $\hat\mu$ | $\|K\|$ | $J^\star(\theta_1=0.1)$ | $c^\star(50\text{ N})$ |
|---|---|---|---|---|---|
| 0.02 | 13.65 | $1.13\times10^{-3}$ | 614 | 0.855 | 4.59 |
| 0.20 | 18.06 | $4.92\times10^{-4}$ | 1176 | 1.620 | 3.73 |
| 0.80 | 25.00 | $1.88\times10^{-4}$ | 2066 | 5.156 | 2.88 |
| 3.20 | 41.83 | $4.16\times10^{-5}$ | 3621 | 24.74 | 1.96 |

Degradación monótona: **cabeza pesada, sistema difícil**. Sin frontera abrupta.

**Longitud del eslabón superior $\ell_3$:**

| $\ell_3$ [m] | $\lambda_{\max}$ | $\hat\mu$ | $\|K\|$ | $c^\star(50\text{ N})$ |
|---|---|---|---|---|
| 0.02 | 42.93 | $3.39\times10^{-5}$ | 1728 | 1.95 |
| 0.10 | 22.00 | $1.95\times10^{-4}$ | 1622 | 3.01 |
| 0.333 | 18.06 | $4.92\times10^{-4}$ | 1176 | 3.73 |
| 1.00 | 17.39 | $5.85\times10^{-4}$ | 858 | 4.22 |
| 2.00 | 17.26 | $5.99\times10^{-4}$ | 741 | 4.45 |

Resultado con buena lectura física: **el eslabón superior corto es el peligroso**.
Un eslabón corto y masivo cae rapidísimo ($\lambda_{\max} = 42.9$ para
$\ell_3 = 2$ cm) y el margen PBH se hunde. Alargarlo *facilita* el control, que es
el conocido efecto de "la escoba larga es más fácil de equilibrar en la palma"
—ahora cuantificado, y ahora en la punta de un triple.

**Masa del carro $M$ (el hallazgo interesante):**

| $M$ [kg] | $\lambda_{\max}$ | $\hat\mu$ | $\|K\|$ | $c^\star(50\text{ N})$ |
|---|---|---|---|---|
| 0.05 | 21.41 | $1.412\times10^{-3}$ | 256 | 2.88 |
| **0.10** | 20.15 | $\mathbf{1.643\times10^{-3}}$ | 311 | 3.07 |
| 0.25 | 18.90 | $1.340\times10^{-3}$ | 467 | 3.36 |
| 1.00 | 18.06 | $4.916\times10^{-4}$ | 1176 | 3.73 |
| 5.00 | 17.81 | $1.068\times10^{-4}$ | 4778 | 3.93 |
| 20.00 | 17.77 | $2.703\times10^{-5}$ | 18043 | 4.00 |

Aquí **el margen PBH no es monótono: tiene un máximo interior alrededor de
$M \approx 0.1$ kg**. Es exactamente el tipo de resultado que se pedía: una
*configuración óptima* que el barrido encuentra y que no es un extremo del
dominio. La interpretación es un compromiso: un carro muy pesado se mueve poco
ante la misma fuerza (autoridad de control baja, $\|K\|$ crece a 18043); uno muy
liviano se mueve mucho pero el acoplamiento con los eslabones se degrada. En el
medio hay un óptimo.

Nótese además que $c^\star$ y $\hat\mu$ **apuntan en direcciones opuestas** en
este barrido: $c^\star$ crece monótonamente con $M$ mientras $\hat\mu$ tiene un
máximo. No hay una única "mejor" configuración; depende de qué se optimice. Eso
también es un resultado, y da material para la discusión del informe.

### 3.4 Diseño de los barridos

**Barridos 1D** (tabla + curva). Un parámetro contra todas las métricas:
$m_3$, $\ell_3$, $M$, $b$, $g$, $R$, $u_{\max}$, $h$.

**Barridos 2D** (mapa de calor + curva de frontera). Los cuatro que más aportan:

| Mapa | Ejes | Métrica | Qué muestra |
|---|---|---|---|
| A | $(\ell_3/\ell_1,\ m_3/m_1)$ | $\hat\mu$ | Región estructuralmente sana |
| B | $(u_{\max},\ \theta_0)$ | éxito/fracaso | **La frontera de operabilidad práctica** |
| C | $(M,\ m_{\text{total}})$ | $\hat\mu$ y $\|K\|$ | El óptimo interior de 3.3 |
| D | $(R,\ \theta_0)$ | $\theta_{\max}$ | Compromiso agresividad-robustez |

El mapa **B** es el entregable estrella: un plano $(u_{\max}, \theta_0)$ pintado
en dos colores con una curva separando "se recupera" de "se cae", y sobre él,
superpuesta, la cota elipsoidal $\theta_0 = \sqrt{c^\star/P_{33}}$ del punto (e).
Ver la cota analítica quedar **por dentro** de la frontera numérica —conservadora
pero válida— es la mejor diapositiva posible para cerrar la presentación.

**Costo computacional.** Un mapa de $40\times40$ con bisección de 10 pasos son
16000 integraciones. Estimando 5 ms cada una, unos 80 s por mapa con un hilo.
Recomendaciones: `Threads.@threads` sobre las filas, `saveat=0.01`, `save_everystep=false`
cuando solo interesa el estado final, y un `callback` de terminación temprana que
aborte en cuanto $\max_j|\theta_j| > \pi/2$ (la mayoría de los fracasos se
detectan en las primeras décimas de segundo). Con eso baja a menos de 20 s por
mapa. **Cachear los resultados** en un `.jld2` o un `.csv` para no recomputar al
regenerar figuras.

### 3.5 API propuesta

```
src/metrics.jl     modulo Metrics    Metricas de viabilidad (genericas en n)
  pbh_controllability_margin(A, B)
  pbh_observability_margin(A, C)
  controllability_condition(A, B)
  optimal_cost(P, x0)
  nonsaturation_level(K, P, R, B, umax)      -> c*
  max_axis_angle(P, c, idx)                  -> sqrt(c / P[idx,idx])
  discretize_zoh(A, B, h)                    -> (Ad, Bd)  via Van Loan
  is_schur(M)

src/sweep.jl       modulo Sweep      Barridos y fronteras
  viability(build, params; Q, R, umax, x0)   -> NamedTuple con todas las metricas
  sweep_1d(build, base, sym, values; ...)    -> tabla
  sweep_2d(build, base, sym1, vals1, sym2, vals2; ...)
  max_recoverable_angle(build, params, K; umax, ...)
  min_actuator_force(build, params, K; theta0, ...)
  max_sampling_period(A, B, K)
  print_sweep_summary(result)

docs/Entrega_2/atlas/make_atlas_figs.jl      Genera los mapas A-D
```

`Metrics` y `Sweep` deben ser **genéricos en $n$**, igual que `Controller`: así
los mismos barridos corren sobre las Configuraciones I y II y producen la tabla
comparativa de tres columnas que es el corazón del nuevo informe.

---

## 4. Catálogo de módulos para robustecer

Cada módulo es independiente. La columna "valor" es mi lectura de cuánto suma al
trabajo frente a lo que cuesta.

| ID | Módulo | Esfuerzo | Valor | Depende de |
|---|---|---|---|---|
| **M1** | Modelo genérico de $N$ eslabones (`ModelNLink`) | Medio | Alto | — |
| **M2** | Configuración III completa (triple + pipeline + animación) | Medio | Alto | M1 |
| **M3** | Métricas de viabilidad (`Metrics`) | Bajo | **Muy alto** | — |
| **M4** | Atlas de operabilidad (`Sweep` + mapas) | Alto | **Muy alto** | M3 |
| **M5** | Observador de Luenberger y LQG | Medio | **Muy alto** | — |
| **M6** | Región de atracción por Lyapunov | Medio | Alto | M3 |
| **M7** | Pruebas de regresión (`test/runtests.jl`) | Bajo | Alto | M1 |
| **M8** | Verificación simbólica de los jacobianos | Bajo | Medio | — |
| **M9** | Control digital y período de muestreo | Bajo | Medio | M3 |
| **M10** | Robustez paramétrica | Medio | Medio | M3 |
| **M11** | LQR con acción integral | Bajo | Medio | — |
| **M12** | Swing-up energético | Alto | Medio | M2 |

A continuación, los que necesitan explicación (M1-M4 ya están desarrollados en
las secciones 2 y 3).

### M5. Observador de Luenberger y LQG

**La grieta que cierra.** Es la más importante del catálogo. El proyecto
*demuestra* que $\operatorname{rank}\mathcal{O} = n$ y luego controla con estado
completo, que en la práctica no se mide. Un observador cierra el argumento.

**La matemática.** Se estima el estado con

$$\dot{\hat{\mathbf{x}}} = A\hat{\mathbf{x}} + Bu + L\left(\mathbf{y} - C\hat{\mathbf{x}}\right)$$

y el error $\mathbf{e} = \mathbf{x}-\hat{\mathbf{x}}$ obedece
$\dot{\mathbf{e}} = (A - LC)\mathbf{e}$. Diseñar $L$ para que $A-LC$ sea Hurwitz
es **el problema dual** de diseñar $K$: por la dualidad de Kalman,

$$L = \operatorname{lqr}(A^\top,\ C^\top,\ Q_o,\ R_o)^\top$$

de modo que **se reutiliza `design_lqr` sin escribir un solo algoritmo nuevo**.
Ese es el punto elegante: la dualidad convierte un problema nuevo en el que ya
está resuelto.

El sistema aumentado en lazo cerrado, en coordenadas $(\mathbf{x},\mathbf{e})$,
tiene matriz **triangular por bloques**:

$$\begin{pmatrix}\dot{\mathbf{x}}\\ \dot{\mathbf{e}}\end{pmatrix} =
\begin{pmatrix}A - BK & BK\\ 0 & A - LC\end{pmatrix}
\begin{pmatrix}\mathbf{x}\\ \mathbf{e}\end{pmatrix}$$

y por tanto $\sigma = \sigma(A-BK)\cup\sigma(A-LC)$: el **principio de
separación**, que sale de una observación de álgebra lineal de dos líneas (el
espectro de una triangular por bloques es la unión de los espectros diagonales).
Es un resultado precioso y encaja perfecto en el marco del curso.

**El barrido asociado.** ¿Qué subconjuntos de sensores mantienen la
observabilidad? Para el triple hay $2^4-1 = 15$ combinaciones de
$\{x,\theta_1,\theta_2,\theta_3\}$. Calcular para cada una el rango y el **margen
PBH de observabilidad** $\min_\lambda\sigma_{\min}\!\begin{pmatrix}A-\lambda I\\ C\end{pmatrix}$
produce una tabla que responde: *¿cuál es el juego mínimo de sensores?* Es una
pregunta de ingeniería con respuesta de álgebra lineal, y es exactamente el tipo
de resultado que luce en una defensa.

**Extensión natural (LQG).** Si se interpretan $Q_o$ y $R_o$ como covarianzas de
ruido de proceso y de medición, $L$ es la ganancia del **filtro de Kalman** y el
conjunto $K+L$ es un controlador LQG. Cuesta un párrafo adicional y sube el nivel
del informe.

### M6. Región de atracción por Lyapunov

**Qué agrega.** Convierte "el LQR aguanta hasta 0.40 rad, según lo que probamos"
en "el LQR **garantiza** convergencia desde cualquier $\mathbf{x}_0$ con
$\mathbf{x}_0^\top P\mathbf{x}_0 \le c$", que es una afirmación matemática, no
empírica.

**La matemática.** Tres capas, de más fácil a más difícil:

1. **Capa de no saturación** (exacta, ya desarrollada en 3.2(e)):
   $c^\star = u_{\max}^2/(KP^{-1}K^\top)$. Cuesta una línea de código.
2. **Capa no lineal.** Escribiendo $\dot{\mathbf{x}} = (A-BK)\mathbf{x} + \mathbf{r}(\mathbf{x})$
   con $\|\mathbf{r}(\mathbf{x})\|\le \gamma\|\mathbf{x}\|^2$ en una bola, se
   obtiene $\dot V \le -\lambda_{\min}(Q+K^\top RK)\|\mathbf{x}\|^2 + 2\gamma\|P\|\,\|\mathbf{x}\|^3$,
   que es negativa mientras $\|\mathbf{x}\| < \lambda_{\min}(Q+K^\top RK)/(2\gamma\|P\|)$.
   Estimando $\gamma$ numéricamente (muestreando el residuo en una malla) se
   obtiene un $c$ menor pero **riguroso también para el no lineal**.
3. **Capa numérica.** La bisección de 3.2(f) da la frontera real.

Presentar las tres juntas —cota analítica conservadora, cota analítica ajustada,
frontera numérica— en una sola figura es un resultado de mucha calidad: se ve
cómo la teoría acota y dónde la simulación aporta.

**Alternativa si el tiempo aprieta:** hacer solo la capa 1 y la 3, y decir
explícitamente que la 2 queda propuesta. Sigue siendo bueno.

### M7. Pruebas de regresión

**Por qué importa.** Los jacobianos analíticos de `linearization.jl` están
escritos a mano. Verifiqué que son correctos, pero nada en el repositorio lo
comprueba de forma reproducible, y con la Configuración III se van a escribir
más. Un `test/runtests.jl` convierte esa verificación en algo permanente.

**Qué probar** (todo con `@test` de la stdlib `Test`, sin dependencias nuevas
salvo `ForwardDiff` si se acepta):

1. **Jacobiano analítico vs. numérico**: para cada configuración,
   $\|A_{\text{analítico}} - \partial f/\partial\mathbf{x}\| < 10^{-8}$ usando
   diferenciación automática (`ForwardDiff.jacobian`) sobre `nonlinear_eom_x!`
   evaluada en el equilibrio. Es la prueba más valiosa del conjunto.
2. **`ModelNLink` reproduce `Model` y `ModelDouble`**: dos implementaciones
   independientes de la misma física deben dar la misma $\dot{\mathbf x}$ para
   estados aleatorios. Esta es la prueba que justifica el diseño aditivo de 2.9.
3. **Valores del informe**: los espectros $\{+4.21,\dots\}$, $\{+8.57,\dots\}$ y
   las ganancias $K$ publicadas, con tolerancia de $10^{-3}$. Blinda los cuatro
   documentos LaTeX contra cambios accidentales.
4. **`solve_care` vs. `MatrixEquations.arec`**: la implementación propia contra la
   biblioteca. `MatrixEquations` ya está en `Project.toml` y no se está usando.
5. **Propiedades algebraicas**: $P = P^\top$, $P \succ 0$, el residuo de la CARE
   $\|A^\top P + PA - PBR^{-1}B^\top P + Q\| < 10^{-8}$, y Cayley-Hamilton
   ($\|\varphi_A(A)\| \approx 0$) como guiño al contenido del curso.

**Costo:** bajo, y desproporcionadamente rentable.

### M8. Verificación simbólica de los jacobianos

> **Nota posterior.** M8 no se implementó, y `Symbolics` acabó retirándose del
> `Project.toml` precisamente por no usarse. Si se retoma este módulo hay que
> volver a añadirlo. El resto de la sección se conserva tal como se planeó.

`Symbolics` estaba en `Project.toml` y **no se usaba en ningún archivo del
repositorio**. Un script que derive simbólicamente $\partial f/\partial\mathbf{x}$
a partir de las EOM y lo compare con las expresiones escritas a mano cerraría ese
cabo suelto, y produciría de paso las fórmulas de $A$ y $B$ de la Configuración
III en forma cerrada para el apéndice del informe. Es una tarde de trabajo con
buen retorno documental.

### M9. Control digital y período de muestreo

Ya desarrollado en 3.2(f). Se implementa entero con `exp` de matrices:
$A_d = e^{Ah}$ y $B_d$ por el truco de Van Loan, criterio de Schur en vez de
Hurwitz, bisección sobre $h$. Conecta con la sección 6.1 de la guía maestra
(la exponencial matricial como operador solución) y da un resultado de
ingeniería inmediato: **el triple exige muestrear por encima de 31 Hz**.

Extensión opcional: rediseñar el LQR **directamente en discreto** (DARE en lugar
de CARE) y comparar. La DARE se resuelve con la matriz **simpléctica** en vez de
la hamiltoniana, lo que da un paralelo estructural muy bonito con lo ya hecho.

### M10. Robustez paramétrica

Diseñar $K$ con los parámetros nominales y evaluar el lazo cerrado sobre una
planta **perturbada** ($m_j$, $\ell_j$, $M$ desviados un $\pm p\%$). Métricas:

- Máximo $p$ para el que $A_{\text{real}} - B_{\text{real}}K$ sigue siendo Hurwitz.
- Barrido 2D: error en $m_3$ contra error en $\ell_3$, coloreando estabilidad.

Responde "¿qué tan bien hay que conocer el sistema?", que es la objeción práctica
más común a un diseño basado en modelo. Se puede hacer barato reutilizando
`Sweep` de M4.

### M11. LQR con acción integral

Aumentar el estado con $\xi = \int (x - x_{\text{ref}})\,dt$:

$$\tilde A = \begin{pmatrix}A & \mathbf{0}\\ \mathbf{c}_x^\top & 0\end{pmatrix},\qquad
\tilde B = \begin{pmatrix}B\\ 0\end{pmatrix}$$

y diseñar LQR sobre $(\tilde A,\tilde B)$. Permite **seguir consignas de posición
del carro** en vez de solo regular a cero: mover el carro un metro manteniendo
los tres eslabones erguidos es una demo espectacular y cuesta poco código.
Requiere verificar que el par aumentado siga siendo controlable, lo cual es una
condición de rango sobre una matriz orlada: más álgebra lineal.

### M12. Swing-up energético

Levantar el péndulo desde abajo con control basado en energía y conmutar al LQR
al entrar en la región de atracción (que M6 cuantifica). Espectacular en la
animación, pero **se sale del marco de álgebra lineal** y para el triple es un
problema de investigación, no de curso. Sugerencia: hacerlo **solo para la
Configuración I** si se quiere el efecto visual, y presentarlo como apéndice o
como cierre de la presentación, dejando claro que es no lineal.

---

## 5. Rutas recomendadas

### Ruta mínima — "cerrar las grietas" (aprox. 1 semana)

`M3` + `M5` + `M7`

No agrega configuración nueva. Cierra las tres grietas más señalables: mide la
controlabilidad de forma continua, construye el observador que faltaba y blinda
los resultados con pruebas. Es la ruta de **mayor valor por unidad de esfuerzo**
si el objetivo es que el trabajo aguante una defensa exigente.

### Ruta media — "la tercera configuración" (aprox. 2 semanas) — **recomendada**

`M1` + `M2` + `M3` + `M4` + `M7`

Agrega el péndulo triple, las métricas y el atlas de operabilidad completo, con
pruebas de regresión. Es la ruta que cumple literalmente lo pedido: nueva
configuración más barridos que identifican óptimos y fronteras. Produce una
segunda entrega autónoma, con resultados nuevos y figuras nuevas.

### Ruta ambiciosa — "el trabajo completo" (aprox. 3-4 semanas)

Ruta media + `M5` + `M6` + `M9`

Añade observador, región de atracción rigurosa y control digital. El resultado es
un trabajo que cubre modelado, análisis, control por realimentación de estado,
estimación de estado, garantías de dominio de validez e implementación digital:
esencialmente un curso completo de control lineal, sostenido sobre álgebra
lineal. Si hay tiempo y ganas, es lo que yo haría.

### Lo que yo priorizaría si solo pudiera escoger dos cosas

`M3` (métricas) y `M5` (observador). El primero porque convierte una afirmación
binaria y frágil en una medida continua y defendible, y porque habilita todo lo
demás. El segundo porque es la grieta más visible del trabajo actual y porque su
solución —la dualidad— es el argumento de álgebra lineal más bonito que queda sin
usar en el proyecto.

---

## 6. Prompts listos para usar

### 6.0 Preámbulo común

Anteponer este bloque a cualquiera de los prompts siguientes:

```text
Contexto: proyecto-pendulos, trabajo de Algebra Lineal Aplicada (UNAL Medellin).
Lee primero docs/guia_maestra.md (secciones 3 y 8) y docs/Entrega_2/plan_de_trabajo.md
para entender la arquitectura y las convenciones.

Convenciones obligatorias:
- Responder en espanol.
- Comentarios de codigo en espanol SIN TILDES NI ENIES, coherentes con
  src/model_double.jl. Sin emojis ni simbolos decorativos.
- Nombres de variables y funciones en ingles.
- Toda notacion matematica en LaTeX.
- Angulos medidos desde la vertical superior; theta = 0 es el equilibrio erguido.
- Orden del estado INTERCALADO: [pos, vel, theta1, omega1, theta2, omega2, ...].
- Encabezado de archivo con banner "# ====" describiendo proposito y variables
  de estado, igual que los modulos existentes.
- No duplicar fisica ni algoritmos: reutilizar los modulos de src/.
- No tocar src/model_simple.jl ni src/model_double.jl.
- Commits en espanol con formato "tipo: descripcion".

Antes de escribir codigo, explicame el plan y espera confirmacion.
```

### 6.1 Prompt M1 — Modelo genérico de $N$ eslabones

```text
Implementa src/model_nlink.jl, modulo ModelNLink: modelo generico de N eslabones
articulados en serie sobre un carro, con barras de inercia arbitraria.

Parametros (Base.@kwdef struct SystemParamsNLink):
  M::Float64            masa del carro [kg]
  m::Vector{Float64}    masas de los eslabones [kg]
  l::Vector{Float64}    longitudes de los eslabones [m]
  a::Vector{Float64}    distancia de la articulacion al centro de masa [m]
  Il::Vector{Float64}   momentos de inercia respecto al CM [kg m^2]
  g::Float64 = 9.81
  b::Float64 = 0.0      friccion viscosa del carro [N s/m]
  c::Vector{Float64}    friccion en las articulaciones (relativa), por defecto ceros

Constructores de conveniencia:
  uniform_rods(M, m, l; g, b)   barras uniformes: a = l/2, Il = (1/12) m l^2
  point_masses(M, m, l; g, b)   masas puntuales:  a = l,   Il = 0

Fisica (deducida por Euler-Lagrange, ver docs/Entrega_2/plan_de_trabajo.md seccion 2):
  beta_j  = m_j a_j + l_j * sum_{i>j} m_i
  J_j     = Il_j + m_j a_j^2 + l_j^2 * sum_{i>j} m_i
  Gam_jk  = l_min(j,k) * beta_max(j,k)

Matriz de masa M(q) con q = [pos, theta_1, ..., theta_N]:
  M[1,1]     = M + sum(m)
  M[1,j+1]   = beta_j cos(theta_j)
  M[j+1,j+1] = J_j
  M[j+1,k+1] = Gam_jk cos(theta_j - theta_k)   para j != k

Lado derecho:
  fila del carro:     F - b*vel + sum_j beta_j sin(theta_j) omega_j^2
  fila del eslabon j: g beta_j sin(theta_j)
                      - sum_{k != j} Gam_jk sin(theta_j - theta_k) omega_k^2
                      - par de friccion articular

Resolver con cholesky(Symmetric(Mq)) \ rhs, igual que model_double.jl.

Exportar: SystemParamsNLink, uniform_rods, point_masses, nonlinear_eom_nlink!,
closed_loop_eom_nlink!, state_derivative_nlink.

Criterio de aceptacion (verificalo ejecutando, no lo asumas):
1. Con N=1, uniform_rods(1.0, [0.3], [1.0]; b=0.1), la derivada del estado debe
   coincidir con Model.state_derivative en al menos 20 estados aleatorios,
   tolerancia 1e-10.
2. Con N=2, point_masses(1.0, [0.3,0.3], [0.5,0.5]; b=0.0), idem contra
   ModelDouble.state_derivative_double.
Reporta el error maximo obtenido en cada caso.
```

### 6.2 Prompt M2 — Configuración III completa

```text
Implementa la Configuracion III (pendulo triple sobre carro, estado en R^8)
siguiendo docs/Entrega_2/plan_de_trabajo.md seccion 2.

1. src/model_triple.jl (modulo ModelTriple): capa delgada sobre ModelNLink que
   expone la API del proyecto: SystemParamsTriple, default_params_triple(),
   nonlinear_eom_triple!, closed_loop_eom_triple!, state_derivative_triple.
   Defaults: M=1.0, m1=m2=m3=0.2 kg, l1=l2=l3=1/3 m, barras uniformes,
   g=9.81, b=0.1, sin friccion articular.

2. En src/linearization.jl agrega linearize_system_triple(params) que devuelva
   un StateSpaceModel de dimension 8, con jacobiano ANALITICO (no numerico):
     M0 constante = M(q) evaluada en theta = 0
     G0 = g * diag(0, beta_1, beta_2, beta_3)
     F0 = matriz de friccion
     Acc = M0 \ G0,  Acv = -(M0 \ F0),  Bc = M0 \ e1
   y ensamblar A y B en orden INTERCALADO.
   C mide pos, theta1, theta2, theta3 (4x8). D = zeros(4,1).
   Mantener el estilo de comentarios de linearize_system_double.

3. src/animation_triple.jl (modulo AnimationTriple, GLMakie), generalizando
   link_positions de animation_double.jl a tres eslabones.

4. main_triple.jl: pipeline igual a main_double.jl (parametros, simulacion libre,
   linealizacion y analisis, LQR, lazo cerrado, graficas, animacion).
   Q = diagm([1,0,10,0,10,0,10,0]), R = 0.1, saturacion 150 N (NO 50 N: con
   norm(K) ~ 1176 el actuador de 50 N es insuficiente).
   Figuras: figures/09_triple_comparativa_lqr.png y
   figures/10_triple_animacion_lqr.mp4.

Criterio de aceptacion (verificalo ejecutando y reporta los valores obtenidos):
- Espectro de lazo abierto:
  {+18.0613, +9.7740, +4.3777, 0, -0.0625, -4.3990, -9.7813, -18.0647}
  (tolerancia 1e-3). Tres modos inestables.
- rank(Ctrb) = 8/8 y rank(Obsv) = 8/8.
- B = [0, 0.9455, 0, -3.6, 0, 0.9818, 0, -0.3273]' (tolerancia 1e-3).
- K del LQR = (-3.16, -6.17, -317.43, -4.37, 911.95, 37.73, -667.37, -61.42),
  con norm(K) = 1176.0 (tolerancia 1e-2 relativa).
- Polos de lazo cerrado: -0.84+-0.79i, -4.57+-2.32i, -10.02+-1.73i, -18.14+-0.97i.
- Verificacion cruzada: A debe coincidir con la jacobiana numerica de
  ModelTriple.nonlinear_eom_triple! en el equilibrio (usa diferencias finitas
  centradas si no quieres agregar ForwardDiff).

Si algun valor no coincide, NO ajustes el criterio: hay un error en el modelo.
```

### 6.3 Prompt M3 — Métricas de viabilidad

```text
Implementa src/metrics.jl, modulo Metrics: metricas para cuantificar que tan
controlable es un sistema, mas alla del criterio binario de rango. Todo generico
en la dimension del estado, igual que Controller.

Funciones (con docstring que explique la matematica, en el estilo de
src/controller.jl):

  pbh_controllability_margin(A, B)
      min sobre lambda en spec(A) de sigma_min([A - lambda I  B]).
      Version cuantitativa del test de Popov-Belevitch-Hautus. Por el teorema
      de Eckart-Young es la distancia en norma 2 a perder el rango.

  pbh_controllability_margin_normalized(A, B)
      lo anterior dividido por opnorm([A B]), para comparar entre sistemas.

  pbh_observability_margin(A, C) y su version normalizada
      analogo con [A - lambda I; C].

  controllability_condition(A, B)
      sigma_max / sigma_min de la matriz de Kalman. Documentar en el docstring
      que crece ~1e4 por eslabon y que para N=5 alcanza 1/eps: por eso el rango
      deja de ser confiable.

  optimal_cost(P, x0) = x0' P x0

  nonsaturation_level(K, P, R, B, umax)
      c* = umax^2 / (K P^-1 K'). Como K = R^-1 B' P, equivale a
      umax^2 R^2 / (B' P B) para R escalar; usar esta forma para no invertir P.
      Es el mayor conjunto de nivel de V = x'Px contenido en |K x| <= umax.

  max_axis_angle(P, c, idx) = sqrt(c / P[idx, idx])
      maximo valor de la componente idx sobre el eje, dentro del elipsoide.

  discretize_zoh(A, B, h)
      Ad = exp(A h), Bd por el truco de Van Loan: exponencial de la matriz
      aumentada [A B; 0 0] de tamano (n+1)x(n+1).

  is_schur(M; tol=1.0) -> maximum(abs.(eigvals(M))) < tol

  print_metrics_summary(A, B, C, K, P, R; umax) que imprima todo con el mismo
  formato de banner y @printf de print_analysis.

Criterio de aceptacion (ejecutalo sobre las tres configuraciones y reporta la
tabla):
- Margen PBH normalizado: simple ~1.72e-2, doble ~2.56e-3, triple ~4.92e-4.
- cond(Ctrb): simple ~3.6e1, doble ~1.6e4, triple ~2.2e8.
- c* con umax = 50 N: simple ~22.94, doble ~8.99, triple ~3.73.
- Angulo garantizado sin saturar sqrt(c*/P[3,3]): simple ~0.962 rad,
  doble ~0.234 rad, triple ~0.152 rad.
```

### 6.4 Prompt M4 — Atlas de operabilidad

```text
Implementa src/sweep.jl (modulo Sweep) y
docs/Entrega_2/atlas/make_atlas_figs.jl, siguiendo docs/Entrega_2/plan_de_trabajo.md
seccion 3.

A) Modulo Sweep, generico en n:

  viability(build, params; Q, R, umax, x0) -> NamedTuple con:
      lambda_max, tau_fall, rank_ctrb, cond_ctrb, pbh_margin, pbh_normalized,
      K, norm_K, J_star, c_star, theta_guaranteed, poles_cl

  sweep_1d(build, base_params, sym, values; kwargs...)
      barre un parametro y devuelve un vector de NamedTuple.

  sweep_2d(build, base_params, sym1, vals1, sym2, vals2; kwargs...)
      devuelve una matriz. Paraleliza con Threads.@threads sobre las filas.

  max_recoverable_angle(build, params, K; umax, x_rail=1.5, T=10.0, iters=12)
      biseccion sobre theta0. Criterio de exito, sobre el modelo NO LINEAL con
      saturacion:
        max_j |theta_j(T)| < 0.02 rad   Y
        max_t |x(t)| <= x_rail          Y
        max_{t,j} |theta_j(t)| < pi/2
      Reutiliza el umbral 0.02 rad de settling_time en make_report_figs.jl.
      Usa un callback de terminacion temprana cuando algun |theta_j| > pi/2.

  min_actuator_force(build, params, K; theta0, ...)  biseccion sobre umax.
  max_sampling_period(A, B, K; ...)                  biseccion sobre h, criterio
      de Schur sobre Ad - Bd K usando Metrics.discretize_zoh.

  print_sweep_summary(result) con el formato de banner del proyecto.

B) Script make_atlas_figs.jl (CairoMakie, paleta de docs/Entrega_1/presentacion/
make_slide_figs.jl: MARINO, PETROLEO, GRISVERDE, SALVIA, LIMON). Cuatro mapas:

  Mapa A: (l3/l1, m3/m1) -> margen PBH normalizado, escala log.
  Mapa B: (umax, theta0) -> exito/fracaso, con la curva de frontera Y
          superpuesta la cota elipsoidal theta0 = sqrt(c*/P[3,3]).
          Este es el mapa principal.
  Mapa C: (M, masa total del pendulo) -> margen PBH y norm(K), dos paneles.
  Mapa D: (R, theta0) -> theta_max.

  Guardar en docs/Entrega_2/atlas/figs/ como atlas_a.png ... atlas_d.png.
  Cachear los resultados numericos en un CSV para no recomputar.
  Imprimir por consola las tablas 1D que van al informe.

Resultados que el barrido DEBE reproducir (ya calculados, sirven de control):
- Escalamiento con N (longitud total 1 m, masa total 0.6 kg, barras uniformes):
    N=1: lambda_max=4.51,  pbh_norm=1.31e-2, cond=4.8e1,  norm(K)=52.9
    N=2: lambda_max=10.59, pbh_norm=1.79e-3, cond=6.2e4,  norm(K)=248.7
    N=3: lambda_max=18.06, pbh_norm=4.92e-4, cond=2.2e8,  norm(K)=1176.0
    N=4: lambda_max=26.54, pbh_norm=1.90e-4, cond=1.7e12, norm(K)=5266.8
    N=5: lambda_max=35.64, pbh_norm=9.03e-5, cond=3.3e16, norm(K)=23221.9
- Barrido de la masa del carro M en el triple: el margen PBH normalizado NO es
  monotono, tiene un maximo interior cerca de M = 0.1 kg (valor ~1.64e-3).
  Confirma este optimo con resolucion fina en M entre 0.05 y 0.5 kg.
- Periodo de muestreo maximo con los parametros POR DEFECTO de cada
  configuracion: simple ~197 ms, doble ~77 ms, triple ~32 ms.

Comenta en la salida cual es la configuracion optima encontrada y donde estan
las fronteras de imposibilidad, distinguiendo imposibilidad ESTRUCTURAL
(margen PBH colapsado) de imposibilidad PRACTICA (saturacion, riel, muestreo).
```

### 6.5 Prompt M5 — Observador de Luenberger y LQG

```text
Implementa src/observer.jl, modulo Observer: estimacion de estado por observador
de Luenberger, explotando la dualidad de Kalman para reutilizar Controller.

  design_observer(A, C, Qo, Ro)
      L = design_lqr(A', C', Qo, Ro).K'   (dualidad: no reimplementar Riccati)
      Devuelve (L=L, eigenvalues_obs=eigvals(A - L*C)).

  design_observer_poles(A, C, desired_poles)
      dual de design_pole_placement sobre (A', C').

  augmented_closed_loop(A, B, C, K, L)
      Matriz [A-BK  BK; 0  A-LC] en coordenadas (x, e). Documentar en el
      docstring que su triangularidad por bloques demuestra el PRINCIPIO DE
      SEPARACION: spec = spec(A-BK) union spec(A-LC).

  observer_eom!(dx, x, p, t)
      Dinamica conjunta planta no lineal + observador lineal. El estado tiene
      2n componentes: las n primeras son el estado real (integrado con la EOM
      no lineal del modelo correspondiente, pasado en p.eom!), las n siguientes
      el estimado. El control es u = -K xhat (NO -K x), con saturacion opcional.
      La medicion es y = C x (opcionalmente con ruido gaussiano p.sigma).

  sensor_subset_analysis(A, C_full, labels)
      Para cada uno de los 2^p - 1 subconjuntos no vacios de filas de C_full,
      reporta rank(Obsv) y el margen PBH de observabilidad
      (usa Metrics.pbh_observability_margin). Devuelve una tabla ordenada por
      margen decreciente. Para el triple, C_full mide pos, theta1, theta2,
      theta3: son 15 combinaciones.

Ademas, en main_triple.jl (o en un script aparte) agrega una demostracion:
simula el triple en lazo cerrado con observador, arrancando el estimador con
error inicial (xhat0 = 0 mientras x0 tiene theta1 = 0.1), y grafica el error de
estimacion e(t) junto con la respuesta. Guarda como
figures/11_triple_observador.png.

Criterio de aceptacion:
- Verificar numericamente el principio de separacion: los eigenvalores de la
  matriz aumentada deben ser la union de spec(A-BK) y spec(A-LC), tolerancia 1e-8.
- Reportar la tabla de subconjuntos de sensores e indicar cual es el juego
  MINIMO de sensores que mantiene la observabilidad con margen aceptable.
```

### 6.6 Prompt M6 — Región de atracción

```text
Implementa src/roa.jl, modulo ROA: estimacion de la region de atraccion del
lazo cerrado, en tres capas de rigor decreciente en garantia y creciente en
ajuste. Ver docs/Entrega_2/plan_de_trabajo.md seccion M6.

  ellipsoid_nonsaturation(K, P, R, B, umax)
      Capa 1 (exacta para el sistema lineal): c* = umax^2 / (K P^-1 K').
      Reutiliza Metrics.nonsaturation_level.

  ellipsoid_nonlinear(A, B, K, P, Q, R, eom!, params; radius, samples)
      Capa 2: estima gamma tal que ||r(x)|| <= gamma ||x||^2 muestreando el
      residuo r(x) = f_nolineal(x) - (A - B K) x sobre una malla o muestreo
      aleatorio en una bola de radio dado. Devuelve el mayor c tal que
      Vdot < 0 en el elipsoide, usando
        Vdot <= -lambda_min(Q + K'RK) ||x||^2 + 2 gamma ||P|| ||x||^3.
      Documentar claramente que gamma es una ESTIMACION muestral, no una cota
      demostrada, y por tanto el resultado es una garantia condicionada.

  roa_bisection(build, params, K; umax, direction, ...)
      Capa 3: frontera numerica real por biseccion sobre el modelo no lineal.
      Reutiliza Sweep.max_recoverable_angle.

  compare_roa_estimates(...) que devuelva las tres cotas juntas.

Figura: en el plano (theta1, omega1) con el resto del estado en cero, dibujar
las tres regiones superpuestas (elipse de capa 1, elipse de capa 2, frontera
numerica de capa 3) mas algunas trayectorias del modelo no lineal que converjan
y otras que no. Guardar en docs/Entrega_2/atlas/figs/roa_triple.png.

Hazlo para las tres configuraciones y reporta la tabla comparativa. Con
umax = 50 N los valores de la capa 1 sobre el eje theta1 ya estan calculados:
simple 0.962 rad, doble 0.234 rad, triple 0.152 rad. Deben reproducirse.
```

### 6.7 Prompt M7 — Pruebas de regresión

```text
Crea test/runtests.jl usando la stdlib Test, con @testset anidados. Agrega Test
(y ForwardDiff si lo usas) a Project.toml en la seccion [extras]/[targets] o
como dependencia, lo que sea mas limpio.

@testset "Modelos":
  - ModelNLink con N=1 (barra uniforme) reproduce Model.state_derivative en
    20 estados aleatorios, tolerancia 1e-10.
  - ModelNLink con N=2 (masas puntuales) reproduce
    ModelDouble.state_derivative_double, idem.
  - La matriz de masa es simetrica definida positiva en 50 configuraciones
    aleatorias de angulos, para las tres configuraciones.

@testset "Linealizacion":
  - Para cada configuracion, la A analitica coincide con la jacobiana numerica
    de la EOM no lineal en el equilibrio (ForwardDiff.jacobian o diferencias
    finitas centradas), tolerancia 1e-6.
  - Idem para B respecto a la entrada.
  - Estructura del espectro: sin friccion, los eigenvalores vienen en pares
    +-lambda mas un cero doble (ver seccion 2.6 del documento).

@testset "Valores publicados":
  - Espectro simple {+4.2105, 0, -0.0769, -4.2266}, tolerancia 1e-3.
  - Espectro doble {+8.5726, +4.0941, 0, 0, -4.0941, -8.5726}.
  - Espectro triple {+18.0613, +9.7740, +4.3777, 0, -0.0625, -4.3990, -9.7813,
    -18.0647}.
  - K simple (-3.16, -4.69, -45.39, -10.93), tolerancia 1e-2.
  - K doble (3.16, 5.82, -191.55, -10.99, 228.32, 36.14).
  - K triple (-3.16, -6.17, -317.43, -4.37, 911.95, 37.73, -667.37, -61.42).
  - Ackermann sobre el simple con polos {-1,-2,-3,-4} da
    K = (-1.75, -3.75, -39.01, -9.60) y reproduce esos polos.

@testset "Controller":
  - solve_care coincide con MatrixEquations.arec (ya esta en Project.toml y no
    se usa), tolerancia 1e-8.
  - P simetrica y definida positiva.
  - Residuo de la CARE: norm(A'P + PA - PBR^-1B'P + Q) < 1e-8.
  - Cayley-Hamilton: el polinomio caracteristico de A evaluado en A es ~0.
  - design_pole_placement reproduce exactamente los polos pedidos.

@testset "Metricas" (si M3 esta implementado):
  - Margen PBH normalizado: simple ~1.72e-2, doble ~2.56e-3, triple ~4.92e-4,
    tolerancia 5% relativa.
  - c* con umax=50: simple ~22.94, doble ~8.99, triple ~3.73.

Que corra con: julia --project=. -e "using Pkg; Pkg.test()".
Documenta en README.md como ejecutarlo.
```

### 6.8 Prompt M9 — Control digital

```text
Extiende el proyecto con analisis de control digital, siguiendo
docs/Entrega_2/plan_de_trabajo.md seccion M9. Reutiliza Metrics.discretize_zoh.

En src/metrics.jl (o en src/digital.jl si prefieres separarlo):

  max_sampling_period(A, B, K; hmin=1e-4, hmax=1.0, iters=40)
      Biseccion geometrica sobre h. Criterio: Ad - Bd K debe ser Schur.

  design_lqr_discrete(Ad, Bd, Qd, Rd)
      LQR discreto resolviendo la DARE por la matriz SIMPLECTICA, en paralelo
      estructural con solve_care y su matriz hamiltoniana. Documentar la
      analogia en el docstring: la CARE usa el hamiltoniano, la DARE el
      simplectico; ambos seleccionan el subespacio invariante estable (Re<0
      en continuo, |z|<1 en discreto).

Script o seccion de main_triple.jl que produzca:
  - Tabla de h_max para las tres configuraciones.
  - Figura comparando la respuesta continua contra la muestreada a
    h = 0.5 h_max, h = 0.9 h_max y h = 1.1 h_max (esta ultima debe divergir).
  - Mapa de polos discretos dentro y fuera del circulo unitario.

Criterio de aceptacion: con los parametros POR DEFECTO de cada configuracion,
h_max debe dar aproximadamente simple 197 ms, doble 77 ms, triple 32 ms
(tolerancia 10%).
Comenta la regla practica observada: h_max * lambda_max vale 0.83, 0.66 y 0.58
respectivamente, es decir, el margen se estrecha al crecer la complejidad.
```

---

## 7. Impacto en los documentos LaTeX

Si se sigue la ruta media o la ambiciosa, esto es lo que hay que tocar. Conviene
tenerlo presente desde el principio para no dejarlo todo al final.

### Informe técnico (`docs/Entrega_1/resumen_tecnico/resumen_tecnico.tex`)

| Sección | Cambio |
|---|---|
| 2. Formalización | Nueva subsección: Configuración III. Preferible: reescribir la deducción en forma **genérica en $N$** y presentar I, II y III como instancias. Es más corto y más elegante que tres deducciones. |
| 3. Herramientas | **Nueva subsección: el test PBH y la distancia a la incontrolabilidad.** Es contenido nuevo de álgebra lineal y justifica la extensión. |
| 3. Herramientas | Si se hace M5: subsección sobre dualidad y principio de separación (la demostración de la triangularidad por bloques es de tres líneas). |
| 5. Resultados | Tercera columna en todas las tablas. Nueva subsección con el atlas de operabilidad y sus mapas. |
| 5. Resultados | Tabla de degradación con $N$ (la de 3.2(c)): es el resultado más contundente. |
| 6. Conclusiones | Reforzar: el mismo álgebra escala, pero **el condicionamiento no**; hay que cambiar de herramienta, no de teoría. |
| Apéndice B | Matrices $\mathbf{M}_0$, $A$, $B$ de la Configuración III (están en 2.7). |
| Apéndice C | Actualizar el árbol de archivos. |

### Resumen ejecutivo

Añadir la columna de la Configuración III a `tab:exec` y una o dos frases sobre
el atlas. Cuidado con el límite de 5 páginas: probablemente haya que recortar
algo. Sugerencia: comprimir la descripción de I y II, que ya están cubiertas.

### Presentación

Frames nuevos sugeridos (manteniendo los 20 minutos, hay que quitar algo):

1. "Tres configuraciones" — sustituye al frame de dos, con diagrama TikZ de tres.
2. "Cuando el rango deja de servir" — la tabla de $\operatorname{cond}(\mathcal{C})$
   contra $N$. Es el frame con más impacto de todo lo nuevo.
3. "El margen PBH" — la definición y el mapa de calor A.
4. "La frontera de operabilidad" — el mapa B con la cota elipsoidal superpuesta.
5. Si se hace M5: "El estado que no se mide" — dualidad y separación.

Recorte sugerido para hacer espacio: comprimir los frames de Cayley-Hamilton y
polinomio minimal en uno solo, y unificar los dos pseudocódigos.

### Guía maestra

Agregar una sección 15 que recoja la segunda entrega, o marcar este documento
como su continuación. Actualizar la sección 9 (resultados) con la tercera
columna y la 3.1 (árbol de archivos).

---

## 8. Checklist de la segunda entrega

**Antes de empezar**

- [ ] Decidir la ruta (sección 5) y los módulos.
- [ ] Confirmar la decisión arquitectónica de 2.9: genérico aditivo o triple autónomo.
- [ ] Crear una rama: `git checkout -b segunda-entrega`.
- [ ] Verificar que Julia está en el `PATH` (hoy no lo está; vive en `~/.julia/juliaup/julia-1.12.5+0.x64.w64.mingw32`).
- [ ] Correr `main_simple.jl` y `main_double.jl` y confirmar que la base sigue reproduciendo los números del informe.

**Durante**

- [ ] Cada módulo nuevo verifica sus criterios de aceptación **ejecutando**, no asumiendo.
- [ ] Ningún cambio en `model_simple.jl`, `model_double.jl` ni en las figuras de la primera entrega.
- [ ] Comentarios en español sin tildes; nada de emojis.
- [ ] Un commit por módulo, en español, formato `tipo: descripcion`.

**Al cerrar**

- [ ] `Pkg.test()` en verde.
- [ ] Regenerar las figuras del informe y de la presentación.
- [ ] Recompilar los tres PDFs y verificar que no falte ninguna figura.
- [ ] Actualizar `README.md` (estructura, resultados esperados, cómo correr los tests y el atlas).
- [ ] Actualizar `docs/guia_maestra.md` (secciones 3.1, 9 y 10).
- [ ] Registrar en el informe **qué se midió** y **qué se dejó propuesto**: la honestidad sobre el alcance vale más que inflar resultados.

---

## Apéndice: resumen de los números de referencia

Todo lo que sigue está calculado y debe reproducirse. Si la implementación no
coincide, el error está en la implementación.

**Configuración III, parámetros por defecto** ($M=1.0$; $m_j=0.2$; $\ell_j=1/3$;
barras uniformes; $g=9.81$; $b=0.1$; $c_j=0$):

| Cantidad | Valor |
|---|---|
| $\beta$ | $(0.16667,\ 0.1,\ 0.03333)$ |
| $J$ | $(0.05185,\ 0.02963,\ 0.00741)$ |
| $\sigma(A)$ | $\{+18.0613,\ +9.7740,\ +4.3777,\ 0,\ -0.0625,\ -4.3990,\ -9.7813,\ -18.0647\}$ |
| Modos inestables | 3 |
| $\operatorname{rank}\mathcal{C}$, $\operatorname{rank}\mathcal{O}$ | $8/8$, $8/8$ |
| $\operatorname{cond}(\mathcal{C})$ | $2.157\times10^{8}$ |
| $\hat\mu$ (PBH normalizado) | $4.916\times10^{-4}$ |
| $B$ | $(0,\ 0.9455,\ 0,\ -3.6,\ 0,\ 0.9818,\ 0,\ -0.3273)^\top$ |
| $K$ (LQR, $Q=\operatorname{diag}(1,0,10,0,10,0,10,0)$, $R=0.1$) | $(-3.16,\ -6.17,\ -317.43,\ -4.37,\ 911.95,\ 37.73,\ -667.37,\ -61.42)$ |
| $\|K\|_2$ | $1176.0$ |
| $\sigma(A-BK)$ | $\{-0.84\pm0.79i,\ -4.57\pm2.32i,\ -10.02\pm1.73i,\ -18.14\pm0.97i\}$ |
| $P_{33}$ | $161.97$ |
| $c^\star$ ($u_{\max}=50$ N) | $3.7291$ |
| $\theta_0$ garantizado sin saturar | $0.152$ rad $\approx 8.7°$ |
| $h_{\max}$ (ZOH) | $32.0$ ms |

**Comparativa de las tres configuraciones** (parámetros por defecto de cada una):

| Cantidad | Simple | Doble | Triple |
|---|---|---|---|
| $n$ | 4 | 6 | 8 |
| Modos inestables | 1 | 2 | 3 |
| $\lambda_{\max}$ | 4.21 | 8.57 | 18.06 |
| $\tau_{\text{caída}}$ [ms] | 237.5 | 116.7 | 55.4 |
| $\operatorname{cond}(\mathcal{C})$ | $3.6\times10^{1}$ | $1.6\times10^{4}$ | $2.2\times10^{8}$ |
| $\hat\mu$ | $1.72\times10^{-2}$ | $2.56\times10^{-3}$ | $4.92\times10^{-4}$ |
| $\|K\|_2$ | 47.0 | 300.5 | 1176.0 |
| $c^\star$ ($u_{\max}=50$ N) | 22.94 | 8.99 | 3.73 |
| $\theta_0$ sin saturar | $55.1°$ | $13.4°$ | $8.7°$ |
| $h_{\max}$ [ms] | 197.4 | 77.3 | 32.0 |
| $h_{\max}\cdot\lambda_{\max}$ | 0.83 | 0.66 | 0.58 |

**Escalamiento con $N$** (familia normalizada: longitud total 1 m, masa total
0.6 kg, barras uniformes, $M=1$, $b=0.1$). Ojo: esta familia **no** coincide con
los parámetros por defecto de I y II, que tienen otras masas; es una serie
aparte, construida para que la comparación sea limpia.

| $N$ | $n$ | $\lambda_{\max}$ | $\tau$ [ms] | $\hat\mu$ | $\operatorname{cond}\mathcal{C}$ | $\|K\|_2$ |
|---|---|---|---|---|---|---|
| 1 | 4 | 4.51 | 222 | $1.31\times10^{-2}$ | $4.8\times10^{1}$ | 52.9 |
| 2 | 6 | 10.59 | 94 | $1.79\times10^{-3}$ | $6.2\times10^{4}$ | 248.7 |
| 3 | 8 | 18.06 | 55 | $4.92\times10^{-4}$ | $2.2\times10^{8}$ | 1176.0 |
| 4 | 10 | 26.54 | 38 | $1.90\times10^{-4}$ | $1.7\times10^{12}$ | 5266.8 |
| 5 | 12 | 35.64 | 28 | $9.03\times10^{-5}$ | $3.3\times10^{16}$ | 23221.9 |
