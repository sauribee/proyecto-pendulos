# Fundamentos de Álgebra Lineal para el Proyecto del Péndulo Invertido

## Resumen Conceptual Exhaustivo — Edición Corregida

**Referencia principal:** Olver, P. J. & Shakiban, C. — *Applied Linear Algebra*, 2nd Edition. Springer, 2018.
**Referencia complementaria:** Ogata, K. — *Ingeniería de Control Moderna*, 5ª Ed. Pearson, 2010.

---

## 1. Sistemas de Ecuaciones Lineales y Eliminación Gaussiana

**[Olver & Shakiban, Capítulo 1: Linear Algebraic Systems]**

### 1.1 El problema fundamental

Todo sistema de ecuaciones lineales se expresa en la forma matricial

$$
A\mathbf{x} = \mathbf{b},
$$

donde $A \in \mathbb{R}^{m \times n}$ es la matriz de coeficientes, $\mathbf{x} \in \mathbb{R}^n$ el vector de incógnitas y $\mathbf{b} \in \mathbb{R}^m$ el vector de términos independientes.

### 1.2 Eliminación gaussiana y factorización LU

**[Olver, §1.3: Gaussian Elimination — Regular Case]**

La eliminación gaussiana transforma la matriz aumentada $[A | \mathbf{b}]$ en forma escalonada mediante operaciones elementales por filas. Estas operaciones se codifican en **matrices elementales** $E_k$ (§1.3, subsección "Elementary Matrices"), de modo que:

$$
E_k \cdots E_2 E_1 A = U
$$

donde $U$ es triangular superior. Definiendo $L = E_1^{-1} E_2^{-1} \cdots E_k^{-1}$ (triangular inferior unitriangular, es decir, con unos en la diagonal), se obtiene la **factorización LU** (§1.3, subsección "The LU Factorization"):

$$
A = LU
$$

Olver denomina **regular** a una matriz cuadrada que admite factorización LU sin intercambio de filas. Cuando se requieren intercambios, se obtiene la **factorización LU permutada** $PA = LU$ (§1.4, subsección "The Permuted LU Factorization").

Para matrices simétricas, la factorización toma la forma especial $A = LDL^T$ (§1.6, subsección "Factorization of Symmetric Matrices"), donde $D$ es diagonal con los pivotes como entradas.

**Relevancia para el proyecto:** Resolver $A\mathbf{x} = \mathbf{b}$ es la operación atómica que subyace al cálculo de ganancias de control, inversión de matrices de controlabilidad, y resolución numérica de la ecuación de Riccati.

### 1.3 Matrices invertibles

**[Olver, §1.5: Matrix Inverses]**

**Definición.** Una matriz cuadrada $A \in \mathbb{R}^{n \times n}$ es **invertible** (o no singular) si existe $A^{-1} \in \mathbb{R}^{n \times n}$ tal que $AA^{-1} = A^{-1}A = I_n$.

Las siguientes condiciones son equivalentes (este resultado se va construyendo a lo largo de los capítulos 1 y 2, y se consolida con los valores propios en el capítulo 8):

1. $A$ es invertible.
2. $\det(A) \neq 0$ (§1.9).
3. $\text{rango}(A) = n$ (§2.4–2.5).
4. El sistema $A\mathbf{x} = \mathbf{b}$ tiene solución única para todo $\mathbf{b}$.
5. $\ker(A) = \{\mathbf{0}\}$ (§2.5).
6. Todos los valores propios de $A$ son no nulos (§8.2).

### 1.4 Determinantes

**[Olver, §1.9: Determinants]**

Olver presenta el determinante mediante la fórmula combinatoria general que involucra permutaciones (fórmula (1.87) en el texto). Las propiedades fundamentales incluyen:

$$
\det(AB) = \det(A)\det(B), \quad \det(A^T) = \det(A), \quad \det(A^{-1}) = \frac{1}{\det(A)}
$$

**Relevancia para el proyecto:** El determinante aparece en el cálculo del polinomio característico $p_A(\lambda) = \det(A - \lambda I)$ para hallar los valores propios del sistema, y su no-nulidad es condición necesaria para la invertibilidad de la matriz de controlabilidad.

---

## 2. Espacios Vectoriales, Bases y Dimensión

**[Olver & Shakiban, Capítulo 2: Vector Spaces and Bases]**

### 2.1 Espacio vectorial

**[Olver, §2.1: Real Vector Spaces]**

**Definición (Olver, §2.1).** Un **espacio vectorial real** $V$ es un conjunto dotado de suma vectorial y multiplicación por escalar que satisfacen los ocho axiomas estándar (clausura, asociatividad, conmutatividad, existencia de neutro aditivo y opuesto, distributividad, etc.).

Los espacios relevantes para el proyecto son $\mathbb{R}^n$ (espacio de estados del péndulo, con $n=4$), $\mathbb{R}^{m \times n}$ (espacio de matrices), y los espacios de funciones continuas $C^0[a,b]$ (relevantes para los fundamentos teóricos del producto interno).

### 2.2 Subespacios

**[Olver, §2.2: Subspaces]**

**Definición.** Un subconjunto $W \subseteq V$ es un **subespacio** si es cerrado bajo combinaciones lineales: para todo $\mathbf{u}, \mathbf{v} \in W$ y $\alpha, \beta \in \mathbb{R}$, se cumple $\alpha\mathbf{u} + \beta\mathbf{v} \in W$.

### 2.3 Generación, independencia lineal, base y dimensión

**[Olver, §2.3: Span and Linear Independence; §2.4: Basis and Dimension]**

**Definición (§2.3).** Un conjunto $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$ es **linealmente independiente** si:

$$
c_1 \mathbf{v}_1 + c_2 \mathbf{v}_2 + \cdots + c_k \mathbf{v}_k = \mathbf{0} \implies c_1 = c_2 = \cdots = c_k = 0
$$

**Definición (§2.4).** Una **base** de $V$ es un conjunto linealmente independiente que genera $V$. La **dimensión** $\dim(V)$ es la cardinalidad de cualquier base.

### 2.4 Los subespacios fundamentales de una matriz

**[Olver, §2.5: The Fundamental Matrix Subspaces]**

Dada $A \in \mathbb{R}^{m \times n}$, Olver define los cuatro subespacios fundamentales en §2.5:

- **Kernel (núcleo):** $\ker(A) = \{\mathbf{x} \in \mathbb{R}^n : A\mathbf{x} = \mathbf{0}\}$ (§2.5, subsección "Kernel and Image")
- **Imagen:** $\text{img}(A) = \{A\mathbf{x} : \mathbf{x} \in \mathbb{R}^n\}$ (§2.5, subsección "Kernel and Image")
- **Coimagen:** $\text{coimg}(A) = \text{img}(A^T)$ (§2.5, subsección "Adjoint Systems, Cokernel, and Coimage")
- **Cokernel:** $\text{coker}(A) = \ker(A^T)$ (§2.5, subsección "Adjoint Systems, Cokernel, and Coimage")

> **Nota terminológica:** Olver usa los términos *kernel* e *image* (no "espacio nulo" y "espacio columna"), y *cokernel* y *coimage* para los subespacios de $A^T$.

**Teorema Fundamental del Álgebra Lineal (Olver, §2.5, subsección homónima).**

$$
\dim(\ker A) + \dim(\text{img}\, A) = n
$$

y los cuatro subespacios satisfacen relaciones de complementariedad ortogonal (resultado que se completa en §4.4 con la ortogonalidad).

### 2.5 Rango matricial

**[Olver, §2.4–2.5]**

El **rango** de $A$ es $\text{rango}(A) = \dim(\text{img}\, A)$, equivalentemente, el número de pivotes en la forma escalonada. Se cumple $\text{rango}(A) = \text{rango}(A^T)$.

**Relevancia para el proyecto:** La controlabilidad del sistema se verifica comprobando que $\text{rango}(\mathcal{C}) = n$. Análogamente, la observabilidad exige que $\ker(\mathcal{O}) = \{\mathbf{0}\}$, es decir, $\text{rango}(\mathcal{O}) = n$.

---

## 3. Productos Internos, Normas y Matrices Definidas Positivas

**[Olver & Shakiban, Capítulo 3: Inner Products and Norms]**

### 3.1 Producto interno

**[Olver, §3.1: Inner Products]**

**Definición 3.1 (Olver).** Un **producto interno** en un espacio vectorial real $V$ es una función $\langle \cdot, \cdot \rangle: V \times V \to \mathbb{R}$ que satisface, para todo $\mathbf{u}, \mathbf{v}, \mathbf{w} \in V$ y escalares $c, d \in \mathbb{R}$:

(i) **Bilinealidad:**
$$\langle c\mathbf{u} + d\mathbf{v},\, \mathbf{w} \rangle = c\langle \mathbf{u}, \mathbf{w} \rangle + d\langle \mathbf{v}, \mathbf{w} \rangle$$
$$\langle \mathbf{u},\, c\mathbf{v} + d\mathbf{w} \rangle = c\langle \mathbf{u}, \mathbf{v} \rangle + d\langle \mathbf{u}, \mathbf{w} \rangle$$

(ii) **Simetría:** $\langle \mathbf{v}, \mathbf{w} \rangle = \langle \mathbf{w}, \mathbf{v} \rangle$

(iii) **Positividad:** $\langle \mathbf{v}, \mathbf{v} \rangle > 0$ cuando $\mathbf{v} \neq \mathbf{0}$, y $\langle \mathbf{0}, \mathbf{0} \rangle = 0$.

La **norma inducida** es $\|\mathbf{v}\| = \sqrt{\langle \mathbf{v}, \mathbf{v} \rangle}$ (ecuación (3.7) de Olver).

En $\mathbb{R}^n$, el producto punto estándar es $\mathbf{v} \cdot \mathbf{w} = \mathbf{v}^T \mathbf{w}$. Olver también desarrolla productos internos en espacios de funciones: el producto interno $L^2$ es $\langle f, g \rangle = \int_a^b f(x)g(x)\,dx$ (Ejemplo 3.4).

### 3.2 Desigualdades fundamentales

**[Olver, §3.2: Inequalities]**

**Teorema 3.5 (Desigualdad de Cauchy–Schwarz).** Para todo producto interno:

$$
|\langle \mathbf{v}, \mathbf{w} \rangle| \leq \|\mathbf{v}\| \cdot \|\mathbf{w}\|
$$

con igualdad si y solo si $\mathbf{v}$ y $\mathbf{w}$ son paralelos.

La demostración de Olver (p. 138) es puramente algebraica: se expande $0 \leq \|\mathbf{v} + t\mathbf{w}\|^2$ como un polinomio cuadrático en $t$ y se evalúa en su mínimo $t = -\langle \mathbf{v}, \mathbf{w} \rangle / \|\mathbf{w}\|^2$.

**Teorema 3.9 (Desigualdad triangular).** La norma asociada a un producto interno satisface:

$$
\|\mathbf{v} + \mathbf{w}\| \leq \|\mathbf{v}\| + \|\mathbf{w}\|
$$

Olver demuestra que la desigualdad triangular es consecuencia directa de Cauchy–Schwarz.

**Ortogonalidad (§3.2, subsección "Orthogonal Vectors").** Dos vectores son **ortogonales** si $\langle \mathbf{v}, \mathbf{w} \rangle = 0$ (Definición 3.6).

### 3.3 Normas generales

**[Olver, §3.3: Norms]**

**Definición 3.12 (Olver).** Una **norma** en un espacio vectorial $V$ asigna un número real no negativo $\|\mathbf{v}\|$ a cada $\mathbf{v} \in V$, satisfaciendo: (i) Positividad, (ii) Homogeneidad: $\|c\mathbf{v}\| = |c|\|\mathbf{v}\|$, (iii) Desigualdad triangular.

Olver presenta las normas $\ell^1$, $\ell^2$ (euclídea), y $\ell^\infty$ en $\mathbb{R}^n$ (ecuaciones (3.26)–(3.28)), así como las correspondientes normas $L^p$ en espacios de funciones. También desarrolla **normas matriciales naturales** (Teorema 3.20, Definición 3.23) y la propiedad multiplicativa $\|AB\| \leq \|A\|\|B\|$ (Teorema 3.21, desigualdad (3.42)).

### 3.4 Matrices definidas positivas

**[Olver, §3.4: Positive Definite Matrices]**

**Motivación.** Olver llega a las matrices definidas positivas buscando la forma más general de un producto interno en $\mathbb{R}^n$. Escribiendo $\mathbf{x}$ y $\mathbf{y}$ en la base canónica y usando bilinealidad:

$$
\langle \mathbf{x}, \mathbf{y} \rangle = \sum_{i,j=1}^n k_{ij}\, x_i\, y_j = \mathbf{x}^T K \mathbf{y}
$$

donde $k_{ij} = \langle \mathbf{e}_i, \mathbf{e}_j \rangle$ (ecuación (3.47) de Olver). La simetría del producto interno exige $K^T = K$, y la positividad exige $\mathbf{x}^T K \mathbf{x} > 0$ para todo $\mathbf{x} \neq \mathbf{0}$.

**Definición 3.26 (Olver).** Una matriz $K \in \mathbb{R}^{n \times n}$ es **definida positiva** si es simétrica ($K^T = K$) y satisface:

$$
\mathbf{x}^T K \mathbf{x} > 0 \quad \text{para todo } \mathbf{0} \neq \mathbf{x} \in \mathbb{R}^n
$$

Se escribe $K > 0$. La advertencia de Olver es crucial: "$K > 0$ no significa que todas las entradas de $K$ sean positivas".

**Teorema 3.27 (Olver).** Todo producto interno en $\mathbb{R}^n$ tiene la forma $\langle \mathbf{x}, \mathbf{y} \rangle = \mathbf{x}^T K \mathbf{y}$, donde $K$ es una matriz simétrica definida positiva.

**Proposición 3.31 (Olver).** Si $K$ es definida positiva, entonces $K$ es no singular.

La demostración es inmediata: si existiera $\mathbf{x} \neq \mathbf{0}$ con $K\mathbf{x} = \mathbf{0}$, entonces $\mathbf{x}^T K \mathbf{x} = \mathbf{x}^T \mathbf{0} = 0$, contradiciendo la positividad.

**Formas cuadráticas (§3.4, ecuación (3.52)).** Dada una matriz simétrica $K$, la **forma cuadrática** asociada es:

$$
q(\mathbf{x}) = \mathbf{x}^T K \mathbf{x} = \sum_{i,j=1}^n k_{ij}\, x_i\, x_j
$$

El carácter de $q$ (definida positiva, semidefinida, indefinida) clasifica a la matriz.

**Matrices semidefinidas positivas e indefinidas.** $K \geq 0$ (semidefinida positiva) si $\mathbf{x}^T K \mathbf{x} \geq 0$ para todo $\mathbf{x}$; puede admitir **direcciones nulas** $\mathbf{z} \neq \mathbf{0}$ con $q(\mathbf{z}) = 0$ (ecuación (3.56) de Olver). Una forma cuadrática es **indefinida** si toma valores tanto positivos como negativos.

### 3.5 Criterio de pivotes positivos (completar el cuadrado)

**[Olver, §3.5: Completing the Square]**

Olver demuestra que el proceso de completar el cuadrado para una forma cuadrática es algebraicamente idéntico a la factorización $K = LDL^T$ de la eliminación gaussiana (ecuación (3.73)–(3.74)).

**Teorema 3.43 (Olver).** Una matriz simétrica es definida positiva si y solo si es regular (admite factorización $LDL^T$ sin intercambio de filas) y tiene **todos los pivotes positivos**.

Equivalentemente, $K > 0$ si y solo si $K = LDL^T$ con $L$ triangular inferior unitriangular y $D$ diagonal con todas las entradas positivas.

**Factorización de Cholesky (§3.5, subsección homónima).** Si $K > 0$, podemos escribir $K = R^T R$ donde $R = D^{1/2}L^T$ es triangular superior con entradas diagonales positivas. Esta es la factorización de Cholesky (p. 171).

### 3.6 Matrices de Gram

**[Olver, §3.4, subsección "Gram Matrices"]**

**Definición 3.33 (Olver).** Sea $V$ un espacio con producto interno y $\mathbf{v}_1, \ldots, \mathbf{v}_n \in V$. La **matriz de Gram** asociada es la matriz $n \times n$ con entradas $k_{ij} = \langle \mathbf{v}_i, \mathbf{v}_j \rangle$.

**Teorema 3.34 (Olver).** Toda matriz de Gram es semidefinida positiva. La matriz de Gram es definida positiva si y solo si $\mathbf{v}_1, \ldots, \mathbf{v}_n$ son linealmente independientes.

La demostración es elegante: $q(\mathbf{x}) = \mathbf{x}^T K \mathbf{x} = \langle \sum x_i \mathbf{v}_i, \sum x_j \mathbf{v}_j \rangle = \|\mathbf{v}\|^2 \geq 0$.

En el caso del producto punto euclídeo, si $A = (\mathbf{v}_1 \cdots \mathbf{v}_n)$, entonces la matriz de Gram es simplemente $K = A^T A$ (ecuación (3.62) de Olver).

**Relevancia para el proyecto:**

- Las matrices de peso $Q$ y $R$ en el funcional de costo LQR deben ser simétricas: $Q \succeq 0$ (semidefinida positiva), $R \succ 0$ (definida positiva).
- La solución $P$ de la ecuación de Riccati es simétrica y definida positiva.
- El funcional $J = \int_0^\infty (\mathbf{x}^T Q \mathbf{x} + \mathbf{u}^T R \mathbf{u})\,dt$ está bien definido y es no negativo precisamente porque $Q \succeq 0$ y $R \succ 0$.

---

## 4. Ortogonalidad

**[Olver & Shakiban, Capítulo 4: Orthogonality]**

### 4.1 Bases ortonormales

**[Olver, §4.1: Orthogonal and Orthonormal Bases]**

Un conjunto $\{\mathbf{q}_1, \ldots, \mathbf{q}_k\}$ es **ortonormal** si $\langle \mathbf{q}_i, \mathbf{q}_j \rangle = \delta_{ij}$.

La ventaja computacional clave (§4.1, subsección "Computations in Orthogonal Bases"): las coordenadas de un vector $\mathbf{v}$ en una base ortonormal se calculan por simple producto interno $c_i = \langle \mathbf{v}, \mathbf{q}_i \rangle$, sin necesidad de resolver un sistema lineal.

### 4.2 Proceso de Gram–Schmidt

**[Olver, §4.2: The Gram–Schmidt Process]**

Dado un conjunto linealmente independiente $\{\mathbf{v}_1, \ldots, \mathbf{v}_k\}$, el proceso produce un conjunto ortonormal $\{\mathbf{q}_1, \ldots, \mathbf{q}_k\}$:

$$
\tilde{\mathbf{q}}_j = \mathbf{v}_j - \sum_{i=1}^{j-1} \langle \mathbf{v}_j, \mathbf{q}_i \rangle \, \mathbf{q}_i, \qquad \mathbf{q}_j = \frac{\tilde{\mathbf{q}}_j}{\|\tilde{\mathbf{q}}_j\|}
$$

Olver también presenta la versión de **Gram–Schmidt modificado** (§4.2, subsección "Modifications of the Gram–Schmidt Process") que es numéricamente más estable.

### 4.3 Matrices ortogonales y factorización QR

**[Olver, §4.3: Orthogonal Matrices]**

**Definición 4.18 (Olver, mencionada en §4.3).** $Q \in \mathbb{R}^{n \times n}$ es **ortogonal** si $Q^{-1} = Q^T$, equivalentemente, si sus columnas forman una base ortonormal de $\mathbb{R}^n$.

Propiedades: preserva normas ($\|Q\mathbf{x}\| = \|\mathbf{x}\|$), preserva productos internos, $\det(Q) = \pm 1$.

**Factorización QR (§4.3, subsección "The QR Factorization").** El proceso de Gram–Schmidt aplicado a las columnas de $A$ produce $A = QR$, donde $Q$ tiene columnas ortonormales y $R$ es triangular superior.

### 4.4 Proyecciones ortogonales y la alternativa de Fredholm

**[Olver, §4.4: Orthogonal Projections and Orthogonal Subspaces]**

Aquí Olver completa la teoría de los subespacios fundamentales mostrando que los cuatro subespacios de una matriz vienen en pares ortogonales (§4.4, subsección "Orthogonality of the Fundamental Matrix Subspaces and the Fredholm Alternative"):

$$
(\text{img}\, A)^\perp = \text{coker}\, A, \qquad (\ker A)^\perp = \text{coimg}\, A
$$

**Relevancia para el proyecto:** La factorización QR se emplea en el algoritmo QR iterativo para calcular valores propios (§9.5) y en métodos numéricamente estables para resolver problemas de mínimos cuadrados.

---

## 5. Minimización y Mínimos Cuadrados

**[Olver & Shakiban, Capítulo 5: Minimization and Least Squares]**

### 5.1 Minimización de funciones cuadráticas

**[Olver, §5.2: Minimization of Quadratic Functions]**

**Teorema 5.2 (Olver).** La función cuadrática $p(\mathbf{x}) = \frac{1}{2}\mathbf{x}^T K \mathbf{x} - \mathbf{f}^T \mathbf{x} + c$ tiene un mínimo único si y solo si $K$ es definida positiva, y el minimizador satisface $K\mathbf{x} = \mathbf{f}$.

Olver enfatiza: "la minimización de funciones cuadráticas es un problema de álgebra lineal" (p. 243).

### 5.2 Mínimos cuadrados

**[Olver, §5.3–5.4: The Closest Point; Least Squares]**

Cuando el sistema $A\mathbf{x} = \mathbf{b}$ no tiene solución exacta, el vector $\hat{\mathbf{x}}$ que minimiza $\|A\mathbf{x} - \mathbf{b}\|^2$ satisface las **ecuaciones normales** (Olver, §5.4, "normal equations"):

$$
A^T A \hat{\mathbf{x}} = A^T \mathbf{b}
$$

Si las columnas de $A$ son linealmente independientes, $A^T A$ es invertible (y de hecho definida positiva, al ser una matriz de Gram por la ecuación (3.62)) y la solución es $\hat{\mathbf{x}} = (A^T A)^{-1} A^T \mathbf{b}$.

**Relevancia para el proyecto:** La estructura $A^T A$ aparece en la teoría de optimización detrás del LQR. Si se usaran datos experimentales para estimar parámetros del péndulo, esto conduciría a un problema de mínimos cuadrados.

---

## 6. Linealidad: Transformaciones Lineales y Cambio de Base

**[Olver & Shakiban, Capítulo 7: Linearity]**

### 6.1 Funciones y transformaciones lineales

**[Olver, §7.1: Linear Functions; §7.2: Linear Transformations]**

**Definición (§7.1).** Una **función lineal** $L: V \to W$ satisface $L(\alpha \mathbf{u} + \beta \mathbf{v}) = \alpha L(\mathbf{u}) + \beta L(\mathbf{v})$.

Toda transformación lineal $L: \mathbb{R}^n \to \mathbb{R}^m$ se representa por una matriz $A$: $L(\mathbf{x}) = A\mathbf{x}$. Las columnas de $A$ son las imágenes de los vectores de la base canónica.

### 6.2 Cambio de base y matrices semejantes

**[Olver, §7.2, subsección "Change of Basis"]**

Si $S$ es la matriz de cambio de base, la representación de la transformación lineal en la nueva base es:

$$
\tilde{A} = S^{-1} A S
$$

Matrices $A$ y $\tilde{A}$ son **semejantes** y comparten: valores propios, determinante, traza, rango, polinomio característico.

**Relevancia para el proyecto:** La representación $\dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u}$ depende de la elección de variables de estado. Un cambio de base $\mathbf{x} = S\mathbf{z}$ transforma el sistema a $\dot{\mathbf{z}} = S^{-1}AS\,\mathbf{z} + S^{-1}B\,\mathbf{u}$, preservando las propiedades dinámicas esenciales.

### 6.3 Operadores autoadjuntos y definidos positivos

**[Olver, §7.5: Adjoints, Positive Definite Operators, and Minimization Principles]**

Olver generaliza la noción de simetría al concepto de **operador autoadjunto** ($L^* = L$) y de **operador definido positivo** (§7.5, subsección "Self-Adjoint and Positive Definite Linear Functions"), mostrando que los resultados de §3.4 para matrices se extienden a operadores lineales generales en espacios con producto interno.

---

## 7. Valores Propios y Vectores Propios

**[Olver & Shakiban, Capítulo 8: Eigenvalues and Singular Values]**

### 7.1 Motivación dinámica

**[Olver, §8.1: Linear Dynamical Systems]**

Olver motiva los valores propios a través de la ecuación diferencial escalar $\dot{u} = au$ con solución $u(t) = ce^{at}$ (§8.1, subsección "Scalar Ordinary Differential Equations"), y su extensión a sistemas $\dot{\mathbf{u}} = A\mathbf{u}$ buscando soluciones de la forma $\mathbf{u}(t) = e^{\lambda t}\mathbf{v}$ (§8.1, subsección "First Order Dynamical Systems").

### 7.2 Definición

**[Olver, §8.2: Eigenvalues and Eigenvectors]**

Sea $A \in \mathbb{R}^{n \times n}$. Un escalar $\lambda$ (posiblemente complejo) es un **valor propio** de $A$ si existe un vector no nulo $\mathbf{v} \neq \mathbf{0}$ tal que $A\mathbf{v} = \lambda \mathbf{v}$. El vector $\mathbf{v}$ es el **vector propio** asociado.

### 7.3 Polinomio característico

**[Olver, §8.2, subsección "Basic Properties of Eigenvalues"]**

Los valores propios son las raíces del **polinomio característico** (ecuación (8.17) de Olver):

$$
p_A(\lambda) = \det(A - \lambda I)
$$

**Teorema 8.11 (Olver).** Una matriz $n \times n$ posee al menos uno y a lo sumo $n$ valores propios complejos distintos.

**Proposición 8.13 (Olver).** La suma de los valores propios iguala la traza, y el producto iguala el determinante:

$$
\lambda_1 + \lambda_2 + \cdots + \lambda_n = \text{tr}(A), \qquad \lambda_1 \lambda_2 \cdots \lambda_n = \det(A)
$$

(ecuaciones (8.25) y (8.26) de Olver).

### 7.4 Ejemplo concreto del proyecto

Para la matriz del sistema linealizado del péndulo invertido (Ogata, Ejemplo 10-5):

$$
A = \begin{pmatrix} 0 & 1 & 0 & 0 \\ 20.601 & 0 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ -0.4905 & 0 & 0 & 0 \end{pmatrix}
$$

El polinomio característico es $p(\lambda) = \lambda^4 - 20.601\lambda^2$. Las raíces son $\lambda_1 \approx +4.539$, $\lambda_2 \approx -4.539$, $\lambda_{3,4} = 0$. La existencia de $\lambda_1 > 0$ confirma la **inestabilidad** del punto de equilibrio (por el Teorema 10.16 que veremos en §10).

### 7.5 Multiplicidad algebraica y geométrica

**[Olver, §8.2–8.3]**

- **Multiplicidad algebraica** $m_a(\lambda)$: multiplicidad de $\lambda$ como raíz de $p_A(\lambda)$.
- **Multiplicidad geométrica**: $\dim(\ker(A - \lambda I))$.

Siempre se tiene $1 \leq m_g(\lambda) \leq m_a(\lambda)$.

### 7.6 Diagonalización y matrices completas

**[Olver, §8.3: Eigenvector Bases]**

Olver usa el término **completa** (*complete*) para una matriz cuyos vectores propios (complejos) forman una base de $\mathbb{C}^n$. Esto es equivalente a la **diagonalización**:

$$
A = S \Lambda S^{-1}, \qquad \Lambda = \text{diag}(\lambda_1, \ldots, \lambda_n)
$$

(§8.3, subsección "Diagonalization"). Las columnas de $S$ son los vectores propios.

**Criterio:** $A$ es completa (diagonalizable) si y solo si $m_g(\lambda_i) = m_a(\lambda_i)$ para todo valor propio. En particular, si $A$ tiene $n$ valores propios distintos, es completa.

**Importancia computacional:** Si $A = S\Lambda S^{-1}$, entonces $A^k = S\Lambda^k S^{-1}$ y $e^{At} = Se^{\Lambda t}S^{-1}$.

### 7.7 Teorema de Cayley–Hamilton

**[Olver, Capítulo 8; mencionado en el índice p. 420]**

Toda matriz satisface su propio polinomio característico: $p_A(A) = 0$.

**Relevancia para el proyecto:** Es la base teórica del método de Ackermann para colocación de polos.

---

## 8. Matrices Simétricas y el Teorema Espectral

**[Olver, §8.5: Eigenvalues of Symmetric Matrices]**

### 8.1 El Teorema Espectral

**Teorema 8.32 (Olver, §8.5).** Sea $A = A^T$ una matriz simétrica real $n \times n$. Entonces:

(a) Todos los valores propios de $A$ son **reales**.

(b) Vectores propios correspondientes a valores propios **distintos** son **ortogonales**.

(c) Existe una **base ortonormal** de $\mathbb{R}^n$ formada por $n$ vectores propios de $A$.

En particular, toda matriz simétrica real es completa y diagonalizable, con **diagonalización ortogonal**:

$$
A = Q\Lambda Q^T
$$

donde $Q$ es ortogonal ($Q^T Q = I$) y $\Lambda = \text{diag}(\lambda_1, \ldots, \lambda_n)$.

### 8.2 Conexión con positividad definida

**Teorema 8.35 (Olver, §8.5).** Una matriz simétrica $K = K^T$ es definida positiva si y solo si todos sus valores propios son estrictamente positivos.

La demostración de Olver usa la base ortonormal de vectores propios del Teorema 8.32: si $\mathbf{x} = c_1\mathbf{u}_1 + \cdots + c_n\mathbf{u}_n$, entonces $\mathbf{x}^T K\mathbf{x} = \lambda_1 c_1^2 + \cdots + \lambda_n c_n^2$, que es positivo para todo $\mathbf{x} \neq \mathbf{0}$ si y solo si todos los $\lambda_i > 0$.

**Resumen de criterios equivalentes para $K > 0$** (combinando §3.4, §3.5 y §8.5):

1. $\mathbf{x}^T K \mathbf{x} > 0$ para todo $\mathbf{x} \neq \mathbf{0}$ (Definición 3.26).
2. Todos los pivotes de $K$ son positivos (Teorema 3.43).
3. Todos los valores propios de $K$ son positivos (Teorema 8.35).
4. $K = LDL^T$ con $D$ diagonal positiva (Teorema 3.43 / §3.5).
5. $K = R^T R$ con $R$ invertible — factorización de Cholesky (§3.5).
6. $K$ define un producto interno: $\langle \mathbf{x}, \mathbf{y} \rangle_K = \mathbf{x}^T K \mathbf{y}$ (Teorema 3.27).

---

## 9. Descomposición en Valores Singulares (SVD)

**[Olver & Shakiban, §8.7: Singular Values]**

### 9.1 Valores singulares

Los **valores singulares** de $A$ son las raíces cuadradas positivas de los valores propios no nulos de la matriz de Gram $K = A^T A$ (que es simétrica y semidefinida positiva). Los vectores propios de $K$ son los **vectores singulares (derechos)** de $A$.

### 9.2 La factorización SVD

**Teorema 8.63 (Olver).** Toda matriz real no nula $A \in \mathbb{R}^{m \times n}$ de rango $r > 0$ admite la factorización:

$$
A = P \Sigma Q^T
$$

donde $P \in \mathbb{R}^{m \times r}$ tiene columnas ortonormales ($P^T P = I$), $\Sigma = \text{diag}(\sigma_1, \ldots, \sigma_r)$ con $\sigma_1 \geq \cdots \geq \sigma_r > 0$, y $Q \in \mathbb{R}^{n \times r}$ tiene columnas ortonormales ($Q^T Q = I$).

La demostración de Olver (p. 455–456) construye los vectores singulares izquierdos como $\mathbf{p}_i = A\mathbf{q}_i / \sigma_i$.

### 9.3 Norma matricial euclídea y número de condición

**Teorema 8.71 (Olver).** La norma matricial euclídea de $A$ iguala su valor singular dominante: $\|A\|_2 = \sigma_1$.

El **número de condición** $\kappa(A) = \sigma_1 / \sigma_r$ (§8.7, subsección "Condition Number and Rank") cuantifica la sensibilidad numérica del sistema.

### 9.4 Pseudoinversa

**Definición 8.67 (Olver).** La **pseudoinversa** de $A = P\Sigma Q^T$ es $A^+ = Q\Sigma^{-1}P^T$.

Si $A$ tiene columnas linealmente independientes: $A^+ = (A^T A)^{-1}A^T$ (Lema 8.68).

**Relevancia para el proyecto:** La SVD proporciona una medida de la "calidad" de la controlabilidad. Si la matriz de controlabilidad $\mathcal{C}$ tiene un valor singular muy pequeño, el sistema es *casi no controlable* en cierta dirección, lo que requiere mucha energía de control.

---

## 10. Iteración y Radio Espectral

**[Olver & Shakiban, Capítulo 9: Iteration]**

### 10.1 Sistemas iterativos lineales

**[Olver, §9.1: Linear Iterative Systems]**

El sistema iterativo $\mathbf{u}^{(k+1)} = T\mathbf{u}^{(k)}$ tiene solución $\mathbf{u}^{(k)} = T^k \mathbf{a}$.

### 10.2 Radio espectral y estabilidad discreta

**[Olver, §9.2: Stability]**

**Definición 9.13 (Olver).** El **radio espectral** es $\rho(T) = \max_i |\lambda_i|$.

**Teorema 9.14 (Olver).** La matriz $T$ es **convergente** ($T^k \to 0$) si y solo si $\rho(T) < 1$.

**Teorema 9.12 (Olver).** El equilibrio cero del sistema iterativo $\mathbf{u}^{(k+1)} = T\mathbf{u}^{(k)}$ es globalmente asintóticamente estable si y solo si todos los valores propios satisfacen $|\lambda_j| < 1$ (están dentro del **círculo unitario**).

**Relevancia para el proyecto:** Si el controlador del péndulo se implementa digitalmente, la estabilidad del sistema discretizado requiere que los polos del sistema en lazo cerrado estén dentro del círculo unitario.

---

## 11. Sistemas Dinámicos Lineales y Exponencial Matricial

**[Olver & Shakiban, Capítulo 10: Dynamics]**

### 11.1 Soluciones propias

**[Olver, §10.1: Basic Solution Techniques]**

El sistema $\dot{\mathbf{u}} = A\mathbf{u}$ tiene soluciones de la forma $\mathbf{u}(t) = e^{\lambda t}\mathbf{v}$ cuando $\lambda$ es valor propio y $\mathbf{v}$ vector propio de $A$.

Si $A$ es completa con $A = S\Lambda S^{-1}$, la **solución general** es (§10.1, subsección "Complete Systems"):

$$
\mathbf{u}(t) = \sum_{i=1}^n c_i\, e^{\lambda_i t}\, \mathbf{v}_i
$$

donde los $c_i$ se determinan por la condición inicial: $\mathbf{c} = S^{-1}\mathbf{u}_0$.

Para matrices incompletas, Olver desarrolla las soluciones mediante **cadenas de Jordan** (§10.1, subsección "The General Case"), con soluciones que involucran potencias de $t$ multiplicando exponenciales.

### 11.2 Exponencial matricial

**[Olver, §10.4: Matrix Exponentials]**

**Definición.** La exponencial matricial se define como la solución del problema de valor inicial matricial $\frac{d}{dt}e^{tA} = Ae^{tA}$, $e^{0} = I$ (ecuación (10.39) de Olver), y tiene la serie de potencias (ecuación (10.47)):

$$
e^{tA} = \sum_{k=0}^{\infty} \frac{(tA)^k}{k!} = I + tA + \frac{t^2 A^2}{2!} + \frac{t^3 A^3}{3!} + \cdots
$$

**Propiedades fundamentales (Olver, §10.4):**

1. $e^{0} = I$
2. $\frac{d}{dt} e^{tA} = Ae^{tA} = e^{tA}A$
3. $(e^{tA})^{-1} = e^{-tA}$ (ecuación (10.45))
4. Si $AB = BA$: $e^{t(A+B)} = e^{tA}e^{tB}$ (pero **no** en general; ecuación (10.46))
5. $\det(e^{tA}) = e^{t\,\text{tr}(A)}$ (Lema 10.28)

**Cálculo práctico.** Si $A = S\Lambda S^{-1}$ (completa): $e^{tA} = S\,\text{diag}(e^{\lambda_1 t}, \ldots, e^{\lambda_n t})\,S^{-1}$ (Ejercicio 10.4.25 de Olver). Para matrices incompletas, se usa la forma canónica de Jordan (Ejercicio 10.4.24).

### 11.3 Estabilidad asintótica

**[Olver, §10.2: Stability of Linear Systems]**

**Definición 10.14 (Olver).** El equilibrio $\mathbf{u}^*$ del sistema $\dot{\mathbf{u}} = \mathbf{f}(\mathbf{u})$ es:

- **Estable** si soluciones cercanas permanecen cercanas.
- **Asintóticamente estable** si es estable y además $\mathbf{u}(t) \to \mathbf{u}^*$ cuando $t \to \infty$.

**Teorema 10.16 (Olver, §10.2 — Criterio fundamental de estabilidad asintótica).** El sistema $\dot{\mathbf{u}} = A\mathbf{u}$ tiene equilibrio cero **asintóticamente estable** si y solo si:

$$
\text{Re}(\lambda_i) < 0 \quad \text{para todo valor propio } \lambda_i \text{ de } A
$$

Es decir, todos los valores propios deben estar en el **semiplano izquierdo** del plano complejo.

Si algún valor propio tiene $\text{Re}(\lambda) > 0$, el equilibrio es **inestable**.

**Teorema 10.19 (Olver).** El equilibrio es (meramente) **estable** si y solo si todos los $\text{Re}(\lambda_i) \leq 0$ y, además, todo valor propio en el eje imaginario ($\text{Re}(\lambda) = 0$) es **completo** (multiplicidad geométrica = algebraica).

### 11.4 Retratos de fase en 2D

**[Olver, §10.3: Two-Dimensional Systems]**

Para un sistema $2 \times 2$, el comportamiento cualitativo depende de la traza $\tau = \text{tr}(A)$ y el determinante $\delta = \det(A)$.

**Proposición 10.22 (Olver).** El sistema planar es: (i) asintóticamente estable sii $\delta > 0$ y $\tau < 0$; (ii) estable sii $\delta \geq 0$, $\tau \leq 0$, y si $\delta = \tau = 0$ entonces $A = O$.

Los tipos de equilibrio (nodo estable/inestable, foco, centro, silla, estrella) se clasifican según los valores propios, resumidos en la Figura 10.4 de Olver (p. 591).

### 11.5 Sistema con entrada (variación de parámetros)

**[Olver, §10.4, subsección "Inhomogeneous Linear Systems"]**

Para $\dot{\mathbf{u}} = A\mathbf{u} + \mathbf{g}(t)$, la solución es (ecuación derivada en p. 605–606):

$$
\mathbf{u}(t) = e^{(t-t_0)A}\mathbf{u}_0 + \int_{t_0}^t e^{(t-\tau)A}\mathbf{g}(\tau)\,d\tau
$$

En el contexto del control, $\mathbf{g}(t) = B\mathbf{u}(t)$, y esta fórmula es la piedra angular de la teoría de control lineal.

---

## 12. Subespacios Invariantes, Linealización y Teorema de la Variedad Central

**[Olver, §10.4, subsección "Invariant Subspaces and Linear Dynamical Systems"]**

### 12.1 Subespacios estable, central e inestable

**Definición 10.32 (Olver).** Para una matriz real $A$, se definen:

- **Subespacio estable** $S$: generado por las partes reales e imaginarias de los vectores propios/cadenas de Jordan con $\text{Re}(\lambda) < 0$.
- **Subespacio central** $C$: asociado a $\text{Re}(\lambda) = 0$.
- **Subespacio inestable** $U$: asociado a $\text{Re}(\lambda) > 0$.

Estos tres subespacios son complementarios: $S \oplus C \oplus U = \mathbb{R}^n$ (Olver, p. 604).

### 12.2 Linealización en un punto de equilibrio

**[Olver, §10.4, p. 604–605]**

Para un sistema no lineal $\dot{\mathbf{u}} = \mathbf{f}(\mathbf{u})$ con equilibrio $\mathbf{u}_0$ ($\mathbf{f}(\mathbf{u}_0) = \mathbf{0}$), Olver define la **linealización** como la matriz Jacobiana $A = \mathbf{f}'(\mathbf{u}_0) = (\partial f_i / \partial u_j)$ evaluada en el equilibrio. La expansión de Taylor de primer orden da el sistema linealizado $\dot{\boldsymbol{\delta}} \approx A\boldsymbol{\delta}$ para perturbaciones $\boldsymbol{\delta} = \mathbf{u} - \mathbf{u}_0$.

### 12.3 Teorema de la Variedad Central

**[Olver, §10.4, p. 604–605]**

Olver menciona el **Teorema de la Variedad Central** (*Center Manifold Theorem*), un resultado celebrado de la dinámica no lineal: en una vecindad de un equilibrio, el sistema no lineal admite tres variedades invariantes (estable, central e inestable) que son tangentes a los correspondientes subespacios de la linealización. En particular, las soluciones en la variedad estable convergen al equilibrio a una tasa exponencial determinada por los valores propios de $A$.

**Aplicación al péndulo invertido:** Las ecuaciones no lineales del péndulo se linealizan en el equilibrio $(\theta, \dot{\theta}, x, \dot{x}) = (0, 0, 0, 0)$ usando $\sin\theta \approx \theta$, $\cos\theta \approx 1$, produciendo la matriz $A$ del espacio de estados. La existencia de un valor propio positivo ($\lambda_1 \approx 4.539$) coloca al equilibrio en el régimen inestable, justificando la necesidad de control por retroalimentación.

---

## 13. Controlabilidad y Observabilidad

**[Ogata, Capítulos 9–10; fundamentación en álgebra lineal de Olver, Caps. 2 y 8]**

### 13.1 Controlabilidad

**Criterio de Kalman (Ogata, §9-6).** El par $(A, B)$ es controlable si y solo si la **matriz de controlabilidad**

$$
\mathcal{C} = \begin{pmatrix} B & AB & A^2B & \cdots & A^{n-1}B \end{pmatrix}
$$

tiene rango completo: $\text{rango}(\mathcal{C}) = n$.

**Interpretación en álgebra lineal:** La condición exige que las columnas de $\mathcal{C}$ generen todo $\mathbb{R}^n$ (es decir, $\text{img}(\mathcal{C}) = \mathbb{R}^n$ en la terminología de Olver). Por el Teorema de Cayley–Hamilton, potencias $A^k$ con $k \geq n$ son combinaciones lineales de $I, A, \ldots, A^{n-1}$, por lo que bastan las primeras $n$ potencias.

### 13.2 Observabilidad

**Criterio de Kalman (Ogata, §9-7).** El par $(A, C)$ es observable si y solo si $\text{rango}(\mathcal{O}) = n$, donde $\mathcal{O} = (C;\, CA;\, \ldots;\, CA^{n-1})^T$.

**Dualidad:** $(A, B)$ es controlable $\iff$ $(A^T, B^T)$ es observable.

---

## 14. Ecuación Algebraica de Riccati y Control Óptimo LQR

**[Ogata, Cap. 10; fundamentación en Olver, Caps. 3, 5, 8, 10]**

### 14.1 Formulación del problema

El regulador cuadrático lineal minimiza:

$$
J = \int_0^\infty \left(\mathbf{x}^T Q\mathbf{x} + \mathbf{u}^T R\mathbf{u}\right)dt
$$

sujeto a $\dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u}$, con $\mathbf{u} = -K\mathbf{x}$.

Las matrices de peso satisfacen: $Q \succeq 0$ (semidefinida positiva — §3.4 de Olver), $R \succ 0$ (definida positiva). El funcional $J$ representa una forma cuadrática en el sentido de la ecuación (3.52) de Olver, integrada sobre el tiempo.

### 14.2 La ecuación algebraica de Riccati (ARE)

Si $(A, B)$ es controlable y $(A, C_Q)$ es observable (donde $Q = C_Q^T C_Q$), existe una única solución $P$ simétrica definida positiva de:

$$
A^T P + PA - PBR^{-1}B^T P + Q = 0
$$

y la ganancia óptima es $K = R^{-1}B^T P$.

### 14.3 Estructura algebraica

La ARE es una **ecuación matricial cuadrática** (el término $PBR^{-1}B^T P$ es cuadrático en $P$).

Los valores propios del sistema en lazo cerrado $A - BK$ coinciden con los $n$ valores propios con parte real negativa del **Hamiltoniano**:

$$
H = \begin{pmatrix} A & -BR^{-1}B^T \\ -Q & -A^T \end{pmatrix} \in \mathbb{R}^{2n \times 2n}
$$

El Hamiltoniano tiene la propiedad de simetría espectral: si $\lambda$ es valor propio, también lo es $-\lambda$.

### 14.4 Función de Lyapunov

La función $V(\mathbf{x}) = \mathbf{x}^T P\mathbf{x}$ (con $P \succ 0$) es una **función de Lyapunov**. Dado que $P$ es definida positiva (Teorema 8.35 de Olver: todos sus valores propios son positivos):

$$
V(\mathbf{x}) > 0 \text{ para } \mathbf{x} \neq \mathbf{0}, \quad V(\mathbf{0}) = 0
$$

y $\dot{V}(\mathbf{x}) = -\mathbf{x}^T(Q + K^T RK)\mathbf{x} \leq 0$, lo que garantiza estabilidad asintótica.

**Este es el punto de convergencia de todos los conceptos:** la forma cuadrática $\mathbf{x}^T P\mathbf{x}$ requiere $P > 0$ (§3.4), la condición $\dot{V} < 0$ depende de los valores propios (§8.2, §10.2), y la existencia de $P$ requiere controlabilidad (§13) verificada vía rango (§2.5).

---

## 15. Flujos gradiente y la conexión energética

**[Olver, §10.2, ecuaciones (10.19)–(10.22)]**

Olver define los **flujos gradiente** $\dot{\mathbf{u}} = -K\mathbf{u}$ donde $K > 0$ es definida positiva. Los valores propios de $-K$ son reales y negativos (por el Teorema 8.35), garantizando estabilidad asintótica (por el Teorema 10.16). El campo vectorial $-K\mathbf{u}$ es el negativo del gradiente de la función cuadrática $q(\mathbf{u}) = \frac{1}{2}\mathbf{u}^T K\mathbf{u}$ (ecuación (10.21)), de modo que las soluciones descienden $q$ lo más rápido posible.

**Relevancia para el proyecto:** El funcional de costo LQR y la función de Lyapunov $V(\mathbf{x}) = \mathbf{x}^T P\mathbf{x}$ funcionan como una energía generalizada. La estabilización del péndulo se entiende como un proceso de minimización energética gobernado por álgebra lineal.

---

## 16. Tabla de Conexiones: Mapa Conceptual

| Concepto de Álgebra Lineal | Referencia Olver | Aplicación en el Proyecto |
| --- | --- | --- |
| Eliminación gaussiana, LU, $LDL^T$ | §1.3, §1.6 | Resolver sistemas, factorizar $Q$, $R$, $P$ |
| Espacios vectoriales, bases, dimensión | §2.1–2.4 | Estado vive en $\mathbb{R}^4$; estructura del espacio de control |
| Kernel, imagen, rango | §2.5 | Test de Kalman: $\text{rango}(\mathcal{C}) = n$ |
| Producto interno (Def. 3.1), norma | §3.1, §3.3 | Funcional de costo LQR, norma del error |
| Cauchy–Schwarz (Thm 3.5), desigualdad triangular (Thm 3.9) | §3.2 | Cotas en análisis de convergencia y robustez |
| Matrices definidas positivas (Def. 3.26, Thm 3.43) | §3.4–3.5 | $Q \succeq 0$, $R \succ 0$, $P \succ 0$ en LQR; Lyapunov $V > 0$ |
| Matrices de Gram (Def. 3.33, Thm 3.34) | §3.4 | $A^T A$ en ecuaciones normales y SVD |
| Formas cuadráticas | §3.4 | $\mathbf{x}^T Q\mathbf{x}$, $\mathbf{u}^T R\mathbf{u}$ en el funcional de costo |
| Factorización de Cholesky | §3.5 | Resolver ARE numéricamente |
| Bases ortonormales, Gram–Schmidt, QR | §4.1–4.3 | Algoritmos numéricos para eigenvalores (§9.5) |
| Proyecciones ortogonales, Fredholm | §4.4 | Complementariedad de subespacios fundamentales |
| Minimización cuadrática (Thm 5.2), mínimos cuadrados | §5.2–5.4 | Optimización LQR; identificación de parámetros |
| Transformaciones lineales, cambio de base | §7.1–7.2 | Cambio de variables de estado |
| Operadores autoadjuntos y definidos positivos | §7.5 | Generalización teórica de la simetría de $Q$, $R$, $P$ |
| Valores propios (§8.2), polinomio característico | §8.2 | Polos del sistema, estabilidad |
| Diagonalización / matrices completas | §8.3 | Cálculo de $e^{At}$, desacople modal |
| Subespacios invariantes | §8.4 | Descomposición estable/central/inestable |
| Teorema Espectral (Thm 8.32), eigenvalores simétricos | §8.5 | Diagonalización ortogonal de $Q$, $R$, $P$ |
| Positividad definida vía eigenvalores (Thm 8.35) | §8.5 | Verificación de $P > 0$ en LQR |
| Jordan canonical form | §8.6 | Soluciones generales de sistemas incompletos |
| SVD (Thm 8.63), pseudoinversa, norma espectral | §8.7 | Robustez numérica, condición de controlabilidad |
| Radio espectral (Def. 9.13), convergencia discreta (Thm 9.14) | §9.2 | Estabilidad de controladores digitales |
| Algoritmo QR, método de la potencia | §9.5 | Cálculo numérico de eigenvalores en la práctica |
| Soluciones propias, plano fase, existencia y unicidad | §10.1 | Dinámica del sistema libre |
| Estabilidad asintótica (Thm 10.16, Thm 10.19) | §10.2 | Criterio de diseño del controlador |
| Retratos de fase 2D (Prop. 10.22) | §10.3 | Visualización del comportamiento dinámico |
| Exponencial matricial, serie (ec. (10.47)) | §10.4 | Solución $\mathbf{u}(t) = e^{At}\mathbf{u}_0$ |
| Subespacios invariantes dinámicos, variedad central | §10.4 | Linealización, justificación de $\sin\theta \approx \theta$ |
| Flujo gradiente | §10.2 | Interpretación energética del LQR |
| Dinámica de estructuras, modos de vibración | §10.5 | Conexión con modos propios del péndulo |

---

## 17. Referencias

1. **Olver, P. J. & Shakiban, C.** (2018). *Applied Linear Algebra*, 2nd Ed. Springer.
   - Capítulo 1 (§1.1–1.9): Eliminación, factorizaciones, inversas, determinantes.
   - Capítulo 2 (§2.1–2.5): Espacios vectoriales, subespacios fundamentales.
   - Capítulo 3 (§3.1–3.5): Producto interno (§3.1), Cauchy–Schwarz (§3.2), normas (§3.3), **matrices definidas positivas** (§3.4), completar el cuadrado y Cholesky (§3.5).
   - Capítulo 4 (§4.1–4.4): Ortogonalidad, Gram–Schmidt, QR, proyecciones.
   - Capítulo 5 (§5.1–5.4): Minimización cuadrática, mínimos cuadrados.
   - Capítulo 7 (§7.1–7.5): Transformaciones lineales, cambio de base, adjuntos.
   - Capítulo 8 (§8.1–8.8): Valores propios (§8.2), diagonalización (§8.3), **Teorema Espectral** (§8.5, Thm 8.32), **positividad definida vía eigenvalores** (§8.5, Thm 8.35), SVD (§8.7, Thm 8.63).
   - Capítulo 9 (§9.1–9.2): Iteración, **radio espectral** (Def. 9.13, Thm 9.14).
   - Capítulo 10 (§10.1–10.5): Dinámica, **estabilidad** (Thm 10.16), retratos de fase (§10.3), **exponencial matricial** (§10.4), subespacios invariantes y variedad central (§10.4), sistemas forzados (§10.4).

2. **Ogata, K.** (2010). *Ingeniería de Control Moderna*, 5ª Ed. Pearson.
   - Capítulo 2: Espacio de estados, funciones de transferencia.
   - Capítulo 3: Modelado del péndulo invertido.
   - Capítulo 9: Controlabilidad y observabilidad.
   - Capítulo 10: Diseño LQR, ecuación de Riccati.
