### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 99b4c781-70d4-4e78-8512-e7f6b4fc6cf4
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	
	using LinearAlgebra
	using Printf
	using DifferentialEquations
	using CairoMakie
	
	# PlutoUI nos da sliders, checkboxes, etc.
	using PlutoUI
end

# ╔═╡ db1608db-421a-4171-97e6-7062a27ebf3b
begin
	include(joinpath(@__DIR__, "..", "src", "model.jl"))
	include(joinpath(@__DIR__, "..", "src", "linearization.jl"))
	include(joinpath(@__DIR__, "..", "src", "controller.jl"))
	
	using .Model
	using .Linearization
	using .Controller
	
	md"Módulos cargados correctamente"
end

# ╔═╡ be61f640-b844-420e-9f9d-af6ecaa5203a
md"""
# Péndulo Invertido — Exploración Interactiva

**Proyecto de Álgebra Lineal Aplicada**

Este notebook te permite explorar el sistema carro-péndulo de forma interactiva.  
Cada vez que modificas un parámetro (un slider, un valor), **todas las celdas dependientes se recalculan automáticamente**.

---
"""

# ╔═╡ 0955a7cb-8aac-470d-a4fe-9b7050342699
md"""
## 1. Módulos del proyecto

Cargamos los tres módulos que definimos en `src/`:
- **Model** → parámetros físicos y ecuaciones de movimiento no lineales
- **Linearization** → Jacobiano, matrices A/B/C/D, eigenvalores, Kalman
- **Controller** → LQR (Riccati), pole placement (Ackermann)
"""

# ╔═╡ 4cd10adc-e79b-4dad-95e2-9c2cc7a40833
md"""
## 2. Parámetros del sistema

Usa los sliders para modificar los parámetros físicos y observar cómo cambia **todo** el análisis en tiempo real.

| Parámetro | Símbolo | Descripción |
|-----------|---------|-------------|
| Masa del carro | ``M`` | Masa de la base móvil |
| Masa del péndulo | ``m`` | Masa de la barra |
| Longitud total | ``L_{bar}`` | Longitud de la barra completa |
| Gravedad | ``g`` | Aceleración gravitacional |
| Fricción | ``b`` | Coeficiente de fricción viscosa del carro |
"""

# ╔═╡ 54b11144-5dfd-4561-8580-4800c147d52e
md"**Masa del carro M [kg]:**"

# ╔═╡ 0565789e-462d-49a8-8534-413709c5b091
@bind M_val PlutoUI.Slider(0.2:0.1:5.0, default=1.0, show_value=true)

# ╔═╡ 93b07fd7-ef7a-43d1-a0fd-51bdd218d055
md"**Masa del péndulo m [kg]:**"

# ╔═╡ 12d58a06-f53a-4de0-a6fa-e1f1f14edf46
@bind m_val PlutoUI.Slider(0.05:0.05:2.0, default=0.3, show_value=true)

# ╔═╡ ea2d982a-30dd-4d28-9665-2f43d851613b
md"**Longitud total de la barra Lbar [m]:**"

# ╔═╡ 0a37ad6f-e44b-4706-bb70-48401d650700
@bind Lbar PlutoUI.Slider(0.2:0.1:3.0, default=1.0, show_value=true)

# ╔═╡ 5138621e-392c-4218-9f5e-9b59b0707cbb
md"**Gravedad g [m/s²]:**"

# ╔═╡ dd76336c-5655-4e1b-9302-525ec8640eeb
@bind g_val PlutoUI.Slider(1.0:0.1:20.0, default=9.81, show_value=true)

# ╔═╡ dc717647-fd75-4c58-97b3-8a423b78dfcf
md"**Fricción b [N·s/m]:**"

# ╔═╡ c691cfef-181c-4c2b-bc2a-fbdda4623861
@bind b_val PlutoUI.Slider(0.0:0.05:2.0, default=0.1, show_value=true)

# ╔═╡ 71f6c6e0-a119-43f6-9a4e-3a7ad1bcf0b0
begin
	# Calcular parámetros derivados
	L_val = Lbar / 2    # distancia al centro de masa
	I_val = (1/12) * m_val * Lbar^2   # inercia de barra uniforme
	
	params = SystemParams(M=M_val, m=m_val, L=L_val, g=g_val, b=b_val, I=I_val)
	
	md"""
	### Parámetros derivados
	
	- Distancia al centro de masa: **L = Lbar/2 = $(round(L_val, digits=3)) m**
	- Momento de inercia (barra uniforme): **I = (1/12)·m·Lbar² = $(round(I_val, digits=5)) kg·m²**
	"""
end

# ╔═╡ 4d23db82-2560-4a28-b84d-fe6395055ae5
md"""
## 3. Linealización — Matrices del espacio de estados

Al linealizar las ecuaciones no lineales de Euler-Lagrange alrededor del equilibrio superior (``\theta = 0``), obtenemos el sistema:

$$\dot{\mathbf{x}} = A\mathbf{x} + B\mathbf{u}, \quad \mathbf{y} = C\mathbf{x} + D\mathbf{u}$$

donde el estado es ``\mathbf{x} = [x,\; \dot{x},\; \theta,\; \dot{\theta}]^T``.
"""

# ╔═╡ 4250622d-0830-4333-9f16-d7198160c5a3
begin
	ss = linearize_system(params)
	md"Sistema linealizado"
end

# ╔═╡ b1e47212-c25b-4d83-934f-0fee7f773174
md"### Matriz A (dinámica del sistema, 4x4):"

# ╔═╡ 9d0a35e5-5768-4d2f-a5aa-1814b53e8721
round.(ss.A, digits=4)

# ╔═╡ 9a603b09-b6ba-4e73-8d51-695015cae9d8
md"### Matriz B (entrada de control, 4x1):"

# ╔═╡ 375e2a1f-27d8-466f-a343-d81c95ce2a9a
round.(ss.B, digits=4)

# ╔═╡ 4426d097-d814-40ee-aa8b-f75d8ee0789d
md"### Matriz C (salidas medibles, 2x4):"

# ╔═╡ 18152fa9-97fe-4400-8bd5-8eec10276d92
ss.C

# ╔═╡ 329d269b-1534-4460-86d8-1e9eebb9bc16
md"""
## 4. Análisis de estabilidad — Eigenvalores de A

Los eigenvalores de ``A`` determinan la estabilidad del sistema en lazo abierto.  
Si **algún** eigenvalor tiene parte real positiva → el sistema es **inestable** (necesita control).

Recordar (Olver, Teorema 10.16): *estabilidad asintótica si y solo si* ``\text{Re}(\lambda_i) < 0 \;\forall\, i``.
"""

# ╔═╡ 027d6af7-96a0-4f7c-8d57-dafd92fbb930
begin
	eig_info = map(enumerate(ss.eigenvalues)) do (i, λ)
		rp = round(real(λ), digits=4)
		ip = round(imag(λ), digits=4)
		estab = rp > 1e-10 ? "INESTABLE" : (rp < -1e-10 ? "Estable" : "Marginal")
		val_str = abs(ip) < 1e-10 ? "$rp" : "$rp + $(ip)i"
		"| λ_$i | $val_str | $estab |"
	end
	
	any_unstable = any(real.(ss.eigenvalues) .> 1e-10)
	status = any_unstable ? "**SISTEMA INESTABLE** — requiere control activo" : "**Sistema estable**"
	
	Markdown.parse("""
| Eigenvalor | Valor | Estabilidad |
|------------|-------|-------------|
$(join(eig_info, "\n"))

**Diagnóstico:** $status
""")
end

# ╔═╡ cca01fd6-94af-48cc-ae8d-84a8ca36551c
md"""
## 5. Controlabilidad y Observabilidad — Criterio de Kalman

**Controlabilidad** (Ogata §9-6): ``\text{rank}[\mathbf{B} \;\; A\mathbf{B} \;\; A^2\mathbf{B} \;\; A^3\mathbf{B}] = n``

**Observabilidad** (Ogata §9-7): ``\text{rank}[\mathbf{C}^T \;\; A^T\mathbf{C}^T \;\; (A^T)^2\mathbf{C}^T \;\; (A^T)^3\mathbf{C}^T] = n``

Si ambos rangos son ``n = 4``, podemos diseñar un controlador y un observador.
"""

# ╔═╡ 86fd0c62-d8f6-4811-b7d3-b52ba62eb769
begin
	ctrl = check_controllability(ss)
	obs = check_observability(ss)
	
	ctrl_status = ctrl.is_controllable ? "CONTROLABLE" : "NO CONTROLABLE"
	obs_status = obs.is_observable ? "OBSERVABLE" : "NO OBSERVABLE"
	
	Markdown.parse("""
| Propiedad | Rango | Requerido | Resultado |
|-----------|-------|-----------|-----------|
| Controlabilidad | $(ctrl.rank) | $(ctrl.required_rank) | $ctrl_status |
| Observabilidad | $(obs.rank) | $(obs.required_rank) | $obs_status |
""")
end

# ╔═╡ d26e24fe-cfa6-4c42-80d7-27175c4e6e43
md"### Matriz de controlabilidad ``\\mathcal{C} = [B \\;\\; AB \\;\\; A^2B \\;\\; A^3B]``:"

# ╔═╡ aea2bad2-82fa-4dc6-879e-1edbaa4f7e4d
round.(ctrl.matrix, digits=4)

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000001
md"""
## 6. Diseño del controlador LQR

El LQR minimiza el funcional de costo:

$$J = \int_0^\infty \left(\mathbf{x}^T Q \mathbf{x} + \mathbf{u}^T R \mathbf{u}\right) dt$$

- ``Q`` **penaliza las desviaciones del estado** (¿qué tan lejos del equilibrio?)
- ``R`` **penaliza el esfuerzo de control** (¿cuánta fuerza usa el actuador?)

Ajusta los pesos y observa cómo cambia la ganancia ``K`` y los polos en lazo cerrado:
"""

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000002
md"**Peso en posición del carro (Q₁₁):**"

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000003
@bind q_x PlutoUI.Slider(0.1:0.1:50.0, default=1.0, show_value=true)

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000004
md"**Peso en ángulo del péndulo (Q₃₃):**"

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000005
@bind q_theta PlutoUI.Slider(0.1:0.5:100.0, default=10.0, show_value=true)

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000006
md"**Peso en esfuerzo de control (R):**"

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000007
@bind r_val PlutoUI.Slider(0.01:0.01:5.0, default=0.1, show_value=true)

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000008
begin
	Q_mat = diagm([q_x, 0.0, q_theta, 0.0])
	R_mat = reshape([r_val], 1, 1)
	
	lqr = design_lqr(ss.A, ss.B, Q_mat, R_mat)
	
	labels_K = ["x", "ẋ", "θ", "θ̇"]
	k_info = ["|  K_$(labels_K[i])  |  $(round(lqr.K[i], digits=4))  |" for i in 1:4]
	
	cl_info = map(enumerate(lqr.eigenvalues_cl)) do (i, λ)
		rp = round(real(λ), digits=4)
		ip = round(imag(λ), digits=4)
		val_str = abs(ip) < 1e-10 ? "$rp" : "$rp + $(ip)i"
		"| λ_$i | $val_str | |"
	end
	
	all_stable = all(real.(lqr.eigenvalues_cl) .< 0)
	cl_status = all_stable ? "**TODOS LOS POLOS ESTABLES** — el control funciona" : "**Sistema aún inestable**"
	
	Markdown.parse("""
### Ganancia de retroalimentación ``K`` (ley de control: ``u = -K\\mathbf{x}``)

| Componente | Valor |
|------------|-------|
$(join(k_info, "\n"))

### Eigenvalores en lazo cerrado ``(A - BK)``

| Eigenvalor | Valor | Estado |
|------------|-------|--------|
$(join(cl_info, "\n"))

**Diagnóstico:** $cl_status
	""")
end

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000009
md"### Solución de la ecuación de Riccati P (norma):"

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000010
round(norm(lqr.P), digits=4)

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000011
md"""
## 7. Simulación comparativa

Simulamos el sistema con una condición inicial ``\theta_0`` y comparamos la respuesta **sin control** (caída libre) vs. **con LQR** (estabilización).
"""

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000012
md"**Ángulo inicial θ₀ [rad]:**"

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000013
@bind θ₀ PlutoUI.Slider(0.01:0.01:0.5, default=0.15, show_value=true)

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000014
begin
	x0 = [0.0, 0.0, θ₀, 0.0]
	
	# ── Simulación sin control ──
	prob_free = ODEProblem(nonlinear_eom!, x0, (0.0, 3.0),
		(params=params, F=0.0))
	sol_free = solve(prob_free, Tsit5(), saveat=0.02)
	
	# ── Simulación con LQR ──
	prob_ctrl = ODEProblem(closed_loop_eom!, x0, (0.0, 10.0),
		(params=params, K=lqr.K, saturate=50.0))
	sol_ctrl = solve(prob_ctrl, Tsit5(), saveat=0.02)
	
	θ_final_free = round(sol_free.u[end][3], digits=3)
	θ_final_ctrl = round(sol_ctrl.u[end][3], digits=6)
	
	Markdown.parse("""
### Resultados de la simulación

| Escenario | θ₀ [rad] | θ final [rad] | Tiempo [s] |
|-----------|----------|---------------|------------|
| Sin control | $θ₀ | $θ_final_free | 3.0 |
| Con LQR | $θ₀ | $θ_final_ctrl | 10.0 |

**Interpretación:** Sin control, el péndulo cae (θ diverge). Con LQR, el ángulo converge a 0 (péndulo estabilizado en la vertical).
	""")
end

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000001
md"""
## 8. Visualización gráfica — Evolución temporal del sistema

Las siguientes gráficas muestran la evolución de las **cuatro variables de estado** ``\mathbf{x}(t) = [x,\, \dot{x},\, \theta,\, \dot{\theta}]^T`` a lo largo del tiempo, comparando el comportamiento **sin control** (caída libre, en rojo) contra la **estabilización con LQR** (en azul).

Estas curvas son la respuesta del sistema al resolver numéricamente las EDOs no lineales con `DifferentialEquations.jl`, renderizadas con `CairoMakie`.
"""

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000002
begin
	# ─────────────────────────────────────────────────────────
	# Tema visual global para todas las gráficas del notebook
	# ─────────────────────────────────────────────────────────
	const COL_CTRL   = RGBf(0.173, 0.447, 0.710)    # azul intenso
	const COL_LIBRE  = RGBf(0.839, 0.153, 0.157)    # rojo
	const COL_FUERZA = RGBf(0.122, 0.467, 0.706)    # azul medio  
	const COL_REGION = RGBAf(0.173, 0.447, 0.710, 0.08)
	const COL_GRID   = RGBAf(0.0, 0.0, 0.0, 0.08)
	
	const TEMA_PENDULO = Theme(
		fontsize = 13,
		Axis = (
			xgridcolor      = COL_GRID,
			ygridcolor      = COL_GRID,
			xgridwidth       = 0.5,
			ygridwidth       = 0.5,
			spinewidth       = 0.8,
			xticksize        = 4,
			yticksize        = 4,
			xlabelpadding    = 4,
			ylabelpadding    = 4,
			titlefont        = :bold,
			titlegap         = 8,
			xlabelfont       = :regular,
			ylabelfont       = :regular,
		),
		Lines = (linewidth = 2.0,),
		Legend = (
			framevisible = false,
			padding      = (6, 6, 4, 4),
			labelsize    = 11,
			patchsize    = (20, 10),
		),
	)
	
	md"Tema visual configurado (CairoMakie)"
end

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000003
begin
	# ────────────────────────────────────────────
	# Extraer trayectorias para graficar
	# ────────────────────────────────────────────
	t_free = sol_free.t
	x_free  = [u[1] for u in sol_free.u]
	v_free  = [u[2] for u in sol_free.u]
	θ_free  = [u[3] for u in sol_free.u]
	ω_free  = [u[4] for u in sol_free.u]
	
	t_ctrl = sol_ctrl.t
	x_ctrl  = [u[1] for u in sol_ctrl.u]
	v_ctrl  = [u[2] for u in sol_ctrl.u]
	θ_ctrl  = [u[3] for u in sol_ctrl.u]
	ω_ctrl  = [u[4] for u in sol_ctrl.u]
	
	# Fuerza de control F(t)
	F_ctrl_vec = map(sol_ctrl.u) do s
		F_raw = -dot(lqr.K[1,:], s)
		clamp(F_raw, -50.0, 50.0)
	end
	
	md"Trayectorias extraídas: $(length(t_free)) puntos (libre), $(length(t_ctrl)) puntos (controlada)"
end

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000010
md"""
### 8.1 Variables de estado — Comparación libre vs. controlado

Panel de cuatro gráficas mostrando cada componente del vector de estado ``\mathbf{x}(t)``.  
La respuesta libre (rojo) solo cubre 3 segundos porque el péndulo cae rápidamente; la respuesta controlada (azul) se extiende 10 segundos para apreciar la convergencia.
"""

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000004
begin
	with_theme(TEMA_PENDULO) do

	fig_states = Figure(size = (820, 700))
	
	# ── Eje 1: Posición del carro x(t) ──
	ax1 = Axis(fig_states[1, 1],
		ylabel = "x  [m]",
		title  = "Posición del carro",
		xticklabelsvisible = false)
	lines!(ax1, t_ctrl, x_ctrl, color = COL_CTRL, label = "Con LQR")
	lines!(ax1, t_free, x_free, color = COL_LIBRE, linestyle = :dash, label = "Sin control")
	hlines!(ax1, [0.0], color = :gray60, linewidth = 0.6, linestyle = :dot)
	axislegend(ax1, position = :rt)
	
	# ── Eje 2: Velocidad del carro ẋ(t) ──
	ax2 = Axis(fig_states[1, 2],
		ylabel = "ẋ  [m/s]",
		title  = "Velocidad del carro",
		xticklabelsvisible = false)
	lines!(ax2, t_ctrl, v_ctrl, color = COL_CTRL)
	lines!(ax2, t_free, v_free, color = COL_LIBRE, linestyle = :dash)
	hlines!(ax2, [0.0], color = :gray60, linewidth = 0.6, linestyle = :dot)
	
	# ── Eje 3: Ángulo del péndulo θ(t) ──
	ax3 = Axis(fig_states[2, 1],
		xlabel = "t  [s]",
		ylabel = "θ  [rad]",
		title  = "Ángulo del péndulo")
	lines!(ax3, t_ctrl, θ_ctrl, color = COL_CTRL)
	lines!(ax3, t_free, θ_free, color = COL_LIBRE, linestyle = :dash)
	hlines!(ax3, [0.0], color = :gray60, linewidth = 0.6, linestyle = :dot)
	# Banda visual: zona de linealización válida (|θ| < ~0.25 rad)
	band!(ax3, [t_ctrl[1], t_ctrl[end]], [-0.25, -0.25], [0.25, 0.25],
		  color = COL_REGION)
	
	# ── Eje 4: Velocidad angular θ̇(t) ──
	ax4 = Axis(fig_states[2, 2],
		xlabel = "t  [s]",
		ylabel = "θ̇  [rad/s]",
		title  = "Velocidad angular")
	lines!(ax4, t_ctrl, ω_ctrl, color = COL_CTRL)
	lines!(ax4, t_free, ω_free, color = COL_LIBRE, linestyle = :dash)
	hlines!(ax4, [0.0], color = :gray60, linewidth = 0.6, linestyle = :dot)
	
	# Título superior
	Label(fig_states[0, :],
		"Evolución temporal del estado  x(t) = [x, ẋ, θ, θ̇]ᵀ",
		fontsize = 15, font = :bold, padding = (0, 0, 0, 4))
	
	rowgap!(fig_states.layout, 10)
	colgap!(fig_states.layout, 14)
	
	fig_states
	end
end

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000011
md"""
### 8.2 Señal de control — Fuerza del actuador ``F(t)``

La fuerza ``F(t) = -K\mathbf{x}(t)`` es la acción que el controlador LQR aplica al carro.  
Observa el pico inicial (máximo esfuerzo para detener la caída) y la posterior convergencia a cero.  
La línea punteada marca el **límite de saturación** del actuador (±50 N).
"""

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000005
begin
	with_theme(TEMA_PENDULO) do

	fig_force = Figure(size = (720, 320))
	
	ax_f = Axis(fig_force[1, 1],
		xlabel = "t  [s]",
		ylabel = "F  [N]",
		title  = "Señal de control  u(t) = −Kx(t)")
	
	# Banda de saturación
	band!(ax_f, [t_ctrl[1], t_ctrl[end]], [-50.0, -50.0], [50.0, 50.0],
		  color = RGBAf(0.95, 0.85, 0.75, 0.25))
	hlines!(ax_f, [50.0, -50.0], color = :orange, linewidth = 1.0, linestyle = :dash,
			label = "Saturación ±50 N")
	
	# Curva de fuerza
	lines!(ax_f, t_ctrl, F_ctrl_vec, color = COL_FUERZA, linewidth = 2.2,
		   label = "F(t)")
	hlines!(ax_f, [0.0], color = :gray60, linewidth = 0.6, linestyle = :dot)
	
	# Pico máximo
	F_max = maximum(abs.(F_ctrl_vec))
	idx_max = argmax(abs.(F_ctrl_vec))
	scatter!(ax_f, [t_ctrl[idx_max]], [F_ctrl_vec[idx_max]],
			 color = COL_FUERZA, markersize = 8, strokecolor = :white, strokewidth = 1)
	text!(ax_f, t_ctrl[idx_max] + 0.15, F_ctrl_vec[idx_max],
		  text = @sprintf("|F|_max = %.1f N", F_max),
		  fontsize = 10, color = :gray30)
	
	axislegend(ax_f, position = :rt)
	
	fig_force
	end
end

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000012
md"""
### 8.3 Retrato de fase — Plano ``(\theta,\; \dot{\theta})``

El retrato de fase muestra la **trayectoria del sistema en el espacio de estados** (proyectado en las coordenadas angulares). El origen ``(0, 0)`` es el punto de equilibrio inestable (péndulo vertical hacia arriba).

- **Rojo** → trayectoria libre: se aleja del origen (inestable)  
- **Azul** → trayectoria controlada: espiral convergente hacia el origen (estabilizado)
"""

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000006
begin
	with_theme(TEMA_PENDULO) do

	fig_phase = Figure(size = (560, 500))
	
	ax_p = Axis(fig_phase[1, 1],
		xlabel = "θ  [rad]",
		ylabel = "θ̇  [rad/s]",
		title  = "Retrato de fase  (θ, θ̇)",
		aspect = DataAspect())
	
	# Trayectoria libre
	lines!(ax_p, θ_free, ω_free, color = COL_LIBRE,
		   linewidth = 1.8, linestyle = :dash, label = "Sin control")
	scatter!(ax_p, [θ_free[1]], [ω_free[1]],
			 color = COL_LIBRE, markersize = 10,
			 marker = :circle, strokecolor = :white, strokewidth = 1)
	
	# Trayectoria controlada
	lines!(ax_p, θ_ctrl, ω_ctrl, color = COL_CTRL,
		   linewidth = 1.8, label = "Con LQR")
	scatter!(ax_p, [θ_ctrl[1]], [ω_ctrl[1]],
			 color = COL_CTRL, markersize = 10,
			 marker = :circle, strokecolor = :white, strokewidth = 1)
	
	# Punto de equilibrio
	scatter!(ax_p, [0.0], [0.0], color = :black, markersize = 10,
			 marker = :xcross, label = "Equilibrio (0, 0)")
	
	axislegend(ax_p, position = :lt)
	
	fig_phase
	end
end

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000013
md"""
### 8.4 Energía del sistema — Cinética, potencial y total

La energía total ``E(t) = T + V`` permite verificar la **coherencia física** de la simulación. Sin control, la energía crece (el péndulo gana energía cinética al caer). Con control, el actuador extrae energía del sistema hasta que ``E \to E_{\text{eq}}``.

- ``T = \frac{1}{2}(M + m)\dot{x}^2 + \frac{1}{2}(I + mL^2)\dot{\theta}^2 + mL\dot{x}\dot{\theta}\cos\theta``
- ``V = mgL\cos\theta``  (con ``V_{\max} = mgL`` en la vertical)
"""

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000007
begin
	# ── Función de energía para el modelo con inercia ──
	function calc_energy(sol, sp)
		M_e, m_e, L_e, g_e, I_e = sp.M, sp.m, sp.L, sp.g, sp.I
		map(sol.u) do s
			_, v, θ, ω = s
			T = 0.5*(M_e + m_e)*v^2 + 0.5*(I_e + m_e*L_e^2)*ω^2 + m_e*L_e*v*ω*cos(θ)
			V = m_e*g_e*L_e*cos(θ)
			(T = T, V = V, E = T + V)
		end
	end
	
	E_ctrl = calc_energy(sol_ctrl, params)
	E_free = calc_energy(sol_free, params)
	
	with_theme(TEMA_PENDULO) do

	fig_energy = Figure(size = (820, 340))
	
	# Panel izquierdo: libre
	ax_e1 = Axis(fig_energy[1, 1],
		xlabel = "t  [s]", ylabel = "Energía  [J]",
		title  = "Sin control (caída libre)")
	
	T_f = [e.T for e in E_free]
	V_f = [e.V for e in E_free]
	E_f = [e.E for e in E_free]
	lines!(ax_e1, t_free, T_f, color = RGBf(0.90, 0.55, 0.13), label = "Cinética T")
	lines!(ax_e1, t_free, V_f, color = RGBf(0.20, 0.63, 0.17), label = "Potencial V")
	lines!(ax_e1, t_free, E_f, color = :gray25, linewidth = 2.5,
		   linestyle = :dash, label = "Total E")
	axislegend(ax_e1, position = :lt)
	
	# Panel derecho: controlado
	ax_e2 = Axis(fig_energy[1, 2],
		xlabel = "t  [s]", ylabel = "Energía  [J]",
		title  = "Con LQR (estabilización)")
	
	T_c = [e.T for e in E_ctrl]
	V_c = [e.V for e in E_ctrl]
	E_c = [e.E for e in E_ctrl]
	lines!(ax_e2, t_ctrl, T_c, color = RGBf(0.90, 0.55, 0.13), label = "Cinética T")
	lines!(ax_e2, t_ctrl, V_c, color = RGBf(0.20, 0.63, 0.17), label = "Potencial V")
	lines!(ax_e2, t_ctrl, E_c, color = :gray25, linewidth = 2.5,
		   linestyle = :dash, label = "Total E")
	axislegend(ax_e2, position = :rt)
	
	Label(fig_energy[0, :],
		"Balance energético del sistema",
		fontsize = 14, font = :bold, padding = (0, 0, 0, 4))
	
	colgap!(fig_energy.layout, 18)
	
	fig_energy
	end
end

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000014
md"""
### 8.5 Eigenvalores en el plano complejo — Lazo abierto vs. cerrado

Visualización de la **ubicación espectral** de los eigenvalores de ``A`` (lazo abierto, sistema inestable) y de ``A - BK`` (lazo cerrado, estabilizado por LQR). Todo eigenvalor a la izquierda del eje imaginario contribuye a la estabilidad asintótica (Olver, Thm. 10.16).
"""

# ╔═╡ c3d4e5f6-2222-4000-a000-000000000008
begin
	with_theme(TEMA_PENDULO) do

	fig_eig = Figure(size = (620, 400))
	
	ax_eig = Axis(fig_eig[1, 1],
		xlabel = "Re(λ)",
		ylabel = "Im(λ)",
		title  = "Eigenvalores — Lazo abierto vs. cerrado")
	
	# Eje imaginario (frontera de estabilidad)
	vlines!(ax_eig, [0.0], color = :gray50, linewidth = 1.2, linestyle = :dash)
	hlines!(ax_eig, [0.0], color = :gray80, linewidth = 0.5)
	
	# Sombrear semiplano izquierdo (zona estable)
	re_lims = let
		all_re = vcat(real.(ss.eigenvalues), real.(lqr.eigenvalues_cl))
		margin = max(abs(minimum(all_re)), abs(maximum(all_re))) * 0.4
		(minimum(all_re) - margin, maximum(all_re) + margin)
	end
	im_lims = let
		all_im = vcat(imag.(ss.eigenvalues), imag.(lqr.eigenvalues_cl))
		margin = max(maximum(abs.(all_im)) * 0.5, 0.5)
		(-margin - 0.3, margin + 0.3)
	end
	band!(ax_eig, [re_lims[1], 0.0], [im_lims[1], im_lims[1]],
		  [im_lims[2], im_lims[2]], color = RGBAf(0.2, 0.7, 0.3, 0.06))
	text!(ax_eig, re_lims[1] * 0.5, im_lims[2] * 0.85,
		  text = "Re(λ) < 0\nEstable", fontsize = 10, color = :gray50, align = (:center, :top))
	
	# Eigenvalores lazo abierto
	scatter!(ax_eig,
		real.(ss.eigenvalues), imag.(ss.eigenvalues),
		color = COL_LIBRE, markersize = 14, marker = :circle,
		strokecolor = :white, strokewidth = 1.5,
		label = "Lazo abierto (A)")
	
	# Eigenvalores lazo cerrado
	scatter!(ax_eig,
		real.(lqr.eigenvalues_cl), imag.(lqr.eigenvalues_cl),
		color = COL_CTRL, markersize = 14, marker = :utriangle,
		strokecolor = :white, strokewidth = 1.5,
		label = "Lazo cerrado (A − BK)")
	
	axislegend(ax_eig, position = :lb)
	
	fig_eig
	end
end

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000015
md"""
## 9. Datos de la trayectoria controlada

Aquí puedes inspeccionar los valores numéricos de la simulación controlada:
"""

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000016
begin
	# Tabla resumida de la trayectoria controlada (cada 0.5 s)
	t_sample = 0.0:0.5:10.0
	traj_rows = map(t_sample) do t
		s = sol_ctrl(t)
		F_ctrl_t = -dot(lqr.K[1,:], s)
		F_sat = clamp(F_ctrl_t, -50.0, 50.0)
		"| $(round(t, digits=1)) | $(round(s[1], digits=4)) | $(round(s[2], digits=4)) | $(round(rad2deg(s[3]), digits=3)) | $(round(s[4], digits=4)) | $(round(F_sat, digits=2)) |"
	end
	
	Markdown.parse("""
| t [s] | x [m] | ẋ [m/s] | θ [°] | θ̇ [rad/s] | F [N] |
|-------|-------|---------|-------|-----------|-------|
$(join(traj_rows, "\n"))
	""")
end

# ╔═╡ a1b2c3d4-1111-4000-8000-000000000017
md"""
---

## Ejercicios sugeridos

1. **Efecto de la masa:** Aumenta `m` (masa del péndulo) a 1.0 kg. ¿Qué pasa con el eigenvalor inestable? ¿Crece o decrece?

2. **Efecto de la longitud:** Prueba `Lbar = 0.3 m` (péndulo corto) vs. `Lbar = 2.0 m` (péndulo largo). ¿Cuál es más difícil de controlar?

3. **Sintonización del LQR:** Con los valores por defecto, sube `Q₃₃` (peso en θ) a 50. ¿Qué le pasa a la ganancia `K_θ`? ¿Y si subes `R` a 2.0?

4. **Ángulo inicial grande:** Sube `θ₀` a 0.4 rad (~23°). ¿El LQR aún logra estabilizar? Recuerda que el controlador es **lineal** pero el modelo es **no lineal** — ¿cuándo falla la aproximación?

5. **Gravedad lunar:** Pon `g = 1.62`. ¿El péndulo es más fácil o más difícil de controlar en la Luna?
"""

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000001
md"""
---

## 10. Animación del sistema carro-péndulo

Esta sección implementa la **visualización mecánica** del sistema: ya no solo gráficas de variables de estado, sino el carro moviéndose sobre un riel y la barra del péndulo oscilando.

El pipeline es:
1. **Un frame estático** — función `draw_pendulum_frame!` que dibuja el sistema dado un estado `[x, ẋ, θ, θ̇]`
2. **Slider de tiempo** — recorre la trayectoria pre-calculada frame a frame
3. **Comparación lado a lado** — caída libre vs. LQR estabilizado
4. **Exportación a GIF** — `CairoMakie.record()` genera un archivo animado

Todo funciona sobre **CairoMakie** (CPU), sin necesidad de GPU ni ventanas externas.
"""

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000002
begin
	"""
	    draw_pendulum_frame!(ax, x_cart, θ, L_draw; kwargs...)
	
	Dibuja un frame completo del sistema carro-péndulo invertido.
	
	Componentes visuales:
	  - Riel (línea punteada gris)
	  - Carro (rectángulo con ruedas)
	  - Barra rígida (línea gruesa desde pivote hasta masa)
	  - Masa puntual (círculo en el extremo)
	  - Pivote (punto en el centro del carro)
	
	Convención geométrica:
	  - θ = 0 → péndulo vertical hacia ARRIBA (equilibrio inestable)
	  - x_tip = x_cart + L·sin(θ)
	  - y_tip = L·cos(θ)
	"""
	function draw_pendulum_frame!(ax, x_cart, θ, L_draw;
		cart_color  = COL_CTRL,
		bob_color   = RGBf(0.85, 0.20, 0.20),
		rail_span   = 3.0,
		cart_w      = 0.40,
		cart_h      = 0.20,
		wheel_r     = 0.05,
		bob_size    = 18,
		rod_width   = 3.5,
		show_angle  = false,
		time_label  = nothing)

		empty!(ax)   # limpiar el eje para redibujar desde cero

		# ── Riel / suelo ──
		lines!(ax, [-rail_span, rail_span], [0.0, 0.0],
			color = :gray70, linewidth = 1.5, linestyle = :dash)

		# ── Sombreado del suelo ──
		band!(ax, [-rail_span, rail_span], [-0.35, -0.35], [0.0, 0.0],
			color = RGBAf(0.0, 0.0, 0.0, 0.03))

		# ── Carro (cuerpo) ──
		poly!(ax,
			Rect(x_cart - cart_w/2, -cart_h/2, cart_w, cart_h),
			color      = (cart_color, 0.85),
			strokecolor = :gray20,
			strokewidth = 1.5)

		# ── Ruedas ──
		for dx in [-0.12, 0.12]
			scatter!(ax, [x_cart + dx], [-cart_h/2 - 0.02],
				markersize  = 10,
				color       = :gray45,
				strokecolor = :gray20,
				strokewidth = 0.8)
		end

		# ── Coordenadas del extremo del péndulo ──
		x_tip = x_cart + L_draw * sin(θ)
		y_tip = L_draw * cos(θ)

		# ── Barra rígida (rod) ──
		lines!(ax, [x_cart, x_tip], [0.0, y_tip],
			color     = :gray25,
			linewidth = rod_width)

		# ── Masa del péndulo (bob) ──
		scatter!(ax, [x_tip], [y_tip],
			markersize  = bob_size,
			color       = bob_color,
			strokecolor = :gray20,
			strokewidth = 1.5)

		# ── Pivote ──
		scatter!(ax, [x_cart], [0.0],
			markersize = 6,
			color      = :gray30)

		# ── Indicador de ángulo (arco) ──
		if show_angle && abs(θ) > 0.01
			n_arc = 30
			r_arc = L_draw * 0.25
			θ_range = range(0, θ, length = n_arc)
			arc_x = x_cart .+ r_arc .* sin.(θ_range)
			arc_y = r_arc .* cos.(θ_range)
			lines!(ax, arc_x, arc_y,
				color = RGBAf(0.8, 0.3, 0.1, 0.6), linewidth = 1.5)
			# Etiqueta del ángulo
			mid_θ = θ / 2
			text!(ax, x_cart + r_arc * 1.4 * sin(mid_θ), r_arc * 1.4 * cos(mid_θ),
				text     = @sprintf("%.1f°", rad2deg(θ)),
				fontsize = 10,
				color    = RGBf(0.8, 0.3, 0.1),
				align    = (:center, :center))
		end

		# ── Etiqueta de tiempo (opcional) ──
		if !isnothing(time_label)
			text!(ax, rail_span - 0.1, L_draw + 0.25,
				text     = time_label,
				fontsize = 12,
				color    = :gray40,
				align    = (:right, :top))
		end
	end

	md"Función `draw_pendulum_frame!` definida — dibuja un frame completo del sistema"
end

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000003
md"""
### 10.1 Frame estático — Exploración manual del estado

Usa el **slider de tiempo** para recorrer la trayectoria controlada frame a frame. Observa cómo el carro se desplaza y la barra oscila hasta estabilizarse en la vertical.
"""

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000004
md"**Tiempo t [s] — trayectoria controlada (LQR):**"

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000005
@bind t_frame PlutoUI.Slider(
	range(0.0, stop=10.0, length=501), default=0.0, show_value=true)

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000006
begin
	with_theme(TEMA_PENDULO) do

	# ── Estado actual interpolado de la solución ──
	state_now = sol_ctrl(t_frame)
	x_now, v_now, θ_now, ω_now = state_now
	F_now = clamp(-dot(lqr.K[1,:], state_now), -50.0, 50.0)

	# ── Longitud visual de la barra ──
	L_draw = Lbar   # longitud total de la barra para dibujo

	fig_frame = Figure(size = (820, 480))

	# ── Panel principal: animación del carro-péndulo ──
	ax_pend = Axis(fig_frame[1, 1:2],
		title      = "Péndulo invertido — Trayectoria LQR",
		xlabel     = "x  [m]",
		ylabel     = "y  [m]",
		aspect     = DataAspect(),
		limits     = (-2.5, 2.5, -0.4, L_draw + 0.4))

	draw_pendulum_frame!(ax_pend, x_now, θ_now, L_draw,
		cart_color = COL_CTRL,
		show_angle = true,
		time_label = @sprintf("t = %.2f s", t_frame))

	# ── Panel inferior: mini-gráfica de θ(t) con cursor temporal ──
	ax_mini = Axis(fig_frame[2, 1],
		xlabel = "t  [s]",
		ylabel = "θ  [rad]",
		height = 120,
		title  = "Ángulo θ(t)")
	lines!(ax_mini, t_ctrl, θ_ctrl, color = COL_CTRL, linewidth = 1.5)
	hlines!(ax_mini, [0.0], color = :gray60, linewidth = 0.5, linestyle = :dot)
	vlines!(ax_mini, [t_frame], color = :orange, linewidth = 1.8)
	scatter!(ax_mini, [t_frame], [θ_now],
		markersize = 8, color = :orange, strokecolor = :white, strokewidth = 1)

	# ── Panel inferior derecho: mini-gráfica de F(t) ──
	ax_force_mini = Axis(fig_frame[2, 2],
		xlabel = "t  [s]",
		ylabel = "F  [N]",
		height = 120,
		title  = "Fuerza de control F(t)")
	lines!(ax_force_mini, t_ctrl, F_ctrl_vec, color = COL_FUERZA, linewidth = 1.5)
	hlines!(ax_force_mini, [0.0], color = :gray60, linewidth = 0.5, linestyle = :dot)
	vlines!(ax_force_mini, [t_frame], color = :orange, linewidth = 1.8)
	scatter!(ax_force_mini, [t_frame], [F_now],
		markersize = 8, color = :orange, strokecolor = :white, strokewidth = 1)

	# ── Indicadores numéricos ──
	Label(fig_frame[0, :],
		@sprintf("x = %.3f m   |   ẋ = %.3f m/s   |   θ = %.2f°   |   ω = %.3f rad/s   |   F = %.1f N",
			x_now, v_now, rad2deg(θ_now), ω_now, F_now),
		fontsize = 12, color = :gray35, padding = (0, 0, 0, 4))

	rowgap!(fig_frame.layout, 12)
	colgap!(fig_frame.layout, 14)

	fig_frame
	end
end

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000007
md"""
### 10.2 Comparación animada — Sin control vs. LQR

Esta es la visualización más importante del proyecto: **el mismo péndulo bajo las mismas condiciones iniciales**, pero a la izquierda sin control (caída libre) y a la derecha estabilizado por el LQR.

Desliza el slider de tiempo y observa el contraste dramático:
- **Izquierda (rojo):** el péndulo cae irremediablemente — la simulación solo dura 3 s
- **Derecha (azul):** el carro se desplaza activamente para mantener la barra en equilibrio
"""

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000008
md"**Tiempo t [s] — comparación simultánea:**"

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000009
@bind t_comp PlutoUI.Slider(
	range(0.0, stop=10.0, length=501), default=0.0, show_value=true)

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000010
begin
	with_theme(TEMA_PENDULO) do

	L_draw_comp = Lbar
	t_max_free  = sol_free.t[end]   # límite de la simulación libre (~3 s)

	fig_comp = Figure(size = (900, 520))

	# ═══════════════════════════════════════════════
	# Panel izquierdo: SIN CONTROL (caída libre)
	# ═══════════════════════════════════════════════
	ax_libre = Axis(fig_comp[1, 1],
		title  = "Sin control (caída libre)",
		xlabel = "x  [m]",
		ylabel = "y  [m]",
		aspect = DataAspect(),
		limits = (-2.5, 2.5, -0.5, L_draw_comp + 0.5))

	if t_comp ≤ t_max_free
		s_free = sol_free(t_comp)
		draw_pendulum_frame!(ax_libre, s_free[1], s_free[3], L_draw_comp,
			cart_color = COL_LIBRE,
			bob_color  = RGBf(0.85, 0.20, 0.20),
			show_angle = true,
			time_label = @sprintf("t = %.2f s", t_comp))
	else
		# El péndulo ya cayó — mostrar el último estado con overlay
		s_last = sol_free(t_max_free)
		draw_pendulum_frame!(ax_libre, s_last[1], s_last[3], L_draw_comp,
			cart_color = RGBf(0.6, 0.4, 0.4),
			bob_color  = RGBf(0.6, 0.35, 0.35),
			time_label = @sprintf("t = %.2f s  (sim. terminó a t=%.1f)", t_comp, t_max_free))
		# Mensaje de caída
		text!(ax_libre, 0.0, L_draw_comp * 0.5,
			text     = "Péndulo caído\n(simulación terminada)",
			fontsize = 14,
			color    = COL_LIBRE,
			align    = (:center, :center))
	end

	# ═══════════════════════════════════════════════
	# Panel derecho: CON LQR (estabilización)
	# ═══════════════════════════════════════════════
	ax_ctrl_anim = Axis(fig_comp[1, 2],
		title  = "Con LQR (estabilización)",
		xlabel = "x  [m]",
		ylabel = "y  [m]",
		aspect = DataAspect(),
		limits = (-2.5, 2.5, -0.5, L_draw_comp + 0.5))

	s_ctrl = sol_ctrl(t_comp)
	draw_pendulum_frame!(ax_ctrl_anim, s_ctrl[1], s_ctrl[3], L_draw_comp,
		cart_color = COL_CTRL,
		bob_color  = RGBf(0.20, 0.50, 0.85),
		show_angle = true,
		time_label = @sprintf("t = %.2f s", t_comp))

	# ═══════════════════════════════════════════════
	# Mini-gráficas inferiores: θ(t) de ambos con cursor
	# ═══════════════════════════════════════════════
	ax_theta_comp = Axis(fig_comp[2, 1:2],
		xlabel = "t  [s]",
		ylabel = "θ  [rad]",
		height = 130,
		title  = "Comparación de θ(t)")
	lines!(ax_theta_comp, t_free, θ_free,
		color = COL_LIBRE, linewidth = 1.5, linestyle = :dash, label = "Sin control")
	lines!(ax_theta_comp, t_ctrl, θ_ctrl,
		color = COL_CTRL, linewidth = 1.5, label = "Con LQR")
	hlines!(ax_theta_comp, [0.0], color = :gray60, linewidth = 0.5, linestyle = :dot)
	vlines!(ax_theta_comp, [t_comp], color = :orange, linewidth = 1.8)
	axislegend(ax_theta_comp, position = :rt)

	# ── Título y layout ──
	Label(fig_comp[0, :],
		"Comparación: Caída libre vs. Estabilización LQR",
		fontsize = 15, font = :bold, padding = (0, 0, 0, 4))

	rowgap!(fig_comp.layout, 10)
	colgap!(fig_comp.layout, 16)

	fig_comp
	end
end

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000011
md"""
### 10.3 Exportar animación como GIF

Genera un archivo `.gif` animado de la comparación lado a lado usando `CairoMakie.record()`. Cada frame se renderiza en CPU y se ensambla en el archivo final.

Marca la casilla para iniciar la exportación (puede tomar 15–30 segundos):
"""

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000012
@bind do_export PlutoUI.CheckBox(default=false)

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000013
begin
	if do_export
		# ── Configuración del GIF ──
		fps_gif      = 30
		t_gif_end    = 6.0       # solo los primeros 6 s (después ya estabilizó)
		t_gif_range  = range(0.0, stop=t_gif_end, step=1/fps_gif)
		n_frames     = length(t_gif_range)
		L_draw_gif   = Lbar
		t_max_f      = sol_free.t[end]

		# ── Directorio de salida ──
		fig_dir = joinpath(@__DIR__, "..", "figures")
		mkpath(fig_dir)
		gif_path = joinpath(fig_dir, "04_comparacion_libre_vs_lqr.gif")

		# ── Figura para el GIF ──
		fig_gif = Figure(size = (900, 480), figure_padding = 10)

		ax_g1 = Axis(fig_gif[1, 1],
			title  = "Sin control",
			xlabel = "x [m]", ylabel = "y [m]",
			aspect = DataAspect(),
			limits = (-2.5, 2.5, -0.5, L_draw_gif + 0.5))

		ax_g2 = Axis(fig_gif[1, 2],
			title  = "Con LQR",
			xlabel = "x [m]", ylabel = "y [m]",
			aspect = DataAspect(),
			limits = (-2.5, 2.5, -0.5, L_draw_gif + 0.5))

		Label(fig_gif[0, :],
			"Péndulo Invertido — Caída libre vs. LQR  (θ₀ = $(round(rad2deg(θ₀), digits=1))°)",
			fontsize = 14, font = :bold)

		colgap!(fig_gif.layout, 16)

		# ── Generar el GIF ──
		record(fig_gif, gif_path, eachindex(t_gif_range); framerate = fps_gif) do i
			t = t_gif_range[i]

			# Panel izquierdo: sin control
			if t ≤ t_max_f
				sf = sol_free(t)
				draw_pendulum_frame!(ax_g1, sf[1], sf[3], L_draw_gif,
					cart_color = COL_LIBRE,
					bob_color  = RGBf(0.85, 0.20, 0.20),
					time_label = @sprintf("t = %.2f s", t))
			else
				sf = sol_free(t_max_f)
				draw_pendulum_frame!(ax_g1, sf[1], sf[3], L_draw_gif,
					cart_color = RGBf(0.6, 0.4, 0.4),
					bob_color  = RGBf(0.6, 0.35, 0.35),
					time_label = @sprintf("t = %.2f s", t))
				text!(ax_g1, 0.0, L_draw_gif * 0.5,
					text = "Péndulo caído", fontsize = 13, color = COL_LIBRE,
					align = (:center, :center))
			end

			# Panel derecho: con LQR
			sc = sol_ctrl(t)
			draw_pendulum_frame!(ax_g2, sc[1], sc[3], L_draw_gif,
				cart_color = COL_CTRL,
				bob_color  = RGBf(0.20, 0.50, 0.85),
				time_label = @sprintf("t = %.2f s", t))
		end

		Markdown.parse("""
**GIF exportado exitosamente**

- Ruta: `$gif_path`
- Frames: $n_frames
- Duración: $(t_gif_end) s a $(fps_gif) fps
- Resolución: 900 x 480 px

*Desmarca la casilla para evitar re-generar el GIF cada vez que cambies un parámetro.*
""")
	else
		md"Marca la casilla de arriba para exportar el GIF de la comparación."
	end
end

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000014
md"""
### 10.4 Exportar video MP4 (mayor calidad)

Si necesitas un video en formato MP4 (mejor compresión, más compatibilidad para presentaciones), puedes ejecutar este fragmento. MP4 requiere que `ffmpeg` esté instalado en tu sistema.
"""

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000015
@bind do_export_mp4 PlutoUI.CheckBox(default=false)

# ╔═╡ d4e5f6a7-3333-4000-b000-000000000016
begin
	if do_export_mp4
		fps_mp4      = 30
		t_mp4_end    = 8.0
		t_mp4_range  = range(0.0, stop=t_mp4_end, step=1/fps_mp4)
		n_frames_mp4 = length(t_mp4_range)
		L_draw_mp4   = Lbar
		t_max_free_mp4 = sol_free.t[end]

		mp4_dir = joinpath(@__DIR__, "..", "figures")
		mkpath(mp4_dir)
		mp4_path = joinpath(mp4_dir, "05_comparacion_libre_vs_lqr.mp4")

		fig_mp4 = Figure(size = (960, 500), figure_padding = 10)

		ax_m1 = Axis(fig_mp4[1, 1],
			title  = "Sin control (caída libre)",
			xlabel = "x [m]", ylabel = "y [m]",
			aspect = DataAspect(),
			limits = (-2.5, 2.5, -0.5, L_draw_mp4 + 0.5))

		ax_m2 = Axis(fig_mp4[1, 2],
			title  = "Con LQR (estabilización)",
			xlabel = "x [m]", ylabel = "y [m]",
			aspect = DataAspect(),
			limits = (-2.5, 2.5, -0.5, L_draw_mp4 + 0.5))

		Label(fig_mp4[0, :],
			"Péndulo Invertido — Comparación de control  (θ₀ = $(round(rad2deg(θ₀), digits=1))°)",
			fontsize = 14, font = :bold)

		colgap!(fig_mp4.layout, 16)

		record(fig_mp4, mp4_path, eachindex(t_mp4_range); framerate = fps_mp4) do i
			t = t_mp4_range[i]

			if t ≤ t_max_free_mp4
				sf = sol_free(t)
				draw_pendulum_frame!(ax_m1, sf[1], sf[3], L_draw_mp4,
					cart_color = COL_LIBRE,
					bob_color  = RGBf(0.85, 0.20, 0.20),
					time_label = @sprintf("t = %.2f s", t))
			else
				sf = sol_free(t_max_free_mp4)
				draw_pendulum_frame!(ax_m1, sf[1], sf[3], L_draw_mp4,
					cart_color = RGBf(0.6, 0.4, 0.4),
					bob_color  = RGBf(0.6, 0.35, 0.35),
					time_label = @sprintf("t = %.2f s", t))
				text!(ax_m1, 0.0, L_draw_mp4 * 0.5,
					text = "Péndulo caído", fontsize = 13, color = COL_LIBRE,
					align = (:center, :center))
			end

			sc = sol_ctrl(t)
			draw_pendulum_frame!(ax_m2, sc[1], sc[3], L_draw_mp4,
				cart_color = COL_CTRL,
				bob_color  = RGBf(0.20, 0.50, 0.85),
				time_label = @sprintf("t = %.2f s", t))
		end

		Markdown.parse("""
**MP4 exportado exitosamente**

- Ruta: `$mp4_path`
- Frames: $n_frames_mp4
- Duración: $(t_mp4_end) s a $(fps_mp4) fps

*Nota: MP4 requiere `ffmpeg`. Si ves un error, instálalo con `sudo apt install ffmpeg` (Linux) o `brew install ffmpeg` (macOS).*
""")
	else
		md"Marca la casilla para exportar el video MP4."
	end
end

# ╔═╡ 487aaede-2a86-11f1-acf5-215541c49e3a


# ╔═╡ Cell order:
# ╠═99b4c781-70d4-4e78-8512-e7f6b4fc6cf4
# ╟─be61f640-b844-420e-9f9d-af6ecaa5203a
# ╟─0955a7cb-8aac-470d-a4fe-9b7050342699
# ╠═db1608db-421a-4171-97e6-7062a27ebf3b
# ╟─4cd10adc-e79b-4dad-95e2-9c2cc7a40833
# ╟─54b11144-5dfd-4561-8580-4800c147d52e
# ╠═0565789e-462d-49a8-8534-413709c5b091
# ╟─93b07fd7-ef7a-43d1-a0fd-51bdd218d055
# ╠═12d58a06-f53a-4de0-a6fa-e1f1f14edf46
# ╟─ea2d982a-30dd-4d28-9665-2f43d851613b
# ╠═0a37ad6f-e44b-4706-bb70-48401d650700
# ╟─5138621e-392c-4218-9f5e-9b59b0707cbb
# ╠═dd76336c-5655-4e1b-9302-525ec8640eeb
# ╟─dc717647-fd75-4c58-97b3-8a423b78dfcf
# ╠═c691cfef-181c-4c2b-bc2a-fbdda4623861
# ╠═71f6c6e0-a119-43f6-9a4e-3a7ad1bcf0b0
# ╟─4d23db82-2560-4a28-b84d-fe6395055ae5
# ╠═4250622d-0830-4333-9f16-d7198160c5a3
# ╟─b1e47212-c25b-4d83-934f-0fee7f773174
# ╠═9d0a35e5-5768-4d2f-a5aa-1814b53e8721
# ╟─9a603b09-b6ba-4e73-8d51-695015cae9d8
# ╠═375e2a1f-27d8-466f-a343-d81c95ce2a9a
# ╟─4426d097-d814-40ee-aa8b-f75d8ee0789d
# ╠═18152fa9-97fe-4400-8bd5-8eec10276d92
# ╟─329d269b-1534-4460-86d8-1e9eebb9bc16
# ╠═027d6af7-96a0-4f7c-8d57-dafd92fbb930
# ╟─cca01fd6-94af-48cc-ae8d-84a8ca36551c
# ╠═86fd0c62-d8f6-4811-b7d3-b52ba62eb769
# ╟─d26e24fe-cfa6-4c42-80d7-27175c4e6e43
# ╠═aea2bad2-82fa-4dc6-879e-1edbaa4f7e4d
# ╟─a1b2c3d4-1111-4000-8000-000000000001
# ╟─a1b2c3d4-1111-4000-8000-000000000002
# ╠═a1b2c3d4-1111-4000-8000-000000000003
# ╟─a1b2c3d4-1111-4000-8000-000000000004
# ╠═a1b2c3d4-1111-4000-8000-000000000005
# ╟─a1b2c3d4-1111-4000-8000-000000000006
# ╠═a1b2c3d4-1111-4000-8000-000000000007
# ╠═a1b2c3d4-1111-4000-8000-000000000008
# ╟─a1b2c3d4-1111-4000-8000-000000000009
# ╠═a1b2c3d4-1111-4000-8000-000000000010
# ╟─a1b2c3d4-1111-4000-8000-000000000011
# ╟─a1b2c3d4-1111-4000-8000-000000000012
# ╠═a1b2c3d4-1111-4000-8000-000000000013
# ╠═a1b2c3d4-1111-4000-8000-000000000014
# ╟─c3d4e5f6-2222-4000-a000-000000000001
# ╠═c3d4e5f6-2222-4000-a000-000000000002
# ╠═c3d4e5f6-2222-4000-a000-000000000003
# ╟─c3d4e5f6-2222-4000-a000-000000000010
# ╠═c3d4e5f6-2222-4000-a000-000000000004
# ╟─c3d4e5f6-2222-4000-a000-000000000011
# ╠═c3d4e5f6-2222-4000-a000-000000000005
# ╟─c3d4e5f6-2222-4000-a000-000000000012
# ╠═c3d4e5f6-2222-4000-a000-000000000006
# ╟─c3d4e5f6-2222-4000-a000-000000000013
# ╠═c3d4e5f6-2222-4000-a000-000000000007
# ╟─c3d4e5f6-2222-4000-a000-000000000014
# ╠═c3d4e5f6-2222-4000-a000-000000000008
# ╟─a1b2c3d4-1111-4000-8000-000000000015
# ╠═a1b2c3d4-1111-4000-8000-000000000016
# ╟─a1b2c3d4-1111-4000-8000-000000000017
# ╟─d4e5f6a7-3333-4000-b000-000000000001
# ╠═d4e5f6a7-3333-4000-b000-000000000002
# ╟─d4e5f6a7-3333-4000-b000-000000000003
# ╟─d4e5f6a7-3333-4000-b000-000000000004
# ╠═d4e5f6a7-3333-4000-b000-000000000005
# ╠═d4e5f6a7-3333-4000-b000-000000000006
# ╟─d4e5f6a7-3333-4000-b000-000000000007
# ╟─d4e5f6a7-3333-4000-b000-000000000008
# ╠═d4e5f6a7-3333-4000-b000-000000000009
# ╠═d4e5f6a7-3333-4000-b000-000000000010
# ╟─d4e5f6a7-3333-4000-b000-000000000011
# ╠═d4e5f6a7-3333-4000-b000-000000000012
# ╠═d4e5f6a7-3333-4000-b000-000000000013
# ╟─d4e5f6a7-3333-4000-b000-000000000014
# ╠═d4e5f6a7-3333-4000-b000-000000000015
# ╠═d4e5f6a7-3333-4000-b000-000000000016
# ╟─487aaede-2a86-11f1-acf5-215541c49e3a
