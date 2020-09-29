### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ fd505a66-0284-11eb-18f2-bd1700f88af0
md"# Workshop 9"

# ╔═╡ 52d73384-027a-11eb-1df0-bbed88425c06
@enum Methods begin
	natural
	adjusted
	clamped
	paraterm
	notaknot
end

# ╔═╡ 4fe4b6f8-02a0-11eb-0d56-5d561455730e
md"## 1. Code: Splines"

# ╔═╡ 466438fa-0265-11eb-133c-4d7f9bfe66cb
function cubic_splines(inter_points, method, slopes=nothing)
	n = length(inter_points)
	xs = [x[1] for x in inter_points]
	ys = [y[2] for y in inter_points]
	
	
	matrix = zeros(n, n)
	RHS = zeros(1, n)'

	
	for i in 2:(n-1)
		RHS[i] = 6.0 * (((ys[i+1] - ys[i]) / (xs[i+1] - xs[i])) - ((ys[i] - ys[i-1]) / (xs[i] - xs[i-1])))
	
		for j in 1:n
			if j == i-1
				matrix[i, j] += (xs[i] - xs[i-1])
			elseif j == i
				matrix[i, j] += 2.0 * (xs[i+1] - xs[i-1])
			elseif j == i+1
				matrix[i, j] += (xs[i+1] - xs[i])
			end
		end
	end
	
	# Natrual Spline Boundary Conditions
	if method == natural 
		matrix[1, 1] = 1.0
		RHS[1] = 0.0
		matrix[n, n] = 1.0
		RHS[n] = 0.0
	elseif method == adjusted
		matrix[1, 1] = 1.0
		RHS[1] = slopes[1]
		matrix[n, n] = 1.0
		RHS[n] = slopes[2]
	elseif method == clamped
		matrix[1, 1] = 2.0
		matrix[1, 2] = 1.0
		RHS[1] = 6.0 * ((ys[2] - ys[1]) / (xs[2] - xs[1]) - slopes[1])
		matrix[n, n-1] = 2.0
		matrix[n, n] = 1.0
		RHS[n] = 6.0 * (slopes[2] - (ys[n] - ys[n-1]) / (xs[n] - xs[n-1]))
	elseif method == paraterm
		matrix[1, 1] = 1.0
		matrix[1, 2] = -1.0
		RHS[1] = 0.0
		matrix[n, n-1] = -1.0
		matrix[n, n] = 1.0
		RHS[n] = 0.0
	elseif method == notaknot
		matrix[1, 1] = xs[3] - xs[2]
		matrix[1, 2] = xs[1] - xs[3]
		matrix[1, 3] = xs[2] - xs[1]
		RHS[1] = 0.0
		matrix[n, n-2] = xs[n] - xs[n-1]
		matrix[n, n-1] = xs[n-2] - xs[n]
		matrix[n, n]   = xs[n-1] - xs[n-2]
		RHS[n] = 0.0
	else
		println("ERROR: Unknown Boundry Condition Method")
		return nothing
	end
	
	
	
	solution = matrix \ RHS
	
	outputs = []
	points = 20
	for i in 1:(n-1)
		dx = xs[i+1] - xs[i]
		
		for j in 2:points
			b = float(j) / float(points)
			a = 1-b
			x = xs[i] + b*dx
			
			c = ((a^3 - a) * dx^2) / 6.0
			d = ((b^3 - b) * dx^2) / 6.0
			
			y = a*ys[i] + b*ys[i+1] + c*solution[i] + d*solution[i+1]
			
			push!(outputs, (x,y))
		end
	end
	return outputs
end

# ╔═╡ ed67cc88-029d-11eb-0c98-8fceab493fd7
md"### Function $f(x)$ and It's Derivatives 
$f(x) = x^3 - 5x^2 + 6x -1$
$f'(x) = 3x^2 - 10x +6$
$f''(x) = 6x - 10$"

# ╔═╡ ea0b4254-0278-11eb-3f56-2b4f0f42ab8f
function f(x)
	(x.^3) .- (5 .*x.^2) .+ (6 .*x) .- 1
end

# ╔═╡ b03a9824-029e-11eb-1199-5dbc03aa75be
function f_prime(x)
	(3 .* x.^2) .- (10 .* x) .+ 6
end

# ╔═╡ d760634a-029e-11eb-23eb-097475fa236a
function f_double_prime(x)
	(6 .*x) .- 10
end

# ╔═╡ 01e53634-0279-11eb-305a-850dc557fdd4
function fx_values_1to6()
	points = []
	for i in 1:6
		push!(points, (i, f(i)))
	end
	return points
end

# ╔═╡ e2be020a-029f-11eb-1701-d908a28f6fc8
md"## 3. Validation"

# ╔═╡ 35364848-02a1-11eb-124b-a14c5bc1e3f7
md"Here we can see that our cubic spline, $y_{spline}(x)$, completely overlaps our oringial function $f(x)$ when using the curvature adjusted boundary conditions."

# ╔═╡ 3eca0b2e-0279-11eb-1094-3d64078f9a70
spline_adjusted = cubic_splines(fx_values_1to6(), adjusted, (f_double_prime(1), f_double_prime(6)))

# ╔═╡ c72c98fc-0278-11eb-29bb-c502157ec6cc
begin 
	if length(size(spline_adjusted)) == 1
		using Plots
		plotly()
		plot([p[1] for p in spline_adjusted], [p[2] for p in spline_adjusted], label="y_spline(x)")
		plot!(f, 0, 6, label="f(x)")
	end
end

# ╔═╡ 80b1158e-029f-11eb-1ead-ed01fae0ced9
md"## 4. Invalidation?"

# ╔═╡ df509b0e-02a2-11eb-3db0-59a8debaffd9
md"From the following two plots, $y_{spline}(x)$ vs. $f(x)$ and $y_{spline}(x) - f(x)$ vs. $x$, we can see that all methods are relatively close to the original function. However, from the error plot it is clear that the natural, clamped, and parabolically terminated boundry conditions do not perform as well as the curvature-adjusted and not-a-knot boundary conditions."

# ╔═╡ 01cac6da-0299-11eb-16cb-3d878730be0d
spline_natural = cubic_splines(fx_values_1to6(), natural)

# ╔═╡ f7298d78-0299-11eb-286e-cb9a696f7c86
spline_clamped = cubic_splines(fx_values_1to6(), clamped, (f_prime(1), f_prime(6)))

# ╔═╡ 411e3d42-029b-11eb-0072-6bfc707542e5
spline_paraterm = cubic_splines(fx_values_1to6(), paraterm)

# ╔═╡ 49c84f8c-029b-11eb-2d46-995485ff671d
spline_notaknot = cubic_splines(fx_values_1to6(), notaknot)

# ╔═╡ 487d62b4-02a2-11eb-1815-2b292c34895f
begin 
	if length(size(spline_adjusted)) == 1
		plotly()
		plot( [p[2] for p in spline_natural], f([p[1] for p in spline_natural]), label="Natural")
		plot!([p[2] for p in spline_adjusted], f([p[1] for p in spline_adjusted]), label="Adjusted")
		plot!([p[2] for p in spline_clamped], f([p[1] for p in spline_clamped]), label="Clamped")
		plot!([p[2] for p in spline_paraterm], f([p[1] for p in spline_paraterm]), label="ParaTerm")
		plot!([p[2] for p in spline_notaknot], f([p[1] for p in spline_notaknot]), label="NotaKnot")
		title!("y_spline(x) vs. f(x)")
		xaxis!("f(x)")
		yaxis!("y_spline(x)")
	end
end

# ╔═╡ 27d494f2-029c-11eb-282d-1f99d853604d
begin 
	if length(size(spline_adjusted)) == 1
		plotly()
		plot( [p[1] for p in spline_natural], [p[2] for p in spline_natural] .- f([p[1] for p in spline_natural]), label="Natural")
		plot!([p[1] for p in spline_adjusted], [p[2] for p in spline_adjusted] .- f([p[1] for p in spline_adjusted]), label="Adjusted")
		plot!([p[1] for p in spline_clamped], [p[2] for p in spline_clamped] .- f([p[1] for p in spline_clamped]), label="Clamped")
		plot!([p[1] for p in spline_paraterm], [p[2] for p in spline_paraterm] .- f([p[1] for p in spline_paraterm]), label="ParaTerm")
		plot!([p[1] for p in spline_notaknot], [p[2] for p in spline_notaknot] .- f([p[1] for p in spline_notaknot]), label="NotaKnot")
		title!("Error")
		xaxis!("x")
		yaxis!("y_spline(x) - f(x)")
	end
end

# ╔═╡ 90b43690-026e-11eb-0c1d-0dcdc65b929f
md"### ---- Testing ----"

# ╔═╡ cb08577a-026b-11eb-02af-c330d801b36f
A = [1 1 0; -1 1 0; 0 0 -0.5]

# ╔═╡ 5ef43404-026c-11eb-296d-3b1fb87305f1
B = [1 1 1]'

# ╔═╡ c55cece6-026e-11eb-2339-0b96e48aac3a
x = B \ A

# ╔═╡ Cell order:
# ╟─fd505a66-0284-11eb-18f2-bd1700f88af0
# ╠═52d73384-027a-11eb-1df0-bbed88425c06
# ╟─4fe4b6f8-02a0-11eb-0d56-5d561455730e
# ╠═466438fa-0265-11eb-133c-4d7f9bfe66cb
# ╟─ed67cc88-029d-11eb-0c98-8fceab493fd7
# ╟─ea0b4254-0278-11eb-3f56-2b4f0f42ab8f
# ╟─b03a9824-029e-11eb-1199-5dbc03aa75be
# ╟─d760634a-029e-11eb-23eb-097475fa236a
# ╟─01e53634-0279-11eb-305a-850dc557fdd4
# ╟─e2be020a-029f-11eb-1701-d908a28f6fc8
# ╟─35364848-02a1-11eb-124b-a14c5bc1e3f7
# ╠═3eca0b2e-0279-11eb-1094-3d64078f9a70
# ╠═c72c98fc-0278-11eb-29bb-c502157ec6cc
# ╟─80b1158e-029f-11eb-1ead-ed01fae0ced9
# ╟─df509b0e-02a2-11eb-3db0-59a8debaffd9
# ╠═01cac6da-0299-11eb-16cb-3d878730be0d
# ╠═f7298d78-0299-11eb-286e-cb9a696f7c86
# ╠═411e3d42-029b-11eb-0072-6bfc707542e5
# ╠═49c84f8c-029b-11eb-2d46-995485ff671d
# ╠═487d62b4-02a2-11eb-1815-2b292c34895f
# ╠═27d494f2-029c-11eb-282d-1f99d853604d
# ╟─90b43690-026e-11eb-0c1d-0dcdc65b929f
# ╠═cb08577a-026b-11eb-02af-c330d801b36f
# ╠═5ef43404-026c-11eb-296d-3b1fb87305f1
# ╠═c55cece6-026e-11eb-2339-0b96e48aac3a
