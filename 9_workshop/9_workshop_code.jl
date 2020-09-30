
@enum Methods begin
	natural
	adjusted
	clamped
	paraterm
	notaknot
end


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


function f(x)
	(x.^3) .- (5 .*x.^2) .+ (6 .*x) .- 1
end

function f_prime(x)
	(3 .* x.^2) .- (10 .* x) .+ 6
end

function f_double_prime(x)
	(6 .*x) .- 10
end

function fx_values_1to6()
	points = []
	for i in 1:6
		push!(points, (i, f(i)))
	end
	return points
end


# 3. Validation

spline_adjusted = cubic_splines(fx_values_1to6(), adjusted, (f_double_prime(1), f_double_prime(6)))

begin
	if length(size(spline_adjusted)) == 1
		using Plots
		plotly()
		plot([p[1] for p in spline_adjusted], [p[2] for p in spline_adjusted], label="y_spline(x)")
		plot!(f, 0, 6, label="f(x)")
	end
end


# 4. Invalidation

spline_natural = cubic_splines(fx_values_1to6(), natural)

spline_clamped = cubic_splines(fx_values_1to6(), clamped, (f_prime(1), f_prime(6)))

spline_paraterm = cubic_splines(fx_values_1to6(), paraterm)

spline_notaknot = cubic_splines(fx_values_1to6(), notaknot)


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
