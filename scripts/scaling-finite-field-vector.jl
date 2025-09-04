using LinearAlgebra
using Distances
using CairoMakie

n = 4
m = 10
prime = 160003
L = floor((prime - 1) / m)

lower_bound_p = 1 + (2 * m)^n
println("Lower bound on p according to proof: $lower_bound_p")

vect_init = rand(0:prime-1, n) # random vector in (Z/pZ)^n
println("Initial vector: $vect_init")

vect = copy(vect_init)

nframes = 250
f = Figure(resolution = (400, 400))
ax = Axis(f[1, 1], title = "Scaling finite field vector", aspect = 1)
poly!(ax, Circle(Point2f(0.0, 0.0), 1.0),
      color=:transparent, strokecolor=:blue, strokewidth=2)

# scatter points for evolving vector
circular_reps = [Point2f(cos(2π * α/prime), sin(2π * α/prime)) for α in vect]
sc = scatter!(ax, [pt[1] for pt in circular_reps],
                 [pt[2] for pt in circular_reps],
                 color=:red, markersize=15)

found_opt_scaler = Ref(false)
opt_mu = Ref{Union{Int,Nothing}}(nothing)

record(f, "figs/videos/optimal-scaler.mp4", 1:nframes) do frame
    μ = found_opt_scaler[] ? opt_mu[] : rand(1:prime-1)

    vect .= (vect_init .* μ) .% prime
    println("Frame $frame: scaling by μ = $μ -> $vect")

    if !found_opt_scaler[] && all((vect .<= L) .| (vect .>= prime - L))
        found_opt_scaler[] = true
        opt_mu[] = μ
        println("Found optimal scaler μ = $(opt_mu[]) at frame $frame")
    end

    circular_reps .= [Point2f(cos(2π * α/prime), sin(2π * α/prime)) for α in vect]
    sc[1] = [pt[1] for pt in circular_reps]
    sc[2] = [pt[2] for pt in circular_reps]
end

