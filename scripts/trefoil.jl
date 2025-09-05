include("../src/datasets.jl")
include("../src/circular-coordinates.jl")
include("../src/utils.jl")
using .Datasets
using .Utils
using .CircularCoordinates

using MultivariateStats, LinearAlgebra, Statistics
using Ripserer, Distances
using CairoMakie
using LaTeXStrings

savefig = false
n_points = 2500

X, θ = Datasets.torus_knot(p=2, q=3, n_samples=n_points, R=5.0, r=2.0, noise_var=0.1, type_noise="gaussian")

# === 2D scatter plot of the point cloud colored by true angle ===
fig_2d = Figure(size = (400, 400), fontsize = 14, figure_padding = 25)
ax_2d = Axis(fig_2d[1, 1]; aspect = 1)
scatter!(ax_2d, X[:, 1], X[:, 2];
    markersize = 4,
    color = θ,
    colormap = :plasma,
    transparency = true
)
hidexdecorations!(ax_2d; ticks=true, ticklabels=true, grid=true)
hideydecorations!(ax_2d; ticks=true, ticklabels=true, grid=true)
display(fig_2d)

# === Persistence Diagram ===
p = 47

distance_matrix = pairwise(Euclidean(), X')
rips = Ripserer.ripserer(distance_matrix; dim_max = 1, modulus = p, verbose = true)

fig_ph = Utils.plot_persistence_diagram(rips; infinity = 5.0, palette = (:blue, :orange, :green))
display(fig_ph)

# === Circular Coordinates ===
th = 3.5

point_cloud = [Tuple(X[i, :]) for i in 1:size(X, 1)]
θ_cc = CircularCoordinates.get_circular_coordinates(point_cloud, p; threshold = 3.5, verbose = true)



scales = [1, 2]
fig = Figure(size = (800, 400), fontsize = 14, figure_padding = 25)
for (i, scale) in enumerate(scales)
    θ_scaled = (θ_cc .* scale) .% 1.0
    ax = Axis(fig[1, i]; aspect = 1)
    scatter!(ax, X[:, 1], X[:, 2];
        markersize = 4,
        color = θ_scaled,
        colormap = :plasma,
        transparency = true
    )
    hidexdecorations!(ax; ticks=true, ticklabels=true, grid=true)
    hideydecorations!(ax; ticks=true, ticklabels=true, grid=true)
end
display(fig)

if savefig
    save("figs/torus-knot-circular-coordinates.pdf", fig)
end