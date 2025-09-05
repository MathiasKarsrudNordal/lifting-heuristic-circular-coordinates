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
n_points = 1000
embedding_dim = 300

# Generate data
θ, X_embedded = Datasets.circle(n_points, embedding_dim)
X = row_tuples = [Tuple(X_embedded[i, :]) for i in 1:size(X_embedded, 1)]

## Print first and last value of θ to verify range
println("First angle: ", θ[1])
println("Last angle: ", θ[end])


## ============================== PCA Projection ==============================
# Apply PCA
model = fit(PCA, X_embedded'; maxoutdim=2)
X_pca = transform(model, X_embedded')'

# Remove points with x value > 2 (only first column)
X_pca[abs.(X_pca[:, 1]) .> 2, 1] .= 0

fig_pca = Figure(size = (600, 600), fontsize = 14, figure_padding = 75)
ax = Axis(fig_pca[1, 1];
    title = latexstring("PCA Projection of \$S^1 \\subseteq \\mathbb{R}^{$embedding_dim}\$"),
    titlesize = 28,   # << increase title font size
    aspect = 1,
)

scatter!(ax, X_pca[:, 1], X_pca[:, 2];
    markersize = 10,
    color = θ,
    colormap = :plasma,
    transparency = true
)

hidexdecorations!(ax; ticks=true, ticklabels=true, grid=true)
hideydecorations!(ax; ticks=true, ticklabels=true, grid=true)

fig_pca
display(fig_pca)

if savefig
    save("figs/pca-high-dim-circle.pdf", fig_pca)
end

## ============================== Persistence Diagram ======================
p = 7
thresh = 20

DistMatrix = pairwise(Euclidean(), X_embedded')  # (n_points × n_points)
rips = Ripserer.ripserer(DistMatrix; dim_max=1, modulus=p, threshold=thresh, verbose=true)

fig_ph = Utils.plot_persistence_diagram(rips; infinity=thresh, palette=(:blue, :orange, :green))
display(fig_ph)

if savefig
    save("figs/persistence-diagram-high-dim-circle.pdf", fig_ph)
end

## ============================== Circular Coordinates ======================
th = 17.3
θ_cc = CircularCoordinates.get_circular_coordinates(X, p; threshold = th, verbose = true)

scales = [1, 5]
fig = Figure(size = (1000, 1000), fontsize = 14, figure_padding = 75)
for (i, scale) in enumerate(scales)
    θ_scaled = (θ_cc .* scale) .% 1.0

    ax_pca = Axis(fig[i, 1]; aspect = 1)
    scatter!(ax_pca, X_pca[:, 1], X_pca[:, 2];
        markersize = 10,
        color = θ_scaled,
        colormap = :plasma,
        transparency = true
    )
    hidexdecorations!(ax_pca; ticks=true, ticklabels=true, grid=true)
    hideydecorations!(ax_pca; ticks=true, ticklabels=true, grid=true)

    ax_cc_vs_orig = Axis(fig[i, 2]; aspect = 1)
    scatter!(ax_cc_vs_orig, θ, θ_scaled;
        markersize = 6,
        color = θ_scaled,
        colormap = :plasma,
        transparency = true
    )
    hidexdecorations!(ax_cc_vs_orig; ticks=true, ticklabels=true, grid=true)
    hideydecorations!(ax_cc_vs_orig; ticks=true, ticklabels=true, grid=true)
end
display(fig)

if savefig
    save("figs/different-order-circular-coordinates-high-dim-circle.pdf", fig)
end