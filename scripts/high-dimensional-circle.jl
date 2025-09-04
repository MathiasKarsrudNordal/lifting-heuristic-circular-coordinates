include("../src/datasets.jl")
using .Datasets
using MultivariateStats, LinearAlgebra, Statistics
using CairoMakie
using LaTeXStrings

savefig = false
n_points = 1000
embedding_dim = 300

# Generate data
θ, X_embedded = circle(n_points, embedding_dim)

## ============================== PCA Projection ==============================
# Apply PCA
model = fit(PCA, X_embedded'; maxoutdim=2)
X_pca = transform(model, X_embedded')'

# Remove points with x value > 2 (only first column)
X_pca[abs.(X_pca[:, 1]) .> 2, 1] .= 0

println("Explained variance ratio = ", principalvars(model) ./ tvar(model))

fig_pca = Figure(resolution = (600, 600), fontsize = 14, figure_padding = 75)
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

## ============================== Circular Coordinates ==============================