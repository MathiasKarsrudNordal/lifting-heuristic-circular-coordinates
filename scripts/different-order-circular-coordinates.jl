include("../src/datasets.jl")
include("../src/circular-coordinates.jl")
using .Datasets
using .CircularCoordinates

using CairoMakie

savefig = true

prime = 47
n_points = 500

X, θ = Datasets.annulus(n_points, 1.0, 1.1)

scale_factors = [1, 2, 5, 10]
θ_cc_array = CircularCoordinates.get_circular_coordinates(X, prime; multpl = scale_factors, verbose = true)

fig_cc = Figure(size = (1000, 500), fontsize = 14, figure_padding = 25)

for (i, θs) in enumerate(θ_cc_array)
    ax = Axis(fig_cc[1, i]; aspect = 1)
    
    scatter!(ax, first.(X), last.(X);
        markersize = 4,
        color = θs,
        colormap = :plasma,
        transparency = true
    )
    
    hidexdecorations!(ax; ticks=true, ticklabels=true, grid=true)
    hideydecorations!(ax; ticks=true, ticklabels=true, grid=true)
end
display(fig_cc)

if savefig
    save("figs/different-order-circular-coordinates.pdf", fig_cc)
end
