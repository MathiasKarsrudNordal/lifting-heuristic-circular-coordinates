include("../src/datasets.jl")
include("../src/circular-coordinates.jl")
using .Datasets
using .CircularCoordinates

using CairoMakie

savefig = true

prime = 47
n_points = 500

X, θ = Datasets.annulus(n_points, 1.0, 1.1)

θ_cc = CircularCoordinates.get_circular_coordinates(X, prime)

scale_factors = [1, 2, 5, 10]
θ_scaled = [(θ_cc .* scale_factor) .% 1.0 for scale_factor in scale_factors]

fig_cc = Figure(size = (1000, 500), fontsize = 14, figure_padding = 25)

for (i, θs) in enumerate(θ_scaled)
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

if savefig
    save("figs/different-order-circular-coordinates.pdf", fig_cc)
end
