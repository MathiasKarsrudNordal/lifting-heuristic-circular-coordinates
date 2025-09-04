include("../src/datasets.jl")
using .Datasets
using Distances
using Plots

N = 1000
data = Datasets.annulus(N, 1, 2, (0, 0))

x_coords = getindex.(data, 1)
y_coords = getindex.(data, 2)

data_matrix = hcat(x_coords, y_coords)'

distance_matrix = Distances.pairwise(Euclidean(), data_matrix);

plt = Plots.scatter(data; label="data", markersize=2, aspect_ratio=1)
display(plt)