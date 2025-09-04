module Datasets

using LinearAlgebra
using Distances

export annulus
export circle

function annulus(n::Int, r1::Float64=1.0, r2::Float64=2.0, offset::Tuple{Float64, Float64}=(0.0, 0.0))
    θ = 2π .* rand(n)

    r = sqrt.(rand(n) .* (r2^2 - r1^2) .+ r1^2)
    x = r .* cos.(θ) .+ offset[1]
    y = r .* sin.(θ) .+ offset[2]
    return [(x[i], y[i]) for i in 1:n]
end

function circle(n_points::Int, embedding_dim::Int)
    embedding_dim < 2 && throw(ArgumentError("embedding_dim must be at least 2"))

    θ = range(0, 2π; length=n_points+1)[1:end-1]  # avoid duplicate point at 2π
    θ_norm = θ ./ (2π)

    X_embedded = zeros(n_points, embedding_dim)

    # Use pairs of dimensions to create circles with increasing frequency
    for (i, freq) in enumerate(1:div(embedding_dim, 2))
        X_embedded[:, 2i-1] = cos.(freq .* θ)
        X_embedded[:, 2i] = sin.(freq .* θ)
    end

    # If embedding_dim is odd, fill the last dimension with
    if isodd(embedding_dim)
        frequency = div(embedding_dim, 2) + 1
        X_embedded[:, end] = cos.(frequency .* θ)
    end

    println("Generated circle with $n_points points in $embedding_dim dimensions.")
    println("Returning angles of shape $(size(θ_norm)) and point cloud of shape $(size(X_embedded)).")

    return θ_norm, X_embedded
end

end # module Datasets