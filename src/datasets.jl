module Datasets

using LinearAlgebra
using Distances
using Random
using Distributions

export annulus, circle, torus_knot

function annulus(n::Int, r1::Float64=1.0, r2::Float64=2.0, offset::Tuple{Float64, Float64}=(0.0, 0.0))
    θ = 2π .* rand(n)

    r = sqrt.(rand(n) .* (r2^2 - r1^2) .+ r1^2)
    x = r .* cos.(θ) .+ offset[1]
    y = r .* sin.(θ) .+ offset[2]
    return [(x[i], y[i]) for i in 1:n], θ ./ (2π)
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

function torus_knot(; p=2, q=3, n_samples=2500, R=5.0, r=2.0, noise_var=0.2, type_noise="gaussian") 
    Random.seed!(0)  # reproducibility

    # uniform theta samples
    θ = rand(Uniform(0, 2π), n_samples)

    if type_noise == "gaussian"
        noise = Normal(0, sqrt(noise_var))
        noise_x = rand(noise, n_samples)
        noise_y = rand(noise, n_samples)
        noise_z = rand(noise, n_samples)
    elseif type_noise == "uniform"
        noise = Uniform(-sqrt(3*noise_var), sqrt(3*noise_var))
        noise_x = rand(noise, n_samples)
        noise_y = rand(noise, n_samples)
        noise_z = rand(noise, n_samples)
    else
        error("Unknown noise type: $type_noise. Use 'gaussian' or 'uniform'.")
    end
    
    # parametric equations of the knot
    x = (R .+ r .* cos.(q .* θ)) .* cos.(p .* θ) .+ noise_x
    y = (R .+ r .* cos.(q .* θ)) .* sin.(p .* θ) .+ noise_y
    z = r .* sin.(q .* θ) .+ noise_z

    X = hcat(x, y, z)
    return X, θ ./ (2π)
end

end # module Datasets