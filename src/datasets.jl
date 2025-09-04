module Datasets

using LinearAlgebra
using Distances

export annulus

function annulus(n, r1=1, r2=2, offset=(0.0, 0.0))
    θ = 2π .* rand(n)

    r = sqrt.(rand(n) .* (r2^2 - r1^2) .+ r1^2)
    x = r .* cos.(θ) .+ offset[1]
    y = r .* sin.(θ) .+ offset[2]
    return [(x[i], y[i]) for i in 1:n]
end

end # module Datasets