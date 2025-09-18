
using SparseArrays

"""
    build_complex_and_boundary_operator(vertices, edges)

Constructs the simplicial complex up to dimension 2 from a graph and computes
the boundary operator matrix D = [∂₁ ∂₂].

A 2-simplex (triangle) is included if and only if its three bounding edges
are all present in the `edges` list.

# Arguments
- `vertices`: An iterable collection of vertex identifiers (e.g., `1:100`).
- `edges`: A `Vector` of 2-element `Tuples` or `Vectors`, e.g., `[(1, 2), (2, 3), ...]`,
  representing the 1-simplices (edges) of the graph.

# Returns
- `D::SparseMatrixCSC{Int, Int}`: The boundary operator matrix.
- `simplex_indices::Dict{Symbol, Any}`: A dictionary containing mappings from
  each simplex to its corresponding index range in the matrix `D`.
  Keys are `:vertices`, `:edges`, and `:triangles`.
"""
function build_complex_and_boundary_operator(vertices, edges)
    # === 1. CANONICALIZE INPUT AND FIND 2-SIMPLICES (TRIANGLES) ===

    # Ensure all edges are in a canonical format (v₁, v₂) where v₁ < v₂
    canonical_edges = Set(sort(collect(e)) |> Tuple for e in edges)

    # Build an adjacency structure for fast neighbor lookups (O(|E|))
    adj = Dict{Int, Set{Int}}(v => Set{Int}() for v in vertices)
    for (u, v) in canonical_edges
        push!(adj[u], v)
        push!(adj[v], u)
    end

    # Find all unique triangles (2-simplices) efficiently (O(Σ deg(v)²))
    # A triangle {u, v, w} exists if (u,v), (v,w), and (u,w) are all edges.
    triangles = Set{NTuple{3, Int}}()
    for u in vertices
        # Iterate over pairs of neighbors of u
        neighbors_of_u = sort(collect(adj[u]))
        for i in 1:length(neighbors_of_u)
            for j in (i+1):length(neighbors_of_u)
                v = neighbors_of_u[i]
                w = neighbors_of_u[j]

                # If v and w are connected, then {u, v, w} is a triangle.
                # We already know (u,v) and (u,w) are edges.
                if w in adj[v] # This is the correct way to check for membership in a Set

                    # Add in sorted order to ensure canonical representation
                    push!(triangles, tuple(sort([u, v, w])...))
                end
            end
        end
    end

    # === 2. CREATE STABLE ORDERING AND INDEX MAPPINGS ===

    # Sort simplices to have a deterministic order for matrix construction
    sorted_vertices = sort(collect(vertices))
    sorted_edges = sort(collect(canonical_edges))
    sorted_triangles = sort(collect(triangles))

    n_verts = length(sorted_vertices)
    n_edges = length(sorted_edges)
    n_triangles = length(sorted_triangles)

    # Dictionaries to map each simplex to its column index in the matrix D
    edge_to_idx = Dict(edge => i + n_verts for (i, edge) in enumerate(sorted_edges))
    triangle_to_idx = Dict(tri => i + n_verts + n_edges for (i, tri) in enumerate(sorted_triangles))

    # === 3. CONSTRUCT THE BOUNDARY OPERATOR MATRIX D ===

    total_dim = n_verts + n_edges + n_triangles
    
    # Use I, J, V format for efficient sparse matrix creation
    I, J, V = Int[], Int[], Int[]

    # Part 1: ∂₁ (Edges -> Vertices)
    # For an edge (u, v), its boundary is ∂(u,v) = v - u
    for (edge, col_idx) in edge_to_idx
        u, v = edge # edge is canonical, so u < v
        # Add -1 for the starting vertex row
        push!(I, u); push!(J, col_idx); push!(V, 1)
        # Add +1 for the ending vertex row
        push!(I, v); push!(J, col_idx); push!(V, -1)
    end

    # Part 2: ∂₂ (Triangles -> Edges)
    # For a triangle (u, v, w), ∂(u,v,w) = (v,w) - (u,w) + (u,v)
    for (triangle, col_idx) in triangle_to_idx
        u, v, w = triangle # triangle is canonical, u < v < w
        
        # Get the row indices for the three bounding edges
        row_idx_vw = edge_to_idx[(v, w)]
        row_idx_uw = edge_to_idx[(u, w)]
        row_idx_uv = edge_to_idx[(u, v)]
        
        # Populate the matrix according to the boundary formula
        push!(I, row_idx_vw); push!(J, col_idx); push!(V, 1)
        push!(I, row_idx_uw); push!(J, col_idx); push!(V, -1)
        push!(I, row_idx_uv); push!(J, col_idx); push!(V, 1)
    end

    D = sparse(I, J, V, total_dim, total_dim)

    # === 4. PREPARE RETURN VALUES ===
    simplex_indices = Dict(
        :vertices => 1:n_verts,
        :edges => (n_verts + 1):(n_verts + n_edges),
        :triangles => (n_verts + n_edges + 1):total_dim
    )

    return D, simplex_indices, n_verts + n_edges
end


function base_change(D::SparseMatrixCSC{Int64, Int64}, cocycle, field, multpl::Int = 1)
    α = spzeros(Int64,size(D,2))
    
    counter = 0
    for col in 1:size(D,2)
        #index = reverse(sort(unique(findnz(D[:,findnz(D[:,col])[1]])[1])))
        index = reverse(sort(findnz(D[:,col])[1]))
        if size(index,1) == 0
            continue
        end
        if size(index, 1) == 3
            break
        end
        
        for c ∈ cocycle
            if Ripserer.vertices(c) == (index[1], index[2])#, index[3])
                counter += 1
                coeff = Int64(coefficient(c) * multpl)
                #println(coeff)
                if coeff <= (field-1)/2
                    α[col,1] = coeff
                else
                    α[col,1] = coeff - field
                end
                #α[col,1] = coeff
            end  
            
        end
    end
    if length(cocycle) != counter
        println("ERRROR: length(cocycle) = ", length(cocycle),", counter = ", counter)
    end
    return α
end


function verify_coefficients_indexed(simps_ho::Vector{Tuple{Int64, Int64}}, coeff_ho::Vector{Int64}, prime::Int)
    # Ensure the number of simplices matches the number of coefficients.
    @assert length(simps_ho) == length(coeff_ho) "The number of simplices must match the number of coefficients."

    # Step 1: Calculate the list of coefficients and term count for each vertex.
    # vertex_coeffs_list stores a list of contributing coefficients for each vertex.
    vertex_coeffs_list = Dict{Int, Vector{Int}}()
    vertex_edge_counts = Dict{Int, Int}()

    for i in 1:length(simps_ho)
        v_start, v_end = simps_ho[i]
        coeff = coeff_ho[i]

        # For the boundary map ∂[v_start, v_end] = v_end - v_start:
        # - Add the coefficient to the end vertex's list.
        # - Add the negative of the coefficient to the start vertex's list.
        
        # Initialize lists if they don't exist
        if !haskey(vertex_coeffs_list, v_end)
            vertex_coeffs_list[v_end] = Int[]
        end
        if !haskey(vertex_coeffs_list, v_start)
            vertex_coeffs_list[v_start] = Int[]
        end

        push!(vertex_coeffs_list[v_end], coeff)
        push!(vertex_coeffs_list[v_start], -coeff)

        # Increment the edge count for both vertices.
        vertex_edge_counts[v_end] = get(vertex_edge_counts, v_end, 0) + 1
        vertex_edge_counts[v_start] = get(vertex_edge_counts, v_start, 0) + 1
    end

    # Create the dictionary of total coefficients from the lists for the return value
    vertex_coeffs_total = Dict{Int, Int}()
    for (vertex, coeffs_list) in vertex_coeffs_list
        vertex_coeffs_total[vertex] = sum(coeffs_list)
    end

    # Step 2: Verify the condition for each vertex.
    results = Dict{Int, Bool}()
    #println("--- Verification Results ---")
    #println("Prime (p): $prime\n")
    
    # Get all unique vertices involved to ensure we check all of them.
    all_vertices = sort(collect(union(keys(vertex_coeffs_list), keys(vertex_edge_counts))))

    for vertex in all_vertices
        coeffs_list = mod.(get(vertex_coeffs_list, vertex, Int[]), prime)
        n = get(vertex_edge_counts, vertex, 1) # Get the number of incident edges.
        
        # Ensure n is not zero to avoid division by zero error.
        if n == 0
            results[vertex] = isempty(coeffs_list) || (mod(sum(coeffs_list), prime) == 0)
            continue
        end

        # Calculate f = floor(p / n)
        f = floor(Int, prime / n)

        # Check if EVERY coefficient in the list lies in the allowed set.
        is_valid = true # Assume valid until a failure is found.
        individual_checks = String[]

        for c in coeffs_list
            c_mod_p = mod(c, prime)
            current_coeff_is_valid = (c_mod_p <= f) || (c_mod_p >= prime - f)
            push!(individual_checks, "$c mod $prime = $c_mod_p => $(current_coeff_is_valid ? "✓" : "✗")")
            if !current_coeff_is_valid
                is_valid = false # If any coefficient is invalid, the entire vertex is invalid.
            end
        end
        results[vertex] = is_valid
        if is_valid == false
        # Print detailed results for clarity
        println("Vertex: $vertex")
        println("  - Coefficients: {$(join(coeffs_list, ","))}")
        println("  - Incident Edges (n): $n")
        println("  - Floor(p/n): $f")
        println("  - Allowed Set: {0..$f} U {$(prime-f)..$(prime-1)}")
        #println("  - Individual Checks: [$(join(individual_checks, ", "))]")
        println("  - Is Valid?: $is_valid\n")
        end
    end

    #return results, vertex_coeffs_total
end