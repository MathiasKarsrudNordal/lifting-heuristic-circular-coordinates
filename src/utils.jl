module Utils

using CairoMakie
using GeometryBasics: Point2f
using LaTeXStrings

export plot_persistence_diagram

function plot_persistence_diagram(result; infinity=Inf, palette=Makie.wong_colors())
    fig = Figure(resolution = (600, 600), fontsize = 14, figure_padding = 75)
    ax = Axis(fig[1, 1];
        titlesize = 28,   # << increase title font size
        title = L"\text{Persistence Diagram}",
        xlabelsize = 24,
        xlabel = L"\text{Birth}",
        ylabelsize = 24,
        ylabel = L"\text{Death}",
        aspect = 1
    )

    # Gather global min/max over all dims for square limits
    births_all = Float64[]
    deaths_all = Float64[]
    for diag in result
        for p in diag
            push!(births_all, p.birth)
            push!(deaths_all, isinf(p.death) ? infinity : p.death)
        end
    end
    if isempty(births_all)
        lo, hi = 0.0, 1.0
    else
        lo = min(minimum(births_all), minimum(deaths_all))
        hi = max(maximum(births_all), maximum(deaths_all))
        # Expand limits a bit for better visualization
        lo -= 0.1 * hi
        hi += 0.1 * hi
    end

    # Plot each homology dimension with its own color + legend label
    for (dim, diag) in enumerate(result)   # dim=1→H₀, dim=2→H₁, ...
        isempty(diag) && continue
        births = [p.birth for p in diag]
        deaths = [isinf(p.death) ? infinity : p.death for p in diag]
        scatter!(ax, births, deaths;
            markersize = 15,
            color = palette[(dim - 1) % length(palette) + 1],
            label = latexstring("H_{$(dim - 1)}")
        )
    end

    # Highlight the last H1 point if it exists
    if length(result) ≥ 2 && !isempty(result[2])
        last_point = result[2][end]
        b = last_point.birth
        d = isinf(last_point.death) ? infinity : last_point.death
        scatter!(ax, [b], [d];
            markersize = 30,
            color = :transparent,
            strokecolor = :black,
            strokewidth = 2,
        )
    end

    # Diagonal y = x (use a vector of Point2f to avoid the tuple-of-tuples error)
    lines!(ax, [Point2f(lo, lo), Point2f(hi, hi)];
        color = :gray, linestyle = :dash)
    lines!(ax, [Point2f(lo, infinity), Point2f(hi, infinity)];
        color = :gray, linestyle = :dash)

    # Square limits and legend
    xlims!(ax, lo, hi)
    ylims!(ax, lo, hi)
    axislegend(ax, position = :rb, labelsize=24)

    ax.xticksvisible = false       # hide the ticks
    ax.xticklabelsvisible = false  # hide the tick labels
    ax.yticksvisible = false
    ax.yticklabelsvisible = false
    return fig
end

end # module Utils


