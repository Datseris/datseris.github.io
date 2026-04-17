using Pkg; Pkg.activate(@__DIR__)
ENV["BGCOLOR"] = "#1e1e20"
ENV["AXISCOLOR"] = :white
ENV["COLORSCHEME"] = "JuliaDynamicsLight"
using DynamicalSystems, CairoMakie, MakieForProjects
set_theme!(make_theme())

colors = COLORS

using Statistics: mean, std
display_timeseries(X::StateSpaceSet, a...; kw...) = display_timeseries(collect(columns(X)), a...; kw...)

function display_timeseries(X, idxs = 1:length(X); decay = 1e-4, kw...)
    fig = Figure(size = (1920, 1080), figure_padding = 0)
    ax = Axis(fig[1,1])
    t = 0:length(first(X))-1
    e = @. exp(-decay*t)
    X = [e .* (x .- mean(x)) for x in X]
    colors = [
        "#855DE4",
        lighten("#855DE4", 0.8),
        "#149eaa",
        lighten("#149eaa", 1.8),
        "#B7BEF1",
        lighten("#B7BEF1", 0.8),
    ]
    for (i, d) in enumerate(idxs)
        x = X[i]
        lines!(ax, t, x; linestyle = :solid, linewidth = 3, color = colors[d], alpha = 0.25, kw...)
    end
    # normalize limits
    ma = maximum(maximum(x) for x in X)
    mi = minimum(minimum(x) for x in X)
    d = ma-mi
    xlims!(ax, 0, t[end])
    ylims!(ax, mi - d/10, ma + 2d)
    hidedecorations!(ax)
    hidespines!(ax)
    display(fig)
    return fig
end

# %% bistable roessler

function roessler_rule(u, p, t)
    x, y, z = u
    a, b, c = p
    dx = -y - z
    dy = x + a*y
    dz = b + z*(x - c)
    return SVector{3}(dx, dy, dz)
end

p0 = [0.29, 0.14, 4.52]
u0s = [[-1.25, -0.72, -0.1], [0.72, 1.28, 0.21]]
# run the above to get states on the attractors
u0s = SVector{3, Float64}[[-3.1167416198371725, -9.384166559175227, 0.016333574889970465], [5.813245732922099, -8.655686826767814, 0.08606051196041611]]
# push!(u0s, [2, 0.1, 5.0])
ds = CoupledODEs(roessler_rule, u0s[1], p0)

X = []

for u0 in u0s
    y, t = trajectory(ds, 100, u0; Ttr = 100, Δt = 0.01)
    push!(X, y[:, 1], y[:, 2])
end

fig = display_timeseries(X; decay = 5e-4)

save(joinpath(@__DIR__, "roessler.png"), fig)

# %% chaotic reversals

ds = Systems.gissinger()
u0s = [[1, 1, 5.0],  [2.0, -2.0, 2.0], [-2.0, -2.0, 4.0],  [-2.0, 2.0, 4.0],]

X = []

for u0 in u0s
    y, t = trajectory(ds, 1000, u0; Ttr = 100, Δt = 0.1)
    push!(X, y[:, 1],  )
end

fig = display_timeseries(X; decay = 4e-4, linewidth = 2)

save(joinpath(@__DIR__, "gissinger.png"), fig)

# %% oscillatory arma
# generate autoregressive process
using Random
rng = Random.MersenneTwister(77163)
η = randn(rng, 5000)
s = ones(5000)
X = []
for _ in 1:3
    η = randn(5000)
    for n in 4:5000
        s[n] = 1.625s[n-1] - 0.284s[n-2] - 0.355s[n-3] + η[n] - 0.96η[n-1]
    end
    push!(X, copy(s))
end

fig = display_timeseries(X; decay = 5e-4)

save(joinpath(@__DIR__, "arma.png"), fig)


# %% Lorenz96

ds = Systems.lorenz96(6; F = 31.0)

X, t = trajectory(ds, 10; Ttr = 156, Δt = 0.01)
X = collect(columns(X))

fig = display_timeseries(X; decay = 3e-3, linewidth = 2, alpha  = 0.5)

save(joinpath(@__DIR__, "lorenz96.png"), fig)




# %% Ceres

using DelimitedFiles

X = readdlm(joinpath(@__DIR__, "11.csv"))

fig = display_timeseries(eachcol(X), [1, 3]; decay = 8e-3, linewidth = 3)

save(joinpath(@__DIR__, "ceres.png"), fig)

