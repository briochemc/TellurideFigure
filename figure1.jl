function make_uniform_DFe_and_FeT(;n = 10000, max_FeT = 100, max_DFe = 3)
    DFe = max_DFe * rand(n)
    FeT = max_FeT * rand(n)
    ivalid = findall(DFe .≤ FeT)
    DFe, FeT = DFe[ivalid], FeT[ivalid]
    return DFe, FeT
end

function make_lognormal(;n = 10000, m = 1, v = 1)
    # convert arithmetic mean and variance to log mean and log std
    μ = log(m / sqrt(1 + v / m^2))
    σ = sqrt(log(1 + v / m^2))
    return exp.(μ .+ σ * randn(n))
end

function make_lognormal_DFe_and_FeT(;n = 10000, m_FeT = 50  , v_FeT = 800,
                                               m_DFe =  1.5, v_DFe =  0.75)
    FeT = make_lognormal(n = n, m = m_FeT, v = v_FeT)
    DFe = make_lognormal(n = n, m = m_DFe, v = v_DFe)
    ivalid = findall(DFe .≤ FeT)
    DFe, FeT = DFe[ivalid], FeT[ivalid]
    return DFe, FeT
end

function make_lognormal_DFe_and_PFe(;n = 10000, m_PFe = 50  , v_PFe = 800,
                                                 m_DFe =  1.5, v_DFe =  0.75)
    PFe = make_lognormal(n = n, m = m_PFe, v = v_PFe)
    DFe = make_lognormal(n = n, m = m_DFe, v = v_DFe)
    return DFe, PFe
end

function logreg(logx, logy)
    # If logy = a logx + b
    # then logy = A [a b]ᵀ with A = [logx ones(length(logx))]
    A = [logx ones(length(logx))]
    # so
    ab = A \ logy
    a, b = ab[1], ab[2]
    return a, b
end

function myloglog(logx, logy)
    x_lim = (floor(minimum(logx)), ceil(maximum(logx)))
    y_lim = (floor(minimum(logy)), ceil(maximum(logy)))

    # make nice labels
    y_ticks = collect(y_lim[1]:y_lim[2])
    no_y_tick_labels = ["" for y in y_ticks]
    y_tick_labels = ["10^{$(convert(Int,y))}" for y in y_ticks]
    x_ticks = collect(x_lim[1]:x_lim[2])
    no_x_tick_labels = ["" for x in x_ticks]
    x_tick_labels = ["10^{$(convert(Int,x))}" for x in x_ticks]

    # main scatter plot
    p21 = scatter(logx, logy,
        marker=:o, markeralpha = 0.1, markersize = 1.0,
        xaxis=("FeT"), yaxis=("DFe / FeT"), linewidth=0,
        xlim = x_lim, ylim = y_lim, legend=nothing,
        xticks = (x_ticks, x_tick_labels), yticks = (y_ticks, y_tick_labels))

    a, b = logreg(logx, logy)
    p21 = plot!(p21, collect(x_lim), a * collect(x_lim) .+ b)

    # side and top histograms
    p22 = histogram(logy, bin = 100, orientation=:horizontal,
        yticks = (y_ticks, no_y_tick_labels), ylim = y_lim, legend=nothing,
        xaxis=("# of samples"))
    p11 = histogram(logx, bin = 100,
        xticks = (x_ticks, no_x_tick_labels), xlim = x_lim, legend=nothing,
        yaxis=("# of samples"))

    # empty plot for top right
    p12 = plot(legend=false,axis=false,grid=false,foreground_color_subplot=:white)

    # plot all of them
    return plot(p11, p12, p21, p22, link=:all, layout=(2,2))
end


using Plots
gr(size = (600, 600))
n = 10000 # number of samples

# Plot as N. Meskhidze
DFe, FeT = make_uniform_DFe_and_FeT()
logy, logx = log10.(DFe ./ FeT), log10.(FeT)
plt = myloglog(logx, logy)
savefig("dissolved_Fe_fraction_from_uniform_DFe_and_FeT.svg")

# Plot using lognormal dist.
DFe, FeT = make_lognormal_DFe_and_FeT()
logy, logx = log10.(DFe ./ FeT), log10.(FeT)
plt = myloglog(logx, logy)
savefig("dissolved_Fe_fraction_from_lognormal_DFe_and_FeT.svg")

# Plot using lognormal dist. but such that a ≠ -1
DFe, FeT = make_lognormal_DFe_and_FeT(n=10000, m_DFe = 10, v_DFe = 800)
logy, logx = log10.(DFe ./ FeT), log10.(FeT)
plt = myloglog(logx, logy)
savefig("dissolved_Fe_fraction_from_lognormal_DFe_and_FeT2.svg")

# Plot using lognormal dist. but by generating PFe rather than FeT
DFe, PFe = make_lognormal_DFe_and_PFe()
logy, logx = log10.(DFe ./ (DFe .+ PFe)), log10.(DFe .+ PFe)
plt = myloglog(logx, logy)
savefig("dissolved_Fe_fraction_from_lognormal_DFe_and_PFe.svg")

