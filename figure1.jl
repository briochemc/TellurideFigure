


function make_uniform_FeT_and_DFe(;n = 10000, max_FeT = 100, max_DFe = 3)
    DFe = max_DFe * rand(n)
    FeT = max_FeT * rand(n)
    ivalid = findall(DFe .≤ FeT)
    DFe, FeT = DFe[ivalid], FeT[ivalid]
    return DFe, FeT
end

function make_lognormal_FeT_and_DFe(;n = 10000, m_FeT = 50  , v_FeT = 800,
                                               m_DFe =  1.5, v_DFe =  0.75)
    μ_FeT = log(m_FeT / sqrt(1 + v_FeT / m_FeT^2))
    σ_FeT = sqrt(log(1 + v_FeT / m_FeT^2))
    μ_DFe = log(m_DFe / sqrt(1 + v_DFe / m_DFe^2))
    σ_DFe = sqrt(log(1 + v_DFe / m_DFe^2))
    DFe = exp.(μ_DFe .+ σ_DFe * randn(n))
    FeT = exp.(μ_FeT .+ σ_FeT * randn(n))
    ivalid = findall(DFe .≤ FeT)
    DFe, FeT = DFe[ivalid], FeT[ivalid]
    return DFe, FeT
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

function logreg(logx, logy)
    # If logy = a logx + b
    # then logy = A [a b]ᵀ with A = [logx ones(length(logx))]
    A = [logx ones(length(logx))]
    # so
    ab = A \ logy
    a, b = ab[1], ab[2]
    return a, b
end

using Plots
gr(size = (600, 600))
n = 10000 # number of samples

# Plot as N. Meskhidze
DFe, FeT = make_uniform_FeT_and_DFe()
logy, logx = log10.(DFe ./ FeT), log10.(FeT)
plt = myloglog(logx, logy)
a, b = logreg(logx, logy)

# Plot using lognormal dist.
DFe, FeT = make_lognormal_FeT_and_DFe()
logy, logx = log10.(DFe ./ FeT), log10.(FeT)
plt = myloglog(logx, logy)
a, b = logreg(logx, logy)

# Plot using lognormal dist. but such that a ≠ -1
DFe, FeT = make_lognormal_FeT_and_DFe(n=10000, m_DFe = 10, v_DFe = 800)
logy, logx = log10.(DFe ./ FeT), log10.(FeT)
plt = myloglog(logx, logy)
a, b = logreg(logx, logy)