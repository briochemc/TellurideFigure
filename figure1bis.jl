using StatsPlots
using DataFrames


using Compat, Random, Distributions
Random.seed!(123)

using Plots
μ = 10
σ = 2
dist_FeT = LogNormal(μ, σ)  # Find correct values for μ and σ
sample_FeT = rand(dist_FeT, 10000)
histogram(log10.(sample_FeT), bin = 100)
