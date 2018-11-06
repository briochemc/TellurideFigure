using StatPlots, DataFrames, MultivariateStats, Compat, Random, Distributions, Plots

# Find observed μ and σ for FeT
μ_FeT = log10(100)
σ_FeT = 2
# Find observed μ and σ for DFe
μ_DFe = log10(1)
σ_DFe = 2
# Create distributions
dist_FeT = LogNormal(μ_FeT, σ_FeT)
dist_DFe = LogNormal(μ_DFe, σ_DFe)
# Sample from distributions and put in DataFrame
n = 10000
M = DataFrame()
M[:DFe] = rand(dist_DFe, n)
M[:FeT] = rand(dist_FeT, n)
M[:valid] = M[:DFe] .≤ M[:FeT]
head(M)

