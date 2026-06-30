using LinearAlgebra, Statistics, DataFrames, FixedEffectModels

include("functions.jl")

bigMA = [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0]
w = [-1.2, -0.4, -0.1, 0.6, 1.1, -0.8, 0.2, 1.4]
w_diff = w .- mean(w)

expected = [0.3, -0.4, 0.2]
Δlnl = 0.1 .+ expected[3] .* w_diff .+ expected[1] .* bigMA .+
    expected[2] .* bigMA .* w_diff
Δlnw = zeros(length(w))

β = EstimateTwoPeriodDIDMoments(
    zeros(length(w)),
    Δlnl,
    zeros(length(w)),
    Δlnw,
    bigMA,
    w;
    wage_center = :all
)

@assert maximum(abs.(β .- expected)) < 1e-10
println("Synthetic two-period moment test passed: ", β)
