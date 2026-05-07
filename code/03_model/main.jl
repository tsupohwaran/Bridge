#==================================================#
# Cross-Sea Bridge, Labor Market Power and Welfare      
# Date: December 3, 2025
# Author: Wu Chengjun, Central University of Finance and Economics
# OS: MacOS 26.2
# Version: 1.10.9
#==================================================#

#==================================================#
# Load packages and functions
# !!!Path must be redefined by users!!!
#==================================================#

projPath = "/Users/pohwaran/Doctorate/Paper/Bridge" # Project root path
stataPath = "/Applications/Stata/StataMP.app/Contents/MacOS/stata-mp" # Path to Stata executable
cd(projPath) 

include(projPath * "/code/00_setup/install_julia_pkgs.jl") # Install required packages
include(projPath * "/code/03_model/load_packages.jl") # Load required packages
include(projPath * "/code/03_model/functions.jl") # Load functions

RunStata(projPath, stataPath, "code/00_setup/install_stata_pkgs.do") # Install required Stata packages. Warning: May take some time if all packages are not installed.

#==================================================#
# Data Preparation
# Stata .do files are used for data cleaning and can be directly run from Julia
# .log files are generated to record the output of each .do file
# !!! Since some data are too large, we do not recommend re-running the entire data preparation process !!!
#==================================================#

# Preparing the Chinese industrial enterprises database (CIED)
RunStata(projPath, stataPath, "code/01_data_prep/01_prep_cied.do")

# Preparing the Chinese Tax Survey Database (CTSD)
RunStata(projPath, stataPath, "code/01_data_prep/02_prep_ctsd.do")

# Markdown estimation
RunStata(projPath, stataPath, "code/01_data_prep/04_estimate_markdown.do")

# Market access calculation
RunStata(projPath, stataPath, "code/01_data_prep/03_prep_market_access.do")

# Effect of the cross-sea bridge
RunStata(projPath, stataPath, "code/02_empirical/bridge_effect.do")

#==================================================#
# Load data
#==================================================#

df = DataFrame(load(projPath * "/data/model/processed/firm_qingdao_model.dta"));
l = df[!, :pop] |> x -> Float64.(x) |> x -> reshape(x, Z, J) |> x -> x[:, 1] |> x -> x ./ sum(x); # normalize total population to 1
lⱼ_data = df[!, :employ] |> x -> Float64.(x) |> x -> reshape(x, Z, J) |> x -> x[1, :] |> x -> x ./ sum(x);
d = reshape(Float64.(df[!, :dzj]), Z, J) |> x -> replace(x, 0.0 => 1e-2);
d′ = reshape(Float64.(df[!, :dzj_prime]), Z, J) |> x -> replace(x, 0.0 => 1e-2);
wⱼ_data = reshape(Float64.(df[!, :wage_inital]), Z, J)[1, :] |> x -> x ./ sum(x .* lⱼ_data); # normalize total wage bill to 1

# zⱼ = reshape(Float64.(df[!, :z]), Z, J)[1, :] |> x -> clamp.(x, quantile(x, 0.05), quantile(x, 0.95)) |> x -> x ./ mean(x) # normalize zⱼ to have mean 1 (winsor 5% at both ends)
# zⱼ = reshape(Float64.(df[!, :z]), Z, J)[1, :] |> x -> x ./ mean(x)

#==================================================#
# Calibration: Back out η and θ from β₁ and β₂
#==================================================#

α = 0.4
β_target = [-0.052, 0.085]
x0 = [2.75, 2.0]

result = optimize(
    x -> ObjectiveFunction(x; l, d, d′, wⱼ_data, lⱼ_data, α, β_target, verbose=true),
    x0,
    NelderMead(),
    Optim.Options(show_trace = true, iterations = 100, g_tol = 1e-8)
)

println("\n" * "="^50)
println("OPTIMIZATION RESULTS")
println("="^50)
η_est, θ_est = Optim.minimizer(result)
println("Estimated η: ", η_est)
println("Estimated θ: ", θ_est)
println("Final objective: ", Optim.minimum(result))
println("Converged: ", Optim.converged(result))

# Verify final moments
β_final = ComputeModelMoments([η_est, θ_est]; l, d, d′, wⱼ_data, lⱼ_data, α)
println("\nTarget  β: ", β_target)
println("Model   β: ", β_final)

#==================================================#
# Simulation given η = 1; θ = 5;
#==================================================#

α = 0.4
η = 3.1; # commute elasticity
θ = 2.14; # shape parameter of Fréchet distribution

# Solve the zⱼ from observed wⱼ
vars = (; wⱼ = wⱼ_data, l, d)
params = (; η, θ, α)
zⱼ = SolveZfromW(vars, params);

# Solve the model
vars = (; l, d, zⱼ)
params = (; η, θ, α)

# solve the model for baseline
# Use observed wages as the warm start and avoid the power update here:
# zⱼ was inverted from wⱼ_data, so this keeps the solver on the same equilibrium branch.
wⱼ, π_zj, ε_zj, lⱼ, εⱼ = SolveModel(vars, params; displayGap = true, damp = 0.6, tol = 1e-7, displaySummary = true, power = false, wⱼ_init = wⱼ_data);

# Calculate correlation between solved wages and observed wages
corr_wages = cor(vec(wⱼ), vec(wⱼ_data))
println("Correlation between wⱼ_solved and wⱼ_data: ", corr_wages) # should be very close to 1

# solve the model for counterfactual
vars′ = (; l, d = d′, zⱼ);
wⱼ′, π_zj′, ε_zj′, lⱼ′, εⱼ′ = SolveModel(vars′, params; displayGap = false, damp = 0.7, tol = 1e-7, displaySummary = true, power = false, wⱼ_init = wⱼ);

# Export wⱼ to Excel
filepath = projPath * "/output/tables/wages.xlsx"
isfile(filepath) && rm(filepath)
XLSX.writetable(filepath, DataFrame(wj = vec(wⱼ)))

# Analyze the results
l̂ⱼ = lⱼ′ ./ lⱼ;
dlnlⱼ = l̂ⱼ .- 1;
dlnMA = log.(sum((d - d′) .* l, dims = 1)');
dlnMA = replace(dlnMA, Inf => -8.0, -Inf => -8.0)
density(dlnMA, title="Kernel Density Estimate of dlnMA", xlabel="dlnMA", ylabel="Density", legend=false)

dlnMA = log.(sum((d - d′) .* l, dims = 1)') |> x -> replace(x, -Inf => -8)
bigMA = Float64.(dlnMA .>= -1)

regDF = DataFrame(bigMA = vec(bigMA), dlnl = vec(dlnlⱼ), w = vec(log.(wⱼ)))
regDF.w_treatedmean = fill(mean(regDF[regDF.bigMA .== 1, :w]), nrow(regDF))
regDF.w_diff = regDF.w .- regDF.w_treatedmean
regModel = lm(@formula(dlnl ~ bigMA + bigMA & w_diff + w_diff), regDF)

β₁ = coef(regModel)[2]
β₂ = coef(regModel)[4]

#==================================================#
# Plot
#==================================================#

# density(vec(wⱼ), title="Kernel Density Estimate of wⱼ", xlabel="wⱼ", ylabel="Density", legend=false)
# labor
begin
    sorted_idx = sortperm(vec(wⱼ))
    scatter(vec(dlnMA)[sorted_idx], vec(clamp.(dlnlⱼ, -Inf, 5))[sorted_idx], 
        marker_z = vec(wⱼ)[sorted_idx],
        color = :RdBu,
        colorbar = true,
        colorbar_title = "wⱼ",
        xlabel = "dlnMA", 
        ylabel = "dlnlⱼ",
        title = "Employment Change vs Market Access Change",
        legend = false,
        markersize = 3,
        alpha = 0.6,
        dpi = 1000)
end
savefig(projPath * "/output/figures/model/labor_change.png")

scatter(vec(dlnMA), vec(log.(lⱼ)), 
    marker_z = vec(wⱼ),
    color = :RdBu,
    colorbar = true,
    colorbar_title = "wⱼ",
    xlabel = "log of dlnMA", 
    ylabel = "log of lⱼ",
    title = "Employment vs Market Access Change",
    legend = false,
    markersize = 3,
    alpha = 0.6,
    dpi = 1000)
savefig(projPath * "/output/figures/model/labor.png")

# wages change
ŵⱼ = wⱼ′ ./ wⱼ;
dlnwⱼ = ŵⱼ .- 1;
scatter(vec(dlnMA), vec(clamp.(dlnwⱼ, -Inf, 0.5)), 
    marker_z = vec(wⱼ),
    color = :RdBu,
    colorbar = true,
    colorbar_title = "wⱼ",
    xlabel = "dlnMA", 
    ylabel = "dlnwⱼ",
    title = "Wage Change vs Market Access Change",
    legend = false,
    markersize = 3,
    alpha = 0.6,
    dpi = 1000)
savefig(projPath * "/output/figures/model/wage_change.png")



# markdown
ν = 1 .+ 1 ./ εⱼ
ν′ = 1 .+ 1 ./ εⱼ′
ν̂ = ν′ ./ ν
dlnν = clamp.(ν̂ .- 1, -0.001, 0.001)
dlnν = ν̂ .- 1
scatter(vec(dlnMA), vec(dlnν), 
    marker_z = vec(wⱼ),
    color = :RdBu,
    colorbar = true,
    colorbar_title = "wⱼ",
    xlabel = "log of dlnMA", 
    ylabel = "dlnν",
    title = "Labor Market Power Change vs Market Access Change",
    legend = false,
    markersize = 3,
    alpha = 0.6,
    dpi = 1000)
savefig(projPath * "/output/figures/model/markdown_change.png")

ν = clamp.(ν, 1, 1.3)
scatter(vec(dlnMA), vec(ν), 
    marker_z = vec(wⱼ),
    color = :RdBu,
    colorbar = true,
    colorbar_title = "wⱼ",
    xlabel = "log of dlnMA", 
    ylabel = "νⱼ",
    title = "Labor Market Power vs Market Access Change",
    legend = false,
    markersize = 3,
    alpha = 0.6,
    dpi = 1000)
savefig(projPath * "/output/figures/model/markdown.png")
