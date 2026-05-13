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

# load model sample and normalize the variables
Z = 128; # number of towns
df = DataFrame(load(projPath * "/data/model/processed/firm_qingdao_model_10.dta"));
J = nrow(df) ÷ Z; # number of firms

l = df[!, :pop] |> x -> Float64.(x) |> 
    x -> reshape(x, Z, J) |> 
    x -> x[:, 1] |> 
    x -> x ./ sum(x); # normalize total population to 1
lⱼ_data = Float64.(disallowmissing(df[!, :employ])) |> 
    x -> reshape(x, Z, J) |> 
    x -> x[1, :] |> 
    x -> x ./ sum(x); # normalize total employment to 1
d = reshape(Float64.(df[!, :dzj]), Z, J) |> x -> replace(x, 0.0 => 1e-2);
d′ = reshape(Float64.(df[!, :dzj_prime]), Z, J) |> x -> replace(x, 0.0 => 1e-2);
wⱼ_data = reshape(Float64.(df[!, :wage_inital]), Z, J)[1, :] |> 
    x -> x ./ sum(x .* lⱼ_data); # normalize total wage bill to 1

# load regression sample
df_reg = DataFrame(load(projPath * "/data/model/processed/firm_two_year_reg.dta"));
firm_ids = collect(df[1:Z:end, :id]);
reg_sample_ids = collect(skipmissing(df_reg[!, :id]));
restrict_to_reg_sample = true; # set false to use all model firms for model moments
moment_firm_mask = MomentFirmMask(J; firm_ids, reg_sample_ids, restrict_to_reg_sample);
println("Model moment sample firms: ", sum(moment_firm_mask), "/", J,
    restrict_to_reg_sample ? " matched to firm_two_year_reg.dta" : " full model sample");

#==================================================#
# Calibration: Back out η and θ from β₁ and β₂
#==================================================#

α = 0.4
β_target = [-0.092242, 0.0939565]
η_bounds = [0.25, 8.0]
θ_bounds = [0.25, 8.0]
employment_change = :log
wage_center = :all

# use grid search to find good starting points for the optimization
# η_grid = [1.5, 3.0, 4.5, 6.0, 7.5]
# θ_grid = [0.75, 1.5, 2.5, 3.5, 4.5]
# grid_results = EvaluateCalibrationGrid(;
#     l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
#     η_grid, θ_grid,
#     inner_tol = 2e-5,
#     inner_maxIter = 3000,
#     continuation_steps = 5,
#     employment_change,
#     wage_center,
#     moment_firm_mask,
#     max_abs_moment = 10.0,
#     verbose = true
# )

# mkpath(projPath * "/output/tables")
# CSV.write(projPath * "/output/tables/calibration_grid.csv", grid_results)

# top_grid = first(grid_results, min(4, nrow(grid_results)))
# println("\nTop calibration grid points:")
# show(top_grid, allrows = true, allcols = true)
# println()

# calibration_starts = [[row.η, row.θ] for row in eachrow(top_grid)]
# pushfirst!(calibration_starts, x0)

# based on the grid search results, we choose a good starting point for the optimization
x0 = [1, 4.8]
calibration = CalibrateEtaTheta(;
    l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
    x0,
    starts = [x0],
    lower = [η_bounds[1], θ_bounds[1]],
    upper = [η_bounds[2], θ_bounds[2]],
    iterations = 250,
    x_abstol = 1e-4,
    f_reltol = 1e-8,
    inner_tol = 1e-5,
    inner_maxIter = 3000,
    continuation_steps = 5,
    employment_change,
    wage_center,
    moment_firm_mask,
    max_abs_moment = 10.0,
    verbose = true,
    show_trace = true
)
result = calibration.result

println("\n" * "="^50)
println("OPTIMIZATION RESULTS")
println("="^50)
η_est, θ_est = calibration.parameters

println("Final objective: ", calibration.objective)
println("Converged: ", calibration.converged)
println("Bounds: η ∈ ", η_bounds, ", θ ∈ ", θ_bounds)

# Verify final moments
β_final = ComputeModelMoments([η_est, θ_est]; l, d, d′, wⱼ_data, lⱼ_data, α,
    inner_tol = 1e-5, inner_maxIter = 5000, inner_display = true,
    continuation_steps = 5, employment_change, wage_center,
    moment_firm_mask)
println("\nTarget  β: ", β_target)
println("Model   β: ", β_final)
println("Estimated η: ", η_est)
println("Estimated θ: ", θ_est)

#==================================================#
# Simulation using calibrated η and θ
#==================================================#

α = 0.4
η = η_est; # commute-wage elasticity
θ = θ_est; # 

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

## Calculate correlation between solved wages and observed wages
# corr_wages = cor(vec(wⱼ), vec(wⱼ_data))
# println("Correlation between wⱼ_solved and wⱼ_data: ", corr_wages) # should be very close to 1

# solve the model for counterfactual
vars′ = (; l, d = d′, zⱼ);
wⱼ′, π_zj′, ε_zj′, lⱼ′, εⱼ′ = SolveModel(vars′, params; displayGap = false, damp = 0.6, tol = 1e-9, displaySummary = true, power = true, wⱼ_init = wⱼ);



# Analyze the results
l̂ⱼ = lⱼ′ ./ lⱼ;
dlnlⱼ = log.(max.(lⱼ′, eps(Float64))) .- log.(max.(lⱼ, eps(Float64)));
dMA = sum((d - d′) .* l, dims = 1)' |> x -> replace(x, -Inf => -8);
bigMA = Float64.(dMA .>= 0.5);

regDF = DataFrame(bigMA = vec(bigMA), dlnl = vec(dlnlⱼ), w = vec(log.(wⱼ)));
regDF = regDF[moment_firm_mask, :];
regDF.w_treatedmean = fill(mean(regDF[!, :w]), nrow(regDF));

regDF.w_diff = regDF.w .- regDF.w_treatedmean;
regModel = lm(@formula(dlnl ~ bigMA + bigMA & w_diff + w_diff), regDF)

β = [coef(regModel)[2], coef(regModel)[4]]

density(dMA, title="Kernel Density Estimate of dMA", xlabel="dMA", ylabel="Density", legend=false)

#==================================================#
# Plot
#==================================================#

# density(vec(wⱼ), title="Kernel Density Estimate of wⱼ", xlabel="wⱼ", ylabel="Density", legend=false)
# labor
begin
    sorted_idx = sortperm(vec(wⱼ))
    scatter(vec(dMA)[sorted_idx], vec(clamp.(dlnlⱼ, -Inf, 5))[sorted_idx], 
        marker_z = vec(wⱼ)[sorted_idx],
        color = :RdBu,
        colorbar = true,
        colorbar_title = "wⱼ",
        xlabel = "dMA", 
        ylabel = "dlnlⱼ",
        title = "Employment Change vs Market Access Change",
        legend = false,
        markersize = 3,
        alpha = 0.6,
        dpi = 1000)
end
savefig(projPath * "/output/figures/model/labor_change.png")

scatter(vec(dMA), vec(log.(lⱼ)), 
    marker_z = vec(wⱼ),
    color = :RdBu,
    colorbar = true,
    colorbar_title = "wⱼ",
    xlabel = "log of dMA", 
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
scatter(vec(dMA), vec(clamp.(dlnwⱼ, -Inf, 0.5)), 
    marker_z = vec(wⱼ),
    color = :RdBu,
    colorbar = true,
    colorbar_title = "wⱼ",
    xlabel = "dMA", 
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
scatter(vec(dMA), vec(dlnν), 
    marker_z = vec(wⱼ),
    color = :RdBu,
    colorbar = true,
    colorbar_title = "wⱼ",
    xlabel = "log of dMA", 
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
