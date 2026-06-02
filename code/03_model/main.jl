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
include("load_packages.jl") # Load required packages
include("functions.jl") # Load functions

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
wⱼ_raw = reshape(Float64.(df[!, :wage_inital]), Z, J)[1, :];
lo, hi = quantile(wⱼ_raw, [0.05, 0.95])
wⱼ_data = clamp.(wⱼ_raw, lo, hi) |>
    x -> x ./ sum(x .* lⱼ_data); # normalize total wage bill to 1
aⱼ_init = zeros(J); # initial guess; firm amenities are inverted from observed employment

# load regression sample
df_reg = DataFrame(load(projPath * "/data/model/processed/firm_two_year_reg.dta"));
firm_ids = collect(df[1:Z:end, :id]);
reg_sample_ids = collect(skipmissing(df_reg[!, :id]));
restrict_to_reg_sample = false; # set false to use all model firms for model moments
moment_firm_mask = MomentFirmMask(J; firm_ids, reg_sample_ids, restrict_to_reg_sample);
moment_sample_label = restrict_to_reg_sample ? "matched regression sample" : "full model sample";
println("Model moment sample firms: ", sum(moment_firm_mask), "/", J,
    " in ", moment_sample_label);

#==================================================#
# Calibration: Back out η and θ from four regression moments
#==================================================#

α = 0.4
# From code/02_empirical/calculate_calibration_4_moments.do:
# regression 1: BIG#post, BIG#post#c.demean_lnw0
# regression 2 (BIG == 1): post#c.lndma, post#c.demean_lnw0#c.lndma
β_target = [-0.0912855, 0.1068467, -0.0302357, 0.0512688]
η_bounds = [0.25, 3.0]
θ_bounds = [0.25, 8.0]
employment_change = :log
wage_center = :all

# use grid search to find good starting points for the optimization.
# η is the elasticity of the wage compensation needed for a commuting-time
# change, so large values are hard to interpret economically. The constrained
# grid keeps η ≤ 3. With the corrected empirical moments, the best tested
# plausible-η starts are at high θ, so include that edge explicitly.
η_grid = [1.0, 1.5, 2.0, 2.5, 3.0]
θ_grid = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0]
grid_results = EvaluateCalibrationGrid(;
    l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
    η_grid, θ_grid,
    inner_tol = 2e-5,
    inner_maxIter = 3000,
    continuation_steps = 5,
    employment_change,
    wage_center,
    moment_firm_mask = moment_firm_mask,
    max_abs_moment = 10.0,
    verbose = true
)

mkpath(projPath * "/output/tables")
CSV.write(projPath * "/output/tables/calibration_grid.csv", grid_results)

grid_results = CSV.read(projPath * "/output/tables/calibration_grid.csv", DataFrame)
top_grid = first(grid_results, min(4, nrow(grid_results)))
println("\nTop calibration grid points:")
show(top_grid, allrows = true, allcols = true)
println()

calibration_starts = [[row.η, row.θ] for row in eachrow(top_grid)]

# Add the best plausible-η probe as a fallback start in case it is not among
# the top grid rows after rerunning with slightly different tolerances.
x0 = [3.0, 8.0]
if !any(start -> isapprox(start[1], x0[1]) && isapprox(start[2], x0[2]), calibration_starts)
    push!(calibration_starts, x0)
end
calibration = CalibrateEtaTheta(;
    l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
    aⱼ_init,
    starts = calibration_starts,
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
    moment_firm_mask = moment_firm_mask,
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
η_est, θ_est = [1.0, 1.0]
β_final = ComputeModelMoments([η_est, θ_est]; l, d, d′, wⱼ_data, lⱼ_data, α,
    aⱼ_init,
    inner_tol = 1e-5, inner_maxIter = 5000, inner_display = true,
    continuation_steps = 5, employment_change, wage_center,
    moment_firm_mask = moment_firm_mask)
println("\nTarget  β: ", β_target)
println("Model   β (", moment_sample_label, "): ", β_final)
println("Estimated η: ", η_est)
println("Estimated θ: ", θ_est)

#==================================================#
# Simulation using calibrated η and θ
#==================================================#

α = 0.4
η = 8.0; # commute-wage elasticity
θ = 0.7560388103691049; # 

# Solve firm amenities and productivity from observed employment and wages
vars = (; wⱼ = wⱼ_data, lⱼ = lⱼ_data, l, d)
params = (; η, θ, α)
primitives = SolveFirmPrimitivesFromData(vars, params; aⱼ_init, displaySummary = true);
zⱼ, aⱼ = primitives.zⱼ, primitives.aⱼ;

# Solve the model
vars = (; l, d, zⱼ, aⱼ)
params = (; η, θ, α)

# solve the model for baseline
# Use observed wages as the warm start and avoid the power update here:
# zⱼ was inverted from wⱼ_data, so this keeps the solver on the same equilibrium branch.
wⱼ, π_zj, ε_zj, lⱼ, εⱼ = SolveModel(vars, params; displayGap = true, damp = 0.6, tol = 1e-7, displaySummary = true, power = false, wⱼ_init = wⱼ_data);

## Calculate correlation between solved wages and observed wages
# corr_emp = cor(vec(lⱼ), vec(lⱼ_data))
# corr_wages = cor(vec(wⱼ), vec(wⱼ_data))
# println("Correlation between wⱼ_solved and wⱼ_data: ", corr_wages) # should be very close to 1

# solve the model for counterfactual
vars′ = (; l, d = d′, zⱼ, aⱼ);
wⱼ′, π_zj′, ε_zj′, lⱼ′, εⱼ′ = SolveModel(vars′, params; displayGap = false, damp = 0.6, tol = 1e-9, displaySummary = true, power = true, wⱼ_init = wⱼ);



# Analyze the results
l̂ⱼ = lⱼ′ ./ lⱼ;
dlnlⱼ = log.(max.(lⱼ′, eps(Float64))) .- log.(max.(lⱼ, eps(Float64)));
dMA = vec(sum((d - d′) .* l, dims = 1));
# density(dMA, title="Kernel Density Estimate of dMA", xlabel="dMA", ylabel="Density", legend=false)
bigMA = Float64.(dMA .> 0.5);
ln_dMA = MomentLogDMA(dMA, bigMA)

regDF = DataFrame(
    bigMA = vec(bigMA)[moment_firm_mask],
    dlnl = vec(dlnlⱼ)[moment_firm_mask],
    w = vec(log.(max.(wⱼ, eps(Float64))))[moment_firm_mask],
    dMA = vec(dMA)[moment_firm_mask]
);

regDF.lndma = fill(NaN, nrow(regDF));
treated_lndma = (regDF.bigMA .== 1) .& (regDF.dMA .> 0);
regDF.lndma[treated_lndma] .= log.(regDF.dMA[treated_lndma]);
regDF.w_center = fill(mean(regDF[!, :w]), nrow(regDF));
regDF.w_diff = regDF.w .- regDF.w_center;
regDF.wdiff_bigMA = regDF.w_diff .* regDF.bigMA;

regModelFull = lm(@formula(dlnl ~ w_diff + bigMA + wdiff_bigMA), regDF);
println("Model moment regression 1, full sample (", moment_sample_label, "):")
println(coeftable(regModelFull))

regDF_treated = regDF[(regDF.bigMA .== 1) .& isfinite.(regDF.lndma), :];
regDF_treated.wdiff_lndma = regDF_treated.w_diff .* regDF_treated.lndma;
regModelTreated = lm(@formula(dlnl ~ w_diff + lndma + wdiff_lndma), regDF_treated);
println("Model moment regression 2, BIG == 1 (", moment_sample_label, "):")
println(coeftable(regModelTreated))

β_names_full = coefnames(regModelFull)
β_names_treated = coefnames(regModelTreated)
β = [
    coef(regModelFull)[findfirst(==("bigMA"), β_names_full)],
    coef(regModelFull)[findfirst(==("wdiff_bigMA"), β_names_full)],
    coef(regModelTreated)[findfirst(==("lndma"), β_names_treated)],
    coef(regModelTreated)[findfirst(==("wdiff_lndma"), β_names_treated)]
]
println("Reported model β (", moment_sample_label, "): ", β)

#==================================================#
# Plot
#==================================================#

# density(vec(wⱼ), title="Kernel Density Estimate of wⱼ", xlabel="wⱼ", ylabel="Density", legend=false)
plot_w = vec(log.(max.(wⱼ, eps(Float64))))
plot_w_center = mean(plot_w[moment_firm_mask])
plot_w_color = plot_w .- plot_w_center
w_color_limit = quantile(abs.(plot_w_color[moment_firm_mask]), 0.95)
w_color_limit = w_color_limit > 0 ? w_color_limit : maximum(abs.(plot_w_color))
w_color_limit = max(w_color_limit, eps(Float64))
w_color_clims = (-w_color_limit, w_color_limit)

# labor
begin
    sorted_idx = sortperm(vec(wⱼ))
    scatter(vec(ln_dMA)[sorted_idx], vec(clamp.(dlnlⱼ, -Inf, 5))[sorted_idx], 
        marker_z = vec(wⱼ)[sorted_idx],
        color = :RdBu,
        colorbar = true,
        colorbar_title = "ln wⱼ - mean(ln wⱼ)",
        xlabel = "ln(dMA)", 
        ylabel = "dlnlⱼ",
        title = "Employment Change vs Market Access Change",
        legend = false,
        markersize = 3,
        alpha = 0.6,
        dpi = 1000)
end

begin
    sorted_idx = sortperm(plot_w_color)
    scatter(vec(ln_dMA)[sorted_idx], vec(dlnlⱼ)[sorted_idx], 
        marker_z = clamp.(plot_w_color[sorted_idx], w_color_clims...),
        clims = w_color_clims,
        color = :RdBu,
        colorbar = true,
        colorbar_title = "ln wⱼ - mean",
        xlabel = "ln(dMA)", 
        ylabel = "dlnlⱼ",
        title = "Employment Change vs Market Access Change",
        legend = false,
        markersize = 3,
        alpha = 0.6,
        dpi = 1000)
end
savefig(projPath * "/output/figures/labor_change.png")

scatter(vec(dMA), vec(log.(lⱼ)), 
    marker_z = clamp.(plot_w_color, w_color_clims...),
    clims = w_color_clims,
    color = :RdBu,
    colorbar = true,
    colorbar_title = "ln wⱼ - mean",
    xlabel = "dMA", 
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
    marker_z = clamp.(plot_w_color, w_color_clims...),
    clims = w_color_clims,
    color = :RdBu,
    colorbar = true,
    colorbar_title = "ln wⱼ - mean",
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
    marker_z = clamp.(plot_w_color, w_color_clims...),
    clims = w_color_clims,
    color = :RdBu,
    colorbar = true,
    colorbar_title = "ln wⱼ - mean",
    xlabel = "dMA", 
    ylabel = "dlnν",
    title = "Labor Market Power Change vs Market Access Change",
    legend = false,
    markersize = 3,
    alpha = 0.6,
    dpi = 1000)
savefig(projPath * "/output/figures/model/markdown_change.png")

ν = clamp.(ν, 1, 1.3)
scatter(vec(dMA), vec(ν), 
    marker_z = clamp.(plot_w_color, w_color_clims...),
    clims = w_color_clims,
    color = :RdBu,
    colorbar = true,
    colorbar_title = "ln wⱼ - mean",
    xlabel = "dMA", 
    ylabel = "νⱼ",
    title = "Labor Market Power vs Market Access Change",
    legend = false,
    markersize = 3,
    alpha = 0.6,
    dpi = 1000)
savefig(projPath * "/output/figures/model/markdown.png")
