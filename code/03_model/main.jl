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

# Calibration moments
RunStata(projPath, stataPath, "code/02_empirical/calculate_calibration_moments.do")

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
restrict_to_reg_sample = true; # set false to use all model firms for model moments
moment_firm_mask = MomentFirmMask(J; firm_ids, reg_sample_ids, restrict_to_reg_sample);
moment_sample_label = restrict_to_reg_sample ? "matched regression sample" : "full model sample";
println("Model moment sample firms: ", sum(moment_firm_mask), "/", J,
    " in ", moment_sample_label);

firm_ind = df[1:Z:end, :ind_code2];

# `df` is a firm-by-origin-town commute matrix, so its `town` column is not the
# firm-location town used by Stata's cluster(town2#ind).
df_reg_full = DataFrame(load(projPath * "/data/regression/processed/regression_qingdao_07_20.dta"));
firm_town_source = df_reg_full[in.(df_reg_full.year, Ref([2010, 2012])), [:id, :town]];
dropmissing!(firm_town_source, [:id, :town]);
firm_town_by_id = combine(groupby(firm_town_source, :id), :town => first => :firm_town);
firm_town_lookup = Dict(row.id => row.firm_town for row in eachrow(firm_town_by_id));
firm_town = [get(firm_town_lookup, id, missing) for id in firm_ids];
town_ind_cluster = [ismissing(town) || ismissing(ind) ? missing : string(town, "#", ind)
    for (town, ind) in zip(firm_town, firm_ind)];

#==================================================#
# Calibration: Back out η and θ from labor and wage moments with fixed α
#==================================================#

moment_target_path = joinpath(projPath, "output", "tables", "calibration_moments.csv")
moment_order = ["labor_bigMA", "labor_bigMA_wdiff", "wage_bigMA"]
moment_targets = CSV.read(moment_target_path, DataFrame)
moment_lookup = Dict(String(row.moment) => Float64(row.beta) for row in eachrow(moment_targets))
missing_moments = setdiff(moment_order, collect(keys(moment_lookup)))
isempty(missing_moments) || error("Missing calibration moments: " * join(missing_moments, ", "))
β_target = [moment_lookup[moment] for moment in moment_order]
println("Calibration target moments:")
show(DataFrame(moment = moment_order, beta = β_target), allrows = true, allcols = true)
println()

η_bounds = [0.1, 5.0]
θ_bounds = [0.1, 10]
α_fixed = 0.4
# α_bounds = [0.1, 0.9]
employment_change = :log
wage_center = :all

# use grid search to find good starting points for the optimization
η_grid = [0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
θ_grid = [0.1, 0.5, 1.0, 5.0, 10.0]
α_grid = [0.2, 0.4, 0.6, 0.8]
grid_results = EvaluateCalibrationGrid(;
    l, d, d′, wⱼ_data, lⱼ_data, β_target, 
    # α = α_fixed,
    α_grid,
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
CSV.write(
    joinpath(projPath, "output", "tables", "calibration_grid.csv"),
    grid_results;
    bom = true
)

grid_results = CSV.read(projPath * "/output/tables/calibration_grid.csv", DataFrame)
top_grid = first(grid_results, min(4, nrow(grid_results)))
println("\nTop calibration grid points:")
show(top_grid, allrows = true, allcols = true)
println()

calibration_starts = [[row.η, row.θ] for row in eachrow(top_grid)]

calibration = CalibrateEtaThetaAlpha(;
    l, d, d′, wⱼ_data, lⱼ_data, β_target,
    α = α_fixed,
    aⱼ_init,
    starts = [[4.0, 5.0]],
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
η_est, θ_est, α_est = calibration.parameters

println("Final objective: ", calibration.objective)
println("Converged: ", calibration.converged)
println("Bounds: η ∈ ", η_bounds, ", θ ∈ ", θ_bounds, ", fixed α = ", α_est)

# Verify final moments
β_final = ComputeModelMoments([η_est, θ_est, α_est]; l, d, d′, wⱼ_data, lⱼ_data,
    aⱼ_init,
    inner_tol = 1e-5, inner_maxIter = 5000, inner_display = true,
    continuation_steps = 1, employment_change, wage_center,
    moment_firm_mask = moment_firm_mask, displayGap = true, damp_cf = 0.98,
    amenity_maxIter = 10000, amenity_damp = 0.85, amenity_tol = 1e-6
)
println("\nTarget  β: ", β_target)
println("Model   β (", moment_sample_label, "): ", β_final)
println("Estimated η: ", η_est)
println("Estimated θ: ", θ_est)
println("Fixed α: ", α_est)

#==================================================#
# Simulation using calibrated η, θ, and α
#==================================================#

η, θ, α = [3.742411481046953, 4.950962175934991, 0.4]

# Solve firm amenities and productivity from observed employment and wages
vars = (; wⱼ = wⱼ_data, lⱼ = lⱼ_data, l, d)
params = (; η, θ, α)
primitives = SolveFirmPrimitivesFromData(vars, params; aⱼ_init, displaySummary = true, displayGap = true, amenity_damp = 0.7, amenity_tol = 1e-8);
zⱼ, aⱼ = primitives.zⱼ, primitives.aⱼ;

# Solve the model
vars = (; l, d, zⱼ, aⱼ)
params = (; η, θ, α)

# solve the model for baseline
# Use observed wages as the warm start and avoid the power update here:
# zⱼ was inverted from wⱼ_data, so this keeps the solver on the same equilibrium branch.
wⱼ, π_zj, ε_zj, lⱼ, εⱼ = SolveModel(vars, params; displayGap = true, damp = 0.6, tol = 1e-7, displaySummary = true, power = false, wⱼ_init = wⱼ_data);

## Calculate correlation between solved wages and observed wages
println("Corrleation between lⱼ and lⱼ_data:", cor(vec(lⱼ), vec(lⱼ_data)))
println("Corrleation between wⱼ and wⱼ_data:", cor(vec(wⱼ), vec(wⱼ_data)))
println("Corrleation between wⱼ and aⱼ:", cor(vec(wⱼ), vec(aⱼ)))
println("Corrleation between wⱼ and lⱼ:", cor(vec(wⱼ), vec(lⱼ)))
# println("Correlation between wⱼ_solved and wⱼ_data: ", corr_wages) # should be very close to 1

# solve the model for counterfactual
vars′ = (; l, d = d′, zⱼ, aⱼ);
wⱼ′, π_zj′, ε_zj′, lⱼ′, εⱼ′ = SolveModel(vars′, params; displayGap = true, damp = 0.98, tol = 1e-9, displaySummary = true, power = true, wⱼ_init = wⱼ, maxIter = 3000);


# Analyze the results
l̂ⱼ = lⱼ′ ./ lⱼ;
dlnlⱼ = log.(max.(lⱼ′, eps(Float64))) .- log.(max.(lⱼ, eps(Float64)));
dMA = sum((d - d′) .* l, dims = 1)' |> x -> replace(x, -Inf => -8);
# density(dMA, title="Kernel Density Estimate of dMA", xlabel="dMA", ylabel="Density", legend=false)
bigMA = Float64.(dMA .>= 0.5);
ln_dMA = log.(sum((d - d′) .* l, dims = 1)') |> x -> replace(x, -Inf => -8)

lnl = vec(log.(max.(lⱼ, eps(Float64))))
lnl′ = vec(log.(max.(lⱼ′, eps(Float64))))
lnw = vec(log.(max.(wⱼ, eps(Float64))))
lnw′ = vec(log.(max.(wⱼ′, eps(Float64))))
β_report = EstimateTwoPeriodDIDMoments(lnl, lnl′, lnw, lnw′, vec(bigMA), lnw;
    wage_center, keep = moment_firm_mask, cluster = town_ind_cluster,
    return_stats = true)
β = β_report.beta
println("Reported model β (", moment_sample_label, "): ", β)
show(DataFrame(moment = moment_order, beta = β_report.beta, se = β_report.se,
    t = β_report.t, n_clusters = β_report.n_clusters, df = β_report.df),
    allrows = true, allcols = true)
println()

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

# Plot model labor reallocation for matched sample
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

# Plot real-data labor reallocation for regression sample
begin
    reg_sample_id_set = Set(string.(reg_sample_ids))
    reg_year_mask = coalesce.(in.(df_reg_full.year, Ref([2010, 2012])), false)
    reg_id_mask = in.(string.(df_reg_full.id), Ref(reg_sample_id_set))
    reg_labor_raw = df_reg_full[reg_year_mask .& reg_id_mask,
        [:id, :year, :employ, :dma, :wage_total]]
    dropmissing!(reg_labor_raw, [:id, :year, :employ, :dma, :wage_total])
    disallowmissing!(reg_labor_raw, [:id, :year, :employ, :dma, :wage_total])
    reg_labor_raw = reg_labor_raw[
        (Float64.(reg_labor_raw.employ) .> 0) .&
        (Float64.(reg_labor_raw.wage_total) .> 0), :]

    if nrow(reg_labor_raw) == 0
        @warn "Skipping real-data labor reallocation plot: no matched 2010/2012 observations."
    else
        reg_labor_raw[!, :lnemp] = log.(Float64.(reg_labor_raw.employ))
        reg_labor_raw[!, :lnw] =
            log.(Float64.(reg_labor_raw.wage_total) ./ Float64.(reg_labor_raw.employ))

        reg_lnw0 = reg_labor_raw.lnw[reg_labor_raw.year .== 2010]
        reg_wlo, reg_whi = quantile(reg_lnw0, [0.05, 0.95])
        reg_labor_raw[!, :lnw_winsor] = clamp.(reg_labor_raw.lnw, reg_wlo, reg_whi)

        reg_base = select(reg_labor_raw[reg_labor_raw.year .== 2010, :],
            :id, :lnemp => :lnemp_2010, :dma => :dma, :lnw_winsor => :lnw0)
        reg_post = select(reg_labor_raw[reg_labor_raw.year .== 2012, :],
            :id, :lnemp => :lnemp_2012)
        sort!(reg_base, :id); unique!(reg_base, :id)
        sort!(reg_post, :id); unique!(reg_post, :id)

        reg_labor_plot = innerjoin(reg_base, reg_post, on = :id)
        if nrow(reg_labor_plot) > 0
            reg_labor_plot[!, :dlnemp] =
                reg_labor_plot.lnemp_2012 .- reg_labor_plot.lnemp_2010
            reg_labor_plot[!, :dlnemp_clip] = min.(reg_labor_plot.dlnemp, 5.0)
            reg_labor_plot[!, :wdiff] = reg_labor_plot.lnw0 .- mean(reg_labor_plot.lnw0)

            reg_dma = Float64.(reg_labor_plot.dma)
            reg_positive_ln_dma = log.(reg_dma[reg_dma .> 0])
            reg_ln_dma_floor = isempty(reg_positive_ln_dma) ? -8.0 :
                minimum(reg_positive_ln_dma) - 0.1
            reg_labor_plot[!, :ln_dMA] = [dma > 0 ? log(dma) : reg_ln_dma_floor
                for dma in reg_dma]

            reg_sorted_idx = sortperm(reg_labor_plot.lnw0)

            p_real_labor = scatter(reg_labor_plot.ln_dMA[reg_sorted_idx],
                reg_labor_plot.dlnemp_clip[reg_sorted_idx],
                marker_z = reg_labor_plot.wdiff[reg_sorted_idx],
                color = :RdBu,
                colorbar = true,
                colorbar_title = "ln w0 - mean(ln w0)",
                xlabel = "ln(dMA)", 
                ylabel = "dlnlⱼ",
                title = "Employment Change vs Market Access Change",
                legend = false,
                markersize = 3,
                alpha = 0.6,
                dpi = 1000)
            savefig(p_real_labor, projPath * "/output/figures/labor_reallocation_realdata_julia.png")
        end
    end
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

#==================================================#
# Diagnostic: differential wage-slope of Δemployment
# Mirrors the labor moment β₂ (labor_bigMA_wdiff).
# The two-period firm+period FE DiD first-differences to
#   dlnlⱼ = c + γ·w_diff + β₁·treated + β₂·(treated·w_diff) + ε,
# so the treated-minus-control gap in the dlnl-on-w_diff slope IS β₂.
# We plot dlnl vs initial wage deviation, split by treatment, with
# binned means (to cut through the cloud) and within-group OLS lines.
#==================================================#
begin
    mask         = vec(moment_firm_mask)
    wdiff_plot   = plot_w_color[mask]            # ln wⱼ - mean, same centering as the reg
    dlnl_plot    = vec(dlnlⱼ)[mask]
    treated_plot = vec(bigMA)[mask] .== 1

    # quantile-binned group means of y against x
    binmeans = function (x, y; nbins = 20)
        edges = quantile(x, range(0, 1, length = nbins + 1))
        edges[1] -= eps(); edges[end] += eps()
        bin = clamp.(searchsortedlast.(Ref(edges), x), 1, nbins)
        present = [b for b in 1:nbins if any(bin .== b)]
        (xb = [mean(x[bin .== b]) for b in present],
         yb = [mean(y[bin .== b]) for b in present])
    end

    diag_df  = DataFrame(dlnl = dlnl_plot, wdiff = wdiff_plot, treated = treated_plot)
    fit_ctrl = lm(@formula(dlnl ~ wdiff), diag_df[.!diag_df.treated, :])
    fit_trt  = lm(@formula(dlnl ~ wdiff), diag_df[diag_df.treated, :])
    slope_ctrl, slope_trt = coef(fit_ctrl)[2], coef(fit_trt)[2]

    xgrid    = range(minimum(wdiff_plot), maximum(wdiff_plot), length = 100)
    pred(f)  = coef(f)[1] .+ coef(f)[2] .* xgrid
    bc_ctrl  = binmeans(wdiff_plot[.!treated_plot], dlnl_plot[.!treated_plot])
    bc_trt   = binmeans(wdiff_plot[treated_plot],   dlnl_plot[treated_plot])
    n_ctrl   = sum(.!treated_plot)
    n_trt    = sum(treated_plot)

    ctrl_raw_color = "#9ECAE1"
    trt_raw_color  = "#F4A3A8"
    ctrl_color     = "#2166AC"
    trt_color      = "#B2182B"
    p_wdiff = plot(xlabel = "initial wage deviation  (ln wⱼ - mean)",
        ylabel = "dlnlⱼ",
        title = "Employment wage-slope (β₂ = $(round(slope_trt - slope_ctrl, digits = 1)))",
        titlefontsize = 10,
        legend = :topleft, dpi = 1000)
    scatter!(p_wdiff, wdiff_plot[.!treated_plot], dlnl_plot[.!treated_plot],
        label = "control firms (N = $(n_ctrl))",
        color = ctrl_raw_color, seriescolor = ctrl_raw_color,
        markercolor = ctrl_raw_color, markerstrokewidth = 0, markersize = 2)
    scatter!(p_wdiff, wdiff_plot[treated_plot], dlnl_plot[treated_plot],
        label = "treated firms (N = $(n_trt))",
        color = trt_raw_color, seriescolor = trt_raw_color,
        markercolor = trt_raw_color, markerstrokewidth = 0, markersize = 1)
    plot!(p_wdiff, xgrid, pred(fit_ctrl), color = ctrl_color, seriescolor = ctrl_color,
        linecolor = ctrl_color, lw = 2, linestyle = :dash,
        label = "control slope γ = $(round(slope_ctrl, digits = 3))")
    plot!(p_wdiff, xgrid, pred(fit_trt), color = trt_color, seriescolor = trt_color,
        linecolor = trt_color, lw = 2, linestyle = :dash,
        label = "treated slope γ+β₂ = $(round(slope_trt, digits = 3))")
    savefig(p_wdiff, projPath * "/output/figures/labor_wdiff_slope.png")
end



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
