#==================================================#
# Plot wage-ladder poaching GE labor change
#==================================================#

ENV["ETA"] = get(ENV, "ETA", "0.9")
ENV["THETA"] = get(ENV, "THETA", "10")
ENV["SIGMA"] = get(ENV, "SIGMA", "0.25")
ENV["ALPHA"] = get(ENV, "ALPHA", "0.8")
ENV["RHO_GRID"] = get(ENV, "RHO_GRID", get(ENV, "RHO_PLOT", "0.25"))
ENV["KAPPA_GRID"] = get(ENV, "KAPPA_GRID", get(ENV, "KAPPA_PLOT", "0"))
ENV["RUN_INVERSION"] = get(ENV, "RUN_INVERSION", "1")
ENV["FAST_AMENITY"] = get(ENV, "FAST_AMENITY", "1")
ENV["RUN_GE"] = "0"
ENV["OUTPUT_NAME"] = get(ENV, "OUTPUT_NAME", "poaching_group_plot_tmp.csv")

include(joinpath(@__DIR__, "diagnose_poaching_groups.jl"))

rho_plot = parse(Float64, get(ENV, "RHO_PLOT", "0.25"))
kappa_plot = parse(Float64, get(ENV, "KAPPA_PLOT", "0"))
use_continuation = get(ENV, "PLOT_CONTINUATION", "1") in
    ["1", "true", "TRUE", "yes", "YES"]

println("Solving poaching GE for plot: rho=$rho_plot, kappa=$kappa_plot")
if use_continuation && (!iszero(rho_plot) || !iszero(kappa_plot))
    ge0 = solve_poaching_ge(M_zg, d_post, w_j_data, z_j, a_j, firm_sector,
        firm_ind, groups.group_lnw, groups.group_ind, 0.0, 0.0;
        tol = ge_tol,
        maxIter = ge_maxIter,
        damp = ge_damp,
        power = ge_power,
        displayGap = ge_display,
        maxLogStep = ge_log_step)
    ge = solve_poaching_ge(M_zg, d_post, ge0.w, z_j, a_j, firm_sector,
        firm_ind, groups.group_lnw, groups.group_ind, rho_plot, kappa_plot;
        tol = ge_tol,
        maxIter = ge_maxIter,
        damp = ge_damp,
        power = ge_power,
        displayGap = ge_display,
        maxLogStep = ge_log_step)
else
    ge = solve_poaching_ge(M_zg, d_post, w_j_data, z_j, a_j, firm_sector,
        firm_ind, groups.group_lnw, groups.group_ind, rho_plot, kappa_plot;
        tol = ge_tol,
        maxIter = ge_maxIter,
        damp = ge_damp,
        power = ge_power,
        displayGap = ge_display,
        maxLogStep = ge_log_step)
end

plot_w = vec(log.(max.(w_j_data, eps(Float64))))
plot_w_center = mean(plot_w[keep_mask])
plot_w_color = plot_w .- plot_w_center
w_color_limit = quantile(abs.(plot_w_color[keep_mask]), 0.95)
w_color_limit = w_color_limit > 0 ? w_color_limit : maximum(abs.(plot_w_color))
w_color_limit = max(w_color_limit, eps(Float64))
w_color_clims = (-w_color_limit, w_color_limit)

dMA_plot = vec(sum((d .- d_post) .* reshape(l, :, 1), dims = 1))
ln_dMA = replace(log.(dMA_plot), -Inf => -8)
dlnl = log.(max.(ge.L, eps(Float64))) .- log.(max.(l_j_data, eps(Float64)))
sorted_idx = sortperm(plot_w_color)

p = scatter(ln_dMA[sorted_idx], dlnl[sorted_idx],
    marker_z = clamp.(plot_w_color[sorted_idx], w_color_clims...),
    clims = w_color_clims,
    color = :RdBu,
    colorbar = true,
    colorbar_title = "ln wⱼ - mean",
    xlabel = "ln(dMA)",
    ylabel = "dlnlⱼ",
    title = "Poaching GE: Employment Change vs Market Access (rho = $(rho_plot))",
    legend = false,
    markersize = 1.5,
    markerstrokewidth = 0,
    alpha = 0.6,
    size = (1100, 760),
    titlefontsize = 12,
    guidefontsize = 11,
    tickfontsize = 9,
    colorbar_titlefontsize = 11,
    colorbar_tickfontsize = 9,
    margin = 5Plots.mm,
    dpi = 1000)

mkpath(joinpath(projPath, "output", "figures", "model"))
outpath = joinpath(projPath, "output", "figures", "model",
    "labor_change_poaching_rho$(replace(string(rho_plot), "." => "p")).png")
savefig(p, outpath)

println("GE converged: ", ge.converged)
println("GE iterations: ", ge.iterations)
println("GE gap: ", ge.gap)
println("Saved figure: ", outpath)
