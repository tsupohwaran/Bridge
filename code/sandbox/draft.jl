#==================================================#
# Baseline and counterfactual markdown distributions by eta
#==================================================#

projPath = normpath(joinpath(@__DIR__, "..", ".."))
cd(projPath)

include(joinpath(projPath, "code", "03_model", "load_packages.jl"))
include(joinpath(projPath, "code", "03_model", "functions.jl"))

const Z = 128
const eta_grid = [0.1, 0.9, 3.0]
const alpha_fixed = 0.8
const theta_fixed = 10.0
const sigma_fixed = 0.25
const amenity_tol = 1e-6
const amenity_damp = 0.6
const counterfactual_tol = 1e-6

model_path = joinpath(projPath, "data", "model", "processed", "firm_qingdao_model_10.dta")
output_figure_dir = joinpath(projPath, "output", "figures", "sandbox")
output_table_dir = joinpath(projPath, "output", "tables")
mkpath(output_figure_dir)
mkpath(output_table_dir)

println("Loading baseline model sample: ", model_path)
df = DataFrame(load(model_path))
nrow(df) % Z == 0 || error("Model data row count must be divisible by Z = $Z")
J = nrow(df) ÷ Z

l = df[!, :pop] |>
    x -> Float64.(x) |>
    x -> reshape(x, Z, J) |>
    x -> x[:, 1] |>
    x -> x ./ sum(x)

l_j_data = Float64.(disallowmissing(df[!, :employ])) |>
    x -> reshape(x, Z, J) |>
    x -> x[1, :] |>
    x -> x ./ sum(x)

d = reshape(Float64.(df[!, :dzj]), Z, J) |> x -> replace(x, 0.0 => 1e-2)
d_prime = reshape(Float64.(df[!, :dzj_prime]), Z, J) |> x -> replace(x, 0.0 => 1e-2)

w_j_raw = reshape(Float64.(df[!, :wage_inital]), Z, J)[1, :]
lo, hi = quantile(w_j_raw, [0.05, 0.95])
w_j_data = clamp.(w_j_raw, lo, hi) |>
    x -> x ./ sum(x .* l_j_data)

firm_ids = collect(df[1:Z:end, :id])
firm_sector_raw = df[1:Z:end, :ind_agg]
any(ismissing, firm_sector_raw) &&
    error("Model sample must contain nonmissing ind_agg for sigma = $sigma_fixed.")
firm_sector = Int.(firm_sector_raw)
a_j_init = zeros(J)

markdown_frames = DataFrame[]
summary_frames = DataFrame[]
markdown_change_frames = DataFrame[]
change_summary_frames = DataFrame[]
markdown_by_eta = Dict{Float64, Vector{Float64}}()
markdown_hat_by_eta = Dict{Float64, Vector{Float64}}()

for eta in eta_grid
    println("\nSolving baseline and counterfactual markdown for eta = ", eta)
    params = (; η = eta, θ = theta_fixed, σ = sigma_fixed, α = alpha_fixed)
    vars = (; wⱼ = w_j_data, lⱼ = l_j_data, l, d, firm_sector)

    primitives = SolveFirmPrimitivesFromData(vars, params;
        aⱼ_init = a_j_init,
        amenity_tol,
        amenity_damp,
        amenity_maxIter = 20000,
        displaySummary = true,
        displayGap = false,
        returnInfo = true)

    z_j = vec(primitives.zⱼ)
    a_j = vec(primitives.aⱼ)
    epsilon_j = vec(primitives.εⱼ)
    markdown = vec(1 .+ 1 ./ epsilon_j)
    all(isfinite.(markdown)) || error("Non-finite markdown values for eta = $eta")
    markdown_by_eta[eta] = markdown

    cf_vars = (; l, d = d_prime, zⱼ = z_j, aⱼ = a_j, firm_sector)
    cf = SolveModel(cf_vars, params;
        displayGap = false,
        displaySummary = true,
        damp = 0.98,
        tol = counterfactual_tol,
        maxIter = 3000,
        power = true,
        wⱼ_init = w_j_data,
        returnInfo = true)

    if !cf.converged
        @warn "Counterfactual solve did not converge for eta = $eta" iterations = cf.iterations gap = cf.gap
    end

    epsilon_j_prime = vec(cf.εⱼ)
    markdown_prime = vec(1 .+ 1 ./ epsilon_j_prime)
    markdown_hat = markdown_prime ./ markdown
    all(isfinite.(markdown_hat)) || error("Non-finite markdown-hat values for eta = $eta")
    markdown_hat_by_eta[eta] = markdown_hat

    qs = quantile(markdown, [0.05, 0.25, 0.50, 0.75, 0.95])
    push!(summary_frames, DataFrame(
        eta = [eta],
        alpha = [alpha_fixed],
        theta = [theta_fixed],
        sigma = [sigma_fixed],
        n_firms = [J],
        amenity_converged = [primitives.amenity_converged],
        amenity_iterations = [primitives.amenity_iterations],
        amenity_gap = [primitives.amenity_gap],
        mean = [mean(markdown)],
        sd = [std(markdown)],
        p05 = [qs[1]],
        p25 = [qs[2]],
        median = [qs[3]],
        p75 = [qs[4]],
        p95 = [qs[5]],
        min = [minimum(markdown)],
        max = [maximum(markdown)]
    ))

    qs_hat = quantile(markdown_hat, [0.05, 0.25, 0.50, 0.75, 0.95])
    push!(change_summary_frames, DataFrame(
        eta = [eta],
        alpha = [alpha_fixed],
        theta = [theta_fixed],
        sigma = [sigma_fixed],
        n_firms = [J],
        counterfactual_converged = [cf.converged],
        counterfactual_iterations = [cf.iterations],
        counterfactual_gap = [cf.gap],
        mean = [mean(markdown_hat)],
        sd = [std(markdown_hat)],
        p05 = [qs_hat[1]],
        p25 = [qs_hat[2]],
        median = [qs_hat[3]],
        p75 = [qs_hat[4]],
        p95 = [qs_hat[5]],
        min = [minimum(markdown_hat)],
        max = [maximum(markdown_hat)]
    ))

    push!(markdown_frames, DataFrame(
        id = firm_ids,
        eta = fill(eta, J),
        epsilon_j = epsilon_j,
        markdown = markdown
    ))

    push!(markdown_change_frames, DataFrame(
        id = firm_ids,
        eta = fill(eta, J),
        epsilon_j = epsilon_j,
        epsilon_j_prime = epsilon_j_prime,
        markdown = markdown,
        markdown_prime = markdown_prime,
        markdown_hat = markdown_hat,
        dln_markdown = log.(markdown_hat)
    ))
end

markdown_df = vcat(markdown_frames...)
summary_df = vcat(summary_frames...)
markdown_change_df = vcat(markdown_change_frames...)
change_summary_df = vcat(change_summary_frames...)

CSV.write(joinpath(output_table_dir, "sandbox_baseline_markdown_by_eta.csv"),
    markdown_df; bom = true)
CSV.write(joinpath(output_table_dir, "sandbox_baseline_markdown_by_eta_summary.csv"),
    summary_df; bom = true)
CSV.write(joinpath(output_table_dir, "sandbox_markdown_hat_by_eta.csv"),
    markdown_change_df; bom = true)
CSV.write(joinpath(output_table_dir, "sandbox_markdown_hat_by_eta_summary.csv"),
    change_summary_df; bom = true)

println("\nBaseline markdown summary:")
show(summary_df, allrows = true, allcols = true)
println()

println("\nMarkdown-hat summary:")
show(change_summary_df, allrows = true, allcols = true)
println()

palette = [:steelblue, :darkorange, :seagreen]
all_markdown = markdown_df.markdown
bulk_xmax = quantile(all_markdown, 0.995)
p = plot(
    xlabel = "Baseline firm markdown, nu_j = 1 + 1 / epsilon_j",
    ylabel = "Density",
    title = "Baseline Markdown Distribution by eta (99.5% bulk)\n" *
        "alpha = $alpha_fixed, theta = $theta_fixed, sigma = $sigma_fixed",
    legend = :topright,
    grid = :y,
    xlims = (minimum(all_markdown), bulk_xmax),
    titlefontsize = 12,
    dpi = 400
)

for (idx, eta) in enumerate(eta_grid)
    markdown = markdown_by_eta[eta]
    density!(p, markdown;
        label = "eta = $(eta)",
        color = palette[idx],
        linewidth = 2.5)
    vline!(p, [mean(markdown)];
        label = "",
        color = palette[idx],
        linestyle = :dash,
        linewidth = 1.4)
end

output_figure_path = joinpath(output_figure_dir, "baseline_markdown_distribution_by_eta.png")
savefig(p, output_figure_path)
println("\nSaved figure: ", output_figure_path)

all_markdown_hat = markdown_change_df.markdown_hat
hat_lo, hat_hi = quantile(all_markdown_hat, [0.005, 0.995])
p_hat = plot(
    xlabel = "Markdown change, nu_hat_j = nu_prime_j / nu_j",
    ylabel = "Density",
    title = "Distribution of Markdown Change by eta (central 99%)\n" *
        "alpha = $alpha_fixed, theta = $theta_fixed, sigma = $sigma_fixed",
    legend = :topright,
    grid = :y,
    xlims = (hat_lo, hat_hi),
    titlefontsize = 12,
    dpi = 400
)

vline!(p_hat, [1.0];
    label = "no change",
    color = :gray40,
    linestyle = :dot,
    linewidth = 1.5)

for (idx, eta) in enumerate(eta_grid)
    markdown_hat = markdown_hat_by_eta[eta]
    density!(p_hat, markdown_hat;
        label = "eta = $(eta)",
        color = palette[idx],
        linewidth = 2.5)
    vline!(p_hat, [mean(markdown_hat)];
        label = "",
        color = palette[idx],
        linestyle = :dash,
        linewidth = 1.4)
end

output_hat_figure_path = joinpath(output_figure_dir, "markdown_hat_distribution_by_eta.png")
savefig(p_hat, output_hat_figure_path)
println("Saved figure: ", output_hat_figure_path)
