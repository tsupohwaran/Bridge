#==================================================#
# Labor-market-power mode comparison figures
#==================================================#

ENV["GKSwstype"] = "100"

projPath = normpath(joinpath(@__DIR__, "..", ".."))
cd(projPath)

include(joinpath(@__DIR__, "load_packages.jl"))
include(joinpath(@__DIR__, "functions.jl"))

function number_slug(x)
    x_float = Float64(x)
    if isapprox(x_float, round(x_float); atol = 1e-12)
        return string(Int(round(x_float)))
    end
    return replace(string(x_float), "." => "")
end

const Z = 128
const η = 3.0
const θ = 5.0
const σ = 1.0
const α = 0.8
const counterfactual_maxIter = 8000
const wage_center = :all
const moment_order = ["labor_bigMA", "labor_bigMA_wdiff", "wage_bigMA"]
const parameter_suffix = "eta$(number_slug(η))_theta$(number_slug(θ))_sigma$(number_slug(σ))"

model_path = joinpath(projPath, "data", "model", "processed", "firm_qingdao_model_10.dta")
println("Loading model sample: ", model_path)
df = DataFrame(load(model_path))
nrow(df) % Z == 0 || error("Model data row count must be divisible by Z = $Z")
J = nrow(df) ÷ Z

l = df[!, :pop] |>
    x -> Float64.(x) |>
    x -> reshape(x, Z, J) |>
    x -> x[:, 1] |>
    x -> x ./ sum(x)
lⱼ_data = Float64.(disallowmissing(df[!, :employ])) |>
    x -> reshape(x, Z, J) |>
    x -> x[1, :] |>
    x -> x ./ sum(x)
d = reshape(Float64.(df[!, :dzj]), Z, J) |> x -> replace(x, 0.0 => 1e-2)
d′ = reshape(Float64.(df[!, :dzj_prime]), Z, J) |> x -> replace(x, 0.0 => 1e-2)
wⱼ_raw = reshape(Float64.(df[!, :wage_inital]), Z, J)[1, :]
lo, hi = quantile(wⱼ_raw, [0.05, 0.95])
wⱼ_data = clamp.(wⱼ_raw, lo, hi) |>
    x -> x ./ sum(x .* lⱼ_data)
aⱼ_init = zeros(J)

firm_sector_raw = df[1:Z:end, :ind_agg]
any(ismissing, firm_sector_raw) &&
    error("Model sample must contain nonmissing ind_agg for nested-logit sector nests.")
firm_sector = Int.(firm_sector_raw)
firm_ind = firm_sector
firm_ids = collect(df[1:Z:end, :id])

reg_path = joinpath(projPath, "data", "model", "processed", "firm_two_year_reg.dta")
println("Loading model moment sample: ", reg_path)
df_reg = DataFrame(load(reg_path))
reg_sample_ids = collect(skipmissing(df_reg[!, :id]))
restrict_to_reg_sample = true
moment_firm_mask = MomentFirmMask(J; firm_ids, reg_sample_ids, restrict_to_reg_sample)
moment_sample_label = restrict_to_reg_sample ? "matched regression sample" : "full model sample"
println("Model moment sample firms: ", sum(moment_firm_mask), "/", J,
    " in ", moment_sample_label)

# `town` is the worker-origin town in the firm-by-town commute matrix; `firm_town`
# is the coordinate-derived firm-location town used for town-sector clustering.
if "firm_town" in names(df)
    firm_town = collect(df[1:Z:end, :firm_town])
else
    reg_full_path = joinpath(projPath, "data", "regression", "processed",
        "regression_qingdao_07_20.dta")
    println("Loading firm-town clusters: ", reg_full_path)
    df_reg_full = DataFrame(load(reg_full_path))
    firm_town_source = df_reg_full[in.(df_reg_full.year, Ref([2010, 2012])), [:id, :town]]
    dropmissing!(firm_town_source, [:id, :town])
    firm_town_by_id = combine(groupby(firm_town_source, :id), :town => first => :firm_town)
    firm_town_lookup = Dict(row.id => row.firm_town for row in eachrow(firm_town_by_id))
    firm_town = [get(firm_town_lookup, id, missing) for id in firm_ids]
end
town_ind_cluster = [ismissing(town) || ismissing(ind) ? missing : string(town, "#", ind)
    for (town, ind) in zip(firm_town, firm_ind)]

coefficient_tables = DataFrame[]

function log_market_access_change(dMA)
    positive = dMA[dMA .> 0]
    floor_value = isempty(positive) ? -8.0 : minimum(log.(positive)) - 0.1
    return [x > 0 ? log(x) : floor_value for x in dMA]
end

function solve_and_plot_lmp_mode(lmp_mode::Symbol, label::String, slug::String)
    println("\nSolving ", label, " model with η=$η, θ=$θ, σ=$σ, α=$α")
    params = (; η, θ, σ, α, lmp_mode)

    primitive_vars = (; wⱼ = wⱼ_data, lⱼ = lⱼ_data, l, d, firm_sector)
    primitives = SolveFirmPrimitivesFromData(primitive_vars, params;
        aⱼ_init,
        displaySummary = true,
        displayGap = false,
        amenity_damp = 0.6,
        amenity_tol = 1e-6,
        amenity_maxIter = 10000)
    primitives.amenity_converged ||
        @warn "Amenity inversion did not converge for $label"

    wⱼ = vec(Float64.(wⱼ_data))
    lⱼ = vec(primitives.lⱼ)
    zⱼ, aⱼ = primitives.zⱼ, primitives.aⱼ

    vars′ = (; l, d = d′, zⱼ, aⱼ, firm_sector)
    cf = SolveModel(vars′, params;
        displayGap = false,
        damp = 0.98,
        tol = 1e-9,
        displaySummary = true,
        power = true,
        wⱼ_init = wⱼ,
        maxIter = counterfactual_maxIter,
        returnInfo = true)
    cf.converged || @warn "Counterfactual solve did not converge for $label"

    dlnlⱼ = log.(max.(vec(cf.lⱼ), eps(Float64))) .-
        log.(max.(lⱼ, eps(Float64)))
    dMA = vec(l' * (d .- d′))
    ln_dMA = log_market_access_change(dMA)
    bigMA = Float64.(dMA .>= 0.5)

    lnl = vec(log.(max.(lⱼ, eps(Float64))))
    lnl′ = vec(log.(max.(vec(cf.lⱼ), eps(Float64))))
    lnw = vec(log.(max.(wⱼ, eps(Float64))))
    lnw′ = vec(log.(max.(vec(cf.wⱼ), eps(Float64))))
    β_report = EstimateTwoPeriodDIDMoments(lnl, lnl′, lnw, lnw′, bigMA, lnw;
        wage_center, keep = moment_firm_mask, cluster = town_ind_cluster,
        return_stats = true)
    coef_table = DataFrame(
        lmp_mode = fill(String(lmp_mode), length(moment_order)),
        label = fill(label, length(moment_order)),
        moment = moment_order,
        beta = β_report.beta,
        se = β_report.se,
        t = β_report.t,
        n_clusters = β_report.n_clusters,
        df = β_report.df
    )
    push!(coefficient_tables, coef_table)
    println("Regression coefficients for ", label, " (", moment_sample_label, "):")
    show(coef_table, allrows = true, allcols = true)
    println()

    sorted_idx = sortperm(wⱼ)
    p = scatter(ln_dMA[sorted_idx], clamp.(dlnlⱼ, -Inf, 5)[sorted_idx],
        marker_z = wⱼ[sorted_idx],
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

    output_dir = joinpath(projPath, "output", "figures", "model")
    mkpath(output_dir)
    outfile = joinpath(output_dir,
        "labor_reallocation_$(slug)_$(parameter_suffix).png")
    savefig(p, outfile)
    println("Saved figure: ", outfile)

    return (;
        label,
        lmp_mode,
        outfile,
        amenity_converged = primitives.amenity_converged,
        amenity_gap = primitives.amenity_gap,
        counterfactual_converged = cf.converged,
        counterfactual_gap = cf.gap,
        mean_markdown = mean(ModelMarkdown(vec(primitives.εⱼ), params))
    )
end

results = [
    solve_and_plot_lmp_mode(:monopsony, "monopsony markdown", "monopsony"),
    solve_and_plot_lmp_mode(:perfectly_elastic, "perfectly elastic / no markdown",
        "perfectly_elastic_no_markdown")
]

coef_results = vcat(coefficient_tables...)
table_dir = joinpath(projPath, "output", "tables")
mkpath(table_dir)
coef_path = joinpath(table_dir,
    "lmp_mode_comparison_coefficients_$(parameter_suffix).csv")
CSV.write(coef_path, coef_results; bom = true)

println("\nGenerated LMP mode comparison figures:")
show(DataFrame(results), allrows = true, allcols = true)
println()

println("\nSaved regression coefficients: ", coef_path)
show(coef_results, allrows = true, allcols = true)
println()
