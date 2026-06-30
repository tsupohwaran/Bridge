#==================================================#
# Fixed-Primitives Theta Grid Sandbox
#==================================================#

projPath = normpath(joinpath(@__DIR__, "..", ".."))
cd(projPath)

include(joinpath(@__DIR__, "load_packages.jl"))
include(joinpath(@__DIR__, "functions.jl"))

const Z = 128
const η_seed = 0.5
const θ_seed = 5.0
const α_fixed = 0.6
const θ_min = 0.1
const θ_max = 30.0

model_path = joinpath(projPath, "data", "model", "processed", "firm_qingdao_model_10.dta")
reg_path = joinpath(projPath, "data", "model", "processed", "firm_two_year_reg.dta")
moment_target_path = joinpath(projPath, "output", "tables", "calibration_moments.csv")
output_table_path = joinpath(projPath, "output", "tables")

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

println("Loading moment sample ids: ", reg_path)
df_reg = DataFrame(load(reg_path))
firm_ids = collect(df[1:Z:end, :id])
reg_sample_ids = collect(skipmissing(df_reg[!, :id]))
restrict_to_reg_sample = true
moment_firm_mask = MomentFirmMask(J; firm_ids, reg_sample_ids, restrict_to_reg_sample)
moment_sample_label = restrict_to_reg_sample ? "matched regression sample" : "full model sample"
println("Model firms: ", J)
println("Model moment sample firms: ", sum(moment_firm_mask), "/", J,
    " in ", moment_sample_label)

moment_order = ["labor_bigMA", "labor_bigMA_wdiff", "labor_wdiff_post"]
moment_targets = CSV.read(moment_target_path, DataFrame)
moment_lookup = Dict(String(row.moment) => Float64(row.beta) for row in eachrow(moment_targets))
missing_moments = setdiff(moment_order, collect(keys(moment_lookup)))
isempty(missing_moments) || error("Missing calibration moments: " * join(missing_moments, ", "))
β_target = [moment_lookup[moment] for moment in moment_order]
println("Calibration target moments:")
show(DataFrame(moment = moment_order, beta = β_target), allrows=true, allcols=true)
println()

println("\nRecovering fixed primitives with η = $η_seed, θ = $θ_seed, α = $α_fixed")
seed_params = (; η = η_seed, θ = θ_seed, α = α_fixed)
primitives = SolveFirmPrimitivesFromData((; wⱼ = wⱼ_data, lⱼ = lⱼ_data, l, d),
    seed_params;
    aⱼ_init,
    amenity_damp = 0.7,
    amenity_tol = 1e-8,
    amenity_maxIter = 10000,
    displaySummary = true,
    displayGap = false,
    returnInfo = true)
zⱼ, aⱼ = vec(primitives.zⱼ), vec(primitives.aⱼ)

println("Primitive inversion diagnostics:")
println("  converged:  ", primitives.amenity_converged)
println("  iterations: ", primitives.amenity_iterations)
println("  gap:        ", primitives.amenity_gap)
println("  mean(a_j):  ", mean(aⱼ))
println("  mean(z_j):  ", mean(zⱼ))

mkpath(output_table_path)
primitive_df = DataFrame(
    id = firm_ids,
    z_j = zⱼ,
    a_j = aⱼ,
    epsilon_j = vec(primitives.εⱼ),
    l_j_target = vec(primitives.lⱼ_target),
    l_j_model = vec(primitives.lⱼ)
)
CSV.write(joinpath(output_table_path, "sandbox_fixed_primitives.csv"),
    primitive_df; bom=true)

fixed_primitives = (; zⱼ, aⱼ)
employment_change = :log
wage_center = :all

coarse_θ_grid = [0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 15.0, 20.0]
println("\nRunning coarse θ grid...")
coarse_results = EvaluateCalibrationGrid(;
    l, d, d′, wⱼ_data, lⱼ_data, β_target,
    α = α_fixed,
    η_grid = [η_seed],
    θ_grid = coarse_θ_grid,
    fixed_primitives,
    inner_tol = 2e-5,
    inner_maxIter = 3000,
    continuation_steps = 5,
    employment_change,
    wage_center,
    moment_firm_mask,
    max_abs_moment = 10.0,
    verbose = true
)
CSV.write(joinpath(output_table_path, "sandbox_theta_grid_fixed_primitives_coarse.csv"),
    coarse_results; bom=true)

finite_coarse = filter(:objective => isfinite, coarse_results)
nrow(finite_coarse) > 0 || error("No finite objective values found in coarse θ grid")
best_coarse = first(finite_coarse, 1)
bestθ = Float64(best_coarse.θ[1])
half_width = max(1.0, 0.25 * bestθ)
fine_lo = max(θ_min, bestθ - half_width)
fine_hi = min(θ_max, bestθ + half_width)
fine_θ_grid = collect(range(fine_lo, fine_hi; length = 41))

println("\nBest coarse θ: ", bestθ)
println("Running fine θ grid on [", fine_lo, ", ", fine_hi, "]...")
fine_results = EvaluateCalibrationGrid(;
    l, d, d′, wⱼ_data, lⱼ_data, β_target,
    α = α_fixed,
    η_grid = [η_seed],
    θ_grid = fine_θ_grid,
    fixed_primitives,
    inner_tol = 2e-5,
    inner_maxIter = 3000,
    continuation_steps = 5,
    employment_change,
    wage_center,
    moment_firm_mask,
    max_abs_moment = 10.0,
    verbose = true
)
CSV.write(joinpath(output_table_path, "sandbox_theta_grid_fixed_primitives_fine.csv"),
    fine_results; bom=true)

combined_results = vcat(coarse_results, fine_results)
sort!(combined_results, :objective)
CSV.write(joinpath(output_table_path, "sandbox_theta_grid_fixed_primitives_combined.csv"),
    combined_results; bom=true)

println("\nTop 10 fixed-primitives θ candidates:")
show(first(combined_results, min(10, nrow(combined_results))),
    allrows=true, allcols=true)
println()

finite_combined = filter(:objective => isfinite, combined_results)
nrow(finite_combined) > 0 || error("No finite objective values found in combined θ grid")
best = first(finite_combined, 1)
println("\nBest fixed-primitives θ result:")
show(best, allrows=true, allcols=true)
println()
