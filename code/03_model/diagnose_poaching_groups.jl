#==================================================#
# Diagnostic: static baseline-employer-group poaching model
#==================================================#

projPath = normpath(joinpath(@__DIR__, "..", ".."))
cd(projPath)

include(joinpath(@__DIR__, "load_packages.jl"))
include(joinpath(@__DIR__, "functions.jl"))

const Z = parse(Int, get(ENV, "MODEL_Z", "128"))
const eta = parse(Float64, get(ENV, "ETA", "0.9"))
const theta = parse(Float64, get(ENV, "THETA", "10.0"))
const sigma = parse(Float64, get(ENV, "SIGMA", "0.25"))
const alpha = parse(Float64, get(ENV, "ALPHA", "0.8"))
const n_wage_groups = parse(Int, get(ENV, "WAGE_GROUPS", "10"))
const run_inversion = get(ENV, "RUN_INVERSION", "0") in ["1", "true", "TRUE", "yes", "YES"]
const amenity_tol = parse(Float64, get(ENV, "AMENITY_TOL", "1e-5"))
const amenity_maxIter = parse(Int, get(ENV, "AMENITY_MAXITER", "500"))
const amenity_damp_env = get(ENV, "AMENITY_DAMP", "")
const amenity_damp = isempty(amenity_damp_env) ? nothing : parse(Float64, amenity_damp_env)

model_path = joinpath(projPath, "data", "model", "processed", "firm_qingdao_model_10.dta")
reg_path = joinpath(projPath, "data", "model", "processed", "firm_two_year_reg.dta")
output_table_path = joinpath(projPath, "output", "tables")
output_name = get(ENV, "OUTPUT_NAME", "poaching_group_beta_diagnostics.csv")
mkpath(output_table_path)

function parse_grid(s::AbstractString)
    vals = strip.(split(s, ","))
    vals = filter(!isempty, vals)
    return parse.(Float64, vals)
end

const rho_grid = parse_grid(get(ENV, "RHO_GRID", "0,0.25,0.5,1,2,3,5"))
const kappa_grid = parse_grid(get(ENV, "KAPPA_GRID", "0,1"))

function rank_groups(x; n_groups::Integer=10)
    n = length(x)
    order = sortperm(x)
    groups = Vector{Int}(undef, n)
    for (rank, idx) in enumerate(order)
        groups[idx] = min(n_groups, max(1, cld(rank * n_groups, n)))
    end
    return groups
end

function weighted_mean(x, w)
    denom = sum(w)
    denom > 0 || return mean(x)
    return sum(x .* w) / denom
end

function choice_from_logq(log_q_zj, firm_sector, theta, sigma)
    J = size(log_q_zj, 2)
    sector_indices = _firm_sector_indices(firm_sector, J)

    if isnothing(sector_indices)
        isapprox(sigma, 1.0; atol=1e-12) ||
            error("firm_sector is required when sigma differs from 1")
        log_xzj = theta .* log_q_zj
        log_xzj_max = maximum(log_xzj, dims=2)
        log_sum_exp = log_xzj_max .+ log.(sum(exp.(log_xzj .- log_xzj_max), dims=2))
        return exp.(log_xzj .- log_sum_exp)
    end

    x_zj = (theta / sigma) .* log_q_zj
    Z_local = size(log_q_zj, 1)
    S = length(sector_indices)
    log_inclusive_zs = Matrix{Float64}(undef, Z_local, S)
    pi_zj_given_s = similar(x_zj)

    for (s_idx, idx) in enumerate(sector_indices)
        x_zs = @view x_zj[:, idx]
        x_max_zs = maximum(x_zs, dims=2)
        log_sum_zs = x_max_zs .+ log.(sum(exp.(x_zs .- x_max_zs), dims=2))
        log_inclusive_zs[:, s_idx] = vec(sigma .* log_sum_zs)
        pi_zj_given_s[:, idx] = exp.(x_zs .- log_sum_zs)
    end

    inclusive_max = maximum(log_inclusive_zs, dims=2)
    log_denom_z = inclusive_max .+
        log.(sum(exp.(log_inclusive_zs .- inclusive_max), dims=2))
    pi_zs = exp.(log_inclusive_zs .- log_denom_z)

    pi_zj = similar(x_zj)
    for (s_idx, idx) in enumerate(sector_indices)
        pi_zj[:, idx] = pi_zj_given_s[:, idx] .* pi_zs[:, s_idx:s_idx]
    end
    return pi_zj
end

function labor_wdiff_regression(lnl_base, lnl_counterfactual, bigMA, lnw;
    keep, wage_center=:all)

    keep_mask = Bool.(keep)
    big_keep = Float64.(bigMA[keep_mask])
    w_keep = Float64.(lnw[keep_mask])
    y = Float64.(lnl_counterfactual[keep_mask] .- lnl_base[keep_mask])

    if wage_center == :treated
        treated = big_keep .== 1
        w_center = any(treated) ? mean(w_keep[treated]) : mean(w_keep)
    elseif wage_center == :all
        w_center = mean(w_keep)
    else
        error("wage_center must be :all or :treated")
    end

    w_diff = w_keep .- w_center
    X = hcat(ones(length(y)), w_diff, big_keep, big_keep .* w_diff)
    coef_vec = X \ y
    residual = y .- X * coef_vec
    dof = max(length(y) - size(X, 2), 1)
    sigma2 = sum(residual .^ 2) / dof
    vcov = sigma2 .* inv(X' * X)
    se = sqrt.(diag(vcov))

    return (;
        beta_labor_wdiff_post = coef_vec[2],
        beta_labor_wdiff_post_se = se[2],
        beta_labor_bigMA = coef_vec[3],
        beta_labor_bigMA_se = se[3],
        beta_labor_bigMA_wdiff = coef_vec[4],
        beta_labor_bigMA_wdiff_se = se[4],
        treated_wdiff_slope = coef_vec[2] + coef_vec[4],
        n = length(y),
        n_control = sum(big_keep .== 0),
        n_treated = sum(big_keep .== 1),
        wage_center_value = w_center
    )
end

function group_metadata(firm_group, firm_ind, lnw, l_j_data)
    group_values = sort(unique(firm_group))
    group_index = Dict(g => i for (i, g) in enumerate(group_values))
    group_id = [group_index[g] for g in firm_group]
    G = length(group_values)

    group_lnw = Vector{Float64}(undef, G)
    group_ind = Vector{Int}(undef, G)
    group_mass = Vector{Float64}(undef, G)
    for g in 1:G
        idx = findall(==(g), group_id)
        group_lnw[g] = weighted_mean(lnw[idx], l_j_data[idx])
        group_ind[g] = firm_ind[idx[1]]
        group_mass[g] = sum(l_j_data[idx])
    end

    return (; group_id, group_values, group_lnw, group_ind, group_mass)
end

function baseline_group_origin_mass(pi_zj, l, l_j_model, l_j_data, group_id, G)
    gamma_zj = pi_zj .* l ./ vec(l_j_model)'
    M_zg = zeros(size(pi_zj, 1), G)
    for j in axes(pi_zj, 2)
        M_zg[:, group_id[j]] .+= l_j_data[j] .* gamma_zj[:, j]
    end
    return M_zg
end

function post_labor_poaching(M_zg, d_post, w_j, a_j, firm_sector, firm_ind,
    group_lnw, group_ind, rho, kappa)

    lnw = log.(max.(w_j, eps(Float64)))
    log_d = log.(d_post)
    log_base = lnw' .+ a_j' .- eta .* log_d
    L_post = zeros(length(w_j))

    for g in eachindex(group_lnw)
        upward = max.(lnw .- group_lnw[g], 0.0)
        industry_penalty = Float64.(firm_ind .!= group_ind[g])
        group_shift = rho .* upward .- kappa .* industry_penalty
        pi_zj = choice_from_logq(log_base .+ group_shift', firm_sector, theta, sigma)
        L_post .+= vec(sum(pi_zj .* M_zg[:, g], dims=1))
    end

    return L_post
end

println("Loading model sample: ", model_path)
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
d_post = reshape(Float64.(df[!, :dzj_prime]), Z, J) |> x -> replace(x, 0.0 => 1e-2)

w_j_raw = reshape(Float64.(df[!, :wage_inital]), Z, J)[1, :]
lo, hi = quantile(w_j_raw, [0.05, 0.95])
w_j_data = clamp.(w_j_raw, lo, hi) |>
    x -> x ./ sum(x .* l_j_data)
lnw = log.(max.(w_j_data, eps(Float64)))

firm_ids = collect(df[1:Z:end, :id])
firm_ind = Int.(df[1:Z:end, :ind_agg])
firm_sector = firm_ind
wage_group = rank_groups(lnw; n_groups = n_wage_groups)
firm_group = firm_ind .* (n_wage_groups + 1) .+ wage_group
groups = group_metadata(firm_group, firm_ind, lnw, l_j_data)

println("Loading moment sample ids: ", reg_path)
df_reg = DataFrame(load(reg_path))
reg_sample_ids = collect(skipmissing(df_reg[!, :id]))
keep_mask = MomentFirmMask(J; firm_ids, reg_sample_ids, restrict_to_reg_sample = true)

dMA = vec(l' * d .- l' * d_post) |> x -> replace(x, -Inf => -8)
bigMA = Float64.(dMA .>= 0.5)

params = (; η = eta, θ = theta, σ = sigma, α = alpha)
if run_inversion
    println("Recovering baseline amenities for origin shares...")
    primitives = SolveFirmPrimitivesFromData((; wⱼ = w_j_data, lⱼ = l_j_data,
        l, d, firm_sector), params;
        aⱼ_init = zeros(J),
        amenity_tol,
        amenity_maxIter,
        amenity_damp,
        displaySummary = true,
        displayGap = false,
        returnInfo = true)
    a_j = vec(primitives.aⱼ)
    baseline_choice = (; π_zj = primitives.π_zj, lⱼ = primitives.lⱼ)
    baseline_mode = "inverted_amenities"
    amenity_converged = primitives.amenity_converged
    amenity_iterations = primitives.amenity_iterations
    amenity_gap = primitives.amenity_gap
else
    a_j = zeros(J)
    baseline_choice = WorkerChoice(w_j_data, l, d, a_j, params; firm_sector)
    baseline_mode = "observed_wages_zero_amenities"
    amenity_converged = missing
    amenity_iterations = missing
    amenity_gap = missing
end

M_zg = baseline_group_origin_mass(baseline_choice.π_zj, l, baseline_choice.lⱼ,
    l_j_data, groups.group_id, length(groups.group_values))
group_mass_gap = maximum(abs.(vec(sum(M_zg, dims=1)) .- groups.group_mass))
total_mass_gap = abs(sum(M_zg) - sum(l_j_data))

println("Diagnostic parameters: eta=$eta, theta=$theta, sigma=$sigma, alpha=$alpha")
println("rho grid: ", rho_grid)
println("kappa grid: ", kappa_grid)
println("Baseline mode: ", baseline_mode)
println("Employer groups: ", length(groups.group_values), " (ind x wage group)")
println("Kept firms: ", sum(keep_mask), "; controls kept: ",
    sum(keep_mask .& (bigMA .== 0)), "; treated kept: ",
    sum(keep_mask .& (bigMA .== 1)))
println("Group mass gap: ", group_mass_gap, "; total mass gap: ", total_mass_gap)

lnl_base = log.(max.(l_j_data, eps(Float64)))
rows = DataFrame[]

for kappa in kappa_grid
    for rho in rho_grid
        println("Running poaching transition: rho=$rho, kappa=$kappa")
        L_post = post_labor_poaching(M_zg, d_post, w_j_data, a_j, firm_sector,
            firm_ind, groups.group_lnw, groups.group_ind, rho, kappa)
        reg = labor_wdiff_regression(lnl_base, log.(max.(L_post, eps(Float64))),
            bigMA, lnw; keep = keep_mask, wage_center = :all)

        push!(rows, DataFrame(
            model_mode = ["static_group_transition"],
            baseline_mode = [baseline_mode],
            employer_group = ["ind_x_wage_group"],
            eta = [eta],
            theta = [theta],
            sigma = [sigma],
            alpha = [alpha],
            rho = [rho],
            kappa = [kappa],
            n_wage_groups = [n_wage_groups],
            n_employer_groups = [length(groups.group_values)],
            n = [reg.n],
            n_control = [reg.n_control],
            n_treated = [reg.n_treated],
            amenity_converged = [amenity_converged],
            amenity_iterations = [amenity_iterations],
            amenity_gap = [amenity_gap],
            group_mass_gap = [group_mass_gap],
            total_mass_gap = [total_mass_gap],
            beta_labor_wdiff_post = [reg.beta_labor_wdiff_post],
            beta_labor_wdiff_post_se = [reg.beta_labor_wdiff_post_se],
            beta_labor_bigMA = [reg.beta_labor_bigMA],
            beta_labor_bigMA_se = [reg.beta_labor_bigMA_se],
            beta_labor_bigMA_wdiff = [reg.beta_labor_bigMA_wdiff],
            beta_labor_bigMA_wdiff_se = [reg.beta_labor_bigMA_wdiff_se],
            treated_wdiff_slope = [reg.treated_wdiff_slope]
        ))
    end
end

out = vcat(rows...)
out_path = joinpath(output_table_path, output_name)
CSV.write(out_path, out; bom = true)

println("\nKey beta_labor_wdiff_post results:")
show(select(out, :rho, :kappa, :baseline_mode, :beta_labor_wdiff_post,
    :beta_labor_bigMA_wdiff, :treated_wdiff_slope),
    allrows = true, allcols = true)
println("\nWrote poaching diagnostics to: ", out_path)
