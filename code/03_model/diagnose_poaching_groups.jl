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
const run_inversion = get(ENV, "RUN_INVERSION", "1") in ["1", "true", "TRUE", "yes", "YES"]
const run_ge = get(ENV, "RUN_GE", "1") in ["1", "true", "TRUE", "yes", "YES"]
const use_fast_amenity = get(ENV, "FAST_AMENITY", "1") in ["1", "true", "TRUE", "yes", "YES"]
const amenity_tol = parse(Float64, get(ENV, "AMENITY_TOL", "1e-5"))
const amenity_maxIter = parse(Int, get(ENV, "AMENITY_MAXITER", "2000"))
const amenity_damp_env = get(ENV, "AMENITY_DAMP", "")
const amenity_damp = isempty(amenity_damp_env) ? 0.6 : parse(Float64, amenity_damp_env)
const ge_tol = parse(Float64, get(ENV, "GE_TOL", "1e-4"))
const ge_maxIter = parse(Int, get(ENV, "GE_MAXITER", "300"))
const ge_damp = parse(Float64, get(ENV, "GE_DAMP", "0.95"))
const ge_power = get(ENV, "GE_POWER", "0") in ["1", "true", "TRUE", "yes", "YES"]
const ge_display = get(ENV, "DISPLAY_GE_GAP", "0") in ["1", "true", "TRUE", "yes", "YES"]
const ge_log_step = parse(Float64, get(ENV, "GE_LOG_STEP", "0.30"))

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

const rho_grid = parse_grid(get(ENV, "RHO_GRID", "0,0.25,0.5,1"))
const kappa_grid = parse_grid(get(ENV, "KAPPA_GRID", "0"))

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

function choice_from_logq(log_q_zj, firm_sector, theta, sigma; sector_indices=nothing)
    J = size(log_q_zj, 2)
    sector_indices = isnothing(sector_indices) ? _firm_sector_indices(firm_sector, J) :
        sector_indices

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

function choice_and_own_elasticity_from_logq(log_q_zj, firm_sector, theta, sigma,
    lambda_j; sector_indices=nothing)

    J = size(log_q_zj, 2)
    sector_indices = isnothing(sector_indices) ? _firm_sector_indices(firm_sector, J) :
        sector_indices

    if isnothing(sector_indices)
        isapprox(sigma, 1.0; atol=1e-12) ||
            error("firm_sector is required when sigma differs from 1")
        log_xzj = theta .* log_q_zj
        log_xzj_max = maximum(log_xzj, dims=2)
        log_sum_exp = log_xzj_max .+ log.(sum(exp.(log_xzj .- log_xzj_max), dims=2))
        pi_zj = exp.(log_xzj .- log_sum_exp)
        eps_zj = lambda_j' .* theta .* (1 .- pi_zj)
        return pi_zj, eps_zj
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

    eps_zj = lambda_j' .* theta .*
        (1 / sigma .+ (1 - 1 / sigma) .* pi_zj_given_s .- pi_zj)
    return pi_zj, eps_zj
end

function SolveAmenitiesFromEmploymentFast(w_j, l, d, l_j_target, params, firm_sector;
    a_j_init=nothing, tol=1e-6, maxIter=2000, damp=0.6,
    displayGap=false, displaySummary=true)

    (; η, θ, σ) = params
    J = size(d, 2)
    sector_indices = _firm_sector_indices(firm_sector, J)
    target = vec(Float64.(l_j_target))
    log_target = log.(target)
    log_w = log.(max.(vec(Float64.(w_j)), eps(Float64)))
    log_d = log.(d)
    logq = isnothing(a_j_init) ? zeros(J) : θ .* vec(Float64.(a_j_init))
    logq .-= mean(logq)

    iter = 0
    gap = Inf
    pi_zj = zeros(size(d))
    l_model = similar(target)

    while iter < maxIter && gap > tol
        iter += 1
        a_j = logq ./ θ
        pi_zj = choice_from_logq(log_w' .+ a_j' .- η .* log_d,
            firm_sector, θ, σ; sector_indices)
        l_model = vec(sum(pi_zj .* l, dims=1))
        gap = maximum(abs.(l_model .- target))
        gap <= tol && break

        update = log_target .- log.(max.(l_model, eps(Float64)))
        logq_new = logq .+ update
        logq_new .-= mean(logq_new)
        logq = damp .* logq .+ (1 - damp) .* logq_new
        logq .-= mean(logq)

        if displayGap && (iter <= 10 || iter % 50 == 0)
            println("Fast amenity inversion iteration: ", iter, ", Gap: ", gap)
        end
    end

    a_j = logq ./ θ
    a_j .-= mean(a_j)
    choice = WorkerChoice(w_j, l, d, a_j, params;
        log_d=log_d, firm_sector=firm_sector, sector_indices=sector_indices)
    final_gap = maximum(abs.(vec(choice.lⱼ) .- target))
    converged = final_gap <= tol

    if displaySummary
        if converged
            println("Successful fast amenity inversion in ", iter,
                " iterations. Final gap: ", final_gap)
        else
            println("Fast amenity inversion reached max iterations. Final gap: ",
                final_gap)
        end
    end

    return (; aⱼ = a_j, choice..., lⱼ_target = target, converged,
        iterations = iter, gap = final_gap, damp)
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
    sector_indices = _firm_sector_indices(firm_sector, length(w_j))

    for g in eachindex(group_lnw)
        upward = max.(lnw .- group_lnw[g], 0.0)
        industry_penalty = Float64.(firm_ind .!= group_ind[g])
        group_shift = rho .* upward .- kappa .* industry_penalty
        pi_zj = choice_from_logq(log_base .+ group_shift', firm_sector, theta, sigma;
            sector_indices)
        L_post .+= vec(sum(pi_zj .* M_zg[:, g], dims=1))
    end

    return L_post
end

function poaching_labor_and_elasticity(M_zg, d_current, w_j, a_j, firm_sector,
    firm_ind, group_lnw, group_ind, rho, kappa)

    lnw = log.(max.(w_j, eps(Float64)))
    log_d = log.(d_current)
    log_base = lnw' .+ a_j' .- eta .* log_d
    J = length(w_j)
    L = zeros(J)
    eps_num = zeros(J)
    sector_indices = _firm_sector_indices(firm_sector, J)

    for g in eachindex(group_lnw)
        upward = max.(lnw .- group_lnw[g], 0.0)
        lambda_j = 1 .+ rho .* Float64.(lnw .> group_lnw[g])
        industry_penalty = Float64.(firm_ind .!= group_ind[g])
        group_shift = rho .* upward .- kappa .* industry_penalty
        pi_zj, eps_zj = choice_and_own_elasticity_from_logq(
            log_base .+ group_shift', firm_sector, theta, sigma, lambda_j;
            sector_indices)
        mass = M_zg[:, g]
        L .+= vec(sum(pi_zj .* mass, dims=1))
        eps_num .+= vec(sum(pi_zj .* mass .* eps_zj, dims=1))
    end

    epsilon = eps_num ./ max.(L, eps(Float64))
    return (; L, epsilon)
end

function make_poaching_buffers(Z_local, J, S)
    return (;
        x = Matrix{Float64}(undef, Z_local, J),
        pi_given = Matrix{Float64}(undef, Z_local, J),
        log_inclusive = Matrix{Float64}(undef, Z_local, S),
        pi_zs = Matrix{Float64}(undef, Z_local, S),
        L = zeros(J),
        eps_num = zeros(J)
    )
end

function poaching_labor_and_elasticity_fast!(buffers, M_zg, log_d, w_j, a_j,
    sector_indices, firm_ind, group_lnw, group_ind, rho, kappa)

    lnw = log.(max.(w_j, eps(Float64)))
    Z_local = size(log_d, 1)
    J = length(w_j)
    fill!(buffers.L, 0.0)
    fill!(buffers.eps_num, 0.0)

    for g in eachindex(group_lnw)
        upward = max.(lnw .- group_lnw[g], 0.0)
        lambda_j = 1 .+ rho .* Float64.(lnw .> group_lnw[g])
        industry_penalty = Float64.(firm_ind .!= group_ind[g])
        group_shift = rho .* upward .- kappa .* industry_penalty

        buffers.x .= (theta / sigma) .*
            (lnw' .+ a_j' .- eta .* log_d .+ group_shift')

        for (s_idx, idx) in enumerate(sector_indices)
            x_zs = @view buffers.x[:, idx]
            x_max_zs = maximum(x_zs, dims=2)
            log_sum_zs = x_max_zs .+ log.(sum(exp.(x_zs .- x_max_zs), dims=2))
            buffers.log_inclusive[:, s_idx] = vec(sigma .* log_sum_zs)
            buffers.pi_given[:, idx] = exp.(x_zs .- log_sum_zs)
        end

        inclusive_max = maximum(buffers.log_inclusive, dims=2)
        log_denom_z = inclusive_max .+
            log.(sum(exp.(buffers.log_inclusive .- inclusive_max), dims=2))
        buffers.pi_zs .= exp.(buffers.log_inclusive .- log_denom_z)

        mass = M_zg[:, g]
        for (s_idx, idx) in enumerate(sector_indices)
            pi_given = @view buffers.pi_given[:, idx]
            pi_sector = pi_given .* buffers.pi_zs[:, s_idx:s_idx]
            eps_sector = lambda_j[idx]' .* theta .*
                (1 / sigma .+ (1 - 1 / sigma) .* pi_given .- pi_sector)
            buffers.L[idx] .+= vec(sum(pi_sector .* mass, dims=1))
            buffers.eps_num[idx] .+= vec(sum(pi_sector .* mass .* eps_sector, dims=1))
        end
    end

    epsilon = buffers.eps_num ./ max.(buffers.L, eps(Float64))
    return (; L = copy(buffers.L), epsilon)
end

function solve_poaching_ge(M_zg, d_current, w_init, z_j, a_j, firm_sector,
    firm_ind, group_lnw, group_ind, rho, kappa; tol=2e-5, maxIter=300,
    damp=0.90, power=false, displayGap=false, maxLogStep=0.50)

    if iszero(rho) && iszero(kappa)
        origin_mass = vec(sum(M_zg, dims=2))
        sol = SolveModel((; l = origin_mass, d = d_current, zⱼ = z_j,
                aⱼ = a_j, firm_sector),
            (; η = eta, θ = theta, σ = sigma, α = alpha);
            damp, tol, maxIter, power, displayGap,
            displaySummary = false,
            wⱼ_init = w_init,
            returnInfo = true)
        return (; w = vec(sol.wⱼ), L = vec(sol.lⱼ), epsilon = vec(sol.εⱼ),
            converged = sol.converged, iterations = sol.iterations, gap = sol.gap)
    end

    w0 = vec(Float64.(w_init))
    w0 ./= mean(w0)
    logw = log.(max.(w0, eps(Float64)))
    iter = 0
    gap = Inf
    agg = nothing
    w_update = similar(w0)
    sector_indices = _firm_sector_indices(firm_sector, length(w0))
    log_d = log.(d_current)
    buffers = make_poaching_buffers(size(d_current, 1), size(d_current, 2),
        length(sector_indices))

    while iter < maxIter && gap > tol
        iter += 1
        w = exp.(logw)
        w ./= mean(w)
        logw = log.(max.(w, eps(Float64)))
        agg = poaching_labor_and_elasticity_fast!(buffers, M_zg, log_d, w, a_j,
            sector_indices, firm_ind, group_lnw, group_ind, rho, kappa)
        w_update = alpha .* z_j .* agg.L .^ (alpha - 1) .* agg.epsilon ./
            (1 .+ agg.epsilon)
        w_update ./= sum(w_update .* agg.L)
        w_update ./= mean(w_update)
        new_logw = log.(max.(w_update, eps(Float64)))
        delta = new_logw .- logw
        gap = maximum(abs.(delta))
        delta = clamp.(delta, -maxLogStep, maxLogStep)
        logw .+= (1 - damp) .* delta
        logw .-= log(mean(exp.(logw)))

        if displayGap && (iter <= 10 || iter % 25 == 0)
            println("Poaching GE iteration: ", iter, ", Gap: ", gap)
        end
    end

    w = exp.(logw)
    w ./= mean(w)
    agg = poaching_labor_and_elasticity_fast!(buffers, M_zg, log_d, w, a_j,
        sector_indices, firm_ind, group_lnw, group_ind, rho, kappa)
    w_final = alpha .* z_j .* agg.L .^ (alpha - 1) .* agg.epsilon ./
        (1 .+ agg.epsilon)
    w_final ./= sum(w_final .* agg.L)

    return (; w = w_final, L = agg.L, epsilon = agg.epsilon,
        converged = gap <= tol, iterations = iter, gap)
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
    if use_fast_amenity
        primitives = SolveAmenitiesFromEmploymentFast(w_j_data, l, d, l_j_data,
            params, firm_sector;
            a_j_init = zeros(J),
            tol = amenity_tol,
            maxIter = amenity_maxIter,
            damp = amenity_damp,
            displayGap = false,
            displaySummary = true)
        a_j = vec(primitives.aⱼ)
        baseline_choice = (; π_zj = primitives.π_zj, lⱼ = primitives.lⱼ)
        baseline_mode = "fast_inverted_amenities"
        amenity_converged = primitives.converged
        amenity_iterations = primitives.iterations
        amenity_gap = primitives.gap
        eps_j_base = vec(primitives.εⱼ)
        l_j_base_model = vec(primitives.lⱼ)
        z_j = w_j_data .* (1 .+ eps_j_base) ./
            (alpha .* l_j_base_model .^ (alpha - 1) .* eps_j_base)
        z_j ./= mean(z_j)
    else
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
        z_j = vec(primitives.zⱼ)
    end
else
    a_j = zeros(J)
    baseline_choice = WorkerChoice(w_j_data, l, d, a_j, params; firm_sector)
    baseline_mode = "observed_wages_zero_amenities"
    amenity_converged = missing
    amenity_iterations = missing
    amenity_gap = missing
    eps_j_base = vec(baseline_choice.εⱼ)
    l_j_base_model = vec(baseline_choice.lⱼ)
    z_j = w_j_data .* (1 .+ eps_j_base) ./
        (alpha .* l_j_base_model .^ (alpha - 1) .* eps_j_base)
    z_j ./= mean(z_j)
end

M_zg = baseline_group_origin_mass(baseline_choice.π_zj, l, baseline_choice.lⱼ,
    l_j_data, groups.group_id, length(groups.group_values))
group_mass_gap = maximum(abs.(vec(sum(M_zg, dims=1)) .- groups.group_mass))
total_mass_gap = abs(sum(M_zg) - sum(l_j_data))

println("Diagnostic parameters: eta=$eta, theta=$theta, sigma=$sigma, alpha=$alpha")
println("rho grid: ", rho_grid)
println("kappa grid: ", kappa_grid)
println("Baseline mode: ", baseline_mode)
println("GE mode: ", run_ge ? "transition + full GE" : "transition only")
if run_ge
    println("GE solver: tol=$ge_tol, maxIter=$ge_maxIter, damp=$ge_damp, power=$ge_power")
end
println("Employer groups: ", length(groups.group_values), " (ind x wage group)")
println("Kept firms: ", sum(keep_mask), "; controls kept: ",
    sum(keep_mask .& (bigMA .== 0)), "; treated kept: ",
    sum(keep_mask .& (bigMA .== 1)))
println("Group mass gap: ", group_mass_gap, "; total mass gap: ", total_mass_gap)

lnl_base = log.(max.(l_j_data, eps(Float64)))
rows = DataFrame[]

for kappa in kappa_grid
    ge_warm_start = copy(w_j_data)
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
            ge_converged = [missing],
            ge_iterations = [missing],
            ge_gap = [missing],
            beta_labor_wdiff_post = [reg.beta_labor_wdiff_post],
            beta_labor_wdiff_post_se = [reg.beta_labor_wdiff_post_se],
            beta_labor_bigMA = [reg.beta_labor_bigMA],
            beta_labor_bigMA_se = [reg.beta_labor_bigMA_se],
            beta_labor_bigMA_wdiff = [reg.beta_labor_bigMA_wdiff],
            beta_labor_bigMA_wdiff_se = [reg.beta_labor_bigMA_wdiff_se],
            treated_wdiff_slope = [reg.treated_wdiff_slope]
        ))

        if run_ge
            println("Running full GE poaching transition: rho=$rho, kappa=$kappa")
            ge = solve_poaching_ge(M_zg, d_post, ge_warm_start, z_j, a_j, firm_sector,
                firm_ind, groups.group_lnw, groups.group_ind, rho, kappa;
                tol = ge_tol,
                maxIter = ge_maxIter,
                damp = ge_damp,
                power = ge_power,
                displayGap = ge_display,
                maxLogStep = ge_log_step)
            ge_warm_start = ge.w
            reg_ge = labor_wdiff_regression(lnl_base,
                log.(max.(ge.L, eps(Float64))), bigMA, lnw;
                keep = keep_mask, wage_center = :all)

            push!(rows, DataFrame(
                model_mode = ["full_ge_group_transition"],
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
                n = [reg_ge.n],
                n_control = [reg_ge.n_control],
                n_treated = [reg_ge.n_treated],
                amenity_converged = [amenity_converged],
                amenity_iterations = [amenity_iterations],
                amenity_gap = [amenity_gap],
                group_mass_gap = [group_mass_gap],
                total_mass_gap = [total_mass_gap],
                ge_converged = [ge.converged],
                ge_iterations = [ge.iterations],
                ge_gap = [ge.gap],
                beta_labor_wdiff_post = [reg_ge.beta_labor_wdiff_post],
                beta_labor_wdiff_post_se = [reg_ge.beta_labor_wdiff_post_se],
                beta_labor_bigMA = [reg_ge.beta_labor_bigMA],
                beta_labor_bigMA_se = [reg_ge.beta_labor_bigMA_se],
                beta_labor_bigMA_wdiff = [reg_ge.beta_labor_bigMA_wdiff],
                beta_labor_bigMA_wdiff_se = [reg_ge.beta_labor_bigMA_wdiff_se],
                treated_wdiff_slope = [reg_ge.treated_wdiff_slope]
            ))
        end
    end
end

out = vcat(rows...)
out_path = joinpath(output_table_path, output_name)
CSV.write(out_path, out; bom = true)

println("\nKey beta_labor_wdiff_post results:")
show(select(out, :rho, :kappa, :model_mode, :baseline_mode,
    :beta_labor_wdiff_post, :beta_labor_bigMA_wdiff, :treated_wdiff_slope,
    :ge_converged, :ge_gap),
    allrows = true, allcols = true)
println("\nWrote poaching diagnostics to: ", out_path)
