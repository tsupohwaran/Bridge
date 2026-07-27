#==================================================#
# Diagnose commute-profile nests
#
# Constructs proxy labor-pool nests from baseline commute times:
#     h_zj ∝ l_z * d_zj^(-lambda)
# clusters firms by h_.j, then computes
#     corr(w_diff_j, Ehat^T_j | control)
# where Ehat^T_j is treated-firm exposure inside the candidate nest.
#==================================================#

using LinearAlgebra
using Statistics
using Random
using StatFiles
using DataFrames
using CSV

const projPath = get(ENV, "PROJ_PATH", "/Users/pohwaran/Doctorate/Paper/Bridge")
const Z = parse(Int, get(ENV, "N_TOWNS", "128"))
const model_path = joinpath(projPath, "data", "model", "processed", "firm_qingdao_model_10.dta")
const reg_path = joinpath(projPath, "data", "model", "processed", "firm_two_year_reg.dta")
const out_path = joinpath(projPath, "output", "tables", "commute_profile_nest_exposure_diagnostics.csv")

function parse_grid(s::AbstractString, ::Type{T}) where {T}
    vals = strip.(split(s, ","))
    vals = filter(!isempty, vals)
    return parse.(T, vals)
end

const lambda_grid = parse_grid(get(ENV, "LAMBDA_GRID", "0.5,1,2,3"), Float64)
const k_grid = parse_grid(get(ENV, "K_GRID", "5,10,15,20"), Int)
const bigma_threshold = parse(Float64, get(ENV, "BIGMA_THRESHOLD", "0.5"))
const maxiter = parse(Int, get(ENV, "KMEANS_MAXITER", "100"))
const rng_seed = parse(Int, get(ENV, "KMEANS_SEED", "20260712"))

function safe_cor(x, y)
    mask = isfinite.(x) .& isfinite.(y)
    if sum(mask) < 3
        return missing
    end
    xx = Float64.(x[mask])
    yy = Float64.(y[mask])
    if iszero(std(xx)) || iszero(std(yy))
        return missing
    end
    return cor(xx, yy)
end

function moment_mask(firm_ids, reg_path)
    if !isfile(reg_path)
        @warn "Regression sample id file not found; matched-reg-sample correlations will be missing." reg_path
        return trues(length(firm_ids))
    end
    df_reg = DataFrame(load(reg_path))
    reg_ids = Set(collect(skipmissing(df_reg[!, :id])))
    return [!ismissing(id) && id in reg_ids for id in firm_ids]
end

function commute_profile_matrix(l, d, lambda)
    weights = l .* d .^ (-lambda)
    denom = vec(sum(weights, dims=1))
    any(denom .<= 0) && error("Each firm must have positive commute-profile denominator.")
    return weights ./ denom'
end

function squared_distances_to_centers(X, centers)
    J = size(X, 1)
    K = size(centers, 1)
    dist = Matrix{Float64}(undef, J, K)
    for k in 1:K
        c = @view centers[k, :]
        for j in 1:J
            x = @view X[j, :]
            dist[j, k] = sum(abs2, x .- c)
        end
    end
    return dist
end

function initialize_centers(X, K; rng)
    J = size(X, 1)
    K <= J || error("K cannot exceed number of firms.")
    centers = Matrix{Float64}(undef, K, size(X, 2))
    first_idx = rand(rng, 1:J)
    centers[1, :] .= X[first_idx, :]
    min_dist = vec(squared_distances_to_centers(X, centers[1:1, :]))

    for k in 2:K
        total = sum(min_dist)
        next_idx = if total <= eps(Float64)
            rand(rng, 1:J)
        else
            cutoff = rand(rng) * total
            acc = 0.0
            chosen = J
            for j in 1:J
                acc += min_dist[j]
                if acc >= cutoff
                    chosen = j
                    break
                end
            end
            chosen
        end
        centers[k, :] .= X[next_idx, :]
        new_dist = vec(squared_distances_to_centers(X, centers[k:k, :]))
        min_dist .= min.(min_dist, new_dist)
    end
    return centers
end

function kmeans_labels(X, K; maxiter=100, seed=20260712)
    rng = MersenneTwister(seed)
    J, P = size(X)
    centers = initialize_centers(X, K; rng)
    labels = zeros(Int, J)
    old_labels = fill(-1, J)

    for _ in 1:maxiter
        dist = squared_distances_to_centers(X, centers)
        for j in 1:J
            labels[j] = argmin(@view dist[j, :])
        end
        labels == old_labels && break
        old_labels .= labels

        counts = zeros(Int, K)
        fill!(centers, 0.0)
        for j in 1:J
            k = labels[j]
            counts[k] += 1
            centers[k, :] .+= @view X[j, :]
        end
        for k in 1:K
            if counts[k] == 0
                farthest = argmax(vec(minimum(dist, dims=2)))
                centers[k, :] .= X[farthest, :]
                labels[farthest] = k
                counts[k] = 1
            else
                centers[k, :] ./= counts[k]
            end
        end
    end
    return labels
end

function exposure_by_nest(H, nest_labels, treated, dMA)
    Z, J = size(H)
    exposure = zeros(J)
    for nest in unique(nest_labels)
        idx = findall(==(nest), nest_labels)
        treated_idx = idx[treated[idx]]
        isempty(treated_idx) && continue
        treated_profile = zeros(Z)
        for k in treated_idx
            treated_profile .+= H[:, k] .* dMA[k]
        end
        for j in idx
            exposure[j] = dot(@view(H[:, j]), treated_profile)
        end
    end
    return exposure
end

function nest_stats(nest_labels, treated)
    nests = unique(nest_labels)
    sizes = [sum(nest_labels .== n) for n in nests]
    treated_counts = [sum(treated[nest_labels .== n]) for n in nests]
    return (;
        n_nests = length(nests),
        min_nest_size = minimum(sizes),
        median_nest_size = median(sizes),
        max_nest_size = maximum(sizes),
        share_nests_with_treated = mean(treated_counts .> 0)
    )
end

println("Loading model sample: ", model_path)
df = DataFrame(load(model_path))
nrow(df) % Z == 0 || error("Model sample row count must be divisible by Z=$Z.")
J = nrow(df) ÷ Z

l = Float64.(df[!, :pop]) |>
    x -> reshape(x, Z, J) |>
    x -> x[:, 1] |>
    x -> x ./ sum(x)
l_j_data = Float64.(disallowmissing(df[!, :employ])) |>
    x -> reshape(x, Z, J) |>
    x -> x[1, :] |>
    x -> x ./ sum(x)
d = reshape(Float64.(df[!, :dzj]), Z, J) |> x -> replace(x, 0.0 => 1e-2)
d_post = reshape(Float64.(df[!, :dzj_prime]), Z, J) |> x -> replace(x, 0.0 => 1e-2)
w_raw = reshape(Float64.(df[!, :wage_inital]), Z, J)[1, :]
lo, hi = quantile(w_raw, [0.05, 0.95])
w = clamp.(w_raw, lo, hi) |> x -> x ./ sum(x .* l_j_data)
lnw = log.(max.(w, eps(Float64)))
w_diff = lnw .- mean(lnw)
firm_ids = collect(df[1:Z:end, :id])
firm_ind = Int.(df[1:Z:end, :ind_agg])

dMA = vec(sum((d .- d_post) .* reshape(l, :, 1), dims=1))
treated = dMA .>= bigma_threshold
control = .!treated
reg_keep = moment_mask(firm_ids, reg_path)

println("Firms: ", J)
println("Controls: ", sum(control), "; treated: ", sum(treated))
println("Matched regression-sample firms: ", sum(reg_keep))
println("lambda grid: ", lambda_grid)
println("K grid: ", k_grid)

rows = DataFrame[]
for lambda in lambda_grid
    H = commute_profile_matrix(reshape(l, :, 1), d, lambda)
    X = permutedims(sqrt.(H))  # rows are firms; Euclidean distance is Hellinger-style.

    for K in k_grid
        println("Clustering commute profiles: lambda=$lambda, K=$K")
        profile_cluster = kmeans_labels(X, K; maxiter, seed = rng_seed + round(Int, 1000lambda) + K)
        nest_labels = [string(ind, "#", cl) for (ind, cl) in zip(firm_ind, profile_cluster)]
        exposure = exposure_by_nest(H, nest_labels, treated, dMA)
        stats = nest_stats(nest_labels, treated)

        all_control = control
        reg_control = control .& reg_keep
        push!(rows, DataFrame(
            nest_definition = ["ind_agg_x_commute_profile_cluster"],
            lambda = [lambda],
            K_profile_clusters = [K],
            bigma_threshold = [bigma_threshold],
            n_firms = [J],
            n_control = [sum(all_control)],
            n_treated = [sum(treated)],
            n_reg_control = [sum(reg_control)],
            n_reg_firms = [sum(reg_keep)],
            n_nests = [stats.n_nests],
            min_nest_size = [stats.min_nest_size],
            median_nest_size = [stats.median_nest_size],
            max_nest_size = [stats.max_nest_size],
            share_nests_with_treated = [stats.share_nests_with_treated],
            corr_wdiff_treated_exposure_control_all = [
                safe_cor(w_diff[all_control], exposure[all_control])
            ],
            corr_wdiff_treated_exposure_control_reg = [
                safe_cor(w_diff[reg_control], exposure[reg_control])
            ],
            mean_exposure_control = [mean(exposure[all_control])],
            sd_exposure_control = [std(exposure[all_control])],
            mean_exposure_reg_control = [mean(exposure[reg_control])],
            sd_exposure_reg_control = [std(exposure[reg_control])]
        ))
    end
end

out = reduce(vcat, rows)
sort!(out, [:lambda, :K_profile_clusters])
mkpath(dirname(out_path))
CSV.write(out_path, out; bom=true)

println("\nCorrelation diagnostics:")
show(select(out, :lambda, :K_profile_clusters,
    :corr_wdiff_treated_exposure_control_all,
    :corr_wdiff_treated_exposure_control_reg,
    :n_nests, :median_nest_size, :share_nests_with_treated),
    allrows=true, allcols=true)
println("\nWrote diagnostics to: ", out_path)
