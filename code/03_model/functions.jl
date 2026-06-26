#==================================================#
# The main function of the paper
# Date: March 2025
# Author: Wu Chengjun, Central University of Finance and Economics
# OS: MacOS 15.3.2
# Version: 1.10.9
#==================================================#

# For convenience. run Stata dofiles from Julia
function RunStata(projPath::String, stataPath::String, dofileName::String)

    # define paths
    dofile = joinpath(projPath, dofileName)
    wrapper = joinpath(tempdir(), basename(dofileName))

    # write wrapper.do
    open(wrapper, "w") do f
        write(f, """
        global proj_path "$projPath"
        do "$dofile"
        """)
    end

    # run wrapper.do
    try
        run(`$stataPath -b do $wrapper`)
    finally
        # delete wrapper.do
        isfile(wrapper) && rm(wrapper, force=true)
    end
end

# For convenience. This function is used to summarize statistics of an array.
function SumStats(arr, label)
    v = vec(arr)
    println("Summary statistics for ", label, " :")
    println("Mean: ", mean(v))
    println("Median: ", median(v))
    println("Min: ", minimum(v))
    println("Max: ", maximum(v))
    println("Standard Deviation: ", std(v))
    println("25th Percentile: ", quantile(v, 0.25))
    println("75th Percentile: ", quantile(v, 0.75))
    println("--------------------------------------------------")
    return nothing
end

# For convenience. This function is used to sum over the specified dimensions and drop them.
sumsqueeze(A; dims) = dropdims(sum(A, dims=dims), dims=dims)

function WorkerChoice(wⱼ, l, d, aⱼ, params::NamedTuple; log_d=nothing)
    (; η, θ) = params
    J = size(d, 2)
    length(wⱼ) == J || error("wⱼ must have length $J")
    length(aⱼ) == J || error("aⱼ must have length $J")
    any(ismissing, wⱼ) && error("wⱼ cannot contain missing values")
    any(ismissing, aⱼ) && error("aⱼ cannot contain missing values")

    wⱼ = wⱼ isa Vector{Float64} ? wⱼ : vec(Float64.(wⱼ))
    aⱼ = aⱼ isa Vector{Float64} ? aⱼ : vec(Float64.(aⱼ))
    log_d = isnothing(log_d) ? log.(d) : log_d

    # Use log-sum-exp trick for numerical stability
    log_xzj = θ .* (log.(wⱼ') .+ aⱼ' .- η .* log_d)  # Z × J matrix
    log_xzj_max = maximum(log_xzj, dims=2)  # Z × 1
    log_sum_exp = log_xzj_max .+ log.(sum(exp.(log_xzj .- log_xzj_max), dims=2))
    π_zj = exp.(log_xzj .- log_sum_exp)
    
    ε_zj = θ .* (1 .- π_zj)
    lⱼ = sum(π_zj .* l, dims = 1)'
    εⱼ = sumsqueeze(π_zj .* l ./ lⱼ' .* ε_zj, dims = 1)

    return (; π_zj, ε_zj, lⱼ, εⱼ)
end

function SolveAmenitiesFromEmployment(vars::NamedTuple, params::NamedTuple;
    aⱼ_init=nothing, tol=1e-10, maxIter=5000, damp=0.5,
    displayGap=false, displaySummary=false, returnInfo=false)

    (; wⱼ, l, d, lⱼ) = vars
    (; θ) = params
    J = size(d, 2)
    θ > 0 || error("θ must be positive")
    length(lⱼ) == J || error("lⱼ must have length $J")
    any(ismissing, lⱼ) && error("lⱼ cannot contain missing values")

    lⱼ_target = vec(Float64.(lⱼ))
    any(lⱼ_target .<= 0) && error("lⱼ must be strictly positive to invert finite amenities")
    logq = isnothing(aⱼ_init) ? zeros(J) : θ .* vec(Float64.(aⱼ_init))
    length(logq) == J || error("aⱼ_init must have length $J")
    logq .-= mean(logq)
    log_d = log.(d)

    iter = 0
    gap = Inf
    while (iter < maxIter) && (gap > tol)
        iter += 1
        aⱼ = logq ./ θ
        choice = WorkerChoice(wⱼ, l, d, aⱼ, params; log_d)
        lⱼ_model = vec(choice.lⱼ)
        gap = maximum(abs.(lⱼ_model .- lⱼ_target))
        gap <= tol && break

        update = log.(lⱼ_target) .- log.(max.(lⱼ_model, eps(Float64)))
        logq_new = logq .+ update
        logq_new .-= mean(logq_new)
        logq = damp .* logq .+ (1 - damp) .* logq_new
        logq .-= mean(logq)

        if displayGap
            println("Amenity inversion iteration: ", iter, ", Gap: ", gap)
        end
    end

    aⱼ = logq ./ θ
    aⱼ .-= mean(aⱼ)
    choice = WorkerChoice(wⱼ, l, d, aⱼ, params; log_d)
    gap = maximum(abs.(vec(choice.lⱼ) .- lⱼ_target))
    converged = gap <= tol

    if displaySummary
        if converged
            println("Successful amenity inversion in ", iter, " iterations. Final gap: ", gap)
        else
            println("Amenity inversion reached max iterations. Final gap: ", gap)
        end
    end

    info = (; aⱼ, lⱼ_target, choice..., converged, iterations = iter, gap)
    return returnInfo ? info : aⱼ
end

# For convenience. This function is used to check if there are any invalid elements (NaN, Missing, Inf) in the array and
# return the positions of these elements.
function CheckElements(array; position::Bool = false)
    # Check for NaN elements
    nan_mask = isnan.(array)
    nan_count = sum(nan_mask)
    nan_positions = findall(nan_mask)

    # Check for missing elements
    missing_mask = ismissing.(array)
    missing_count = sum(missing_mask)
    missing_positions = findall(missing_mask)

    # Check for Inf elements
    inf_mask = isinf.(array)
    inf_count = sum(inf_mask)
    inf_positions = findall(inf_mask)

    # Print results
    if nan_count > 0
        println("Found $nan_count NaN elements in the array")
        if position
            println("NaN positions: $nan_positions")
        end
    else
        println("No NaN elements found in the array")
    end

    if missing_count > 0
        println("Found $missing_count missing elements in the array")
        if position
            println("Missing positions: $missing_positions")
        end
    else
        println("No missing elements found in the array")
    end

    if inf_count > 0
        println("Found $inf_count Inf elements in the array")
        if position
            println("Inf positions: $inf_positions")
        end
    else
        println("No Inf elements found in the array")
    end

    # Return positions if requested, otherwise return boolean indicating if array is valid
    if position
        return (nan_positions, missing_positions, inf_positions)
    else
        return (nan_count == 0 && missing_count == 0 && inf_count == 0)
    end
end

# Generic function for convergence
function Converge(UpdateRule::Function, init::Array;
    tol=1e-6, maxIter=1e3,
    damp=0.4,
    power=false,
    normalize=nothing,
    displayGap=false,
    displaySummary=false,
    returnInfo=false)

    # Initialize the variables
    iter = 0
    XDiff = 1.0
    X = isnothing(normalize) ? init : normalize(init)
    newX = similar(X)

    # Start the loop
    while (iter < maxIter) && (XDiff > tol)
        iter += 1
        newX = UpdateRule(X)
        if !isnothing(normalize)
            newX = normalize(newX)
        end
        XDiff = maximum(abs.(newX .- X))
        X = power ? damp * X + (1 - damp) * X .* (newX ./ X) .^ 0.5 : damp * X .+ (1 - damp) * newX
        if !isnothing(normalize)
            X = normalize(X)
        end

        if displayGap
            println("Iteration: ", iter, ", Gap: ", XDiff)
        end
    end

    converged = XDiff <= tol

    if displaySummary
        if converged
            println("Successful convergence in ", iter, " iterations. Final gap: ", XDiff)
        else
            println("Maximum number of iterations reached. Final gap: ", XDiff)
        end
    end

    info = (; converged, iterations = iter, gap = XDiff)
    return returnInfo ? (X, info) : X
end

# Solve the equilibrium
function SolveModel(vars::NamedTuple, params::NamedTuple;
    damp=0.8, tol=1e-10, maxIter=1e3, power=false,
    displayGap=false, displaySummary=false,
    wⱼ_init=nothing,
    returnInfo=false)  # optional warm start

    # unpack the data
    (; l, d, zⱼ, aⱼ) = vars
    (; η, θ, α) = params

    # initial guess (use warm start if provided)
    J_local = size(d, 2)
    wⱼ = isnothing(wⱼ_init) ? ones(J_local) : copy(wⱼ_init)
    log_d = log.(d)

    # update rule
    function UpdateRule(wⱼ)
        choice = WorkerChoice(wⱼ, l, d, aⱼ, params; log_d)
        π_zj, ε_zj, lⱼ, εⱼ = choice.π_zj, choice.ε_zj, choice.lⱼ, choice.εⱼ

        wⱼ = α .* zⱼ .* lⱼ .^ (α - 1) .* εⱼ ./ (1 .+ εⱼ)
        wⱼ = wⱼ ./ sum(wⱼ .* lⱼ) # normalize total wage bill to 1
        return wⱼ, π_zj, ε_zj, lⱼ, εⱼ
    end

    # Convergence
    wⱼ, convergence = Converge(x -> UpdateRule(x)[1], wⱼ;
        tol=tol, maxIter=maxIter,
        damp=damp, power=power,
        normalize = x -> x ./ mean(x),
        displayGap=displayGap,
        displaySummary=displaySummary,
        returnInfo=true)

    solution = UpdateRule(wⱼ)
    if returnInfo
        return (;
            wⱼ = solution[1],
            π_zj = solution[2],
            ε_zj = solution[3],
            lⱼ = solution[4],
            εⱼ = solution[5],
            convergence...)
    end

    return solution
end

function SolveFirmPrimitivesFromData(vars::NamedTuple, params::NamedTuple;
    aⱼ_init=nothing, amenity_tol=1e-10, amenity_maxIter=5000, amenity_damp=0.5,
    displayGap=false, displaySummary=false, returnInfo=true)

    (; wⱼ, l, d, lⱼ) = vars
    (; α) = params

    amenity = SolveAmenitiesFromEmployment((; wⱼ, l, d, lⱼ), params;
        aⱼ_init, tol=amenity_tol, maxIter=amenity_maxIter, damp=amenity_damp,
        displayGap, displaySummary, returnInfo=true)

    lⱼ = amenity.lⱼ
    εⱼ = amenity.εⱼ
    zⱼ = wⱼ .* (1 .+ εⱼ) ./ (α .* lⱼ .^ (α - 1) .* εⱼ)
    zⱼ = zⱼ ./ mean(zⱼ) # normalize zⱼ to have mean 1

    solution = (;
        zⱼ,
        aⱼ = amenity.aⱼ,
        π_zj = amenity.π_zj,
        ε_zj = amenity.ε_zj,
        lⱼ,
        εⱼ,
        lⱼ_target = amenity.lⱼ_target,
        amenity_converged = amenity.converged,
        amenity_iterations = amenity.iterations,
        amenity_gap = amenity.gap
    )

    return returnInfo ? solution : (zⱼ, amenity.aⱼ)
end

# Solve zⱼ from observed wⱼ
function SolveZfromW(vars::NamedTuple, params::NamedTuple)

    # unpack the data
    (; wⱼ, l, d, aⱼ) = vars
    (; η, θ, α) = params

    choice = WorkerChoice(wⱼ, l, d, aⱼ, params; log_d = log.(d))
    lⱼ, εⱼ = choice.lⱼ, choice.εⱼ

    zⱼ = wⱼ .* (1 .+ εⱼ) ./ (α .* lⱼ .^ (α - 1) .* εⱼ)
    zⱼ = zⱼ ./ mean(zⱼ) # normalize zⱼ to have mean 1
    return zⱼ
end

function _twoperiod_fe_dataframe(y_base, y_counterfactual, bigMA, w_diff, cluster)
    n_firms = length(y_base)
    firm = repeat(collect(1:n_firms), inner=2)
    period = repeat([0, 1], n_firms)
    y = vec(reshape(vcat(y_base', y_counterfactual'), 2 * n_firms))
    w_diff_long = repeat(w_diff, inner=2)
    big_long = repeat(bigMA, inner=2)
    post = Float64.(period)

    df = DataFrame(
        y = y,
        firm = firm,
        period = period,
        w_diff = w_diff_long,
        big_post = big_long .* post,
        big_wdiff_post = big_long .* w_diff_long .* post,
    )
    if !isnothing(cluster)
        any(ismissing, cluster) && error("cluster cannot contain missing values in the kept sample")
        df[!, :cluster] = repeat(string.(cluster), inner=2)
    end
    return df
end

function _coef_positions(model, names)
    model_names = string.(coefnames(model))
    positions = [findfirst(==(name), model_names) for name in names]
    any(isnothing, positions) && error("Missing coefficient(s): " *
        join(names[isnothing.(positions)], ", "))
    return Int.(positions)
end

function _twoperiod_fixed_effect_model(y_base, y_counterfactual, bigMA, w_diff;
    cluster=nothing, include_wdiff_interaction::Bool=true)

    df = _twoperiod_fe_dataframe(y_base, y_counterfactual, bigMA, w_diff, cluster)
    rhs = include_wdiff_interaction ?
        term(:big_post) + term(:big_wdiff_post) + fe(:firm) + fe(:period) +
            fe(:period)&term(:w_diff) :
        term(:big_post) + fe(:firm) + fe(:period) + fe(:period)&term(:w_diff)
    formula = term(:y) ~ rhs
    model = isnothing(cluster) ? reg(df, formula) : reg(df, formula, Vcov.cluster(:cluster))
    wanted = include_wdiff_interaction ? ["big_post", "big_wdiff_post"] : ["big_post"]
    idx = _coef_positions(model, wanted)
    beta = coef(model)[idx]
    vc = vcov(model)[idx, idx]
    se = sqrt.(diag(vc))
    t = beta ./ se
    n_clusters = isnothing(cluster) ? missing : length(unique(string.(cluster)))
    return (; model, beta, se, t, vcov=vc, n=nobs(model), n_clusters,
        df=dof_residual(model))
end

function EstimateTwoPeriodDIDMoments(lnl_base, lnl_counterfactual,
    lnw_base, lnw_counterfactual, bigMA, w;
    wage_center=:treated, keep=nothing, cluster=nothing, return_stats::Bool=false)

    J = length(w)
    length(lnl_base) == J || error("lnl_base must have length $J")
    length(lnl_counterfactual) == J || error("lnl_counterfactual must have length $J")
    length(lnw_base) == J || error("lnw_base must have length $J")
    length(lnw_counterfactual) == J || error("lnw_counterfactual must have length $J")
    length(bigMA) == J || error("bigMA must have length $J")

    keep_mask = isnothing(keep) ? trues(J) : Bool.(keep)
    length(keep_mask) == J || error("keep must have length $J")
    if !isnothing(cluster)
        length(cluster) == J || error("cluster must have length $J")
    end

    lnl_base_keep = Float64.(vec(lnl_base)[keep_mask])
    lnl_counterfactual_keep = Float64.(vec(lnl_counterfactual)[keep_mask])
    lnw_base_keep = Float64.(vec(lnw_base)[keep_mask])
    lnw_counterfactual_keep = Float64.(vec(lnw_counterfactual)[keep_mask])
    bigMA_keep = Float64.(vec(bigMA)[keep_mask])
    w_keep = Float64.(vec(w)[keep_mask])
    cluster_keep = isnothing(cluster) ? nothing : collect(cluster)[keep_mask]

    invalid = return_stats ?
        (; beta = [Inf, Inf, Inf], se = [Inf, Inf, Inf], t = [Inf, Inf, Inf],
            vcov_labor = fill(Inf, 2, 2), vcov_wage = fill(Inf, 1, 1),
            n = length(lnl_base_keep), n_clusters = fill(missing, 3), df = fill(missing, 3),
            model_labor = nothing, model_wage = nothing) :
        [Inf, Inf, Inf]

    if wage_center == :treated
        treated = bigMA_keep .== 1
        w_center = any(treated) ? mean(w_keep[treated]) : mean(w_keep)
    elseif wage_center == :all
        w_center = mean(w_keep)
    else
        error("wage_center must be :treated or :all")
    end
    w_diff = w_keep .- w_center

    try
        labor = _twoperiod_fixed_effect_model(lnl_base_keep, lnl_counterfactual_keep,
            bigMA_keep, w_diff; cluster=cluster_keep, include_wdiff_interaction=true)
        wage = _twoperiod_fixed_effect_model(lnw_base_keep, lnw_counterfactual_keep,
            bigMA_keep, w_diff; cluster=cluster_keep, include_wdiff_interaction=false)

        beta = [labor.beta[1], labor.beta[2], wage.beta[1]]
        if !return_stats
            return beta
        end

        se = [labor.se[1], labor.se[2], wage.se[1]]
        t = [labor.t[1], labor.t[2], wage.t[1]]
        n_clusters = isnothing(cluster_keep) ? fill(missing, 3) :
            [labor.n_clusters, labor.n_clusters, wage.n_clusters]
        df = [labor.df, labor.df, wage.df]

        return (; beta, se, t, vcov_labor=labor.vcov, vcov_wage=wage.vcov,
            n=labor.n, n_clusters, df, model_labor=labor.model,
            model_wage=wage.model)
    catch err
        if err isa InterruptException
            rethrow(err)
        end
        return invalid
    end
end

function MomentFirmMask(J::Integer; moment_firm_mask=nothing, firm_ids=nothing,
    reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false)

    if !isnothing(moment_firm_mask)
        length(moment_firm_mask) == J || error("moment_firm_mask must have length $J")
        any(ismissing, moment_firm_mask) && error("moment_firm_mask cannot contain missing values")
        keep = Bool.(moment_firm_mask)
    elseif restrict_to_reg_sample
        isnothing(firm_ids) && error("firm_ids is required when restrict_to_reg_sample=true")
        isnothing(reg_sample_ids) && error("reg_sample_ids is required when restrict_to_reg_sample=true")
        length(firm_ids) == J || error("firm_ids must have length $J")

        reg_ids = Set(collect(skipmissing(reg_sample_ids)))
        keep = [!ismissing(id) && id in reg_ids for id in firm_ids]
    else
        keep = trues(J)
    end

    any(keep) || error("Moment sample is empty")
    return keep
end

function _invalid_moments()
    return [Inf, Inf, Inf]
end

function _unpack_calibration_params(params_to_estimate, α_fixed)
    params = Float64.(vec(params_to_estimate))
    if length(params) == 2
        isnothing(α_fixed) && error("α must be supplied when params_to_estimate has length 2")
        η, θ = params
        α = Float64(α_fixed)
    elseif length(params) == 3
        η, θ, α = params
    else
        error("params_to_estimate must have length 2 ([η, θ]) or length 3 ([η, θ, α])")
    end
    return η, θ, α
end

# Compute model moments given parameters
function ComputeModelMoments(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, α=nothing,
    aⱼ_init=nothing, amenity_tol=1e-10, amenity_maxIter=5000, amenity_damp=0.5,
    inner_tol=1e-5, inner_maxIter=3000, inner_display=false, require_convergence=true,
    continuation_steps=5, employment_change=:log, wage_center=:treated,
    moment_firm_mask=nothing, firm_ids=nothing, reg_sample_ids=nothing,
    restrict_to_reg_sample::Bool=false, displayGap = false, damp_cf = nothing)
    η, θ, α = _unpack_calibration_params(params_to_estimate, α)
    continuation_steps = max(1, Int(continuation_steps))
    
    # Ensure parameters are in valid range
    if η <= 0 || θ <= 0 || α <= 0 || α >= 1 ||
        !isfinite(η) || !isfinite(θ) || !isfinite(α)
        return _invalid_moments()
    end
    
    # Higher θ can make the counterfactual fixed point sharper, so use more damping.
    damp_cf = isnothing(damp_cf) ? clamp(0.65 + 0.085 * θ, 0.75, 0.97) : damp_cf
    
    params = (; η, θ, α)
    
    # Step 2a: Invert for firm amenities and productivity using observed
    # employment and wage in the baseline data.
    primitives = SolveFirmPrimitivesFromData((; wⱼ = wⱼ_data, lⱼ = lⱼ_data, l, d), params;
        aⱼ_init, amenity_tol, amenity_maxIter, amenity_damp,
        displaySummary=inner_display, displayGap = displayGap, returnInfo=true)
    zⱼ, aⱼ = primitives.zⱼ, primitives.aⱼ
    if require_convergence && !primitives.amenity_converged
        return _invalid_moments()
    end
    
    # Check for invalid zⱼ
    if any(isnan.(zⱼ)) || any(isinf.(zⱼ)) || any(isnan.(aⱼ)) || any(isinf.(aⱼ))
        return _invalid_moments()
    end
    
    # Step 2b: Use the inverted baseline directly. By construction, the
    # inverted amenities match lⱼ_data and the normalized productivity makes
    # wⱼ_data satisfy the baseline wage FOC, so a separate baseline fixed point
    # solve is redundant on the calibration path.
    wⱼ = vec(Float64.(wⱼ_data))
    lⱼ = vec(Float64.(primitives.lⱼ))
    
    # Check for invalid wages
    if any(isnan.(wⱼ)) || any(isinf.(wⱼ))
        return _invalid_moments()
    end
    
    # Step 2c: Solve counterfactual model along a commuting-time path.
    # This continuation step is much more stable than jumping from d to d′.
    cf = nothing
    wⱼ_cf_init = wⱼ
    log_d = log.(d)
    log_d′ = log.(d′)
    for step in 1:continuation_steps
        path_share = step / continuation_steps
        d_path = exp.((1 - path_share) .* log_d .+ path_share .* log_d′)
        vars_cf = (; l, d = d_path, zⱼ, aⱼ)
        if inner_display && continuation_steps > 1
            println("Counterfactual continuation step ", step, "/", continuation_steps)
        end
        cf = SolveModel(vars_cf, params;
            damp=damp_cf, tol=inner_tol, power=false, maxIter=inner_maxIter,
            wⱼ_init=wⱼ_cf_init, displaySummary=inner_display, displayGap=displayGap, returnInfo=true)

        if require_convergence && !cf.converged
            return _invalid_moments()
        end

        wⱼ_cf_init = cf.wⱼ
    end
    wⱼ′, lⱼ′ = cf.wⱼ, cf.lⱼ

    # Check for invalid counterfactual wages
    if any(isnan.(wⱼ′)) || any(isinf.(wⱼ′))
        return _invalid_moments()
    end
    
    # Step 2d: Compute regression moments from appended two-period firm data.
    # employment_change is kept in the signature for caller compatibility; this
    # DID specification uses log employment levels, matching the lnl outcome.
    dMA = vec(l' * d .- l' * d′) |> x -> replace(x, -Inf => -8)
    bigMA = Float64.(dMA .>= 0.5)
    
    w = vec(log.(max.(wⱼ, eps(Float64))))
    keep = MomentFirmMask(length(w); moment_firm_mask, firm_ids, reg_sample_ids,
        restrict_to_reg_sample)
    lnl = vec(log.(max.(lⱼ, eps(Float64))))
    lnl′ = vec(log.(max.(lⱼ′, eps(Float64))))
    lnw = vec(log.(max.(wⱼ, eps(Float64))))
    lnw′ = vec(log.(max.(wⱼ′, eps(Float64))))

    return EstimateTwoPeriodDIDMoments(lnl, lnl′, lnw, lnw′, bigMA, w;
        wage_center, keep)
end

function ComputeModelMomentsFixedPrimitives(params_to_estimate; l, d, d′, wⱼ_data,
    lⱼ_data, fixed_primitives, α=nothing, inner_tol=1e-5, inner_maxIter=3000,
    inner_display=false, require_convergence=true, continuation_steps=5,
    employment_change=:log, wage_center=:treated, moment_firm_mask=nothing,
    firm_ids=nothing, reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false,
    displayGap=false, damp_base=nothing, damp_cf=nothing)

    η, θ, α = _unpack_calibration_params(params_to_estimate, α)
    continuation_steps = max(1, Int(continuation_steps))

    if η <= 0 || θ <= 0 || α <= 0 || α >= 1 ||
        !isfinite(η) || !isfinite(θ) || !isfinite(α)
        return _invalid_moments()
    end

    J = size(d, 2)
    length(wⱼ_data) == J || error("wⱼ_data must have length $J")
    length(lⱼ_data) == J || error("lⱼ_data must have length $J")
    hasproperty(fixed_primitives, :zⱼ) ||
        error("fixed_primitives must contain zⱼ")
    hasproperty(fixed_primitives, :aⱼ) ||
        error("fixed_primitives must contain aⱼ")

    zⱼ = vec(Float64.(fixed_primitives.zⱼ))
    aⱼ = vec(Float64.(fixed_primitives.aⱼ))
    length(zⱼ) == J || error("fixed_primitives.zⱼ must have length $J")
    length(aⱼ) == J || error("fixed_primitives.aⱼ must have length $J")
    if any(isnan.(zⱼ)) || any(isinf.(zⱼ)) || any(isnan.(aⱼ)) || any(isinf.(aⱼ))
        return _invalid_moments()
    end

    params = (; η, θ, α)
    damp_default = clamp(0.65 + 0.085 * θ, 0.75, 0.97)
    damp_base = isnothing(damp_base) ? damp_default : damp_base
    damp_cf = isnothing(damp_cf) ? damp_default : damp_cf

    base = SolveModel((; l, d, zⱼ, aⱼ), params;
        damp=damp_base, tol=inner_tol, power=false, maxIter=inner_maxIter,
        wⱼ_init=wⱼ_data, displaySummary=inner_display, displayGap=displayGap,
        returnInfo=true)
    if require_convergence && !base.converged
        return _invalid_moments()
    end
    if any(isnan.(base.wⱼ)) || any(isinf.(base.wⱼ))
        return _invalid_moments()
    end

    cf = nothing
    wⱼ_cf_init = base.wⱼ
    log_d = log.(d)
    log_d′ = log.(d′)
    for step in 1:continuation_steps
        path_share = step / continuation_steps
        d_path = exp.((1 - path_share) .* log_d .+ path_share .* log_d′)
        vars_cf = (; l, d = d_path, zⱼ, aⱼ)
        if inner_display && continuation_steps > 1
            println("Counterfactual continuation step ", step, "/", continuation_steps)
        end
        cf = SolveModel(vars_cf, params;
            damp=damp_cf, tol=inner_tol, power=false, maxIter=inner_maxIter,
            wⱼ_init=wⱼ_cf_init, displaySummary=inner_display,
            displayGap=displayGap, returnInfo=true)
        if require_convergence && !cf.converged
            return _invalid_moments()
        end
        wⱼ_cf_init = cf.wⱼ
    end
    if any(isnan.(cf.wⱼ)) || any(isinf.(cf.wⱼ))
        return _invalid_moments()
    end

    # employment_change is kept in the signature for caller compatibility; this
    # DID specification uses log employment levels, matching the lnl outcome.
    dMA = vec(l' * d .- l' * d′) |> x -> replace(x, -Inf => -8)
    bigMA = Float64.(dMA .>= 0.5)

    w = vec(log.(max.(wⱼ_data, eps(Float64))))
    keep = MomentFirmMask(length(w); moment_firm_mask, firm_ids, reg_sample_ids,
        restrict_to_reg_sample)
    lnl = vec(log.(max.(base.lⱼ, eps(Float64))))
    lnl′ = vec(log.(max.(cf.lⱼ, eps(Float64))))
    lnw = vec(log.(max.(base.wⱼ, eps(Float64))))
    lnw′ = vec(log.(max.(cf.wⱼ, eps(Float64))))

    return EstimateTwoPeriodDIDMoments(lnl, lnl′, lnw, lnw′, bigMA, w;
        wage_center, keep)
end

# Objective function for estimation
function ObjectiveFunction(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, β_target, α=nothing,
    aⱼ_init=nothing, amenity_tol=1e-10, amenity_maxIter=5000, amenity_damp=0.5,
    verbose=false, inner_tol=1e-5, inner_maxIter=3000, inner_display=false,
    require_convergence=true, continuation_steps=5, employment_change=:log,
    wage_center=:treated, max_abs_moment=10.0, moment_firm_mask=nothing,
    firm_ids=nothing, reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false)
    length(β_target) == 3 || error("β_target must have length 3")
    η, θ, α_est = _unpack_calibration_params(params_to_estimate, α)
    
    # Return large penalty for invalid parameters
    if η <= 0 || θ <= 0 || α_est <= 0 || α_est >= 1 ||
        !isfinite(η) || !isfinite(θ) || !isfinite(α_est)
        return 1e10
    end
    
    β_model = ComputeModelMoments(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, α,
        aⱼ_init=aⱼ_init, amenity_tol=amenity_tol, amenity_maxIter=amenity_maxIter,
        amenity_damp=amenity_damp, inner_tol=inner_tol, inner_maxIter=inner_maxIter,
        inner_display=inner_display,
        require_convergence=require_convergence, continuation_steps=continuation_steps,
        employment_change=employment_change, wage_center=wage_center,
        moment_firm_mask=moment_firm_mask, firm_ids=firm_ids,
        reg_sample_ids=reg_sample_ids, restrict_to_reg_sample=restrict_to_reg_sample)
    
    # Check for invalid model output
    if any(isnan.(β_model)) || any(isinf.(β_model)) || any(abs.(β_model) .> max_abs_moment)
        return 1e10
    end
    
    # Sum of squared differences (returns scalar, not matrix!)
    diff = β_model - β_target
    obj = sum(diff .^ 2)
    
    if verbose
        println("η=$η, θ=$θ, α=$α_est → β_model=$β_model, obj=$obj")
    end
    
    return obj
end

function EvaluateCalibrationGrid(; l, d, d′, wⱼ_data, lⱼ_data, β_target, α=nothing,
    η_grid, θ_grid, α_grid=nothing, inner_tol=2e-5, inner_maxIter=3000,
    aⱼ_init=nothing, amenity_tol=1e-10, amenity_maxIter=5000, amenity_damp=0.5,
    continuation_steps=5, employment_change=:log, wage_center=:treated,
    max_abs_moment=10.0, verbose=true, moment_firm_mask=nothing,
    firm_ids=nothing, reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false,
    fixed_primitives=nothing)

    length(β_target) == 3 || error("β_target must have length 3")
    if isnothing(α_grid)
        isnothing(α) && error("Either α_grid or fixed α must be supplied")
        α_values = [Float64(α)]
    else
        α_values = Float64.(α_grid)
    end

    grid_points = [(Float64(η), Float64(θ), Float64(α_value))
        for η in η_grid for θ in θ_grid for α_value in α_values]
    total = length(grid_points)
    η_col = Vector{Float64}(undef, total)
    θ_col = Vector{Float64}(undef, total)
    α_col = Vector{Float64}(undef, total)
    β_labor_bigMA_col = Vector{Float64}(undef, total)
    β_labor_bigMA_wdiff_col = Vector{Float64}(undef, total)
    β_wage_bigMA_col = Vector{Float64}(undef, total)
    objective_col = Vector{Float64}(undef, total)
    progress = Ref(0)
    progress_lock = ReentrantLock()

    if verbose
        println("Evaluating ", total, " calibration grid points with ",
            Threads.nthreads(), " Julia thread(s).")
    end

    Threads.@threads :dynamic for idx in eachindex(grid_points)
        η, θ, α_value = grid_points[idx]
        if isnothing(fixed_primitives)
            β_model = ComputeModelMoments([η, θ, α_value]; l, d, d′, wⱼ_data, lⱼ_data,
                inner_tol=inner_tol, inner_maxIter=inner_maxIter,
                inner_display=false, require_convergence=true,
                aⱼ_init=aⱼ_init, amenity_tol=amenity_tol,
                amenity_maxIter=amenity_maxIter, amenity_damp=amenity_damp,
                continuation_steps=continuation_steps,
                employment_change=employment_change, wage_center=wage_center,
                moment_firm_mask=moment_firm_mask, firm_ids=firm_ids,
                reg_sample_ids=reg_sample_ids,
                restrict_to_reg_sample=restrict_to_reg_sample)
        else
            β_model = ComputeModelMomentsFixedPrimitives([η, θ, α_value];
                l, d, d′, wⱼ_data, lⱼ_data, fixed_primitives,
                inner_tol=inner_tol, inner_maxIter=inner_maxIter,
                inner_display=false, require_convergence=true,
                continuation_steps=continuation_steps,
                employment_change=employment_change, wage_center=wage_center,
                moment_firm_mask=moment_firm_mask, firm_ids=firm_ids,
                reg_sample_ids=reg_sample_ids,
                restrict_to_reg_sample=restrict_to_reg_sample)
        end

        if any(isnan.(β_model)) || any(isinf.(β_model)) || any(abs.(β_model) .> max_abs_moment)
            obj = 1e10
        else
            obj = sum((β_model .- β_target) .^ 2)
        end

        η_col[idx] = η
        θ_col[idx] = θ
        α_col[idx] = α_value
        β_labor_bigMA_col[idx] = β_model[1]
        β_labor_bigMA_wdiff_col[idx] = β_model[2]
        β_wage_bigMA_col[idx] = β_model[3]
        objective_col[idx] = obj

        if verbose
            lock(progress_lock)
            try
                progress[] += 1
                println("Grid ", progress[], "/", total,
                    " (thread ", Threads.threadid(), "): η=", η, ", θ=", θ,
                    ", α=", α_value, " → β_model=", β_model, ", obj=", obj)
            finally
                unlock(progress_lock)
            end
        end
    end

    results = DataFrame(
        η = η_col,
        θ = θ_col,
        α = α_col,
        β_labor_bigMA = β_labor_bigMA_col,
        β_labor_bigMA_wdiff = β_labor_bigMA_wdiff_col,
        β_wage_bigMA = β_wage_bigMA_col,
        objective = objective_col
    )
    sort!(results, :objective)
    return results
end

function _box_to_unconstrained(x, lower, upper)
    x_inner = clamp.(Float64.(x), lower .+ eps.(lower), upper .- eps.(upper))
    return log.((x_inner .- lower) ./ (upper .- x_inner))
end

function _unconstrained_to_box(y, lower, upper)
    return lower .+ (upper .- lower) ./ (1 .+ exp.(-Float64.(y)))
end

function CalibrateEtaThetaAlpha(; l, d, d′, wⱼ_data, lⱼ_data, β_target, α=nothing,
    x0=[0.5, 50.0, 0.4], starts=nothing, lower=[0.1, 1.0, 0.1],
    upper=[1.0, 100.0, 0.9], iterations=250, x_abstol=1e-4, f_reltol=1e-8,
    aⱼ_init=nothing, amenity_tol=1e-10, amenity_maxIter=5000, amenity_damp=0.5,
    inner_tol=1e-5, inner_maxIter=3000, continuation_steps=5,
    employment_change=:log, wage_center=:treated, max_abs_moment=10.0,
    verbose=true, show_trace=true, moment_firm_mask=nothing, firm_ids=nothing,
    reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false)

    length(β_target) == 3 || error("β_target must have length 3")
    estimate_alpha = isnothing(α)
    α_fixed = estimate_alpha ? nothing : Float64(α)
    if !estimate_alpha && (α_fixed <= 0 || α_fixed >= 1 || !isfinite(α_fixed))
        error("fixed α must be finite and lie between 0 and 1")
    end

    n_estimated = estimate_alpha ? 3 : 2
    lower = Float64.(lower)
    upper = Float64.(upper)
    if !estimate_alpha && length(lower) == 3
        lower = lower[1:2]
    end
    if !estimate_alpha && length(upper) == 3
        upper = upper[1:2]
    end
    length(lower) == n_estimated ||
        error("lower must have length $n_estimated")
    length(upper) == n_estimated ||
        error("upper must have length $n_estimated")

    function estimated_start_vector(start, label)
        start_vec = Float64.(vec(start))
        if length(start_vec) == n_estimated
            return start_vec
        elseif !estimate_alpha && length(start_vec) == 3
            return start_vec[1:2]
        else
            error("$label must have length $n_estimated" *
                  (estimate_alpha ? "" : " ([η, θ]) or length 3 ([η, θ, α])"))
        end
    end

    full_parameters(x) = estimate_alpha ? x : [x[1], x[2], α_fixed]

    function objective_y(y)
        x = _unconstrained_to_box(y, lower, upper)
        return ObjectiveFunction(x; l, d, d′, wⱼ_data, lⱼ_data, β_target, α=α_fixed,
            aⱼ_init=aⱼ_init, amenity_tol=amenity_tol, amenity_maxIter=amenity_maxIter,
            amenity_damp=amenity_damp, verbose=verbose, inner_tol=inner_tol,
            inner_maxIter=inner_maxIter, inner_display=false, require_convergence=true,
            continuation_steps=continuation_steps, employment_change=employment_change,
            wage_center=wage_center, max_abs_moment=max_abs_moment,
            moment_firm_mask=moment_firm_mask, firm_ids=firm_ids,
            reg_sample_ids=reg_sample_ids,
            restrict_to_reg_sample=restrict_to_reg_sample)
    end

    default_start = isnothing(x0) ? (lower .+ upper) ./ 2 : estimated_start_vector(x0, "x0")
    start_list = isnothing(starts) ? [default_start] : starts
    isempty(start_list) && error("CalibrateEtaThetaAlpha requires at least one starting value")
    runs = NamedTuple[]
    best_run = nothing

    for (start_id, start) in enumerate(start_list)
        start_vec = estimated_start_vector(start, "Each calibration start")
        x_start = clamp.(start_vec, lower .+ 1e-8, upper .- 1e-8)
        y0 = _box_to_unconstrained(x_start, lower, upper)

        if verbose
            println("\nLocal calibration start ", start_id, "/", length(start_list),
                ": x0=", x_start)
        end

        result = optimize(
            objective_y,
            y0,
            NelderMead(),
            Optim.Options(
                show_trace = show_trace,
                iterations = iterations,
                x_abstol = x_abstol,
                f_reltol = f_reltol
            )
        )

        estimated_params = _unconstrained_to_box(Optim.minimizer(result), lower, upper)
        params = full_parameters(estimated_params)
        run = (;
            parameters = params,
            estimated_parameters = estimated_params,
            objective = Optim.minimum(result),
            converged = Optim.converged(result),
            result,
            start = x_start
        )
        push!(runs, run)

        if isnothing(best_run) || run.objective < best_run.objective
            best_run = run
        end
    end

    return (;
        parameters = best_run.parameters,
        objective = best_run.objective,
        converged = best_run.converged,
        result = best_run.result,
        runs,
        estimated_parameters = best_run.estimated_parameters,
        fixed_alpha = α_fixed,
        lower,
        upper
    )
end

function CalibrateEtaTheta(; l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
    x0=[2.75, 2.0], starts=nothing, lower=[0.25, 0.25], upper=[8.0, 8.0],
    iterations=250, x_abstol=1e-4, f_reltol=1e-8,
    aⱼ_init=nothing, amenity_tol=1e-10, amenity_maxIter=5000, amenity_damp=0.5,
    inner_tol=1e-5, inner_maxIter=3000, continuation_steps=5,
    employment_change=:log, wage_center=:treated, max_abs_moment=10.0,
    verbose=true, show_trace=true, moment_firm_mask=nothing, firm_ids=nothing,
    reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false)

    lower = Float64.(lower)
    upper = Float64.(upper)

    function objective_y(y)
        x = _unconstrained_to_box(y, lower, upper)
        return ObjectiveFunction(x; l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
            aⱼ_init=aⱼ_init, amenity_tol=amenity_tol, amenity_maxIter=amenity_maxIter,
            amenity_damp=amenity_damp, verbose=verbose, inner_tol=inner_tol, inner_maxIter=inner_maxIter,
            inner_display=false, require_convergence=true,
            continuation_steps=continuation_steps,
            employment_change=employment_change, wage_center=wage_center,
            max_abs_moment=max_abs_moment, moment_firm_mask=moment_firm_mask,
            firm_ids=firm_ids, reg_sample_ids=reg_sample_ids,
            restrict_to_reg_sample=restrict_to_reg_sample)
    end

    default_start = isnothing(x0) ? (lower .+ upper) ./ 2 : x0
    start_list = isnothing(starts) ? [default_start] : starts
    isempty(start_list) && error("CalibrateEtaTheta requires at least one starting value")
    runs = NamedTuple[]
    best_run = nothing

    for (start_id, start) in enumerate(start_list)
        x_start = clamp.(Float64.(start), lower .+ 1e-8, upper .- 1e-8)
        y0 = _box_to_unconstrained(x_start, lower, upper)

        if verbose
            println("\nLocal calibration start ", start_id, "/", length(start_list),
                ": x0=", x_start)
        end

        result = optimize(
            objective_y,
            y0,
            NelderMead(),
            Optim.Options(
                show_trace = show_trace,
                iterations = iterations,
                x_abstol = x_abstol,
                f_reltol = f_reltol
            )
        )

        params = _unconstrained_to_box(Optim.minimizer(result), lower, upper)
        run = (;
            parameters = params,
            objective = Optim.minimum(result),
            converged = Optim.converged(result),
            result,
            start = x_start
        )
        push!(runs, run)

        if isnothing(best_run) || run.objective < best_run.objective
            best_run = run
        end
    end

    return (;
        parameters = best_run.parameters,
        objective = best_run.objective,
        converged = best_run.converged,
        result = best_run.result,
        runs,
        lower,
        upper
    )
end
