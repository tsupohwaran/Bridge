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
    (; l, d, zⱼ) = vars
    (; η, θ, α) = params

    # initial guess (use warm start if provided)
    J_local = size(d, 2)
    wⱼ = isnothing(wⱼ_init) ? ones(J_local) : copy(wⱼ_init)

    # update rule
    function UpdateRule(wⱼ)
        # Use log-sum-exp trick for numerical stability
        log_xzj = θ .* (log.(wⱼ') .- η .* log.(d))  # Z × J matrix
        log_xzj_max = maximum(log_xzj, dims=2)  # Z × 1
        log_sum_exp = log_xzj_max .+ log.(sum(exp.(log_xzj .- log_xzj_max), dims=2))
        π_zj = exp.(log_xzj .- log_sum_exp)
        
        ε_zj = θ .* (1 .- π_zj)
        lⱼ = sum(π_zj .* l, dims = 1)'
        εⱼ = sumsqueeze(π_zj .* l ./ lⱼ' .* ε_zj, dims = 1)

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

# Solve zⱼ from observed wⱼ
function SolveZfromW(vars::NamedTuple, params::NamedTuple)

    # unpack the data
    (; wⱼ, l, d) = vars
    (; η, θ, α) = params

    # Use log-sum-exp trick for numerical stability
    log_xzj = θ .* (log.(wⱼ') .- η .* log.(d))  # Z × J matrix
    log_xzj_max = maximum(log_xzj, dims=2)  # Z × 1
    log_sum_exp = log_xzj_max .+ log.(sum(exp.(log_xzj .- log_xzj_max), dims=2))
    π_zj = exp.(log_xzj .- log_sum_exp)
    
    ε_zj = θ .* (1 .- π_zj)
    lⱼ = sum(π_zj .* l, dims = 1)'
    εⱼ = sumsqueeze(π_zj .* l ./ lⱼ' .* ε_zj, dims = 1)

    zⱼ = wⱼ .* (1 .+ εⱼ) ./ (α .* lⱼ .^ (α - 1) .* εⱼ)
    zⱼ = zⱼ ./ mean(zⱼ) # normalize zⱼ to have mean 1
    return zⱼ
end

function EstimateTwoPeriodDIDMoments(regDF::DataFrame)
    # Equivalent to:
    # reghdfe lnl i1.bigMA#i1.post c.w_diff#i1.bigMA#i1.post,
    #     absorb(id year c.w_diff#year)
    # for a balanced two-period panel. The code below applies the
    # Frisch-Waugh-Lovell residualization implied by these absorbed effects.
    baseline = sort(regDF[regDF.post .== 0, :], :id)
    counterfactual = sort(regDF[regDF.post .== 1, :], :id)

    if nrow(baseline) != nrow(counterfactual) || !all(baseline.id .== counterfactual.id)
        return [Inf, Inf]
    end

    # Firm FE residualization in a two-period panel is equivalent to taking
    # firm-level changes. This is a transformation of the lnl regression, not a
    # different dependent-variable specification.
    Δlnl = Float64.(counterfactual.lnl .- baseline.lnl)
    bigMA = Float64.(baseline.bigMA)
    w_diff = Float64.(baseline.w_diff)
    X_absorb = hcat(ones(length(Δlnl)), w_diff)
    X_target = hcat(bigMA, w_diff .* bigMA)

    if rank(X_absorb) < size(X_absorb, 2)
        return [Inf, Inf]
    end

    # FWL residualization
    y_resid = Δlnl - X_absorb * (X_absorb \ Δlnl)
    X_resid = X_target - X_absorb * (X_absorb \ X_target)

    if rank(X_resid) < size(X_resid, 2)
        return [Inf, Inf]
    end

    β = X_resid \ y_resid
    return [β[1], β[2]]
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

# Compute model moments given parameters
function ComputeModelMoments(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, α,
    inner_tol=1e-5, inner_maxIter=3000, inner_display=false, require_convergence=true,
    continuation_steps=5, employment_change=:log, wage_center=:treated,
    moment_firm_mask=nothing, firm_ids=nothing, reg_sample_ids=nothing,
    restrict_to_reg_sample::Bool=false)
    η, θ = params_to_estimate
    continuation_steps = max(1, Int(continuation_steps))
    
    # Ensure parameters are in valid range
    if η <= 0 || θ <= 0
        return [Inf, Inf]
    end
    
    # Higher θ can make the fixed point sharper. The baseline is anchored by
    # observed wages, while the bridge counterfactual needs more damping.
    damp_base = clamp(0.45 + 0.04 * θ, 0.55, 0.85)
    damp_cf = clamp(0.65 + 0.085 * θ, 0.75, 0.97)
    
    params = (; η, θ, α)
    
    # Step 2a: Invert for zⱼ given current (η, θ)
    vars = (; wⱼ = wⱼ_data, l, d)
    zⱼ = SolveZfromW(vars, params)
    
    # Check for invalid zⱼ
    if any(isnan.(zⱼ)) || any(isinf.(zⱼ))
        return [Inf, Inf]
    end
    
    # Step 2b: Solve baseline model with incremental continuation for high θ
    vars_base = (; l, d, zⱼ)
    # zⱼ is backed out from wⱼ_data, so use wⱼ_data as the anchor for the baseline solve.
    base = SolveModel(vars_base, params;
        damp=damp_base, tol=inner_tol, power=false, maxIter=inner_maxIter,
        wⱼ_init=wⱼ_data, displaySummary=inner_display, returnInfo=true)
    wⱼ, π_zj, ε_zj, lⱼ, εⱼ = base.wⱼ, base.π_zj, base.ε_zj, base.lⱼ, base.εⱼ

    if require_convergence && !base.converged
        return [Inf, Inf]
    end
    
    # Check for invalid wages
    if any(isnan.(wⱼ)) || any(isinf.(wⱼ))
        return [Inf, Inf]
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
        vars_cf = (; l, d = d_path, zⱼ)
        if inner_display && continuation_steps > 1
            println("Counterfactual continuation step ", step, "/", continuation_steps)
        end
        cf = SolveModel(vars_cf, params;
            damp=damp_cf, tol=inner_tol, power=false, maxIter=inner_maxIter,
            wⱼ_init=wⱼ_cf_init, displaySummary=inner_display, returnInfo=true)

        if require_convergence && !cf.converged
            return [Inf, Inf]
        end

        wⱼ_cf_init = cf.wⱼ
    end
    wⱼ′, π_zj′, ε_zj′, lⱼ′, εⱼ′ = cf.wⱼ, cf.π_zj, cf.ε_zj, cf.lⱼ, cf.εⱼ

    # Check for invalid counterfactual wages
    if any(isnan.(wⱼ′)) || any(isinf.(wⱼ′))
        return [Inf, Inf]
    end
    
    # Step 2d: Compute regression moments from appended two-period firm data.
    # employment_change is kept in the signature for caller compatibility; this
    # DID specification uses log employment levels, matching the lnl outcome.
    dMA = sum((d - d′) .* l, dims = 1)' |> x -> replace(x, -Inf => -8)
    bigMA = Float64.(dMA .>= 0.5)
    
    w = vec(log.(max.(wⱼ, eps(Float64))))
    keep = MomentFirmMask(length(w); moment_firm_mask, firm_ids, reg_sample_ids,
        restrict_to_reg_sample)

    regDF_firm = DataFrame(
        id = collect(1:length(w))[keep],
        bigMA = vec(bigMA)[keep],
        w = w[keep]
    )
    if wage_center == :treated
        treated = regDF_firm.bigMA .== 1
        w_center = any(treated) ? mean(regDF_firm[treated, :w]) : mean(regDF_firm[!, :w])
    elseif wage_center == :all
        w_center = mean(regDF_firm[!, :w])
    else
        error("wage_center must be :treated or :all")
    end
    regDF_firm.w_diff = regDF_firm.w .- w_center

    J = nrow(regDF_firm)
    regDF = vcat(
        DataFrame(
            id = regDF_firm.id,
            year = zeros(Int, J),
            post = zeros(Int, J),
            bigMA = regDF_firm.bigMA,
            w_diff = regDF_firm.w_diff,
            lnl = vec(log.(max.(lⱼ, eps(Float64))))[keep]
        ),
        DataFrame(
            id = regDF_firm.id,
            year = ones(Int, J),
            post = ones(Int, J),
            bigMA = regDF_firm.bigMA,
            w_diff = regDF_firm.w_diff,
            lnl = vec(log.(max.(lⱼ′, eps(Float64))))[keep]
        )
    )

    return EstimateTwoPeriodDIDMoments(regDF)
end

# Objective function for estimation
function ObjectiveFunction(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
    verbose=false, inner_tol=1e-5, inner_maxIter=3000, inner_display=false,
    require_convergence=true, continuation_steps=5, employment_change=:log,
    wage_center=:treated, max_abs_moment=10.0, moment_firm_mask=nothing,
    firm_ids=nothing, reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false)
    η, θ = params_to_estimate
    
    # Return large penalty for invalid parameters
    if η <= 0 || θ <= 0 || !isfinite(η) || !isfinite(θ)
        return 1e10
    end
    
    β_model = ComputeModelMoments(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, α,
        inner_tol=inner_tol, inner_maxIter=inner_maxIter, inner_display=inner_display,
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
        println("η=$η, θ=$θ → β_model=$β_model, obj=$obj")
    end
    
    return obj
end

function EvaluateCalibrationGrid(; l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
    η_grid, θ_grid, inner_tol=2e-5, inner_maxIter=3000,
    continuation_steps=5, employment_change=:log, wage_center=:treated,
    max_abs_moment=10.0, verbose=true, moment_firm_mask=nothing,
    firm_ids=nothing, reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false)

    results = DataFrame(
        η = Float64[],
        θ = Float64[],
        β1 = Float64[],
        β2 = Float64[],
        objective = Float64[]
    )

    total = length(η_grid) * length(θ_grid)
    counter = 0
    for η in η_grid, θ in θ_grid
        counter += 1
        β_model = ComputeModelMoments([η, θ]; l, d, d′, wⱼ_data, lⱼ_data, α,
            inner_tol=inner_tol, inner_maxIter=inner_maxIter,
            inner_display=false, require_convergence=true,
            continuation_steps=continuation_steps,
            employment_change=employment_change, wage_center=wage_center,
            moment_firm_mask=moment_firm_mask, firm_ids=firm_ids,
            reg_sample_ids=reg_sample_ids, restrict_to_reg_sample=restrict_to_reg_sample)

        if any(isnan.(β_model)) || any(isinf.(β_model)) || any(abs.(β_model) .> max_abs_moment)
            obj = 1e10
        else
            obj = sum((β_model .- β_target) .^ 2)
        end

        push!(results, (η, θ, β_model[1], β_model[2], obj))
        if verbose
            println("Grid ", counter, "/", total, ": η=", η, ", θ=", θ,
                " → β_model=", β_model, ", obj=", obj)
        end
    end

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

function CalibrateEtaTheta(; l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
    x0=[2.75, 2.0], starts=nothing, lower=[0.25, 0.25], upper=[8.0, 8.0],
    iterations=250, x_abstol=1e-4, f_reltol=1e-8,
    inner_tol=1e-5, inner_maxIter=3000, continuation_steps=5,
    employment_change=:log, wage_center=:treated, max_abs_moment=10.0,
    verbose=true, show_trace=true, moment_firm_mask=nothing, firm_ids=nothing,
    reg_sample_ids=nothing, restrict_to_reg_sample::Bool=false)

    lower = Float64.(lower)
    upper = Float64.(upper)

    function objective_y(y)
        x = _unconstrained_to_box(y, lower, upper)
        return ObjectiveFunction(x; l, d, d′, wⱼ_data, lⱼ_data, α, β_target,
            verbose=verbose, inner_tol=inner_tol, inner_maxIter=inner_maxIter,
            inner_display=false, require_convergence=true,
            continuation_steps=continuation_steps,
            employment_change=employment_change, wage_center=wage_center,
            max_abs_moment=max_abs_moment, moment_firm_mask=moment_firm_mask,
            firm_ids=firm_ids, reg_sample_ids=reg_sample_ids,
            restrict_to_reg_sample=restrict_to_reg_sample)
    end

    start_list = isnothing(starts) ? [x0] : starts
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
