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
    displaySummary=false)

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

    if displaySummary
        if iter < maxIter
            println("Successful convergence in ", iter, " iterations.")
        else
            println("Maximum number of iterations reached.")
        end
    end

    return X
end

# Solve the equilibrium
function SolveModel(vars::NamedTuple, params::NamedTuple;
    damp=0.8, tol=1e-10, maxIter=1e3, power=false,
    displayGap=false, displaySummary=false,
    wⱼ_init=nothing)  # optional warm start

    # unpack the data
    (; l, d, zⱼ) = vars
    (; η, θ, α) = params

    # initial guess (use warm start if provided)
    wⱼ = isnothing(wⱼ_init) ? ones(J) : copy(wⱼ_init)

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
    wⱼ = Converge(x -> UpdateRule(x)[1], wⱼ;
        tol=tol, maxIter=maxIter,
        damp=damp, power=power,
        normalize = x -> x ./ mean(x),
        displayGap=displayGap,
        displaySummary=displaySummary)

    return UpdateRule(wⱼ)
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

# Compute model moments given parameters
function ComputeModelMoments(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, α)
    η, θ = params_to_estimate
    
    # Ensure parameters are in valid range
    if η <= 0 || θ <= 0
        return [Inf, Inf]
    end
    
    # Adaptive damping: higher θ needs more damping for stability
    damp = clamp(0.5 + 0.03 * θ, 0.6, 0.97)
    
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
    wⱼ, π_zj, ε_zj, lⱼ, εⱼ = SolveModel(vars_base, params;
        damp=damp, tol=1e-7, power=false, maxIter=1000,
        wⱼ_init=wⱼ_data, displaySummary=true)
    
    # Check for invalid wages
    if any(isnan.(wⱼ)) || any(isinf.(wⱼ))
        return [Inf, Inf]
    end
    
    # Step 2c: Solve counterfactual model (use baseline as warm start)
    vars_cf = (; l, d = d′, zⱼ)
    wⱼ′, π_zj′, ε_zj′, lⱼ′, εⱼ′ = SolveModel(vars_cf, params;
        damp=damp, tol=1e-7, power=false, maxIter=1000,
        wⱼ_init=wⱼ, displaySummary=true)
    
    # Check for invalid counterfactual wages
    if any(isnan.(wⱼ′)) || any(isinf.(wⱼ′))
        return [Inf, Inf]
    end
    
    # Step 2d: Compute regression moments
    l̂ⱼ = lⱼ′ ./ lⱼ
    dlnlⱼ = l̂ⱼ .- 1
    dlnMA = log.(sum((d - d′) .* l, dims = 1)') |> x -> replace(x, -Inf => -8)
    bigMA = Float64.(dlnMA .>= -1)
    
    regDF = DataFrame(bigMA = vec(bigMA), dlnl = vec(dlnlⱼ), w = vec(log.(wⱼ)))
    regDF.w_treatedmean = fill(mean(regDF[regDF.bigMA .== 1, :w]), nrow(regDF))
    regDF.w_diff = regDF.w .- regDF.w_treatedmean
    regModel = lm(@formula(dlnl ~ bigMA + bigMA & w_diff + w_diff), regDF)
    
    β₁ = coef(regModel)[2]
    β₂ = coef(regModel)[4]
    
    return [β₁, β₂]
end

# Objective function for estimation
function ObjectiveFunction(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, α, β_target, verbose=false)
    η, θ = params_to_estimate
    
    # Return large penalty for invalid parameters
    if η <= 0 || θ <= 0 || !isfinite(η) || !isfinite(θ)
        return 1e10
    end
    
    β_model = ComputeModelMoments(params_to_estimate; l, d, d′, wⱼ_data, lⱼ_data, α)
    
    # Check for invalid model output
    if any(isnan.(β_model)) || any(isinf.(β_model))
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
