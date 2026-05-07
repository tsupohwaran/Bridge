#=================================================================#
# Prepare data for the paper
# Date: March 2025
# Author: Wu Chengjun, Central University of Finance and Economics
# OS: MacOS 15.3.2
# Version: 1.10.9
#=================================================================#

#==================================================#
# Load packages and functions
# !!!Path must be redefined by users!!!
#==================================================#

projPath = "/Users/pohwaran/Doctorate/Paper/Bridge"
cd(projPath)
include(projPath * "/code/03_model/load_packages.jl")
include(projPath * "/code/03_model/functions.jl")

#==================================================#
# Counterfactual
#==================================================#

# trade cost and tariff
@load "data/model/raw/tariff_imp_sec_exp_00-21.jld2" τ_00 τ_07 τ_17 τ_19
@load "data/model/raw/trade_cost_asym.jld2" κʲ_00 κʲ_07 κʲ_17 κʲ_19

#### 00-07 ####
@load "data/model/raw/model_data_raw_00.jld2" inputData vars params
inputData, vars, params = SolveModel(inputData, vars, params, ones(size(params.τʲ)), params.τʲ; deficit = true, numer = 2);

# Counterfactual: only tariff change
τʲ′ = τ_07;
κ̂ʲ_tariff = (1 .+ τʲ′) ./ (1 .+ params.τʲ); # no change for non-trariff trade cost
_, _, _, changes, check = SolveModel(inputData, vars, params, κ̂ʲ_tariff, τʲ′; deficit = true, numer = 2);

# Counterfactual: both tariff and non-tariff change
κ̂ʲ = κʲ_07 ./ κʲ_00; # total trade cost change
κ̂ʲ[:, 16:end, :] .= 1.0;  # no change for non-trade sectors
_, _, _, changes, check = SolveModel(inputData, vars, params, κ̂ʲ, τʲ′; deficit = true, numer = 2);

# Counterfactual: only non-tariff change
κ̂ʲ_nontariff = κ̂ʲ ./ κ̂ʲ_tariff;
_, _, _, changes, check = SolveModel(inputData, vars, params, κ̂ʲ_nontariff, τʲ′; numer = 2);

changes.dlnOʷ
changes.dlnOᵉᵘ
changes.dlnWʷ
cliparray(changes.dlnIn)
cliparray(changes.dlnw)
cliparray(changes.dlnOₙ)
cliparray(changes.dlnW)
cliparray(changes.dlnpᶠ)
cliparray(changes.dlncʲ[:, 9])
changes.dlnOʲʰ[15, 9]
changes.dlnYʲ[:, 9]
changes.dlnpʲ[:, 9]
changes.dlnw

#==================================================#
# Counterfactual: US trade war 17-19
#==================================================#

@load "data/model/raw/model_data_raw_17.jld2" inputData vars params
inputData, vars, params = SolveModel(inputData, vars, params, ones(size(params.τʲ)), params.τʲ; deficit=true, numer = 2);

# Counterfactual: only tariff change
τʲ′ = τ_19;
κ̂ʲ = (1 .+ τʲ′) ./ (1 .+ params.τʲ); # no change for non-trariff trade cost
_, _, _, changes, check = SolveModel(inputData, vars, params, κ̂ʲ, τʲ′; deficit=true, numer = 2);

# Counterfactual: both tariff and non-tariff change
κ̂ʲ = κʲ_19 ./ κʲ_17; # total trade cost change
κ̂ʲ[:, 16:end, :] .= 1.0;  # no change for non-trade sectors
_, _, _, changes, check = SolveModel(inputData, vars, params, κ̂ʲ, τʲ′; deficit = true, numer = 2);

# Counterfactual: No US trade war
τʲ′[16, :, :] .= params.τʲ[16, :, :];
τʲ′[:, :, 16] .= params.τʲ[:, :, 16];
κ̂ʲ = (1 .+ τʲ′) ./ (1 .+ params.τʲ);
_, _, _, changes, check = SolveModel(inputData, vars, params, κ̂ʲ, τʲ′; numer = 2);

#==================================================#
# Counterfactual: US "recipocal" tariff 
# 47.3% for China
#==================================================#

@load "data/model/processed/model_data_17.jld2" inputData vars params

# Counterfactual: only tariff change
τʲ′ = copy(params.τʲ); # warning !!! "=" is reference
τʲ′[16, :, 9] .= 0.473;
κ̂ʲ_tariff = (1 .+ τʲ′) ./ (1 .+ params.τʲ); # no change for non-trariff trade cost
κ̂ʲ_tariff[:, 16:end, :] .= 1.0;  # no change for non-trade sectors
_, _, _, changes, check = SolveModel(inputData, vars, params, κ̂ʲ_tariff, τʲ′; numer = 2);

