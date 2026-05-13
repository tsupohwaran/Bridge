@show size(df)
@show count(ismissing, df[!, :employ])
@show findfirst(ismissing, df[!, :employ])
