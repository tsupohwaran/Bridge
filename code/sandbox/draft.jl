SumStats(wⱼ_data, "dlnlⱼ")
findall(dlnlⱼ .> 0)


wⱼ_data = wⱼ_raw |>
    x -> x ./ sum(x .* lⱼ_data); # normalize total wage bill to 1