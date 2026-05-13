* corrlation between firm size and markdown
use "taxsurvey_markdown_08_20.dta", clear
gen ln_md = ln(md_TL)
gen ln_wage_total = ln(wage_total)
gen ln_employ = ln(employ_avg)

eststo: reghdfe ln_md ln_wage_total, absorb(id year) vce(cluster id)
eststo: reghdfe ln_md ln_employ, absorb(id year) vce(cluster id)

* drop state owned enterprises
eststo: reghdfe ln_md ln_wage_total if firm_type != "soe", absorb(id year) vce(cluster id)
eststo: reghdfe ln_md ln_employ if firm_type != "soe", absorb(id year) vce(cluster id)

esttab, se label b(3) se(3) star(* 0.1 ** 0.05 *** 0.01) stats(N r2, fmt(3 3) labels("Observations" "R-squared")) title("Correlation between firm size and markdown") varwidth(20)