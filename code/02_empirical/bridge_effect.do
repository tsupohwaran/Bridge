*==================================================*
* Effect of Cross-sea Bridge on Firms
* Date: December 9, 2025
* Author: Wu Chengjun, Central University of Finance and Economics
* OS: MacOS 26.1
* Version: Stata MP 18.0
*==================================================*

* set path
clear all
global proj_path = "/Users/pohwaran/Doctorate/Paper/Bridge"
global firm_data_path = "$proj_path/data"
global temp_path = "$firm_data_path/temp"
global processed_path = "$firm_data_path/processed"
global cied_data_path = "$firm_data_path/raw/cied"
global ctsd_data_path = "$firm_data_path/raw/ctsd"
global geo_data_path = "$proj_path/data/raw/geo"

*==================================================*
* Step 1: Data Preparation for Facts Analysis
*==================================================*

use "$processed_path/regression_qingdao_07_20.dta", clear
merge 1:1 id year using "$processed_path/result_markdown_est_qingdao_07_20.dta", keep(1 3) nogen

/* preserve
    duplicates drop id, force
    count
    kdensity dln_ma, normal xlabel(-8.5(0.5)5.0, angle(45))
restore */

* Dependent variable
gen ln_employ = log(employ)
gen ln_md = log(md_TL)
gen ln_wage = log(wage_total / employ)

* Cross-sea Bridge Shock
gen post = (year > 2011)
gen big_ma = 0 if (dln_ma <= -1) | missing(dln_ma)
replace big_ma = 1 if missing(big_ma)

* control variables
egen ftype = group(firm_type)
egen town_year = group(town year)
gen ln_w0 = log(wage_inital)
gen ln_age = ln(age)

save "$processed_path/regression_qingdao_07_20.dta", replace


* merge wage with no security
use "$processed_path/regression_qingdao_07_20.dta", clear
merge 1:1 id year using "$temp_path/ctsd_qingdao_07_20_nosec.dta", nogen

* merge with exposure variable
use "$processed_path/regression_qingdao_07_20.dta", clear
egen firm_id = group(id)
xtset firm_id year

merge m:1 id using "$processed_path/firm_exposure_qingdao.dta", nogen keep(3)
gen ln_exposure_huangdao = log(exposure_huangdao)
gen ln_exposure_jimo = log(exposure_jimo)

sort firm_id year 
gen wage_hat = ln_wage - L.ln_wage
reg ln_wage ln_exposure_huangdao ln_exposure_jimo if year == 2012 | year == 2013
reg wage_hat dln_ma ln_exposure_huangdao ln_exposure_jimo if year == 2013

*==================================================*
* Step 1: Data Preparation for Facts Analysis
*==================================================*

use "$processed_path/regression_qingdao_07_20.dta", clear
global year_choice="2010 2011 2012 2013"
global control = "i.export_bool#year ftype#year ind_code2#year"

****** Labor Employment Effect ******
* full sample
eststo clear
preserve
    drop if missing(ln_employ, ln_age)
    bys id: gen repeat = _N
    drop if repeat != 5 //keep firms with at least 5 years of data

    sum ln_w0 if (year == 2010) & (big_ma == 1)
    gen ln_w0_mean = r(mean)
    gen ln_w0_diff = (ln_w0 - ln_w0_mean) // 使用不同样本时，这个mean值需要重新计算

    sum ln_w0_diff if (year == 2010) & (big_ma == 1), detail
    sca ln_w0_diff_1 = r(p5)
    sca ln_w0_diff_2 = r(p25)
    sca ln_w0_diff_3 = r(p50)
    sca ln_w0_diff_4 = r(p75)
    sca ln_w0_diff_5 = r(p90)
    sca ln_w0_diff_6 = r(p95)

    * baseline parallel trend test
    reghdfe ln_employ i1.big_ma#i($year_choice).year ///
    ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    coefplot, vertical keep(1.big_ma#*.year) ///
    xlabel(1 "2010" 2 "2011" 3 "2012" 4 "2013") ///
    title("Employment Effect by Year (税调+工企)") ///
    ciopts(lpattern(dash) recast(rcap) msize(medium)) ///
    addplot(line @b @at) ///
    msymbol(circle_hollow) ///
    msize(*0.5) color(gs0) ///
    yline(0, lpattern(dash)) yline(0, lp(dash) lcolor(red)) ///
    graphregion(fcolor(gs16) lcolor(gs16)) ///
    plotregion(lpattern(blank)) scheme(s1mono) ///
    legend(order(1 "95% CI") position(2) ring(0) bplacement(neast))

    * Heterogeneous trend test with wage difference
    reghdfe ln_employ i1.big_ma#i($year_choice).year ///
    c.ln_w0_diff#i1.big_ma#i($year_choice).year ///
    c.ln_w0_diff#i($year_choice).year ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    do "$proj_path/code/01_data_prep/05_plot_wage_trend.do"

    * Heterogeneous trend test with wage difference (Dummy)
    gen high_ln_w0_diff = (ln_w0_diff > 0)
    reghdfe ln_employ i1.big_ma#i($year_choice).year ///
    i1.high_ln_w0_diff#i1.big_ma#i($year_choice).year ///
    i1.high_ln_w0_diff#i($year_choice).year ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    coefplot, vertical keep(1.high_ln_w0_diff#1.big_ma#*.year) ///
    xlabel(1 "2010" 2 "2011" 3 "2012" 4 "2013") ///
    title("Employment Effect by Year (税调+工企)") ///
    ciopts(lpattern(dash) recast(rcap) msize(medium)) ///
    addplot(line @b @at) ///
    msymbol(circle_hollow) ///
    msize(*0.5) color(gs0) ///
    yline(0, lpattern(dash)) yline(0, lp(dash) lcolor(red)) ///
    graphregion(fcolor(gs16) lcolor(gs16)) ///
    plotregion(lpattern(blank)) scheme(s1mono) ///
    legend(order(1 "95% CI") position(2) ring(0) bplacement(neast))

    do "$proj_path/code/01_data_prep/05_plot_wage_trend_dummy.do"

    eststo: reghdfe ln_employ i1.big_ma#i1.post ///
    ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    eststo: reghdfe ln_employ i1.big_ma#i1.post ///
    c.ln_w0_diff#i1.big_ma#i1.post ///
    c.ln_w0_diff#i1.post ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)
restore

* only CTSD sample
preserve
    drop if missing(ln_employ, ln_age)
    bys id: gen repeat = _N
    drop if repeat != 5 //keep firms with at least 5 years of data
    keep if strlen(id) == 9

    sum ln_w0 if (year == 2010) & (big_ma == 1)
    gen ln_w0_mean = r(mean)
    gen ln_w0_diff = (ln_w0 - ln_w0_mean) // 使用不同样本时，这个mean值需要重新计算

    * baseline parallel trend test
    reghdfe ln_employ i1.big_ma#i($year_choice).year ///
    ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    coefplot, vertical keep(1.big_ma#*.year) ///
    xlabel(1 "2010" 2 "2011" 3 "2012" 4 "2013") ///
    title("Employment Effect by Year (税调+工企)") ///
    ciopts(lpattern(dash) recast(rcap) msize(medium)) ///
    addplot(line @b @at) ///
    msymbol(circle_hollow) ///
    msize(*0.5) color(gs0) ///
    yline(0, lpattern(dash)) yline(0, lp(dash) lcolor(red)) ///
    graphregion(fcolor(gs16) lcolor(gs16)) ///
    plotregion(lpattern(blank)) scheme(s1mono) ///
    legend(order(1 "95% CI") position(2) ring(0) bplacement(neast))

    * Heterogeneous trend test with wage difference
    reghdfe ln_employ i1.big_ma#i($year_choice).year ///
    c.ln_w0_diff#i1.big_ma#i($year_choice).year ///
    c.ln_w0_diff#i($year_choice).year ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    do "$proj_path/code/01_data_prep/05_plot_wage_trend.do"

    * Heterogeneous trend test with wage difference (Dummy)
    gen high_ln_w0_diff = (ln_w0_diff > 0)
    reghdfe ln_employ i1.big_ma#i($year_choice).year ///
    i1.high_ln_w0_diff#i1.big_ma#i($year_choice).year ///
    i1.high_ln_w0_diff#i($year_choice).year ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    coefplot, vertical keep(1.high_ln_w0_diff#1.big_ma#*.year) ///
    xlabel(1 "2010" 2 "2011" 3 "2012" 4 "2013") ///
    title("Employment Effect by Year (税调+工企)") ///
    ciopts(lpattern(dash) recast(rcap) msize(medium)) ///
    addplot(line @b @at) ///
    msymbol(circle_hollow) ///
    msize(*0.5) color(gs0) ///
    yline(0, lpattern(dash)) yline(0, lp(dash) lcolor(red)) ///
    graphregion(fcolor(gs16) lcolor(gs16)) ///
    plotregion(lpattern(blank)) scheme(s1mono) ///
    legend(order(1 "95% CI") position(2) ring(0) bplacement(neast))

    do "$proj_path/code/01_data_prep/05_plot_wage_trend_dummy.do"

    eststo: reghdfe ln_employ i1.big_ma#i1.post ///
    ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    eststo: reghdfe ln_employ i1.big_ma#i1.post ///
    c.ln_w0_diff#i1.big_ma#i1.post ///
    c.ln_w0_diff#i1.post ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)
restore

****** Markdown Effect ******
preserve
    drop if missing(ln_md, ln_age)
    bys id: gen repeat = _N
    drop if repeat != 5 //keep firms with at least 5 years of data

    sum ln_w0 if (year == 2010) & (big_ma == 1)
    gen ln_w0_mean = r(mean)
    gen ln_w0_diff = (ln_w0 - ln_w0_mean) // 使用不同样本时，这个mean值需要重新计算

    * baseline parallel trend test
    reghdfe ln_md i1.big_ma#i($year_choice).year ///
    ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    coefplot, vertical keep(1.big_ma#*.year) ///
    xlabel(1 "2010" 2 "2011" 3 "2012" 4 "2013") ///
    title("Markdown Effect by Year (税调+工企)") ///
    ciopts(lpattern(dash) recast(rcap) msize(medium)) ///
    addplot(line @b @at) ///
    msymbol(circle_hollow) ///
    msize(*0.5) color(gs0) ///
    yline(0, lpattern(dash)) yline(0, lp(dash) lcolor(red)) ///
    graphregion(fcolor(gs16) lcolor(gs16)) ///
    plotregion(lpattern(blank)) scheme(s1mono) ///
    legend(order(1 "95% CI") position(2) ring(0) bplacement(neast))

    * Heterogeneous trend test with wage difference
    reghdfe ln_md i1.big_ma#i($year_choice).year ///
    c.ln_w0_diff#i1.big_ma#i($year_choice).year ///
    c.ln_w0_diff#i($year_choice).year ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    do "$proj_path/code/01_data_prep/05_plot_wage_trend.do"

    * Heterogeneous trend test with wage difference (Dummy)
    gen high_ln_w0_diff = (ln_w0_diff > 0)
    reghdfe ln_md i1.big_ma#i($year_choice).year ///
    i1.high_ln_w0_diff#i1.big_ma#i($year_choice).year ///
    i1.high_ln_w0_diff#i($year_choice).year ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    coefplot, vertical keep(1.high_ln_w0_diff#1.big_ma#*.year) ///
    xlabel(1 "2010" 2 "2011" 3 "2012" 4 "2013") ///
    title("Markdown Effect by Year (税调+工企)") ///
    ciopts(lpattern(dash) recast(rcap) msize(medium)) ///
    addplot(line @b @at) ///
    msymbol(circle_hollow) ///
    msize(*0.5) color(gs0) ///
    yline(0, lpattern(dash)) yline(0, lp(dash) lcolor(red)) ///
    graphregion(fcolor(gs16) lcolor(gs16)) ///
    plotregion(lpattern(blank)) scheme(s1mono) ///
    legend(order(1 "95% CI") position(2) ring(0) bplacement(neast))

    do "$proj_path/code/01_data_prep/05_plot_wage_trend_dummy.do"

    eststo: reghdfe ln_md i1.big_ma#i1.post ///
    ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    eststo: reghdfe ln_md i1.big_ma#i1.post ///
    c.ln_w0_diff#i1.big_ma#i1.post ///
    c.ln_w0_diff#i1.post ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)
restore

****** Wage Effect ******
preserve
    drop if missing(ln_wage, ln_age)
    bys id: gen repeat = _N
    drop if repeat != 5 //keep firms with at least 5 years of data

    sum ln_w0 if (year == 2010) & (big_ma == 1)
    gen ln_w0_mean = r(mean)
    gen ln_w0_diff = (ln_w0 - ln_w0_mean) // 使用不同样本时，这个mean值需要重新计算

    * baseline parallel trend test
    reghdfe ln_wage i1.big_ma#i($year_choice).year ///
    ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    coefplot, vertical keep(1.big_ma#*.year) ///
    xlabel(1 "2010" 2 "2011" 3 "2012" 4 "2013") ///
    title("Wage Effect by Year (税调+工企)") ///
    ciopts(lpattern(dash) recast(rcap) msize(medium)) ///
    addplot(line @b @at) ///
    msymbol(circle_hollow) ///
    msize(*0.5) color(gs0) ///
    yline(0, lpattern(dash)) yline(0, lp(dash) lcolor(red)) ///
    graphregion(fcolor(gs16) lcolor(gs16)) ///
    plotregion(lpattern(blank)) scheme(s1mono) ///
    legend(order(1 "95% CI") position(2) ring(0) bplacement(neast))

    * Heterogeneous trend test with wage difference
    reghdfe ln_wage i1.big_ma#i($year_choice).year ///
    c.ln_w0_diff#i1.big_ma#i($year_choice).year ///
    c.ln_w0_diff#i($year_choice).year ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    do "$proj_path/code/01_data_prep/05_plot_wage_trend.do"

    * Heterogeneous trend test with wage difference (Dummy)
    gen high_ln_w0_diff = (ln_w0_diff > 0)
    reghdfe ln_wage i1.big_ma#i($year_choice).year ///
    i1.high_ln_w0_diff#i1.big_ma#i($year_choice).year ///
    i1.high_ln_w0_diff#i($year_choice).year ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    coefplot, vertical keep(1.high_ln_w0_diff#1.big_ma#*.year) ///
    xlabel(1 "2010" 2 "2011" 3 "2012" 4 "2013") ///
    title("Wage Effect by Year (税调+工企)") ///
    ciopts(lpattern(dash) recast(rcap) msize(medium)) ///
    addplot(line @b @at) ///
    msymbol(circle_hollow) ///
    msize(*0.5) color(gs0) ///
    yline(0, lpattern(dash)) yline(0, lp(dash) lcolor(red)) ///
    graphregion(fcolor(gs16) lcolor(gs16)) ///
    plotregion(lpattern(blank)) scheme(s1mono) ///
    legend(order(1 "95% CI") position(2) ring(0) bplacement(neast))

    do "$proj_path/code/01_data_prep/05_plot_wage_trend_dummy.do"

    eststo: reghdfe ln_wage i1.big_ma#i1.post ///
    ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)

    eststo: reghdfe ln_wage i1.big_ma#i1.post ///
    c.ln_w0_diff#i1.big_ma#i1.post ///
    c.ln_w0_diff#i1.post ln_age if (repeat == 5), ///
    a(id year $control) cluster(town_year)
restore

esttab using "$processed_path/bridge_effect_estimates.csv", b(%12.3fc) se ar2 nogap star(* 0.10 ** 0.05 *** 0.01) replace
