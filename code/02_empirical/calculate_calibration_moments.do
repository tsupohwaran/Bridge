* set path
clear all
global proj_path = "/Users/pohwaran/Doctorate/Paper/Bridge"
do "$proj_path/code/00_setup/data_paths.do"


use "$regression_processed_path/regression_qingdao_07_20.dta", clear
keep if year>=2009 & year<=2010
collapse (max) export export_intensity, by(id)
gen export_bool=(export>0)
replace export_intensity=1 if export_intensity>1 & export_intensity<.
gen lnexp_intensity=ln(export_intensity)
replace lnexp_intensity=0 if export_bool==0
save "$regression_temp_path/id_exp.dta", replace

************************************************
use "$regression_processed_path/regression_qingdao_07_20.dta", clear
merge m:1 id using "$regression_temp_path/id_exp.dta", keepusing(export_bool lnexp_intensity) nogenerate

encode firm_type,gen(ftype)
encode town, gen(town2)
gen lnage = ln(age)
drop if dma==.
gen w_ns=wage_total/empl

keep if year>=2009 & year<=2012
gen w=w_ns
gen w2010=w if year==2010
bysort id: egen w0=mean(w2010)
gen lnemp=ln(emp)
// gen lnmd=ln(md_TL)
gen lnw=ln(w)
gen lnw0=ln(w0)
winsor2 lnw0, cuts(5 95) replace

gen BIG=(dma>0.5) //no missing value


gen nempl=(lnemp!=.)
// gen nmd=(lnmd!=.)
gen nw=(lnw!=.)
bysort id: egen Nempl=total(nempl)
// bysort id: egen Nmd=total(nmd)
bysort id: egen Nw=total(nw)


global year_choice="2010 2011 2012"
global ref_year = 2009
global control = "i.export_bool#year c.lnexp_intensity#year  ftype#year ind_code2#year"

gen post=(year>=2011)

egen x=mean(lnw0) if year==2010
egen mlnw0=mean(x)
gen demean_lnw0=lnw0-mlnw0 //demean: lnw0
drop if missing(demean_lnw0)

gen lndma=ln(dma) if BIG==1
su lndma
replace lndma = r(min)-0.1 if BIG==0
su lndma if year == 2010
gen lndma_demean =lndma - r(mean)  //demean: lndma

****************************************************


*******
*moment
*******

reghdfe lnemp i1.BIG#i1.post c.demean_lnw0#i1.BIG#i1.post lnage if (year==2010 | year==2012), a(id year c.demean_lnw0#year $control) cluster(town2#ind)
gen byte in_emp_reg = e(sample)

reghdfe lnw i1.BIG#i1.post lnage if (year==2010 | year==2012), a(id year c.demean_lnw0#year $control) cluster(town2#ind)
gen byte in_w_reg = e(sample)

gen byte in_common_reg = in_emp_reg==1 & in_w_reg==1

reghdfe lnemp i1.BIG#i1.post c.demean_lnw0#i1.BIG#i1.post c.demean_lnw0#i1.post lnage if in_common_reg, a(id year  $control) cluster(town2#ind)
scalar beta_labor_bigMA = _b[1.BIG#1.post]
scalar beta_labor_bigMA_wdiff = _b[c.demean_lnw0#1.BIG#1.post]

reghdfe lnw i1.BIG#i1.post lnage if in_common_reg, a(id year c.demean_lnw0#year $control) cluster(town2#ind)
scalar beta_wage_bigMA = _b[1.BIG#1.post]

cap mkdir "$proj_path/output"
cap mkdir "$output_table_path"
cap mkdir "$proj_path/output/figures"

* Real-data analogue of the Julia labor reallocation cloud.
* x: ln(dMA); y: Delta ln employment; color: initial wage deviation.
preserve
    keep if in_common_reg == 1 & (year == 2010 | year == 2012)

    bysort id: egen lnemp_2010 = max(cond(year == 2010, lnemp, .))
    bysort id: egen lnemp_2012 = max(cond(year == 2012, lnemp, .))
    bysort id: egen dma_plot = max(cond(year == 2010, dma, .))
    bysort id: egen wdiff_plot = max(cond(year == 2010, demean_lnw0, .))

    gen double dlnemp_plot = lnemp_2012 - lnemp_2010
    gen double dlnemp_plot_clip = dlnemp_plot
    replace dlnemp_plot_clip = 5 if dlnemp_plot_clip > 5 & !missing(dlnemp_plot_clip)

    gen double ln_dMA_plot = ln(dma_plot) if dma_plot > 0
    quietly count if !missing(ln_dMA_plot)
    local n_positive_dma = r(N)

    if `n_positive_dma' > 0 {
        quietly summarize ln_dMA_plot, meanonly
        replace ln_dMA_plot = r(min) - 0.1 if missing(ln_dMA_plot) & !missing(dma_plot)

        keep if !missing(dlnemp_plot_clip, ln_dMA_plot, wdiff_plot)
        bysort id: keep if _n == 1

        quietly count
        local n_plot = r(N)

        if `n_plot' > 0 {
            xtile wage_bin = wdiff_plot, nq(7)

            twoway ///
                (scatter dlnemp_plot_clip ln_dMA_plot if wage_bin == 1, ///
                    mcolor("103 0 31") msymbol(O) msize(tiny) mlcolor(none)) ///
                (scatter dlnemp_plot_clip ln_dMA_plot if wage_bin == 2, ///
                    mcolor("178 24 43") msymbol(O) msize(tiny) mlcolor(none)) ///
                (scatter dlnemp_plot_clip ln_dMA_plot if wage_bin == 3, ///
                    mcolor("214 96 77") msymbol(O) msize(tiny) mlcolor(none)) ///
                (scatter dlnemp_plot_clip ln_dMA_plot if wage_bin == 4, ///
                    mcolor("150 150 150") msymbol(O) msize(tiny) mlcolor(none)) ///
                (scatter dlnemp_plot_clip ln_dMA_plot if wage_bin == 5, ///
                    mcolor("146 197 222") msymbol(O) msize(tiny) mlcolor(none)) ///
                (scatter dlnemp_plot_clip ln_dMA_plot if wage_bin == 6, ///
                    mcolor("67 147 195") msymbol(O) msize(tiny) mlcolor(none)) ///
                (scatter dlnemp_plot_clip ln_dMA_plot if wage_bin == 7, ///
                    mcolor("33 102 172") msymbol(O) msize(tiny) mlcolor(none)), ///
                xtitle("ln(dMA)", size(medsmall)) ///
                ytitle("Delta ln employment, 2012 - 2010", size(medsmall)) ///
                title("Employment Change vs Market Access Change", size(medium)) ///
                yline(0, lpattern(dash) lcolor(gs10)) ///
                legend(title("ln w0 - mean", size(small)) ///
                       order(1 "Q1 lowest" 2 "Q2" 3 "Q3" 4 "Q4 middle" ///
                             5 "Q5" 6 "Q6" 7 "Q7 highest") ///
                       position(3) ring(1) cols(1) size(small) region(lcolor(none))) ///
                note("Nonpositive dMA values are placed just left of the positive ln(dMA) support.", size(vsmall)) ///
                graphregion(fcolor(white) lcolor(white)) ///
                plotregion(fcolor(white))

            graph export "$proj_path/output/figures/labor_reallocation_realdata.png", ///
                replace width(6000)
        }
    }
restore

preserve
    keep if in_common_reg == 1 & (year == 2010 | year == 2012)

    bysort id: egen lnemp_2010 = max(cond(year == 2010, lnemp, .))
    bysort id: egen lnemp_2012 = max(cond(year == 2012, lnemp, .))
    bysort id: egen wdiff_plot = max(cond(year == 2010, demean_lnw0, .))
    bysort id: egen treated_plot = max(BIG)

    gen dlnemp_plot = lnemp_2012 - lnemp_2010
    keep if !missing(dlnemp_plot, wdiff_plot, treated_plot)
    bysort id: keep if _n == 1

    count if treated_plot == 0
    local n_ctrl = r(N)
    count if treated_plot == 1
    local n_trt = r(N)

    quietly regress dlnemp_plot c.wdiff_plot if treated_plot == 0
    local slope_ctrl_num = _b[wdiff_plot]
    local slope_ctrl = strtrim(string(`slope_ctrl_num', "%9.3f"))

    quietly regress dlnemp_plot c.wdiff_plot if treated_plot == 1
    local slope_trt_num = _b[wdiff_plot]
    local slope_trt = strtrim(string(`slope_trt_num', "%9.3f"))
    local slope_gap_num = `slope_trt_num' - `slope_ctrl_num'
    local slope_gap = strtrim(string(`slope_gap_num', "%9.3f"))

    twoway ///
        (scatter dlnemp_plot wdiff_plot if treated_plot == 0, ///
            mcolor(navy%40) msymbol(O) msize(vsmall) mlcolor(none)) ///
        (scatter dlnemp_plot wdiff_plot if treated_plot == 1, ///
            mcolor(maroon%40) msymbol(O) msize(vsmall) mlcolor(none)) ///
        (lfit dlnemp_plot wdiff_plot if treated_plot == 0, ///
            lcolor(navy) lpattern(dash) lwidth(medthick)) ///
        (lfit dlnemp_plot wdiff_plot if treated_plot == 1, ///
            lcolor(maroon) lpattern(dash) lwidth(medthick)), ///
        xtitle("initial wage deviation (ln w0 - mean)") ///
        ytitle("Delta ln employment, 2012 - 2010") ///
        title("Employment wage-slope: treated vs control (gap = `slope_gap')") ///
        legend(order(1 "control firms (N = `n_ctrl')" ///
                     2 "treated firms (N = `n_trt')" ///
                     3 "control slope = `slope_ctrl'" ///
                     4 "treated slope = `slope_trt'") ///
               position(11) ring(0) cols(1) region(lcolor(black))) ///
        graphregion(fcolor(white) lcolor(white)) ///
        plotregion(fcolor(white))

    graph export "$proj_path/output/figures/labor_wdiff_slope_regression.png", ///
        replace width(6000)
restore

preserve
    clear
    set obs 3
    gen str24 moment = ""
    gen double beta = .
    replace moment = "labor_bigMA" in 1
    replace beta = scalar(beta_labor_bigMA) in 1
    replace moment = "labor_bigMA_wdiff" in 2
    replace beta = scalar(beta_labor_bigMA_wdiff) in 2
    replace moment = "wage_bigMA" in 3
    replace beta = scalar(beta_wage_bigMA) in 3
    export delimited using "$output_table_path/calibration_moments.csv", replace
restore

keep if in_common_reg==1
duplicates drop id, force
keep id
save "$model_processed_path/firm_two_year_reg.dta", replace
