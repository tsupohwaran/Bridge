*==================================================*
* Markdown estimation
* Date: December 9, 2025
* Author: Wu Chengjun, Central University of Finance and Economics
* OS: MacOS 26.1
* Version: Stata MP 18.0
*==================================================*

*!!! need change to your path !!!*
clear all
global proj_path = "/Users/pohwaran/Doctorate/Paper/Bridge"
do "$proj_path/code/00_setup/data_paths.do"

*==================================================*
* Step 1: markdown estimation
* We encapsulate main estimation code in a function (MdEstTl), where "Tl" means translog production function
*==================================================*

* loop over 12 industries
forvalues ind = 1/12 {
    use "$regression_temp_path/markdown_est_07_20.dta", clear
    keep if ind == `ind'
    if _N == 0 continue
    xtset firm_id year, yearly
    do "$proj_path/code/utils/04_markdown_func.do" // main estimation code
    save "$ctsd_temp_path/md_est_07_20_ind`ind'.dta", replace
}

* append results
use "$ctsd_temp_path/md_est_07_20_ind1.dta", clear
forvalues ind = 2/12 {
    append using "$ctsd_temp_path/md_est_07_20_ind`ind'.dta"
}
gen md_TL = (1 / mu_DLW_TL) * theta_l_tl / alpha_l
drop if md_TL <= 0 | md_TL == .
drop if alpha_m >= 1

* winsorize md_TL by industry-year
egen ind_year = group(ind year)
winsor2 md_TL, cuts(5 95) by(ind_year) suffix(_w)
order id firm_name year md_TL md_TL_w tfp_TL omega_TL
save "$ctsd_processed_path/result_markdown_est_07_20.dta", replace

* calculate weighted mean of markdown by sector
* Step1: weighted mean of markdown by industry-year, weight by employ_avg
* Step2: mean of markdown by industry
use "$ctsd_processed_path/result_markdown_est_07_20.dta", clear
bys ind year: egen markdown_wt_jt = wtmean(md_TL_w), weight(l)
gduplicates drop ind year, force
bys ind: egen markdown_wt_j = wtmean(markdown_wt_jt), weight(employ)
gduplicates drop ind, force
keep ind ind_name markdown_wt_j
sort ind
save "$ctsd_processed_path/result_markdown_est_07_20_sector.dta", replace

* calculate weighted mean of markdown by year
use "$ctsd_processed_path/result_markdown_est_07_20.dta", clear
bys ind year: egen markdown_wt_jt = wtmean(md_TL_w), weight(l)
gduplicates drop ind year, force
bys year: egen markdown_wt_t = wtmean(markdown_wt_jt), weight(employ)
egen markdown_wt = wtmean(markdown_wt_jt), weight(employ)
gduplicates drop year, force
keep year markdown_wt_t markdown_wt

twoway (line markdown_wt_t year), ///
    title("Weighted Mean of Markdown by Year") ///
    xtitle("Year") ///
    xlabel(2008(1)2020) ///
    ytitle("Markdown") ///
    legend(off)

* merge qingdao
use "$ctsd_temp_path/ctsd_qingdao_07_20.dta", clear
keep sdid
gduplicates drop sdid, force
tempfile qingdao_sdid
save `qingdao_sdid', replace

use "$ctsd_processed_path/result_markdown_est_07_20.dta", clear
merge 1:1 sdid using `qingdao_sdid', keep(3) nogen
keep id year md_TL md_TL_w tfp_TL omega_TL
save "$ctsd_processed_path/result_markdown_est_qingdao_07_20.dta", replace
