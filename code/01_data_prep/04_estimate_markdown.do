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
global firm_data_path = "$proj_path/data"
global temp_path = "$firm_data_path/temp"
global processed_path = "$firm_data_path/processed"
global ctsd_data_path = "$firm_data_path/raw/ctsd"
global geo_data_path = "$proj_path/data/raw/geo"

*==================================================*
* Step 1: markdown estimation
* We encapsulate main estimation code in a function (MdEstTl), where "Tl" means translog production function
*==================================================*

* loop over 12 industries
forvalues ind = 1/12 {
    use "$processed_path/markdown_est_07_20.dta", clear
    keep if ind == `ind'
    if _N == 0 continue
    xtset firm_id year, yearly
    do "$proj_path/code/01_data_prep/04_markdown_func.do" // main estimation code
    save "$temp_path/md_est_07_20_ind`ind'.dta", replace
}

* append results
use "$temp_path/md_est_07_20_ind1.dta", clear
forvalues ind = 2/12 {
    append using "$temp_path/md_est_07_20_ind`ind'.dta"
}
gen md_TL = (1 / mu_DLW_TL) * theta_l_tl / alpha_l
drop if md_TL <= 0 | md_TL == .
drop if alpha_m >= 1

* winsorize md_TL by industry-year
egen ind_year = group(ind year)
winsor2 md_TL, cuts(5 95) by(ind_year) suffix(_w)
order id firm_name year md_TL md_TL_w
save "$processed_path/result_markdown_est_07_20.dta", replace

* calculate weighted mean of markdown by sector
* Step1: weighted mean of markdown by industry-year, weight by employ_avg
* Step2: mean of markdown by industry
use "$processed_path/result_markdown_est_07_20.dta", clear
bys ind year: egen markdown_wt_jt = wtmean(md_TL_w), weight(l)
gduplicates drop ind year, force
bys ind: egen markdown_wt_j = wtmean(markdown_wt_jt), weight(employ)
gduplicates drop ind, force
keep ind ind_name markdown_wt_j
sort ind
save "$processed_path/result_markdown_est_07_20_sector.dta", replace 

* calculate weighted mean of markdown by year
use "$processed_path/result_markdown_est_07_20.dta", clear
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
use "$processed_path/result_markdown_est_07_20.dta", clear
merge 1:1 sdid using "$processed_path/ctsd_location_qingdao_07_20.dta", keep(1 3) nogen
bys id: egen 市_mode = mode(市), minmode
replace 市 = 市_mode
drop 市_mode
keep if 市 == "青岛市"
keep id year md_TL md_TL_w 
save "$processed_path/result_markdown_est_qingdao_07_20.dta", replace