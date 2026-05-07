*==================================================*
* Prepare CIED data
* Date: October 6, 2025
* Author: Wu Chengjun, Central University of Finance and Economics
* OS: MacOS 26.1 Beta
* Version: Stata MP 18.0
*==================================================*

*!!! need change to your path !!!*
clear all
global proj_path = "/Users/pohwaran/Doctorate/Paper/Bridge"
global firm_data_path = "$proj_path/data"
global ctsd_data_path = "$firm_data_path/raw/ctsd"
global geo_data_path = "$proj_path/data/raw/geo"
global temp_path = "$firm_data_path/temp"
global processed_path = "$firm_data_path/processed"

*==================================================*
* Step 1: prepare data
* load and rename variables from 2007 to 2020
*==================================================*

* Define common variables to keep across all years
local keepvars "sdid year industry main_business_revenue main_business_cost output_vat employ_avg security_total wage_total depreciation fixed_assets_original_val ownership export"
local keepvars_noexport "sdid year industry main_business_revenue main_business_cost output_vat employ_avg security_total wage_total depreciation fixed_assets_original_val ownership"

forv year = 2007/2020 {
    use "$ctsd_data_path/ctsd_`year'.dta", clear
    
    qui if `year' == 2007 {
        ren (i17 p133 p136 v6 m309 m290 m288 b261 m298 r41 i1) ///
            (open_year main_business_revenue main_business_cost output_vat employ_avg security_total wage_total fixed_assets_original_val depreciation export ownership)
        keep `keepvars' open_year
        replace open_year = . if (open_year < 1500) | (open_year > 2020)
        replace open_year = 1978 if (open_year < 1978)
    }
    qui else if inrange(`year', 2008, 2011) {
        ren (hylb kysjn kysjy lrb_zysr lrb_zycb zzs_xxse qt_npjzgs qt_jtsbfl qt_jtgzjj zcb_mgdzcyz qt_jtzj zzckts_xse djzclx) ///
            (industry open_year open_month main_business_revenue main_business_cost output_vat employ_avg security_total wage_total fixed_assets_original_val depreciation export ownership)
        if `year' >= 2010 {
            replace employ_avg = cond(missing((qt_nczgs + qt_nmzgs) / 2), cond(missing(qt_nczgs), qt_nmzgs, qt_nczgs), (qt_nczgs + qt_nmzgs) / 2)
        }
        replace open_year = ustrregexra(open_year, "[^0-9]", "")
        replace open_month = ustrregexra(open_month, "[^0-9]", "")
        destring open_year, replace force
        destring open_month, replace force
        replace open_year = . if open_year < 1500 | open_year > 2020
        replace open_month = . if open_month < 1 | open_month > 12
        replace open_year = 1978 if open_year < 1978
        keep `keepvars' open_year open_month
    }
    qui else if `year' == 2012 {
        ren (f312 f315 f6 f189 f434 f355 f396 f45 登记注册类型) ///
            (main_business_revenue main_business_cost output_vat security_total wage_total fixed_assets_original_val depreciation export ownership)
        gen employ_avg = cond(missing((f435 + f436) / 2), cond(missing(f435), f436, f435), (f435 + f436) / 2)
        keep `keepvars'
    }
    qui else if `year' == 2013 {
        ren (f315 f318 f7 f195 f437 f358 f399 f46 登记注册类型) ///
            (main_business_revenue main_business_cost output_vat security_total wage_total fixed_assets_original_val depreciation export ownership)
        gen employ_avg = cond(missing((f438 + f439) / 2), cond(missing(f438), f439, f438), (f438 + f439) / 2)
        keep `keepvars'
    }
    qui else if inrange(`year', 2014, 2015) {
        ren (f270 f273 f7 f200 f390 f313 f354 f45 登记注册类型) ///
            (main_business_revenue main_business_cost output_vat security_total wage_total fixed_assets_original_val depreciation export ownership)
        gen employ_avg = cond(missing((f391 + f392) / 2), cond(missing(f391), f392, f391), (f391 + f392) / 2)
        keep `keepvars'
    }
    qui else if `year' == 2016 {
        ren (f273 f276 f7 f203 f398 f316 f357 f45 登记注册类型) ///
            (main_business_revenue main_business_cost output_vat security_total wage_total fixed_assets_original_val depreciation export ownership)
        gen employ_avg = cond(missing((f399 + f400) / 2), cond(missing(f399), f400, f399), (f399 + f400) / 2)
        keep `keepvars'
    }
    qui else if inrange(`year', 2017, 2020) {
        ren (开业成立时间 行业代码 主营业务收入 主营业务成本 增值税_总销项税额 本年已纳各类社会保障性基金 计提工资 年末固定资产原值 本年计提折旧 纳税人登记注册类型代码) ///
            (open_year industry main_business_revenue main_business_cost output_vat security_total wage_total fixed_assets_original_val depreciation ownership)
        replace open_year = . if open_year < 1500 | open_year > 2020
        replace open_year = 1978 if open_year < 1978
        gen employ_avg = cond(missing((年初职工数 + 年末职工数) / 2), cond(missing(年初职工数), 年末职工数, 年初职工数), (年初职工数 + 年末职工数) / 2)
        keep `keepvars_noexport' open_year
        replace open_year = . if open_year < 1500 | open_year > 2020
        replace open_year = 1978 if open_year < 1978
    }
    
    destring ownership, replace force
    save "$temp_path/ctsd_`year'.dta", replace
}

* append data from 2007 to 2020
use "$temp_path/ctsd_2007.dta", replace
forv year = 2008/2020 {
    qui append using "$temp_path/ctsd_`year'.dta"
    di "Appended year `year'"
}
forv year = 2007/2020 {
    erase "$temp_path/ctsd_`year'.dta"
}
merge 1:1 sdid using "$ctsd_data_path/ctsd_basic_info_07_20.dta", nogen
ren (法人代码 企业名称) (id firm_name)
order id firm_name year

save "$temp_path/ctsd_07_20.dta", replace

*==================================================*
* Step 2: define firm identifier, unify open year and sector, keep manufacturing firms
*==================================================*

use "$temp_path/ctsd_07_20.dta", clear

* drop observations with negative or missing key variables
foreach var of varlist fixed_assets_original_val - output_vat {
    replace `var' = 0 if (`var' < 0) | missing(`var')
}
egen all_missing = rowtotal(fixed_assets_original_val wage_total security_total depreciation employ_avg main_business_revenue main_business_cost export output_vat)
drop if all_missing == 0
drop all_missing

* Exact ID Matching: If the IDs of the two years are identical, they are considered to be the same company.
* Exact Name Matching: If the IDs are different (for example, a change of registration number), but the names are identical, they are considered to be the same company.
replace firm_name = subinstr(firm_name, "股份有限", "", .)
replace firm_name = subinstr(firm_name, "集团有限", "", .)
replace firm_name = subinstr(firm_name, "有限责任", "", .)
replace firm_name = subinstr(firm_name, "有限公司", "", .)
replace firm_name = subinstr(firm_name, "有限", "", .)
replace firm_name = subinstr(firm_name, "责任", "", .)
replace firm_name = subinstr(firm_name, "股份", "", .)
replace firm_name = subinstr(firm_name, "公司", "", .)
replace firm_name = subinstr(firm_name, "厂", "", .)
replace firm_name = subinstr(firm_name, " ", "", .)
replace firm_name = subinstr(firm_name, "(集团)", "", .)
replace firm_name = subinstr(firm_name, "（集团）", "", .)
replace firm_name = subinstr(firm_name, "（", "", .)
replace firm_name = subinstr(firm_name, "）", "", .)
replace firm_name = subinstr(firm_name, "(", "", .)
replace firm_name = subinstr(firm_name, ")", "", .)
replace firm_name = subinstr(firm_name, "回族自治区", "", .)
replace firm_name = subinstr(firm_name, "壮族自治区", "", .)
replace firm_name = subinstr(firm_name, "维吾尔自治区", "", .)
replace firm_name = subinstr(firm_name, "自治区", "", .)
replace firm_name = subinstr(firm_name, "省", "", .)
replace firm_name = subinstr(firm_name, "市", "", .)
replace firm_name = subinstr(firm_name, "区", "", .)
replace firm_name = subinstr(firm_name, "县", "", .)
replace firm_name = subinstr(firm_name, "-", "", .)

destring id, generate(id_num) force
replace 行政区划代码 = mod(行政区划代码, 1E6)
tostring 行政区划代码, generate(admin_code_str)
replace id = ustrregexs(1) if ustrregexm(id, "(.*)\(")
replace id = lower(id)
replace id = id + admin_code_str if id_num < 10000
bys firm_name: egen id_mod = mode(id), minmode
replace id = id_mod if firm_name != ""
drop id_mod id_num
drop if (id == "") & (firm_name == "")

bys id: egen firm_name_mod = mode(firm_name), minmode
replace firm_name = firm_name_mod if id != ""
drop firm_name_mod
replace id = firm_name if id == ""

* drop repeated observations in the same year
* keep observations with more available variables
gen zero_missing_count = 0
foreach var of varlist open_year - output_vat {
    capture confirm numeric variable `var'
    if !_rc {
        replace zero_missing_count = zero_missing_count + (`var' == 0 | `var' == .)
    }
}
bys id year: gen n = _N
bys id year: egen zero_missing_count_min = min(zero_missing_count)
drop if (zero_missing_count > zero_missing_count_min) & (n > 1)
gduplicates drop id year, force
drop n zero_missing_count zero_missing_count_min

* unify open year (already set 1978 as the minimum open year in Step 1)
bys id: egen open_year_mod = mode(open_year), minmode
replace open_year = open_year_mod
drop open_year_mod 

* unify sector and keep manufacturing firms (cic2 (2002) between 13 and 43)
* use the most frequent cic2, if tie, use the largest cic2
gen ind_code1 = substr(industry, 1, 1)
gen ind_code4 = substr(industry, 2, 4)
gen ind_code11 = ind_code4 if year >= 2011
merge m:1 ind_code11 using "$processed_path/ind_code_convert_2011_to_2002.dta", keep(1 3) nogen
replace ind_code4 = ind_code02 if ind_code02 != ""
gen ind_code2 = substr(ind_code4, 1, 2)
destring ind_code2, replace
bys id: egen ind_code2_mode = mode(ind_code2), minmode
replace ind_code2 = ind_code2_mode
replace ind_code1 = "C" if ind_code2  >= 13 & ind_code2 <= 42
bys id: egen ind_code1_mode = mode(ind_code1), minmode
replace ind_code1 = ind_code1_mode
drop ind_code2_mode ind_code1_mode ind_code11 ind_code02

* firm ownership
* merge with firm type (soe, poe, foe and other) 
* Source: https://www.stats.gov.cn/sj/tjbz/gjtjbz/202302/t20230213_1902746.html
replace ownership = 200 if ownership == 2
replace ownership = 110 if ownership == float(110.1)
merge m:1 ownership using "$processed_path/FirmType.dta", nogen
bys id: egen firm_type_mode = mode(firm_type), minmode
replace firm_type = firm_type_mode
drop firm_type_mode

sort id year
keep if ind_code2 >= 13 & ind_code2 <= 42
save "$temp_path/ctsd_07_20_Step1.dta", replace

*==================================================*
* Step 2: Gross output (revenue)
* revenue = operating revenue + sales tax
*==================================================*

use "$temp_path/ctsd_07_20_Step1.dta", clear

* use main business revenue + sales tax as revenue
replace main_business_revenue = . if (main_business_revenue <= 0)
replace output_vat = . if (output_vat <= 0)
gen revenue = main_business_revenue + output_vat

* merge output deflator (PPI) using output price index (2007 = 100)
ren ind_code2 cic2
merge m:1 cic2 year using "$processed_path/Index.dta", keep(3) nogen
gen revenue_real = revenue / output_index * 100

save "$temp_path/ctsd_07_20_Step2.dta", replace

*==================================================*
* Step 3: employ_avg and intermidiate inputs
* 1. employment can be directly observed
* 2. intermidiate inputs = operating cost - wage - depreciation
* where wage = wage + security payment
* 
* for security payment:
* 1. we sum up 6 types of social security payments
* 2. for missing security payment, replace it using the reported total social security payment if available
* 3. if both are missing, predict security payment using regression
*==================================================*

use "$temp_path/ctsd_07_20_Step2.dta", clear

* employment
replace employ_avg = . if (employ_avg <= 0)

* security payment
replace security_total = . if (security_total <= 0)
replace wage_total = . if (wage_total <= 0)
replace wage_total = security_total + wage_total // wage = wage + security payment

* intermidiate goods
replace main_business_cost = . if (main_business_cost <= 0)
replace depreciation = . if (depreciation < 0)
gen inter = main_business_cost - wage_total - depreciation
replace inter = . if (inter <= 0)

* intermidiate goods price index is the weighted average of PPI by I-O table
gen inter_real = inter / inter_index * 100

save "$temp_path/ctsd_07_20_Step3.dta", replace

*==================================================*
* Step 4: real capital stock (perpetual inventory method)
* k_t = (1-delta) * k_(t-1) + I_t, delta = 9%
* 1. predict nominal capital stock using first observed capital stock and average growth rate of nominal capital stock
* 2. calculate nominal investment using difference of nominal capital stock in two adjacent years, if missing, use predicted nominal capital stock
* 3. only keep firms entry after 1950
*==================================================*    

use "$temp_path/ctsd_07_20_Step3.dta", clear

* prepare data for capital stock calculation
keep id year open_year fixed_assets_original_val cic2
order id year open_year fixed_assets_original_val cic2
gen n = 1 // index for sample

* first observed nominal capital stock which is used to predict nominal capital stock
replace fixed_assets_original_val = . if fixed_assets_original_val <= 0
bys id: egen start_year_temp = min(year) if fixed_assets_original_val != .
bys id: egen start_year = min(start_year_temp)
drop start_year_temp
save "$temp_path/capital_temp.dta", replace

* fill the unbalanced panel from each firm's open year to 2020
* Note: for codes contain local macro, need run together. run line by line will cause error.
use "$temp_path/capital_temp.dta", clear
preserve
    keep id
    duplicates drop
    tempfile firms
    save `firms'
restore

clear
set obs 43
gen year = 1977 + _n   // 1978 to 2020
tempfile years
save `years'

use `firms', clear
cross using `years' // now we get a balanced panel from 1950 to 2020 for all firms

* merge with main data and drop if year < open_year and before 1985
sort id year
merge 1:1 id year using "$temp_path/capital_temp.dta", nogen
sort id open_year
bys id: replace open_year = open_year[1] if open_year == .
sort id year
drop if (year < open_year) & (open_year != .) // now we have a unbalanced panel from each firm's open year to 2020

* replace start_year and cic2 to filled missing values
ren start_year start_year_temp
bys id: egen start_year = min(start_year_temp)
drop start_year_temp
drop if start_year == .
replace open_year = start_year if open_year == .

ren cic2 cic2_temp
bys id: egen cic2 = min(cic2_temp)
drop cic2_temp

* average growth rate of nominal capital stock
sort id year
bys id: egen end_year_temp = max(year) if fixed_assets_original_val != .
bys id: egen end_year = min(end_year_temp)
drop end_year_temp

gen capital_nomi_end_temp = fixed_assets_original_val if year == end_year
bys id: egen capital_nomi_end = min(capital_nomi_end_temp)
drop capital_nomi_end_temp

gen capital_nomi_start_temp = fixed_assets_original_val if year == start_year
bys id: egen capital_nomi_start = min(capital_nomi_start_temp)
drop capital_nomi_start_temp

gen capital_nomi_growth_av = (capital_nomi_end / capital_nomi_start) ^ (1 / (end_year - start_year)) - 1
drop if capital_nomi_growth_av == 0 // these firms don't have sufficient information to predict capital stock

* predict nominal capital stock using first observed capital stock and average growth rate of nominal capital stock
gen capital_nomi_pred = capital_nomi_start * ((1 + capital_nomi_growth_av) ^ (year - start_year)) // predicted nominal capital stock for all years
drop if capital_nomi_pred == . // these firms don't have sufficient information to predict capital stock

* calculate nominal investment using difference of nominal capital stock in two adjacent years. for missing investment, use predicted nominal capital stock
sort id year
bys id: gen invest_nomi = (fixed_assets_original_val[_n] - fixed_assets_original_val[_n - 1])
bys id: gen invest_nomi_pred = (capital_nomi_pred[_n] - capital_nomi_pred[_n - 1])
replace invest_nomi = invest_nomi_pred if invest_nomi == .

bys id: replace invest_nomi = fixed_assets_original_val if _n == 1
bys id: replace invest_nomi = capital_nomi_pred if _n == 1 & invest_nomi == . // k_0, from predicted capital stock

* merge with capital price index
* because capital price index can only be observed from 1991 to 2019, we replace missing index with PPI from 1978 to 1990, 2020 and RPI from 1950 to 1977
merge m:1 year using "$processed_path/CaptialIndex.dta", nogen keep(3)
replace invest_nomi = invest_nomi / invest_index * 100

* calculate real capital stock using perpetual inventory method
* k_t = (1-delta) * k_(t-1) + I_t, delta = 9%
sort id year
gen capital_real = .
bys id: replace capital_real = invest_nomi if _n == 1
bys id: replace capital_real = (1 - 0.09) * capital_real[_n - 1] + invest_nomi if _n > 1

keep if n == 1 // keep observations the original sample
keep id year capital_real
replace capital_real = . if (capital_real <= 0)
sort id year
merge 1:1 id year using "$temp_path/ctsd_07_20_Step3.dta", keep(1 3) nogen

save "$temp_path/ctsd_07_20_Step4.dta", replace

*==================================================*
* Step 5: final merge and save
*==================================================*

use "$temp_path/ctsd_07_20_Step4.dta", clear

* drop some sample with unreasonable values
drop if missing(revenue_real, capital_real, employ_avg, inter_real)
drop if revenue <= capital_real
drop if revenue <= wage_total
drop if revenue <= inter
drop if revenue < 30
drop if inter < 30

* drop if firms only appear once
bys id: gen n = _N
drop if n == 1
order id year

* log variables
gen q = ln(revenue_real)
gen k = ln(capital_real)
gen l = ln(employ_avg)
gen m = ln(inter_real)

* keep variables we need
* output_index is needed because we will calculate alpha_m using corrected revenue (see line 151-153 in MdEstTl.do)
gen alpha_l = wage_total / revenue
keep id year sdid firm_name q k l m alpha_l inter output_index cic2 wage_total

label var wage_total "工资+社保总额"
label var inter "中间投入总额"
label var output_index "产出价格指数"
label var cic2 "行业代码(2002版)"
label var alpha_l "劳动收入份额"

* merge sectors with little firms
gen ind = 1 if cic2 >= 13 & cic2 <= 16 // 食品饮料
replace ind = 2 if cic2 >= 17 & cic2 <= 19 // 纺织服装
replace ind = 3 if cic2 >= 20 & cic2 <= 21 // 木材家具
replace ind = 4 if cic2 >= 22 & cic2 <= 23 // 造纸印刷
replace ind = 5 if cic2 == 24 // 文教用品
replace ind = 6 if cic2 == 25 // 石油加工
replace ind = 7 if cic2 >= 26 & cic2 <= 30 // 化学医药
replace ind = 8 if cic2 == 31 // 非金属制品
replace ind = 9 if cic2 >= 32 & cic2 <= 34 // 金属制造
replace ind = 10 if cic2 >= 35 & cic2 <= 36 // 机械设备
replace ind = 11 if cic2 == 37 // 运输设备
replace ind = 12 if cic2 >= 39 & cic2 <= 41 // 电气设备
drop if cic2 == 42 // other manufacturing sectors

sort ind year
merge m:1 ind year using "$processed_path/Employ.dta", nogen
ren emp_ employ

egen firm_id = group(id)
xtset firm_id year, yearly
gen year_index = d.year
gen year_index2 = F.year_index
drop if year_index == . & year_index2 == .
gduplicates drop id year, force
sort id year

save "$processed_path/markdown_est_07_20.dta", replace

* calculate number of observations by sector
use "$processed_path/markdown_est_07_20.dta", clear
bys ind: gen n_firm = _N
gduplicates drop ind, force
keep ind ind_name n_firm

* calculate number of observations by sector in qingdao
merge 1:1 id year using "$processed_path/markdown_est_07_20.dta", keep(3) nogen

bys ind: gen n_firm = _N
gduplicates drop ind, force
keep ind ind_name n_firm

*==================================================*
* Step 6: Sample for regression and model
*==================================================*

use "$temp_path/ctsd_07_20_Step3.dta", clear

* keep firms existing in 2010
gen if_year_2010 = (year == 2010)
bys id: egen exist_2010 = max(if_year_2010)
keep if exist_2010 == 1
gduplicates drop id year, force
save "$temp_path/ctsd_07_20_2010exist.dta", replace

* keep firms in Qingdao
use "$processed_path/ctsd_location_qingdao_07_20.dta", clear
keep if 市 == "青岛市"
ren (经度 纬度 县) (longitude latitude county)
keep sdid longitude latitude
tempfile location_qingdao
save `location_qingdao', replace

* merge location info
use "$temp_path/ctsd_07_20_2010exist.dta", clear
merge 1:1 sdid using `location_qingdao', keep(1 3) nogen
replace longitude = . if year > 2010
replace latitude = . if year > 2010
bys id: egen location_year = max(year) if !missing(longitude, latitude)
replace longitude = . if year != location_year
replace latitude = . if year != location_year
bys id: egen longitude_mode = mode(longitude), minmode
replace longitude = longitude_mode
bys id: egen latitude_mode = mode(latitude), minmode
replace latitude = latitude_mode
drop longitude_mode latitude_mode location_year
drop if missing(longitude, latitude)

* keep key variables
drop revenue
ren main_business_revenue revenue
gen wage = wage_total / employ_avg
gen wage_inital_temp_2010 = wage if (year == 2010)
bys id: egen wage_inital_2010 = min(wage_inital_temp_2010)

gen wage_inital_temp_2011 = wage if (year == 2011)
bys id: egen wage_inital_2011 = min(wage_inital_temp_2011)

drop wage_inital_temp_2010 wage_inital_temp_2011

ren employ_avg employ
gen age = year - open_year + 1

gen export_intensity = export / revenue
bys id: egen export_total = sum(export)
gen export_bool = (export_total > 0)

ren cic2 ind_code2

keep id year revenue employ longitude latitude age firm_type export export_bool export_intensity wage_total wage_inital_2010 wage_inital_2011 ind_code2 sdid
drop if missing(longitude, latitude)

* regression analsys, employ_avg, longitude, latitude
save "$temp_path/ctsd_qingdao_07_20.dta", replace

*==================================================*
* contain firms with wage but no security
*==================================================*

use "$temp_path/ctsd_07_20_Step2.dta", clear

* employment
replace employ_avg = . if (employ_avg <= 0)

* intermidiate goods
replace main_business_cost = . if (main_business_cost <= 0)
replace depreciation = . if (depreciation < 0)
gen inter = main_business_cost - wage_total - depreciation
replace inter = . if (inter <= 0)

* intermidiate goods price index is the weighted average of PPI by I-O table
gen inter_real = inter / inter_index * 100

* keep firms existing in 2010
gen if_year_2010 = (year == 2010)
bys id: egen exist_2010 = max(if_year_2010)
keep if exist_2010 == 1
gduplicates drop id year, force
save "$temp_path/ctsd_07_20_2010exist_nosec.dta", replace

* keep firms in Qingdao
use "$processed_path/ctsd_location_qingdao_07_20.dta", clear
keep if 市 == "青岛市"
ren (经度 纬度 县) (longitude latitude county)
keep sdid longitude latitude
save `location_qingdao', replace

* merge location info
use "$temp_path/ctsd_07_20_2010exist_nosec.dta", clear
merge 1:1 sdid using `location_qingdao', keep(1 3) nogen
replace longitude = . if year > 2010
replace latitude = . if year > 2010
bys id: egen location_year = max(year) if !missing(longitude, latitude)
replace longitude = . if year != location_year
replace latitude = . if year != location_year
bys id: egen longitude_mode = mode(longitude), minmode
replace longitude = longitude_mode
bys id: egen latitude_mode = mode(latitude), minmode
replace latitude = latitude_mode
drop longitude_mode latitude_mode location_year
drop if missing(longitude, latitude)

* keep key variables
drop revenue
ren main_business_revenue revenue
gen wage = wage_total / employ_avg
gen wage_inital_temp_2010 = wage if (year == 2010)
bys id: egen wage_inital_2010 = min(wage_inital_temp_2010)
drop wage_inital_temp_2010

gen wage_inital_temp_2011 = wage if (year == 2011)
bys id: egen wage_inital_2011 = min(wage_inital_temp_2011)
drop wage_inital_temp_2011

ren employ_avg employ
gen age = year - open_year + 1

bys id: egen export_total = sum(export)
gen export_bool = (export_total > 0)

ren cic2 ind_code2

keep id year revenue employ longitude latitude age export_intensity firm_type wage_total wage_inital_2010 wage_inital_2011 ind_code2 sdid export export_bool
drop if missing(longitude, latitude)

* regression analsys, employ_avg, longitude, latitude
save "$temp_path/ctsd_qingdao_07_20_nosec.dta", replace