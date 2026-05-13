*==================================================*
* Prepare data for China Industrial Enterprise Database (CIED)
* Date: December 9, 2025
* Author: Wu Chengjun, Central University of Finance and Economics
* OS: MacOS 26.1
* Version: Stata MP 18.0
*==================================================*

* set path
clear all
global proj_path = "/Users/pohwaran/Doctorate/Paper/Bridge"
do "$proj_path/code/00_setup/data_paths.do"

*==================================================*
* step 1: Model Sample Preparation (CTSD + CIED)
*==================================================*

* prepare CTSD data for year 2010
use "$regression_temp_path/match_cied_ctsd_07_14.dta", clear
merge m:1 sdid using "$ctsd_processed_path/ctsd_qingdao_07_20.dta", keep(2 3) nogen
gen index = 1 if !missing(group) // mark matched firms in CTSD data
keep if year == 2010 | year == 2011

* append with CIED data for year 2011
append using "$cied_processed_path/cied_qingdao_11.dta"
drop if missing(wage_inital_2010) & year == 2010
drop if missing(wage_inital_2011) & year == 2011 // drop if missing wage info for both years

bys group year: gen n = _N
drop if (n == 2) & missing(index)
drop n index

tostring group, replace
replace id = group if id == ""

* keep relevant variables for market accessibility calculation
* keep id z longitude latitude
keep id year longitude latitude wage_inital_2010 wage_inital_2011 employ

* merge town info based on firm coordinates
geoinpoly latitude longitude using "$geo_raw_path/shapefiles/town_coord.dta"
ren _ID ID
merge m:1 ID using "$geo_raw_path/shapefiles/town_db.dta", nogen
ren (乡 县) (town county)
drop ID 省 市 treat geom

* export final firm data
sort year id
save "$model_temp_path/firm_qingdao_model.dta", replace

preserve
keep if year == 2010
drop wage_inital_2011
ren wage_inital_2010 wage_inital
save "$model_temp_path/firm_qingdao_model_10.dta", replace
restore

preserve
keep if year == 2011
drop wage_inital_2010
ren wage_inital_2011 wage_inital
save "$model_temp_path/firm_qingdao_model_11.dta", replace
restore

gduplicates drop id, force
keep id longitude latitude
export excel "$geo_temp_path/firm_qingdao_model.xlsx", firstrow(variables) replace

*==================================================*
* Step 2: Regression Sample Preparation (CTSD + CIED)
*==================================================*

* prepare CTSD data for year 2007-2020
use "$regression_temp_path/match_cied_ctsd_07_14.dta", clear
merge m:1 sdid using "$ctsd_temp_path/ctsd_qingdao_07_20.dta", keep(2 3) nogen
drop if missing(export_bool, firm_type, ind_code2)
gen index = 1 if !missing(group)

* prepare CIED data for year 2007-2014
append using "$cied_temp_path/cied_qingdao_07_14.dta"
bys group year: gen n = _N
drop if (n == 2) & missing(index) // drop repeated observations
drop n index
drop if missing(export_bool, firm_type, ind_code2)

tostring group, replace
replace id = group if id == ""
duplicates drop id year, force

* export final firm data
order id year
bys id: gen n = _N
drop if n == 1

* merge town info based on firm coordinates
drop county town
geoinpoly latitude longitude using "$geo_raw_path/shapefiles/town_coord.dta"
ren _ID ID
merge m:1 ID using "$geo_raw_path/shapefiles/town_db.dta", nogen
ren (乡 县) (town county)
drop ID 省 市 treat geom

sort id year
drop 年份 - 企业名称 n
save "$regression_temp_path/firm_qingdao_reg.dta", replace

duplicates drop id, force
sort id
keep id longitude latitude county
export excel using "$geo_temp_path/firm_list_qingdao_reg.xlsx", firstrow(variables) replace

*==================================================*
* Step 3: Government_QingDao.dta Preparation
*==================================================*

use "$geo_raw_path/town_gov_coords_2023.dta", clear
keep if city == "青岛市"
keep town 纬度 经度
ren (纬度 经度) (latitude longitude)
sort town
export excel using "$geo_raw_path/govern_qingdao.xlsx", firstrow(variables) replace

*==================================================*
* Step 4: Population Data Preparation
*==================================================*

import excel "$geo_processed_path/pop_census_qingdao_2010.xlsx", sheet("Sheet2") firstrow clear
bys town: egen pop = sum(I) // 15-64 population
duplicates drop town, force

* split the population of Zhushan subdistrict equally to Yinzhu and Tieshan subdistricts
summarize pop if town == "珠山街道"
local zhushan_pop = r(sum)
local split_pop = `zhushan_pop' / 2
drop if town == "珠山街道"
replace pop = pop + `split_pop' if town == "隐珠街道"
replace pop = pop + `split_pop' if town == "铁山街道"
keep town pop

save "$geo_processed_path/pop_census_qingdao_2010.dta", replace

*==================================================*
* Step 5: Market Access Calculation _python
*==================================================*

python set exec "/opt/miniconda3/envs/bridge/bin/python"
python script "$proj_path/code/00_setup/install_python_pkgs.py"
python script "$proj_path/code/01_data_prep/03_market_access_func.py"

*==================================================*
* Step 6: Calculate market accessibility using QGIS
*==================================================*

foreach suffix in "_model" "_reg" {
    import delimited "$geo_processed_path/commute_time_qingdao`suffix'.csv", clear
    keep firm_id town_id travel_time_min
    ren (firm_id town_id travel_time_min) (id town dzj_prime)
    tempfile market_access
    save `market_access', replace
    import delimited "$geo_processed_path/commute_time_no_bridge_qingdao`suffix'.csv", clear
    keep firm_id town_id travel_time_min
    ren (firm_id town_id travel_time_min) (id town dzj)
    merge 1:1 id town using `market_access', nogen

    * clean town names to match with population data
    replace town = "九水路街道" if town == "世园街道"
    replace town = "九水路街道" if town == "九水街道"
    replace town = "" if town == "双山街道"
    replace town = "大信镇" if town == "大信街道"
    replace town = "洛阳路街道" if town == "开平路街道"
    replace town = "张家楼镇" if town == "张家楼街道"
    replace town = "永安路街道" if town == "沧口街道"
    replace town = "" if town == "湖岛街道"
    replace town = "灵山镇" if town == "灵山街道"
    replace town = "王台镇" if town == "王台街道"
    replace town = "" if town == "胶南街道"
    replace town = "胶莱镇" if town == "胶莱街道"
    replace town = "胶西镇" if town == "胶西街道"
    replace town = "蓝村镇" if town == "蓝村街道"
    replace town = "藏南镇" if town == "藏马镇"
    replace town = "中韩街道" if town == "金家岭街道"
    replace town = "" if town == "金湖路街道"
    replace town = "平度外向型工业加工区" if town == "平度经济开发区"
    merge m:1 town using "$geo_processed_path/pop_census_qingdao_10", keep(3) nogen
    drop if missing(dzj, dzj_prime)

    gduplicates drop id town, force
    bys id: gen n = _N
    keep if n == 128
    drop n

    if "`suffix'" == "_model" {
        preserve
        merge m:1 id using "$model_temp_path/firm_qingdao_model_10.dta", keep(3) nogen
        sort year id town
        save "$model_processed_path/firm_qingdao_model_10.dta", replace
        restore

        merge m:1 id using "$model_temp_path/firm_qingdao_model_11.dta", keep(3) nogen
        sort year id town
        save "$model_processed_path/firm_qingdao_model_11.dta", replace
    }
    else {
        bys id: egen dma = wtmean(dzj - dzj_prime), weight(pop)
        gen dln_ma = ln(dma)
        * replace dln_ma = -10 if missing(dln_ma)
        gduplicates drop id, force
        tempfile market_access
        save `market_access', replace

        use "$regression_temp_path/firm_qingdao_reg.dta", clear
        merge m:1 id using `market_access', keep(3) nogen

        save "$regression_processed_path/regression_qingdao_07_20.dta", replace
    }
}
