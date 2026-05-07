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
* Step 1: append raw CIED data and keep relevant variables
* Source: https://mp.weixin.qq.com/s/xI5q0Vov4OnxgLdLJtqxyg?scene=1&click_id=30
*==================================================*

* Load CIED raw data for Qingdao (1998-2014) and filter by city
forv year = 1998/2014{
    use "$cied_raw_path/cied_`year'.dta", clear
    keep if 市 == "青岛市"
    save "$cied_temp_path/cied_qingdao_`year'_temp.dta", replace
}

* Combine all years into single dataset
use "$cied_temp_path/cied_qingdao_1998_temp.dta", clear
forv year = 1999/2014{
    append using "$cied_temp_path/cied_qingdao_`year'_temp.dta"
}
ren 年份 year
save "$cied_temp_path/cied_qingdao_98_14.dta", replace

* Remove intermediate year files
forv year = 1998/2014{
    erase "$cied_temp_path/cied_qingdao_`year'_temp.dta"
}

* Standardize industry codes, firm founding year, key variables, and firm type
* keep only firms in manufacturing sector (code "C")
use "$cied_temp_path/cied_qingdao_98_14.dta", clear
bys group: egen ind_code1 = mode(行业门类代码), minmode
keep if ind_code1 == "C"

* funding year
replace 开业成立时间年 = . if 开业成立时间年 <= 1000 | 开业成立时间年 > year
bys group: egen open_year = mode(开业成立时间年), minmode
gen age = year - open_year + 1

* export 
ren 其中出口交货值千元 export
replace export = 0 if missing(export)
bys group: egen export_total = sum(export)
gen export_bool = (export_total > 0)

* firm ownership (3 kinds)
gen firm_type = "soe" if (登记注册类型 == "国有与集体联营企业") | (登记注册类型 == "国有企业") | (登记注册类型 == "国有独资公司") | (登记注册类型 == "国有联营企业")
replace firm_type = "foe" if (登记注册类型 == "中外合作经营企业") | (登记注册类型 == "中外合资经营企业") | (登记注册类型 == "合作经营企业（港或澳台资）") | (登记注册类型 == "合资经营企业（港或澳台资）") | (登记注册类型 == "外商投资股份有限公司") | (登记注册类型 == "外资（独资）企业") | (登记注册类型 == "港澳台商投资股份有限公司") | (登记注册类型 == "港澳台独资企业")
replace firm_type = "poe" if firm_type == "私营合伙企业" | (登记注册类型 == "私营有限责任公司") | (登记注册类型 == "私营独资企业") | (登记注册类型 == "私营股份有限公司")
replace firm_type =  "other" if missing(firm_type)
bys group: egen firm_type_freq = mode(firm_type), minmode
replace firm_type = firm_type_freq

* key variables
ren 工业销售产值_当年价格千元 revenue
gen export_intensity = export / revenue
ren 应付工资薪酬总额千元 wage_total
ren 年末从业人员合计人 employ
replace employ = 全部从业人员年平均人数人 if missing(employ) | (employ == 0)
ren 经度 longitude
ren 纬度 latitude
replace wage_total = . if (wage_total <= 0)
replace employ = . if (employ <= 0)
gen wage = wage_total / employ
gen wage_inital_temp = wage if (year == 2011)
bys group: egen wage_inital_2011 = min(wage_inital_temp)

* keep only firm existing in 2011
gen if_year_2011 = (year == 2011)
bys group: egen exist_2011 = max(if_year_2011)
keep if exist_2011 == 1

* industry code 2-digit
gen ind_code2_temp = 行业大类代码 if year == 2011
bys group: egen ind_code2 = mode(行业大类代码), minmode

* keep relevant variables
ren 县 county
keep group year wage_total wage_inital_2011 employ revenue longitude latitude export export_bool export_intensity firm_type age ind_code2 county
drop if missing(longitude, latitude)
keep if year >= 2007

* merge town info based on firm coordinates
geoinpoly latitude longitude using "$geo_raw_path/shapefiles/town_coord.dta"
ren _ID ID
merge m:1 ID using "$geo_raw_path/shapefiles/town_db.dta", nogen
ren 乡 town
drop ID - 市 treat 县

sort group year
save "$cied_temp_path/cied_qingdao_07_14.dta", replace

keep if year == 2011
save "$ctsd_processed_path/firm_cied_qingdao.dta", replace
