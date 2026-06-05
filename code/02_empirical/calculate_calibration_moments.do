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

keep if e(sample)==1
duplicates drop id, force
keep id
save "$model_processed_path/firm_two_year_reg.dta", replace