*==================================================*
* Plot Heterogeneous trend test with wage difference
*==================================================*


preserve


mat T=J(4,2,.)
forvalues i=2010/2013{
	
	local k=`i'-2009
	mat T[`k',1]=_b[i1.big_ma#i`i'.year]
	mat T[`k',2]=_b[i1.big_ma#i`i'.year]+_b[i1.high_ln_w0_diff#i1.big_ma#i`i'.year]
	
	  
}

clear
svmat T
gen year=_n+2009
set obs 5
replace year=2009 if _n==5
forvalues i=1/2{
    replace T`i'=0 if year==2009
}
sort year
scatter T1 year, c(l) || ///
scatter T2 year, c(l) ///
legend(label(1 low_lnw_diff_2010) label(2 high_lnw_diff_2010) )

restore


