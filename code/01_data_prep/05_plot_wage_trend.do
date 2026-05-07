****将ln_w0_diff分位10组
preserve


mat T=J(4,6,.)
forvalues i=2010/2013{
    forvalues j=1/6{
	   local k=`i'-2009
	   mat T[`k',`j']=_b[i1.big_ma#i`i'.year]+_b[i1.big_ma#i`i'.year#c.ln_w0_diff]*ln_w0_diff_`j'
	}
}

clear
svmat T
gen year=_n+2009
set obs 5
replace year=2009 if _n==5
forvalues i=1/6{
    replace T`i'=0 if year==2009
}
sort year
scatter T1 year, c(l) || ///
scatter T2 year, c(l) || ///
scatter T3 year, c(l) || ///
scatter T4 year, c(l) || ///
scatter T5 year, c(l) || ///
scatter T6 year, c(l) ///
legend(label(1 ln_w0_diff_p5) label(2 ln_w0_diff_p25) label(3 ln_w0_diff_p50) label(4 ln_w0_diff_p75) label(5 ln_w0_diff_p90) label(6 ln_w0_diff_p95))

restore
