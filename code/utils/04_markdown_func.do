******************************************************************
* modified version of the code provided by Yeh
******************************************************************

/*
	
	Program: 
	
	- This program performs the GMM-IV estimation procedure as described on p. 2446 - 2449
	- The code:
		- estimates production function coefficients in the spirit of 
		  Ackerberg, Caves and Frazer (WP2006) - ACF2006
		
	- Steps in program:
		- Construct variables for polynomial approximation of control function (based on materials as in Levinsohn and Petrin (ReStud2003)
			- No identification of any coefficients in first stage though
		- Estimate expected output net of productivity based on polynomial approximation
		- GMM-IV procedure to estimate production function coefficients
			- COBB-DOUGLAS PRODUCTION:
				- Markup is identified from coefficient on material inputs
				- Correct gross output share with residual from first stage: measurement error correction 
				- Variation across plants within an industry comes from "corrected" material cost share alone
			- TRANSLOG PRODUCTION:
				1. Markup is identified from coefficient on material inputs
				2. Correct gross output share with residual from first stage: measurement error correction 
				3. Variation across plants within an industry comes from two sources: "corrected" material cost share and varying output-materials elasticity
	
*/

clear mata

*** POLYNOMIAL APPROXIMATION OF REAL OUTPUT ***

qui{
	
	* Create variables for polynominal approximation in first stage 
	* - Total of 4*M + M + 6*M*N generated variables
	
	* Choose degree of polynomial
	local M=3
	local N=3
	
	* Singular components
	* - Total of 4*M terms

	forvalues i=1/`M' {
	
		gen l`i'=l^(`i')
		gen m`i'=m^(`i')
		gen k`i'=k^(`i')

			* Interaction terms
			*  - Interaction of two variables only
			*  - Total of 6*M*N interaction terms
			forvalues j=1/`N' {
		
			gen k`i'l`j'=k^(`i')*l^(`j')
			gen k`i'm`j'=k^(`i')*m^(`j')
			gen l`i'm`j'=l^(`i')*m^(`j')
		
			}
	
	}
	
	* Interaction terms
	*  - Interaction of three variables (same degree only)
	*  - Total of M terms
	gen lkm=l*k*m
	forvalues i=2(1)`M'  {
	
		gen l`i'k`i'm`i'=l`i'*k`i'*m`i'
		
	}
	
}
	
*** INITIALIZATION ***

* GMM-IV procedures require non-trivial initialization of coefficients
* - Use biased OLS coefficients as starting values

* Translog
xi: reg q k l m k1l1 k1m1 l1m1 k2 l2 m2 i.year

mat beta_OLS_TL_10 = J(1,10,0)

mat beta_OLS_TL_10[1,1]=_b[_cons]
mat beta_OLS_TL_10[1,2]=_b[k]
mat beta_OLS_TL_10[1,3]=_b[l]
mat beta_OLS_TL_10[1,4]=_b[m]
mat beta_OLS_TL_10[1,5]=_b[k1l1]
mat beta_OLS_TL_10[1,6]=_b[k1m1]
mat beta_OLS_TL_10[1,7]=_b[l1m1]
mat beta_OLS_TL_10[1,8]=_b[k2]
mat beta_OLS_TL_10[1,9]=_b[l2]
mat beta_OLS_TL_10[1,10]=_b[m2]

mat list beta_OLS_TL_10

local tl1 = _b[_cons]
local tl2 = _b[k]
local tl3 = _b[l]
local tl4 = _b[m]
local tl5 = _b[k1l1]
local tl6 = _b[k1m1]
local tl7 = _b[l1m1]
local tl8 = _b[k2]
local tl9 = _b[l2]
local tl10 = _b[m2]

************************************************************	
*** GMM-IV ESTIMATION OF PRODUCTION FUNCTION COEFFCIENTS ***
************************************************************

*** FIRST STAGE
* Note: - Include year fixed effects to exclude variation due to trends over time
*		- These trends are pushing up gross output over time which can bias the coefficients on all inputs
*	- Note that first stage is similar for both Cobb-Douglas and translog production function estimation

xi: reg q l* m* k* i.year

predict phi
predict epsilon, res
* Note: measurement error is denoted in natural log levels
label var phi "Expected output net of productivity - phi_it"
label var epsilon "Measurement error (first stage)"

gen phi_lag=L.phi
gen k_lag=L.k
gen l_lag=L.l
gen m_lag=L.m
gen k_lag2=k_lag^2
gen l_lag2=l_lag^2
gen m_lag2=m_lag^2

gen kl_lag = k*l_lag
gen km_lag = k*m_lag

gen k_lagl_lag = k_lag*l_lag
gen k_lagm_lag = k_lag*m_lag
gen l_lagm_lag = l_lag*m_lag

* Compute corrected share
* Note: - correct value of gross output from measurement error which is estimated in first stage
*		- material expenditures/value of gross output ratio is based on "corrected" value of gross output

gen lq_c = q - epsilon
gen q_c = exp(lq_c)
gen alpha_m = inter / (q_c * output_index / 100)

sort firm_id year
gen const=1

* Cleaning of data:
*	- Get rid of non-useful/missing variables
drop if k==.
drop if k_lag==.
drop if l==.
drop if l_lag==.
drop if m==.
drop if m_lag==.
drop if phi==.
drop if phi_lag==.

qui count
disp "Sample size equals: `r(N)' observations"

*** SECOND STAGE

******************** MATA START ********************


*** MATA PROGRAMS FOR GMM-IV MINIMIZATION PROBLEM ***

	mata:
	
	beta_OLS_TL=st_matrix("beta_OLS_TL_10")

	void GMM_DLW_TL(todo,betas,crit,g,H)
	{

		real matrix PHI, PHI_LAG, Z, X, X_lag, Y, C
		real matrix OMEGA, OMEGA_lag, OMEGA_lag2, OMEGA_lag3
		real matrix OMEGA_lag_pol, g_b, XI

		PHI=st_data(.,("phi"))
		PHI_LAG=st_data(.,("phi_lag"))
		Z=st_data(.,("const","k","l_lag","m_lag","kl_lag","km_lag","l_lagm_lag","k2","l_lag2","m_lag2"))
		X=st_data(.,("const","k","l","m","k1l1","k1m1","l1m1","k2","l2","m2"))
		X_lag=st_data(.,("const","k_lag","l_lag","m_lag","k_lagl_lag","k_lagm_lag","l_lagm_lag","k_lag2","l_lag2","m_lag2"))
		Y=st_data(.,("q"))
		C=st_data(.,("const"))

		OMEGA=PHI-X*betas'
		OMEGA_lag=PHI_LAG-X_lag*betas'
		OMEGA_lag2=OMEGA_lag:*OMEGA_lag
		OMEGA_lag3=OMEGA_lag2:*OMEGA_lag
		OMEGA_lag_pol=(C,OMEGA_lag,OMEGA_lag2,OMEGA_lag3)
		g_b = invsym(OMEGA_lag_pol'OMEGA_lag_pol)*OMEGA_lag_pol'OMEGA
		XI=OMEGA-OMEGA_lag_pol*g_b
		crit=(Z'XI)'(Z'XI)
		
	}

	void DLW_TRANSLOG()
	{
		transmorphic S
		real rowvector p

		S=optimize_init()
		optimize_init_evaluator(S, &GMM_DLW_TL())
		optimize_init_evaluatortype(S,"d0")
		optimize_init_technique(S, "nm")
		optimize_init_nmsimplexdeltas(S, 0.1)
		optimize_init_which(S,"min")
		optimize_init_params(S,(`tl1',`tl2',`tl3',`tl4',`tl5',`tl6',`tl7',`tl8',`tl9',`tl10'))
		p=optimize(S)
		p
		st_matrix("beta_dlwtl",p)
	}

	end


******************** MATA END ********************

* GMM-IV: Translog
cap program drop dlw_translog
program dlw_translog, rclass
preserve 
sort firm_id year
mata DLW_TRANSLOG()
end

*** Estimates with polynomial approximation for productivity ("control function")

* Translog
* Note: "beta = (const, k, l, m, e, kl, km, ke, lm, le, me, k2, l2, m2, e2)"
dlw_translog

qui  {

	gen beta_c1=beta_dlwtl[1,1]
	gen beta_k1=beta_dlwtl[1,2]
	gen beta_l1=beta_dlwtl[1,3]
	gen beta_m1=beta_dlwtl[1,4]
	
	gen beta_k1l1=beta_dlwtl[1,5]
	gen beta_k1m1=beta_dlwtl[1,6]
	gen beta_l1m1=beta_dlwtl[1,7]
	
		gen beta_k2=beta_dlwtl[1,8]
		gen beta_l2=beta_dlwtl[1,9]
		gen beta_m2=beta_dlwtl[1,10]

		* Translog productivity implied by the DLW/ACF control function.
		gen omega_TL = phi - (beta_c1 + beta_k1*k + beta_l1*l + beta_m1*m ///
			+ beta_k1l1*k1l1 + beta_k1m1*k1m1 + beta_l1m1*l1m1 ///
			+ beta_k2*k2 + beta_l2*l2 + beta_m2*m2)
		gen tfp_TL = exp(omega_TL)
		label var omega_TL "Log productivity from translog production function"
		label var tfp_TL "TFP from translog production function"
	
		* "betam_tl" contains the output elasticity with respect to material inputs ("theta^m = beta_m + beta_km*k + beta_lm* l + beta_me*e + 2*beta_m2*m")
		gen theta_m_tl = beta_m1 + beta_k1m1*k + beta_l1m1*l + 2*beta_m2*m
		gen theta_l_tl = beta_l1 + beta_k1l1*k + beta_l1m1*m + 2*beta_l2*l
	
	* dLW markup
	gen mu_DLW_TL = theta_m_tl / alpha_m

}

