
 * This do file: computes the variance-covariance matrix of the estimated road technology 
 * 		corresponding to column 3 in Table III (OLS hourly density spec)
 * This version: Oct 2020
 * Author: Gabriel Kreindler
 
clear all
pause off
set more off

* Path based on global set in profile.do
	local path "$congestion_pricing_path"
	cd "`path'"

***************	
** Main relationship: density on instant delay
***************	

*** Hour level data
	use "data/coded_road_tech/volumes_h_level", clear

	merge 1:1 hour using "data/coded_road_tech/density_h"
	assert _m==3
	drop _m

	merge 1:1 hour using "data/coded_road_tech/google_maps_h_level"
	keep if _m==3
	drop _m

	tsset hour

	*** do NOT normalize density -- already normalized
	sum density_norm, d 

*** Linear fit
	newey2 delay_gm density_norm, lag(3)
	
	*** save coefficients and variance-covariance matrix
	matrix list e(b)
	svmat e(b), names("temp")
	rename temp1 c_density
	rename temp2 c_const

	matrix list e(V)
	svmat e(V), names("temp")
	rename temp1 vcov_density
	rename temp2 vcov_const

	* save
	preserve
		keep c_* vcov_* 
		keep if _n <= 2
		export delimited using "data/coded_model/road_tech/vcov_density_linear.csv", replace
	restore

	drop c_* vcov_*

*** Power -- Non-linear least square
	nl (delay_gm = {l0=2.0} + {l1=1.5} * (density_norm ^ {gamma=1.0})), vce(hac nwest 3)

	*** save coefficients and variance-covariance matrix
	matrix list e(b)
	svmat e(b), names("temp")
	rename temp1 c_const
	rename temp2 c_density
	rename temp3 c_gamma

	matrix list e(V)
	svmat e(V), names("temp")
	rename temp1 vcov_const
	rename temp2 vcov_density
	rename temp3 vcov_gamma

	* save
	preserve
		keep c_* vcov_* 
		keep if _n <= 3
		export delimited using "data/coded_model/road_tech/vcov_density_pow.csv", replace
	restore
	
